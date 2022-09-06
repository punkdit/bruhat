#!/usr/bin/env python

"""
previous version: higher.py

"""

from time import time
start_time = time()

import z3


class Sort(object):
    def __init__(self, name):
        self.name = name
    def __str__(self):
        return self.name


class Expr(object):
    def __init__(self, theory, key, sort):
        self.theory = theory
        self.key = key
        self.sort = sort

    def equate(self, other):
        self.theory.equate(self.key, other.key)

    def __hash__(self):
        assert 0

    def __eq__(self, other):
        assert 0
        return self.theory.is_equal(self.key, other.key)

    def __mul__(self, other):
        theory = self.theory
        op = theory.get_op("*", self.sort, other.sort)
        term = op(self, other)
        return term

    @staticmethod
    def match(src, tgt, send=None, depth=0):
        if send is None:
            send = {}
        indent = "  "*depth
        #print(indent+"match:", src, "->", tgt, send)
        if src.sort != tgt.sort:
            return None
        if send:
            src = src.substitute(send) # ?
        if isinstance(src, Variable):
            if src is not tgt and tgt.find(src):
                return None
            if src.key in send and send[src.key] is not tgt:
                return None # <------------ return
            send[src.key] = tgt
            #print(indent+"match: Variable")
            return send # <-------------- return
        if isinstance(tgt, Variable):
            return None # <-------------- return
        assert isinstance(src, Term)
        assert isinstance(tgt, Term)
        if src.op != tgt.op:
            return None # fail
        for left, right in zip(src.args, tgt.args):
            send = Expr.match(left, right, send, depth=depth+1) # recurse
            if send is None:
                return None # fail
        return send

    @staticmethod
    def unify(lhs, rhs, send=None, depth=0):
        indent = "  "*depth
        #print(indent+"unify:", lhs, "&", rhs)
        if lhs.sort != rhs.sort:
            return None
        if send is None:
            send = {}
        if send:
            lhs = lhs.substitute(send)
            rhs = rhs.substitute(send)
        if isinstance(rhs, Variable):
            lhs, rhs = rhs, lhs
        if isinstance(lhs, Variable):
            if lhs is not rhs and rhs.find(lhs):
                return None # fail (recursive)
            #print(indent, "<== {%s : %s}"%(lhs, rhs))
            assert lhs.key not in send
            if lhs is not rhs:
                send = dict(send)
                for k,v in list(send.items()):
                    send[k] = v.substitute({lhs.key:rhs})
                send[lhs.key] = rhs
            return send
        assert isinstance(lhs, Term)
        assert isinstance(rhs, Term)
        if lhs.op != rhs.op:
            return None # fail
        for left, right in zip(lhs.args, rhs.args):
            send = Expr.unify(left, right, send, depth=depth+1)
            if send is None:
                return None # fail
        return send



class Variable(Expr):
    def __init__(self, theory, name, sort):
        self.name = name
        self.sort = sort
        key = id(self)
        Expr.__init__(self, theory, key, sort)

    def __str__(self):
        return self.name
        #return "%s:%s"%(self.name, self.sort)
    __repr__ = __str__

    def find(self, v):
        assert isinstance(v, Variable)
        return self is v

    def substitute(self, send):
        return send.get(self.key, self)

    def all_vars(self):
        return {self.key : self}


class Operator(object):
    def __init__(self, theory, name, rsort, argsorts, inline=False, postfix=False):
        self.theory = theory
        self.name = name
        self.rsort = rsort
        self.argsorts = argsorts
        self.inline = inline
        self.postfix = postfix

    def __str__(self):
        return self.name

    def __call__(self, *args):
        theory = self.theory
        argsorts = self.argsorts
        assert len(args) == len(argsorts)
        for (arg, sort) in zip(args, argsorts):
            assert arg.sort == sort
        term = theory.Term(self, *args, 
            sort=self.rsort, inline=self.inline, postfix=self.postfix)
        return term


class Term(Expr):
    def __init__(self, theory, op, *args, sort=None, inline=False, postfix=False):
        assert isinstance(op, Operator)
        self.op = op
        self.args = args
        self.inline = inline
        self.postfix = postfix
        key = (op,) + tuple(arg.key for arg in args)
        Expr.__init__(self, theory, key, sort)

    def __str__(self):
        args = self.args
        op = self.op
        if self.inline:
            assert len(args) == 2
            l, r = args
            s = "(%s%s%s)" % (l, op, r)
        elif self.postfix:
            assert len(args) == 1
            l, = args
            s = "%s.%s"%(l, op)
        elif len(args) == 0:
            s = str(op) # Const
        else:
            s = "%s%s"%(op, args)
        return s
    __repr__ = __str__

    def find(self, v):
        assert isinstance(v, Variable)
        for arg in self.args:
            if arg.find(v):
                return True
        return False

    def substitute(self, send): # hotspot
        op = self.op
        args = [arg.substitute(send) for arg in self.args]
        theory = self.theory
        term = theory.Term(op, *args, 
            sort=self.sort, inline=self.inline, postfix=self.postfix)
        return term

    def all_vars(self):
        vs = {}
        for arg in self.args:
            for v in arg.all_vars().values():
                vs[v.key] = v
        return vs




class Rewrite(object):
    def __init__(self, theory, lhs, rhs):
        assert isinstance(lhs, Expr)
        assert isinstance(rhs, Expr)
        assert theory is lhs.theory
        assert theory is rhs.theory
        assert lhs.sort == rhs.sort
        theory = lhs.theory
        vs = lhs.all_vars() | rhs.all_vars()
        # Make some uniqe Variable's just for this Rewrite:
        fwd = {v.key : theory.Variable("_"+v.name, v.sort) for v in vs.values()}
        # & remember how to get back:
        self.rev = {fwd[v.key].key : v for v in vs.values()}
        lhs = lhs.substitute(fwd)
        rhs = rhs.substitute(fwd)
        self.lhs = lhs
        self.rhs = rhs

    def __str__(self):
        return "%s -> %s" % (self.lhs, self.rhs)

    def match(self, expr, directional=True):
        assert isinstance(expr, Expr)
        lhs, rhs = self.lhs, self.rhs
        #print("Rewrite.match", lhs, "-->", expr)
        send = Expr.match(lhs, expr)
        #print("Rewrite.match =", send)
        if send is not None:
            rhs = rhs.substitute(send)
            rhs = rhs.substitute(self.rev)
            return rhs

        
    



class Theory(object):
    def __init__(self, solver=None):
        if solver is None:
            solver = Solver()
        self.solver = solver
        self.cache = {}
        self.sigs = {}
        self.eqns = []
        self.lookup = {}

    def Variable(self, name, sort):
        expr = Variable(self, name, sort)
        #self.cache[expr.key] = expr # ?
        solver, lookup = self.solver, self.lookup
        v = solver.get_var()
        lookup[expr.key] = v
        return expr

    def Term(self, op, *args, sort=None, inline=False, postfix=False):
        cache = self.cache
        expr = Term(self, op, *args, sort=sort, inline=inline, postfix=postfix)
        if expr.key in cache:
            return cache[expr.key] # <------------- return

        # we made a new Term
        cache[expr.key] = expr
        solver, lookup = self.solver, self.lookup
        op = solver.get_op(op.name, len(args))
        lhs = op(*[lookup[arg.key] for arg in args])
        lookup[expr.key] = lhs

        print("Term:", expr)
        for eqn in self.eqns:
            other = eqn.match(expr)
            if other is None:
                continue
            assert other.key != expr.key, "continue?"
            #print("\t", eqn)
            #print("\t\t", other)
            print("\t", expr, "==", other)
            rhs = lookup[other.key]
            solver.equate(lhs, rhs)

        return expr

    def Const(self, name, sort):
        #return self.Term(name, sort=sort)
        op = self.Operator(name, sort, ())
        return op()

    def Operator(self, name, rsort, argsorts, inline=False, postfix=False):
        sig = (name, tuple(argsorts))
        assert sig not in self.sigs
        op = Operator(self, name, rsort, argsorts, inline=inline, postfix=postfix)
        self.sigs[sig] = op
        return op

    def get_op(self, name, *argsorts):
        sig = (name, argsorts)
        return self.sigs[sig]

    def Rewrite(self, lhs, rhs):
        eq = Rewrite(self, lhs, rhs)
        self.eqns.append(eq)
        return eq

    def Equation(self, lhs, rhs):
        self.Rewrite(lhs, rhs)
        self.Rewrite(rhs, lhs)


class Solver(object):
    def __init__(self):
        pass

    def get_var(self, stem="v"):
        pass

    def get_op(self, name, arity=2):
        op = lambda *args : None
        return op

    def equate(self, lhs, rhs):
        pass

    def is_equal(self, lhs, rhs):
        return False


class Z3Solver(Solver):
    def __init__(self):
        self.solver = z3.Solver()
        self.sort = z3.DeclareSort("THE_SORT")
        self.v_count = 0
        self.oplookup = {}

    def get_var(self, stem="v"):
        name = stem+str(self.v_count)
        self.v_count += 1
        v = z3.Const(name, self.sort)
        #print(v)
        return v

    def get_op(self, name, arity=2):
        if name in self.oplookup:
            f = self.oplookup[name]
            assert f.arity() == arity
        else:
            f = z3.Function(name, [self.sort]*(arity+1))
            self.oplookup[name] = f
        return f

    def equate(self, lhs, rhs):
        print("Z3Solver.equate", lhs, rhs)
        assert len(str(lhs)) < 20
        assert isinstance(lhs, z3.ExprRef)
        assert isinstance(rhs, z3.ExprRef)
        self.solver.add( lhs == rhs )

    def is_equal(self, lhs, rhs):
        #print("Z3Solver.is_equal", lhs, rhs)
        assert isinstance(lhs, z3.ExprRef)
        assert isinstance(rhs, z3.ExprRef)
        solver = self.solver
        solver.push()
        solver.add( lhs != rhs )
        result = solver.check()
        solver.pop()
        return result == z3.unsat



def test_solver():

    solver = Z3Solver()
    a = solver.get_var("a")
    b = solver.get_var("b")
    c = solver.get_var("c")
    d = solver.get_var("d")

    solver.equate(a, b)
    assert solver.is_equal(a, b)
    solver.equate(b, c)
    assert solver.is_equal(a, c)
    assert not solver.is_equal(a, d)

    mul = solver.get_op("*", 2)

    lhs = mul(mul(a, b), c)
    rhs = mul(a, mul(b, c))
    solver.equate( lhs, rhs )

    assert solver.is_equal( mul(lhs, d), mul(rhs, d) )
    assert not solver.is_equal( mul(lhs, d), mul(d, rhs) )


def __test_expr():

    solver = Z3Solver()
    equate = lambda lhs, rhs : lhs.equate(rhs)

    a, b, c, d = [solver.variable(name) for name in 'abcd']

    assert str(a*a) == "(a*a)"

    assert a*(b*c) != (a*b)*c
    assert b*(c*d) != (b*c)*d

    equate( a*(b*c), (a*b)*c )
    assert a*(b*c) == (a*b)*c

    assert b*(c*d) != (b*c)*d

    equate( b*(c*d), (b*c)*d )
    assert b*(c*d) == (b*c)*d


def test_base_theory():

    theory = Theory()
    sort = Sort("sort")

    a, b, c, d, e, f, g = [theory.Variable(name, sort) for name in 'abcdefg']
    one = theory.Const("1", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a is not b
    assert a is not theory.Variable("a", sort)
    assert a*b is a*b

    assert str(a*b) == "(a*b)"

    lhs, rhs = a*b, c*d
    send = Expr.unify(a*b, c*d)
    assert lhs.substitute(send) is rhs.substitute(send)

    lhs, rhs = (a*b)*c, a*(b*c)
    send = Expr.match( (a*b)*c, d*e )
    assert send is None
    send = Expr.match( d*e, (a*b)*c )
    assert send
    assert (d*e).substitute(send) is (a*b)*c

    send = Expr.match( (a*b)*c, (d*e)*(f*g) )
    assert send
    assert ((a*b)*c).substitute(send) is (d*e)*(f*g)


    eqn = Rewrite(theory, (a*b)*c, a*(b*c))

    lhs = (a*b)*(c*d)
    rhs = eqn.match(lhs)
    assert rhs is a*(b*(c*d))

    lhs = one*a
    rhs = one*(one*b)
    assert Expr.match(lhs, rhs)

    eqn = Rewrite(theory, one*a, a)
    lhs = one*(one*a)
    rhs = eqn.match(lhs)
    assert rhs is one*a

    _sort = Sort("_sort")
    _a = theory.Variable("_a", _sort)
    assert Expr.match(a, _a) is None


def test_theory():

    solver = Z3Solver()
    theory = Theory(solver)
    sort = Sort("monoid")

    a, b, c, d = [theory.Variable(name, sort) for name in 'abcd']
    one = theory.Const("1", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a is not b
    assert a is not theory.Variable("a", sort)
    assert a*b is a*b

    assert str(a*b) == "(a*b)"

    theory.Equation( (a*b)*c, a*(b*c) )
    #theory.Rewrite( one*a, a )
    #theory.Rewrite( a*one, a )

    assert a*one == a


if __name__ == "__main__":

    print("\n\n")
    test_solver()
    test_base_theory()
    test_theory()

    #test_category()
    #test_bicategory()
    #test_globular()

    print("OK: ran in %.3f seconds.\n"%(time() - start_time))







