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
    def __repr__(self):
        return "Sort(%r)"%(self.name,)


class Expr(object):
    def __init__(self, theory, key, sort):
        self.theory = theory
        self.key = key
        self.sort = sort

    def equate(self, other):
        assert 0, "use theory.Equation etc."
        self.theory.equate(self, other)

    def __hash__(self):
        assert 0

    def __eq__(self, other):
        assert isinstance(other, Expr), other
        return self.theory.is_equal(self, other)

    def __mul__(self, other):
        assert isinstance(other, Expr), other
        theory = self.theory
        op = theory.get_operator("*", self.sort, other.sort)
        if op is None:
            raise AttributeError("operator %s*%s not found"%(self.sort, other.sort))
        term = op(self, other)
        return term

    def __lshift__(self, other):
        assert isinstance(other, Expr), other
        theory = self.theory
        op = theory.get_operator("<<", self.sort, other.sort)
        if op is None:
            raise AttributeError("operator %s<<%s not found"%(self.sort, other.sort))
        term = op(self, other)
        return term

    def __matmul__(self, other):
        assert isinstance(other, Expr), other
        theory = self.theory
        op = theory.get_operator("@", self.sort, other.sort)
        if op is None:
            raise AttributeError("operator %s@%s not found"%(self.sort, other.sort))
        term = op(self, other)
        return term

    def __getattr__(self, attr):
        theory = self.theory
        op = theory.get_operator(attr, self.sort)
        if op is None:
            raise AttributeError("attr %r not found"%(attr,))
        return op(self)

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
        if isinstance(src, Const):
            if src.key in send and send[src.key] is not tgt:
                return None # <------------ return
            send[src.key] = tgt
            #print(indent+"match: Const")
            return send # <-------------- return
        if isinstance(tgt, Const):
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

#    @staticmethod
#    def unify(lhs, rhs, send=None, depth=0):
#        indent = "  "*depth
#        #print(indent+"unify:", lhs, "&", rhs)
#        if lhs.sort != rhs.sort:
#            return None
#        if send is None:
#            send = {}
#        if send:
#            lhs = lhs.substitute(send)
#            rhs = rhs.substitute(send)
#        if isinstance(rhs, Const):
#            lhs, rhs = rhs, lhs
#        if isinstance(lhs, Const):
#            if lhs is not rhs and rhs.find(lhs):
#                return None # fail (recursive)
#            #print(indent, "<== {%s : %s}"%(lhs, rhs))
#            assert lhs.key not in send
#            if lhs is not rhs:
#                send = dict(send)
#                for k,v in list(send.items()):
#                    send[k] = v.substitute({lhs.key:rhs})
#                send[lhs.key] = rhs
#            return send
#        assert isinstance(lhs, Term)
#        assert isinstance(rhs, Term)
#        if lhs.op != rhs.op:
#            return None # fail
#        for left, right in zip(lhs.args, rhs.args):
#            send = Expr.unify(left, right, send, depth=depth+1)
#            if send is None:
#                return None # fail
#        return send



class Const(Expr):
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
        assert isinstance(v, Const)
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

    def __repr__(self):
        return "%s(%r)"%(self.__class__.__name__, self.name)

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
            s = str(op) # Nullary
        else:
            s = "%s%s"%(op, args)
        return s
    __repr__ = __str__

    def __getitem__(self, idx):
        return self.args[idx]

    def find(self, v):
        assert isinstance(v, Const)
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
    def __init__(self, theory, src, tgt):
        assert isinstance(src, Expr)
        assert isinstance(tgt, Expr)
        assert theory is src.theory
        assert theory is tgt.theory
        assert src.sort == tgt.sort
        self.src = src
        self.tgt = tgt

    def __str__(self):
        return "%s -> %s" % (self.src, self.tgt)

    def match(self, expr, directional=True):
        assert isinstance(expr, Expr)
        src, tgt = self.src, self.tgt
        #print("Rewrite.match", src, "-->", expr)
        send = Expr.match(src, expr)
        #print("Rewrite.match =", send)
        if send is not None:
            tgt = tgt.substitute(send)
            return tgt


class Theory(object):
    DEBUG = False
    def info(self, *msg):
        if self.DEBUG:
            print(*msg)

    def __init__(self, solver=None):
        if solver is None:
            solver = Solver()
        self.solver = solver
        self.cache = {}
        self.sigs = {}
        self.eqns = []
        self.lookup = {}

    def Const(self, name, sort, expr=None):
        if expr is None:
            expr = Const(self, name, sort)
        solver, lookup = self.solver, self.lookup
        v = solver.get_var(name)
        lookup[expr.key] = v
        self.info("Const:", expr)
        return expr

    def apply_rewrite(self, eqn, expr):
        other = eqn.match(expr)
        if other is None:
            return
        assert other.key != expr.key, "continue?"
        self.info("theory.apply_rewrite:", expr, "==", other)
        solver, lookup = self.solver, self.lookup
        lhs, rhs = lookup[expr.key], lookup[other.key]
        #print("solver.equate", lhs, rhs)
        solver.equate(lhs, rhs)

    def Term(self, op, *args, sort=None, inline=False, postfix=False):
        cache = self.cache
        expr = Term(self, op, *args, sort=sort, inline=inline, postfix=postfix)
        if expr.key in cache:
            return cache[expr.key] # <------------- return

        # we made a new Term
        cache[expr.key] = expr
        solver, lookup = self.solver, self.lookup
        op = solver.get_operator(op.name, len(args))
        lhs = op(*[lookup[arg.key] for arg in args])
        lookup[expr.key] = lhs

        #print("Term:", expr)
        for eqn in self.eqns:
            self.apply_rewrite(eqn, expr)

        return expr

    def add_operator(self, op):
        sig = (op.name, tuple(op.argsorts))
        assert sig not in self.sigs
        self.sigs[sig] = op

    def Operator(self, name, rsort, argsorts, inline=False, postfix=False):
        op = Operator(self, name, rsort, argsorts, inline=inline, postfix=postfix)
        if argsorts:
            self.add_operator(op)
        return op

    def get_operator(self, name, *argsorts):
        sig = (name, argsorts)
        return self.sigs.get(sig, None)

    def Nullary(self, name, sort):
        op = self.Operator(name, sort, ())
        expr = op()
        return expr

    def Rewrite(self, lhs, rhs):
        eqn = Rewrite(self, lhs, rhs)
        self.info("theory.Rewrite:", eqn)
        self.eqns.append(eqn)
        solver, lookup = self.solver, self.lookup
        lhs = lookup[lhs.key]
        for expr in list(self.cache.values()):
            self.apply_rewrite(eqn, expr)
        return eqn

    def Equation(self, lhs, rhs):
        self.Rewrite(lhs, rhs)
        self.Rewrite(rhs, lhs)

    def equate(self, lhs, rhs):
        assert 0

    def is_equal(self, lhs, rhs):
        solver, lookup = self.solver, self.lookup
        lhs, rhs = lookup[lhs.key], lookup[rhs.key]
        return solver.is_equal(lhs, rhs)


class Solver(object):
    DEBUG = False
    def info(self, *msg):
        if self.DEBUG:
            print(*msg)

    def __init__(self):
        pass

    def get_var(self, stem="v"):
        pass

    def get_operator(self, name, arity=2):
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

    def get_operator(self, name, arity=2):
        if name in self.oplookup:
            f = self.oplookup[name]
            assert f.arity() == arity
        else:
            f = z3.Function(name, [self.sort]*(arity+1))
            self.oplookup[name] = f
        return f

    def equate(self, lhs, rhs):
        self.info("Z3Solver.equate:", lhs, "==", rhs)
        #assert len(str(lhs)) < 30
        assert isinstance(lhs, z3.ExprRef)
        assert isinstance(rhs, z3.ExprRef)
        self.solver.add( lhs == rhs )

    def is_equal(self, lhs, rhs):
        assert isinstance(lhs, z3.ExprRef)
        assert isinstance(rhs, z3.ExprRef)
        solver = self.solver
        solver.push()
        solver.add( lhs != rhs )
        result = solver.check()
        solver.pop()
        self.info("Z3Solver.is_equal:", lhs, "==", rhs, "?", 
            result==z3.unsat)
        return result == z3.unsat


def distinct(items):
    n = len(items)
    for i in range(n):
      for j in range(i+1, n):
        if items[i] == items[j]:
            return False
    return True
        

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

    mul = solver.get_operator("*", 2)

    lhs = mul(mul(a, b), c)
    rhs = mul(a, mul(b, c))
    solver.equate( lhs, rhs )

    assert solver.is_equal( mul(lhs, d), mul(rhs, d) )
    assert not solver.is_equal( mul(lhs, d), mul(d, rhs) )


def test_base_theory():

    theory = Theory()
    sort = Sort("sort")

    a, b, c, d, e, f, g = [theory.Const(name, sort) for name in 'abcdefg']
    one = theory.Nullary("1", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a is not b
    assert a is not theory.Const("a", sort)
    assert a*b is a*b

    assert str(a*b) == "(a*b)"

    for (src, tgt) in [
        (a, a*b),
        (a*b, b*a),
        (a*b, b*(c*d)),
    ]:
        send = Expr.match( src, tgt )
        assert send
        assert src.substitute(send) is tgt

    send = Expr.match( a*b, a )
    assert send is None

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
    assert rhs is a*(b*(c*d)), rhs

    lhs = one*a
    rhs = one*(one*b)
    assert Expr.match(lhs, rhs)

    eqn = Rewrite(theory, one*a, a)
    lhs = one*(one*a)
    rhs = eqn.match(lhs)
    assert rhs is one*a

    _sort = Sort("_sort")
    _a = theory.Const("_a", _sort)
    assert Expr.match(a, _a) is None


def test_theory():

    solver = Z3Solver()
    #solver.DEBUG = True
    theory = Theory(solver)
    sort = Sort("monoid")

    a, b = [theory.Const(name, sort) for name in 'ab']
    one = theory.Nullary("one", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a*b != b*a

    #theory.DEBUG = True
    theory.Rewrite( one*a, a )
    theory.Rewrite( a*one, a )

    assert a*one == a
    assert a*one*one*a == a*a
    assert a*b != b*a
    for u in [one, a, b, a*b, b*a]:
      for v in [one, a, b, a*b, b*a]:
        assert u*one*v == u*v


def test_monoid_theory():

    solver = Z3Solver()
    theory = Theory(solver)
    sort = Sort("monoid")

    a, b, c, d = [theory.Const(name, sort) for name in 'abcd']
    one = theory.Nullary("one", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a is not b
    assert a is not theory.Const("a", sort)
    assert a*b is a*b

    assert str(a*b) == "(a*b)"

    assert (a*b)*c != a*(b*c)

    theory.Equation( (a*b)*c, a*(b*c) )
    assert (a*b)*c == a*(b*c)

    lhs = ((a*b)*c)*d
    assert lhs == (a*b)*(c*d)
    assert lhs == a*(b*(c*d))
    assert lhs == a*((b*c)*d)
    assert lhs == (a*(b*c))*d

    for u in [a, b, c, d]:
        theory.Rewrite( one*u, u )
        theory.Rewrite( u*one, u )

    assert a*one == a

    assert a*one*one*a == a*a
    assert a*b != b*a

    assert distinct([a, b, c, d])


def test_category_theory():

    solver = Z3Solver()
    theory = Theory(solver)

    cell0 = Sort("cell0") # object's
    cell1 = Sort("cell1") # morphism's

    Operator = theory.Operator
    Const = theory.Const
   
    Operator("identity", cell1, [cell0], postfix=True)
    Operator("*", cell1, [cell1, cell1], inline=True)

    X, Y, Z, W = [Const(name, cell0) for name in "XYZW"]
    f, g, h = [Const(name, cell1) for name in "fgh"]

    assert X.identity is X.identity
    assert X.identity == X.identity
    assert X.identity != Y.identity


    #      h      g      f    
    #  W <--- Z <--- Y <--- X
    # 

    assert g*f != h*g
    theory.Equation( (h*g)*f, h*(g*f) )
    assert g*f != h*g

    def morphism(f, X, Y):
        theory.Rewrite( f*X.identity, f )
        theory.Rewrite( Y.identity*f, f )

    morphism(f, X, Y)
    assert g*f != h*g
    morphism(g, Y, Z)
    assert g*f != h*g
    morphism(h, Z, W)
    assert distinct([X, Y, Z, W])
    assert distinct([X.identity, Y.identity, Z.identity, W.identity])

    assert g*f != h*g

    assert (h*g)*Y.identity == h*g
    assert g*f != h*g



if __name__ == "__main__":

    print("\n\n")
    test_solver()
    test_base_theory()
    test_theory()
    test_monoid_theory()
    test_category_theory()

    print("OK: ran in %.3f seconds.\n"%(time() - start_time))







