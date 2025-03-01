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


# Variable's are universally quantified,
# Const's are existentially quantified.

class Expr(object):
    def __init__(self, theory, key, sort, is_const):
        self.theory = theory
        self.key = key
        self.sort = sort
        self.is_const = is_const

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
        if isinstance(src, Variable):
            if src.key in send and send[src.key] is not tgt:
                return None # <------------ return
            send[src.key] = tgt
            #print(indent+"match: Const")
            return send # <-------------- return
        if isinstance(tgt, Variable):
            return None # <-------------- return
        if src is tgt:
            return send
        if isinstance(src, Const) or isinstance(tgt, Const):
            return None
        assert isinstance(src, Term), (src, tgt)
        assert isinstance(tgt, Term), (src, tgt)
        if src.op != tgt.op:
            return None # fail
        for left, right in zip(src.args, tgt.args):
            send = Expr.match(left, right, send, depth=depth+1) # recurse
            if send is None:
                return None # fail
        return send


class Variable(Expr):
    def __init__(self, theory, name, sort):
        self.name = name
        self.sort = sort
        key = id(self)
        Expr.__init__(self, theory, key, sort, False)

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


class Const(Expr):
    def __init__(self, theory, name, sort):
        self.name = name
        self.sort = sort
        key = id(self)
        Expr.__init__(self, theory, key, sort, True)

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
        is_const = True
        for arg in args:
            if not arg.is_const:
                is_const = False
        Expr.__init__(self, theory, key, sort, is_const)

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

    def substitute(self, send):
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
        assert not src.is_const, 'makes no sense..'
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
        self.cache = {} # cache all constant expressions
        self.sigs = {}
        self.eqns = []
        self.lookup = {}

    def Const(self, name, sort, expr=None):
        if expr is None:
            expr = Const(self, name, sort)
        solver, lookup = self.solver, self.lookup
        v = solver.get_const(name)
        lookup[expr.key] = v
        self.info("Const:", expr)
        return expr

    def Variable(self, name, sort):
        expr = Variable(self, name, sort)
        self.info("Variable:", expr)
        return expr

    def apply_rewrite(self, eqn, expr):
        other = eqn.match(expr)
        if other is None:
            return
        assert other.key != expr.key, "continue?"
        self.info("theory.apply_rewrite:", expr, "==", other)
        self.equate(expr, other)

    def Term(self, op, *args, sort=None, inline=False, postfix=False):
        cache = self.cache
        expr = Term(self, op, *args, sort=sort, inline=inline, postfix=postfix)
        if not expr.is_const:
            return expr # <------------- return

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
            self.apply_rewrite(eqn, expr) # may recurse !

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
        if lhs.is_const:
            assert rhs.is_const
            self.equate(lhs, rhs) # now these constant's are equal
            return
        eqn = Rewrite(self, lhs, rhs)
        self.info("theory.Rewrite:", eqn)
        self.eqns.append(eqn)

        solver, lookup = self.solver, self.lookup
        for expr in list(self.cache.values()):
            assert expr.is_const, expr
            self.apply_rewrite(eqn, expr)
        return eqn

    def Equation(self, lhs, rhs):
        self.Rewrite(lhs, rhs)
        self.Rewrite(rhs, lhs)

    def equate(self, lhs, rhs):
        solver, lookup = self.solver, self.lookup
        lhs, rhs = lookup[lhs.key], lookup[rhs.key]
        solver.equate(lhs, rhs)

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

    def get_const(self, stem="v"):
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

    def get_const(self, stem="v"):
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
    a = solver.get_const("a")
    b = solver.get_const("b")
    c = solver.get_const("c")
    d = solver.get_const("d")

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
    t, u, v, w, x, y, z = [theory.Variable(name, sort) for name in 'tuvwxyz']
    #one = theory.Nullary("1", sort)
    one = theory.Const("one", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a is not b
    assert a is not theory.Const("a", sort)
    assert a*b is a*b

    assert str(a*b) == "(a*b)"

    for (src, tgt) in [
        (u, u*v),
        (u*v, v*u),
        (u*v, v*(w*w)),
        ((u*v)*w, (x*y)*(z*t)),
    ]:
        send = Expr.match( src, tgt )
        assert send
        assert src.substitute(send).key == tgt.key

    send = Expr.match( u*v, v )
    assert send is None

    send = Expr.match( (u*v)*w, u*v )
    assert send is None

    eqn = Rewrite(theory, (u*v)*w, u*(v*w))

    lhs = (a*b)*(c*d)
    rhs = eqn.match(lhs)
    assert rhs is a*(b*(c*d)), rhs

    lhs = one*u
    rhs = one*(one*v)
    assert Expr.match(lhs, rhs)

    eqn = Rewrite(theory, one*u, u)
    lhs = one*(one*a)
    rhs = eqn.match(lhs)
    assert rhs is one*a

    _sort = Sort("_sort")
    _a = theory.Const("_a", _sort)
    assert Expr.match(u, _a) is None


def test_theory():

    solver = Z3Solver()
    #solver.DEBUG = True
    theory = Theory(solver)
    sort = Sort("monoid")

    a, b = [theory.Const(name, sort) for name in 'ab']
    x, y = [theory.Variable(name, sort) for name in 'xy']
    #one = theory.Nullary("one", sort)
    one = theory.Const("one", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a*b != b*a

    #theory.DEBUG = True
    theory.Rewrite( one*x, x )
    theory.Rewrite( x*one, x )

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
    u, v, w = [theory.Variable(name, sort) for name in 'uvw']
    one = theory.Const("one", sort)
    theory.Operator("*", sort, [sort, sort], inline=True)

    assert a is not b
    assert a is not theory.Const("a", sort)
    assert a*b is a*b

    assert str(a*b) == "(a*b)"

    assert (a*b)*c != a*(b*c)

    theory.Equation( (u*v)*w, u*(v*w) )
    assert (a*b)*c == a*(b*c)

    lhs = ((a*b)*c)*d
    assert lhs == (a*b)*(c*d)
    assert lhs == a*(b*(c*d))
    assert lhs == a*((b*c)*d)
    assert lhs == (a*(b*c))*d

    theory.Rewrite(one*u, u)
    theory.Rewrite(u*one, u)

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
    Variable = theory.Variable
   
    Operator("identity", cell1, [cell0], postfix=True)
    Operator("*", cell1, [cell1, cell1], inline=True)

    X, Y, Z, W = [Const(name, cell0) for name in "XYZW"]
    f, g, h = [Const(name, cell1) for name in "fgh"]

    vX, vY, vZ, vW = [Variable(name, cell0) for name in "XYZW"]
    vf, vg, vh = [Variable(name, cell1) for name in "fgh"]

    assert X.identity is X.identity
    assert X.identity == X.identity
    assert X.identity != Y.identity


    #      h      g      f    
    #  W <--- Z <--- Y <--- X
    # 

    assert g*f != h*g
    theory.Equation( (vh*vg)*vf, vh*(vg*vf) )
    assert g*f != h*g

    for cell in [X, Y, Z, W]:
        theory.Rewrite( vf*cell.identity, vf )
        theory.Rewrite( cell.identity*vf, vf )

    assert g*f != h*g
    assert g*f != h*g
    assert distinct([X, Y, Z, W])
    assert distinct([X.identity, Y.identity, Z.identity, W.identity])

    assert g*f != h*g

    assert (h*g)*Y.identity == h*g
    assert g*f != h*g

    assert X.identity * X.identity == X.identity

    # add an isomorphism ... ?
    i = Const("i", cell1)
    j = Const("j", cell1)
    theory.Rewrite(i*j, X.identity)
    theory.Rewrite(j*i, X.identity)


    # Here's the problem. We want to do the following chain:
    #
    #           [1]           [2]                  [3]
    #  (f*i)*j ----> f*(i*j) ----> f*X.identity -------> f
    #  
    # but only z3 knows about [2]. We haven't triggered the
    # Rewrite for [3] so we never tell z3 about it.
    #  

    #assert f*i*j == f*X.identity # uncomment this to tell z3 about [3].
    assert f*i*j == f # FAIL


if __name__ == "__main__":

    print("\n\n")
    test_solver()
    test_base_theory()
    test_theory()
    test_monoid_theory()
    test_category_theory()

    print("OK: ran in %.3f seconds.\n"%(time() - start_time))







