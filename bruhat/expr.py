#!/usr/bin/env python

"""
previous version: higher.py

"""

from time import time
start_time = time()

import z3


class Expr(object):
    def __init__(self, theory, key):
        self.theory = theory
        self.key = key

    def equate(self, other):
        self.theory.equate(self.key, other.key)

    def __hash__(self):
        assert 0

    def __eq__(self, other):
        assert 0
        return self.theory.is_equal(self.key, other.key)

    def __mul__(self, other):
        theory = self.theory
        term = theory.Term("*", self, other, inline=True)
        return term

    @staticmethod
    def unify(lhs, rhs, send={}, depth=0):
        indent = "  "*depth
        #print(indent+"unify:", lhs, "&", rhs)
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
    def __init__(self, theory, name):
        self.name = name
        key = id(self)
        Expr.__init__(self, theory, key)

    def __str__(self):
        return self.name
    __repr__ = __str__

    def find(self, v):
        assert isinstance(v, Variable)
        return self is v

    def substitute(self, send):
        return send.get(self.key, self)

    def all_vars(self):
        return {self.key : self}



class Term(Expr):
    def __init__(self, theory, op, *args, inline=False, postfix=False):
        self.op = op
        self.args = args
        self.inline = inline
        self.postfix = postfix
        key = (op,) + tuple(arg.key for arg in args)
        Expr.__init__(self, theory, key)

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
            s = op # Const
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
        term = theory.Term(op, *args, inline=self.inline, postfix=self.postfix)
        return term

    def all_vars(self):
        vs = {}
        for arg in self.args:
            for v in arg.all_vars().values():
                vs[v.key] = v
        return vs




class Equation(object):
    def __init__(self, lhs, rhs):
        assert isinstance(lhs, Expr)
        assert isinstance(rhs, Expr)
        assert lhs.theory is rhs.theory
        theory = lhs.theory
        vs = lhs.all_vars() | rhs.all_vars()
        # Make some uniqe Variable's just for this Equation:
        fwd = {v.key : theory.Variable(v.name) for v in vs.values()}
        # & remember how to get back:
        self.rev = {fwd[v.key].key : v for v in vs.values()}
        lhs = lhs.substitute(fwd)
        rhs = rhs.substitute(fwd)
        self.lhs = lhs
        self.rhs = rhs

    def match(self, expr, directional=True):
        assert isinstance(expr, Expr)
        lhs, rhs = self.lhs, self.rhs
        send = Expr.unify(lhs, expr)
        if send is not None:
            rhs = rhs.substitute(send)
            rhs = rhs.substitute(self.rev)
            return rhs
        
    



class Theory(object):
    def __init__(self):
        self.cache = {}

    def Variable(self, name):
        expr = Variable(self, name)
        #self.cache[expr.key] = expr
        return expr

    def Term(self, op, *args, inline=False, postfix=False):
        cache = self.cache
        expr = Term(self, op, *args, inline=inline, postfix=postfix)
        if expr.key in cache:
            expr = cache[expr.key]
        else:
            cache[expr.key] = expr
        return expr

    def Const(self, name):
        return self.Term(name)


class SMT(Theory):
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
        #print("SMT.equate", lhs, rhs)
        assert isinstance(lhs, z3.ExprRef)
        assert isinstance(rhs, z3.ExprRef)
        self.solver.add( lhs == rhs )

    def is_equal(self, lhs, rhs):
        #print("SMT.is_equal", lhs, rhs)
        assert isinstance(lhs, z3.ExprRef)
        assert isinstance(rhs, z3.ExprRef)
        solver = self.solver
        solver.push()
        solver.add( lhs != rhs )
        result = solver.check()
        solver.pop()
        return result == z3.unsat

    def variable(self, name):
        return Variable(self, name)

    def operator(self, name, arity=2):
        return Operator(self, name, arity)



def test_solver():

    solver = SMT()
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

    solver = SMT()
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


def test_theory():

    theory = Theory()

    a, b, c, d = [theory.Variable(name) for name in 'abcd']
    one = theory.Const("1")

    assert a is not b
    assert a is not theory.Variable("a")
    assert a*b is a*b

    assert str(a*b) == "(a*b)"

    lhs, rhs = a*b, c*d
    send = Expr.unify(a*b, c*d)

    assert lhs.substitute(send) is rhs.substitute(send)

    eqn = Equation( (a*b)*c, a*(b*c) )

    lhs = (a*b)*(c*d)
    rhs = eqn.match(lhs)
    assert rhs is a*(b*(c*d))

    lhs = one*a
    rhs = one*(one*b)
    assert Expr.unify(lhs, rhs)

    eqn = Equation( one*a, a )
    lhs = one*(one*a)
    rhs = eqn.match(lhs)
    assert rhs is one*a



if __name__ == "__main__":

    print("\n\n")
    test_solver()
    test_theory()

    #test_category()
    #test_bicategory()
    #test_globular()

    print("OK: ran in %.3f seconds.\n"%(time() - start_time))







