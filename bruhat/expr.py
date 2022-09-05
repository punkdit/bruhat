#!/usr/bin/env python

"""
previous version: higher.py

"""

from time import time
start_time = time()

import z3


class Expr(object):
    def __init__(self, solver):
        self.solver = solver

    def equate(self, other):
        self.solver.equate(self, other)

    def __mul__(self, other):
        operator = self.solver.get_op("*", 2)
        expr = Term(operator, (self, other), inline=True)
        return expr


class Variable(Expr):
    def __init__(self, solver, name):
        Expr.__init__(self, solver)
        self.name = name
        self.v = solver.get_var()

    def __str__(self):
        return self.name
    __repr__ = __str__


class Term(Expr):
    def __init__(self, operator, args, inline=False, postfix=False):
        self.operator = operator
        assert operator.arity() == len(args)
        self.args = args
        self.inline = inline
        self.postfix = postfix

    def __str__(self):
        args = self.args
        operator = self.operator
        if self.inline:
            assert len(args) == 2
            l, r = args
            s = "(%s%s%s)" % (l, operator, r)
        elif self.postfix:
            assert len(args) == 1
            l, = args
            s = "%s.%s"%(l, operator)
        else:
            s = "%s%s"%(operator, args)
        return s
    __repr__ = __str__


#class Operator(object):
#    def __init__(self, solver, name, arity=2):
#        self.op = solver.get_op(name, arity)
    



class Solver(object):
    def __init__(self):
        self.solver = z3.Solver()
        self.sort = z3.DeclareSort("THE_SORT")
        self.v_count = 0
        self.oplookup = {}

    def get_var(self, stem="v"):
        name = stem+str(self.v_count)
        self.v_count += 1
        v = z3.Const(name, self.sort)
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
        self.solver.add( lhs == rhs )

    def is_equal(self, lhs, rhs):
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

    solver = Solver()
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


def test_expr():

    solver = Solver()

    a = solver.variable("a")

    print(a*a)


if __name__ == "__main__":

    print("\n\n")
    test_solver()
    test_expr()

    #test_category()
    #test_bicategory()
    #test_globular()

    print("OK: ran in %.3f seconds.\n"%(time() - start_time))







