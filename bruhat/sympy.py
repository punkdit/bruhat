#!/usr/bin/env python3
"""
a simple (faster) replacement for sympy
"""

from time import time
start_time = time()
from functools import reduce
from operator import mul, add

import numpy

from bruhat.argv import argv


class Expr(object):
    @classmethod
    def promote(cls, item):
        if isinstance(item, Expr):
            return item
        if isinstance(item, (float, int, numpy.number)):
            return Const(item)
        if isinstance(item, str):
            return Symbol(item)
        raise ValueError("whats this %r"%type(item))
    def diff(self, v):
        assert isinstance(v, str), repr(v)
        print("diff:", self, v)
        assert 0
    def subs(self, values):
        print("subs", self, values)
        assert 0
    def __repr__(self):
        return self.__str__()
    def __eq__(self, other):
        return str(self) == str(other)
    def __add__(self, other):
        other = Expr.promote(other)
        if self.is_zero():
            return other
        if other.is_zero():
            return self
        if isinstance(self, Add):
            lhs = self.exprs
        else:
            lhs = (self,)
        if isinstance(other, Add):
            rhs = other.exprs
        else:
            rhs = (other,)
        exprs = lhs + rhs
        return Add(*exprs)
    def __radd__(self, other):
        other = Expr.promote(other)
        return other + self
    def __sub__(self, other):
        other = Expr.promote(other)
        #if isinstance(other, Neg):
        #    return self + other.val
        return self + (-1*other)
    def __neg__(self):
        if self.is_zero():
            return self
        return -1*self
    def __rsub__(self, other):
        other = Expr.promote(other)
        return other + -1*self

#    def __sub__(self, other):
#        other = Expr.promote(other)
#        if isinstance(other, Neg):
#            return self + other.val
#        return Sub(self, other)
#    def __rsub__(self, other):
#        other = Expr.promote(other)
#        return Sub(other, self)
#    def __neg__(self):
#        if self.is_zero():
#            return self
#        return Neg(self)
    def __mul__(self, other):
        other = Expr.promote(other)
        if other.is_zero() or self.is_zero():
            return Const(0.)
        if other.is_one():
            return self
        if self.is_one():
            return other
        if isinstance(self, Mul):
            lhs = self.exprs
        else:
            lhs = (self,)
        if isinstance(other, Mul):
            rhs = other.exprs
        else:
            rhs = (other,)
        exprs = lhs + rhs
        return Mul(*exprs)
    def __rmul__(self, other):
        other = Expr.promote(other)
        return other * self
    def is_zero(self):
        return False
    def is_one(self):
        return False

class Const(Expr):
    def __init__(self, val=0.):
        self.val = val
    def __str__(self):
        return str(self.val)
    def is_zero(self):
        return self.val == 0
    def is_one(self):
        return self.val == 1
    def diff(self, v):
        return Const()
    def subs(self, values):
        return self.val
Zero = Const
One = lambda : Const(1)

class Symbol(Expr):
    def __init__(self, name):
        self.name = str(name)
    def __str__(self):
        return self.name
    def diff(self, v):
        if v == self.name:
            return One()
        return Zero()
    def subs(self, values):
        return values[self.name]

#class Neg(Expr):
#    def __init__(self, val):
#        val = self.promote(val)
#        self.val = val
#    def __str__(self):
#        return "(-%s)"%(self.val,)
#    def __neg__(self):
#        return self.val


class MultiOp(Expr):
    op = None
    def __init__(self, *exprs):
        self.exprs = exprs
    def __str__(self):
        #return "(%s %s %s)"%(self.lhs, self.op, self.rhs)
        return "(%s)"%(self.op.join(str(expr) for expr in self.exprs))

class Add(MultiOp):
    op = "+"
    def diff(self, v):
        return reduce(add, [e.diff(v) for e in self.exprs])
    def subs(self, values):
        return reduce(add, [e.subs(values) for e in self.exprs])
    

#class Sub(MultiOp):
#    op = "-"

class Mul(MultiOp):
    op = "*"
    def diff(self, v):
        exprs = self.exprs
        n = len(exprs)
        rows = []
        #print("diff", self, v)
        for i in range(n):
            row = [ (exprs[i].diff(v) if i==j else exprs[j]) for j in range(n) ]
            row = reduce(mul, row)
            #print("\trow:", row)
            rows.append(row)
        result = reduce(add, rows)
        #print("result:", result)
        return result
    def subs(self, values):
        return reduce(mul, [e.subs(values) for e in self.exprs])


def test():
    a, b, c, d, e = [Symbol(c) for c in 'abcde']

    assert a.diff("a") == One()
    assert (a+b).diff("a") == One()
    assert (a*b).diff("a") == b
    x = (a*b*c*d + e*d*c)
    assert x.diff("a") == b*c*d
    

if __name__ == "__main__":

    test()



