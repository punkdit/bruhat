#!/usr/bin/env python3

"""
single variable polynomials
"""

import sys, os


class Poly(object):
    def __init__(self, cs):
        self.cs = {}
        for deg, coeff in cs.items():
            if coeff==0: # strip these out
                continue
            self.cs[deg] = coeff

    def __str__(self):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort(reverse=True)
        terms = []
        for deg in keys:
            coeff = cs[deg]
            assert coeff!=0
            if coeff==1 and deg==1:
                term = "x"
            elif deg==0:
                term = str(coeff)
            elif deg==1:
                term = "%d*x" % (coeff,)
            elif coeff==1:
                term = "x**%d" % deg
            else:
                term = "%d*x**%d" % (coeff, deg)
            terms.append(term)
        s = "+".join(terms)
        return s

    def __repr__(self):
        return "Poly(%s)"%self.cs

    def __eq__(self, other):
        return self.cs == other.cs

    def __ne__(self, other):
        return self.cs != other.cs

    def __hash__(self):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort()
        value = tuple((deg, cs[deg]) for deg in keys)
        return hash(value)

    def __add__(self, other):
        other = Poly.promote(other)
        cs = dict(self.cs)
        for deg, coeff in other.cs.items():
            cs[deg] = cs.get(deg, 0) + coeff
        return Poly(cs)

    def __sub__(self, other):
        other = Poly.promote(other)
        cs = dict(self.cs)
        for deg, coeff in other.cs.items():
            cs[deg] = cs.get(deg, 0) - coeff
        return Poly(cs)

    def __mul__(self, other):
        other = Poly.promote(other)
        cs = {}
        for d0, c0 in self.cs.items():
          for d1, c1 in other.cs.items():
            deg = d0+d1
            cs[deg] = cs.get(deg, 0) + c0*c1
        return Poly(cs)

    def __rmul__(self, value):
        cs = {}
        for deg, coeff in self.cs.items():
            cs[deg] = value*coeff
        return Poly(cs)

    def __radd__(self, value):
        cs = dict(self.cs)
        cs[0] = cs.get(0, 0) + value
        return Poly(cs)

    @classmethod
    def identity(cls):
        return cls({0:1})

    @classmethod
    def zero(cls):
        return cls({})

    @classmethod
    def promote(cls, value):
        if isinstance(value, Poly):
            return value
        return Poly({0:value})

    def __pow__(self, n):
        assert int(n)==n
        assert n>=0
        p = Poly.identity()
        for i in range(n):
            p = self*p
        return p


def test():

    I = Poly({0:1})
    x = Poly({1:1})

    assert (x+I)**5 == x**5+5*x**4+10*x**3+10*x**2+5*x+1


if __name__ == "__main__":

    test()


