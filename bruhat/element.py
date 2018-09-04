#!/usr/bin/env python3

"""
"""

import sys, os

from argv import argv


class Type(object):
    pass


class Ring(Type):

    def neg(self, a):
        zero = self.zero()
        a = self.sub(zero, a)
        return a

    def hash(self):
        return hash(self.__class__)

    def __eq__(self, other):
        assert isinstance(other, Type)
        return self.__class__ == other.__class__ # ?

    def __ne__(self, other):
        assert isinstance(other, Type)
        return self.__class__ != other.__class__ # ?



class Element(object):
    def __init__(self, tp):
        assert isinstance(tp, Type)
        self.tp = tp
        #self.promote = tp.promote # ??

    def __add__(self, other):
        tp = self.tp
        other = tp.promote(other)
        a = tp.add(self, other)
        return a

    def __sub__(self, other):
        tp = self.tp
        other = tp.promote(other)
        a = tp.sub(self, other)
        return a

    def __mul__(self, other):
        tp = self.tp
        other = tp.promote(other)
        a = tp.mul(self, other)
        return a

    def __rmul__(self, value):
        tp = self.tp
        other = tp.promote(value)
        a = tp.mul(other, self)
        return a

    def __radd__(self, value):
        tp = self.tp
        other = tp.promote(value)
        a = tp.add(other, self)
        return a

    def __neg__(self):
        tp = self.tp
        a = tp.neg(self)
        return a

    def __div__(self, other):
        tp = self.tp
        other = tp.promote(other)
        a = tp.div(self, other)
        return a

    def __rdiv__(self, other):
        tp = self.tp
        other = tp.promote(other)
        a = tp.div(other, self)
        return a

    def __pow__(self, n):
        assert int(n)==n
        assert n>=0
        p = self.tp.one()
        for i in range(n):
            p = self*p
        return p



class Integers(Ring):

    # XXX these should be singletons (so we can test equality with "is")

    def __init__(self):
        Ring.__init__(self)

    def add(self, a, b):
        assert a.tp is self # do we need this here?
        assert b.tp is self # do we need this here?
        return RingElement(a.i + b.i, self)
    
    def sub(self, a, b):
        return RingElement(a.i - b.i, self)
    
    def mul(self, a, b):
        return RingElement(a.i * b.i, self)
    
    def neg(self, elem):
        return RingElement(-elem.i, self)
    
    def one(self):
        return RingElement(1, self)

    def zero(self):
        return RingElement(0, self)

    def promote(self, value):
        if isinstance(value, RingElement):
            assert value.tp is self
            return value
        assert int(value)==value
        return RingElement(value, self)


class Field(Ring):

    # XXX these should be singletons (so we can test equality with "is")

    def __init__(self, p):
        Ring.__init__(self)
        self.p = p

    def hash(self):
        return hash((self.__class__, self.p))

    def __eq__(self, other):
        assert self.__class__ is other.__class__
        return self.p == other.p

    def __ne__(self, other):
        assert self.__class__ is other.__class__
        return self.p != other.p

    def add(self, a, b):
        assert a.tp is self
        assert b.tp is self
        p = self.p
        return FieldElement((a.i + b.i)%p, self)
    
    def sub(self, a, b):
        p = self.p
        return FieldElement((a.i - b.i)%p, self)
    
    def mul(self, a, b):
        p = self.p
        return FieldElement((a.i * b.i)%p, self)
    
    def neg(self, elem):
        p = self.p
        return FieldElement((p-elem.i)%p, self)
    
    def inverse(self, elem):
        p = self.p
        assert 0<i<p
        for j in range(1, p):
            if (j*i)%p == 1:
                break
        else:
            assert 0
        return FieldElement(j, self)
    
    def one(self):
        return FieldElement(1, self)

    def zero(self):
        return FieldElement(0, self)

    def promote(self, value):
        if isinstance(value, FieldElement):
            assert value.tp is self
            return value
        assert int(value)==value
        value = value % self.p
        return FieldElement(value, self)


class RingElement(Element): # XXX rename this as Integer ... ??? XXX
    def __init__(self, i, tp):
        Element.__init__(self, tp)
        assert int(i)==i
        self.i = i

    def __eq__(self, other):
        tp = self.tp
        other = tp.promote(other)
        return self.i == other.i

    def __ne__(self, other):
        tp = self.tp
        other = tp.promote(other)
        return self.i != other.i

    def __hash__(self):
        return hash((self.tp, self.i))

    def __str__(self):
        return str(self.i)

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__, self.i)


class FieldElement(RingElement):
    def __init__(self, i, tp):
        Element.__init__(self, tp)
        assert int(i)==i
        assert 0<=i<tp.p
        self.i = i



class PolynomialRing(Ring):

    def __init__(self, base):
        Ring.__init__(self)
        self.base = base

    def zero(self):
        return Polynomial({}, self)

    def one(self):
        one = self.base.one()
        return Polynomial({0:one}, self)

    def promote(self, value):
        if isinstance(value, Polynomial):
            return value
        value = self.base.promote(value)
        return Polynomial({0:value}, self)

    def add(self, a, b):
        cs = dict(a.cs)
        for deg, coeff in b.cs.items():
            cs[deg] = cs.get(deg, 0) + coeff
        return Polynomial(cs, self)

    def sub(self, a, b):
        cs = dict(a.cs)
        for deg, coeff in b.cs.items():
            cs[deg] = cs.get(deg, 0) - coeff
        return Polynomial(cs, self)

    def mul(self, a, b):
        cs = {}
        for d0, c0 in a.cs.items():
          for d1, c1 in b.cs.items():
            deg = d0+d1
            cs[deg] = cs.get(deg, 0) + c0*c1
        return Polynomial(cs, self)



class Polynomial(Element):
    def __init__(self, cs, tp):
        Element.__init__(self, tp)
        self.cs = {}
        for deg, coeff in cs.items():
            if coeff==0: # strip these out
                continue
            coeff = tp.base.promote(coeff)
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
                term = "%s*x" % (coeff,)
            elif coeff==1:
                term = "x**%d" % deg
            else:
                term = "%s*x**%d" % (coeff, deg)
            terms.append(term)
        s = "+".join(terms)
        return s

    def __repr__(self):
        return "Polynomial(%s, %s)"%(self.cs, self.tp)

    def __eq__(self, other):
        assert self.tp == other.tp
        return self.cs == other.cs

    def __ne__(self, other):
        assert self.tp == other.tp
        return self.cs != other.cs

    def __hash__(self):
        cs = self.cs
        keys = list(cs.keys())
        keys.sort()
        value = tuple((deg, cs[deg]) for deg in keys)
        return hash(value)




def test():

    f = Field(5)
    one = f.one()
    zero = f.zero()

    assert one==1
    print(one)

    a = zero
    for i in range(1000):
        a = one + a
        if a==zero:
            break
    assert i==4, i

    # -------------------------

    Z = Integers()
    ring = PolynomialRing(Z)

    one = ring.one()
    x = Polynomial({1:1}, ring)

    assert (x+one)**5 == x**5+5*x**4+10*x**3+10*x**2+5*x+1
    print((x+one)**3)


    # -------------------------

    field = Field(5)
    ring = PolynomialRing(field)

    one = ring.one()
    x = Polynomial({1:1}, ring)

    assert (x+one)**5 == x**5+1




if __name__ == "__main__":

    test()





