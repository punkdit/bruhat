#!/usr/bin/env python3

import math, os


class Number(object):
    def __init__(self, a):
        self.a = a
        self.shape = ()

    def __str__(self):
        return str(self.a)

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__, self.a)

    def __hash__(self):
        return hash(str(self))

    def promote(self, a):
        if isinstance(a, int):
            a = Number(a)
        return a

    def __add__(self, other):
        assert self.__class__ is other.__class__
        return self.a + other.a

    def __sub__(self, other):
        assert self.__class__ is other.__class__
        return self.a - other.a

    def __mul__(self, other):
        assert self.__class__ is other.__class__
        return self.a * other.a

    def __eq__(self, other):
        assert self.__class__ is other.__class__
        return self.a == other.a

    def __ne__(self, other):
        assert self.__class__ is other.__class__
        return self.a != other.a

    def __neg__(self):
        return Number(-self.a)

    def conj(self):
        return self

    def is_real(self):
        return True

    def is_zero(self):
        return self.a == 0

    def get_zero(self):
        return Number(0)


class Double(Number):
    def __init__(self, a, b):
        if isinstance(a, int):
            a = Number(a)
        if isinstance(b, int):
            b = Number(b)
        assert a.shape == b.shape
        self.shape = (a.shape, b.shape)
        assert isinstance(a, Number)
        self.a = a
        self.b = b

    def get_zero(self):
        return Double(self.a.get_zero(), self.b.get_zero())

    def promote(self, a):
        if isinstance(a, int):
            a = Number(a)
        assert isinstance(a, Number)
        if a.shape == self.shape:
            return a
        assert str(a.shape) in str(self.shape)
        a = self.a.promote(a)
        return Double(a, self.b.get_zero())

    def __repr__(self):
        #a, b = self.pair
        return "%s(%s, %s)"%(self.__class__.__name__, self.a, self.b)

    def __str__(self):
        return "(%s, %s)"%(self.a, self.b)

    def __add__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        a = self.a + other.a
        b = self.b + other.b
        return self.__class__(a, b)

    def __sub__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        a = self.a - other.a
        b = self.b - other.b
        return self.__class__(a, b)

    def __mul__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        a, b = self.a, self.b
        c, d = other.a, other.b
        x = self.__class__(a*c - d.conj()*b, d*a + b*c.conj())
        return x

    def __eq__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        return self.a == other.a and self.b == other.b

    def __ne__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        return self.a != other.a or self.b != other.b

    def __neg__(self):
        return self.__class__(-self.a, -self.b)

    def conj(self):
        return self.__class__(self.a.conj(), -self.b)

    def norm2(self):
        return self.conj() * self

    def is_real(self):
        return self.a.is_real() and self.b.is_zero()



def is_commutative(items):
    for a in items:
      for b in items:
        if a*b != b*a:
            return False
    return True


def is_associative(items):
    for a in items:
      for b in items:
        for c in items:
            if a*(b*c) != (a*b)*c:
                return False
    return True


def is_alternative(items):
    for a in items:
      for b in items:
        if a*(b*b) != (a*b)*b:
            return False
    return True


def test():

    x = Number(2)
    y = Double(2, 0)
    assert x==y

    one = Double(1, 0)
    i = Double(0, 1)
    assert i*i == -1

    cplex = [one, i]
    assert is_commutative(cplex)
    assert is_associative(cplex)

    # ----------- double: quaternions

    zero = Double(Double(0, 0), Double(0, 0))
    one = Double(Double(1, 0), Double(0, 0))
    i = Double(Double(0, 1), Double(0, 0))
    j = Double(Double(0, 0), Double(1, 0))
    k = Double(Double(0, 0), Double(0, 1))

    for x in [i, j, k]:
        assert x*x == -1
        for y in [i, j, k]:
            if x==y:
                continue
            assert x*y == -y*x
    assert i*j == -j*i
    assert i*j*k == -1

    quaternions = [one, i, j, k]
    assert not is_commutative(quaternions)
    assert is_associative(quaternions)

    # ----------- double: octonions

    octonions = [
        Double(one, zero), Double(zero, one),
        Double(i, zero), Double(j, zero), Double(k, zero),
        Double(zero, i), Double(zero, j), Double(zero, k)]

    imag = octonions[1:]
    for i in imag:
        assert i*i == -1

    
    assert not is_commutative(octonions)
    assert not is_associative(octonions)
    assert is_alternative(octonions)

    # ----------- double: sedenions

    one = Double(one, zero)
    zero = Double(zero, zero)

    sedenions = [Double(one, zero), Double(zero, one)]
    for i in octonions[1:]:
        sedenions.append(Double(i, zero))
        sedenions.append(Double(zero, i))

    assert not is_commutative(sedenions)
    assert not is_associative(sedenions)
    assert is_alternative(sedenions) # um...

    # try some more sedenions here:
    items = list(sedenions)
    for a in sedenions:
      for b in sedenions:
        items.append(a+b)
    assert not is_alternative(items)

    


if __name__ == "__main__":
    test()



