#!/usr/bin/env python3

"""
Implement the Cayley-Dickson _construction to
get complex numbers, quaternions, octonions, sedenions, etc.

copied from sedenions.py

"""

import math, os
from functools import reduce
from operator import add, mul
from time import time
start_time = time()

from bruhat import element
from bruhat.util import choose, cross
from bruhat.argv import argv


class Number(object):
    def __init__(self, a):
        self.a = a
        self.shape = ()

    def __str__(self):
        return str(self.a)

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__, self.a)

    def __hash__(self):
        assert 0
        print("__hash__", self.__class__.__name__)
        return hash(str(self))

    def promote(self, a):
        if not isinstance(a, Number):
            a = Number(a) # wrap it up
        return a

    def __add__(self, other):
        #assert 0
        assert self.__class__ is other.__class__
        return Number(self.a + other.a)

    def __sub__(self, other):
        #assert 0
        assert self.__class__ is other.__class__
        return Number(self.a - other.a)

    def __mul__(self, other):
        #assert 0
        assert self.__class__ is other.__class__
        return Number(self.a * other.a)

    def __rmul__(self, other):
        return Number(self.a * other)

    def __truediv__(self, other):
        return self.__class__(self.a / other)

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
        if not isinstance(a, Number):
            a = Number(a)
        if not isinstance(b, Number):
            b = Number(b)
        assert a.shape == b.shape
        self.shape = (a.shape, b.shape)
        assert isinstance(a, Number)
        self.a = a
        self.b = b

    def get_zero(self):
        return Double(self.a.get_zero(), self.b.get_zero())

    def promote(self, a):
        if not isinstance(a, Number):
            a = Number(a)
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

    def __rmul__(self, other):
        return self.__class__(self.a * other, self.b * other)

    def __truediv__(self, other):
        return self.__class__(self.a / other, self.b / other)

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

    def __hash__(self):
        return hash(str(self))

    def __neg__(self):
        return self.__class__(-self.a, -self.b)

    def conj(self):
        return self.__class__(self.a.conj(), -self.b)

    def norm2(self):
        return self.conj() * self

    def is_real(self):
        return self.a.is_real() and self.b.is_zero()

    def is_zero(self):
        return self.a.is_zero() and self.b.is_zero()



def is_commutative(items):
    for a in items:
      for b in items:
        if a*b != b*a:
            return False
    return True


def is_anticommutative(items):
    for a in items:
      for b in items:
        if a!=b and a*b != -b*a:
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


def main():

    ring = element.Q
    one, zero = ring.one, ring.zero

    x = Number(2*one)
    y = Double(2*one, zero)
    assert x==y

    # ----------- double: complex --------------------

    one = Double(ring.one, ring.zero)
    i = Double(ring.zero, ring.one)
    zero = Double(ring.zero, ring.zero)
    assert i*i == -1

    cplex = [one, i]
    assert is_commutative(cplex)
    assert is_associative(cplex)

    # ----------- double: quaternions --------------------

    zero, one, i, j, k = [
        Double(zero, zero),
        Double(one, zero),
        Double(i, zero),
        Double(zero, one),
        Double(zero, i),
    ]

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
    assert is_anticommutative(quaternions[1:])
    assert is_associative(quaternions)

    # ----------- double: octonions --------------------

    octonions = [
        Double(one, zero), 
        Double(zero, one),
        Double(i, zero), Double(j, zero), Double(k, zero),
        Double(zero, i), Double(zero, j), Double(zero, k)]

    one = octonions[0]
    imag = octonions[1:]
    for i in imag:
        assert i*i == -1
    i0, i1, i2, i3, i4, i5, i6 = imag
    inf = one

    lhs = inf + i0 + i2 + i6
    rhs = inf + i0 + i1 + i3
    print(lhs * rhs)
    print(i0 + i2 + i3 + i5)

    return

    assert not is_commutative(octonions)
    assert not is_associative(octonions)
    assert is_anticommutative(octonions[1:])
    assert is_alternative(octonions)

    assert one + one == 2*one
    units = [one, -one]
    for i in imag:
        units += [i, -i]
    units = set(units)

    return

    bag = set()
    signss = list(cross([(-1,1)]*4))
    for items in choose(imag, 4):
     for signs in signss:
        ii = [x*i for (x,i) in zip(signs, items)]
        k = reduce(add, ii)
        k = k/2
        bag.add(k)
    print(len(bag))
    bag = list(bag)
    N = len(bag)
    for idx in range(N):
      for jdx in range(idx, N):
        i, j = bag[idx], bag[jdx]
        if i*j != one:
            continue
        if i not in units:
            print(i)
            units.add(i)
        if j not in units:
            print(j)
            units.add(j)
#      for a in units:
#        for b in units:
#          if a is b:
#            continue
#          assert a != b
#          assert hash(a) != hash(b)
#          assert str(a) != str(b)
        
    print()
    print(len(units))

    for a in units:
      for b in units:
        assert a*b in units


if __name__ == "__main__":

    from argv import argv
    name = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%name)
    else:
        fn = eval(name)
        fn()

    print("OK: %.3f seconds\n"%(time() - start_time))






