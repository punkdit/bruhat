#!/usr/bin/env python

"""
PSL(2,Q)

"""

from random import randint

import numpy

from bruhat.argv import argv
if argv.fast:
    from bruhat import _element as element
else:
    from bruhat import element

#from bruhat.chain import Space, Lin

"""
Aut(H) = {
z |->  az + b
       ------
       cz + d
}
= PSL(2, R)

"""

class Mat(object):
    def __init__(self, ring, a, b, c, d):
        self.ring = ring
        pos = (a, b, c, d)
        neg = (-a, -b, -c, -d)

        # we need to make this canonical for hash( )
        zero = ring.zero
        if a == zero and b > zero:
            pass # OK
        elif a == zero and b < zero:
            pos, neg = neg, pos # swap
        elif a == zero:
            assert False, (a, b)
        elif a < zero:
            pos, neg = neg, pos # swap
        self.pos = pos
        self.neg = neg

    @classmethod
    def construct(cls, ring, a, b, c, d):
        a = ring.promote(a)
        b = ring.promote(b)
        c = ring.promote(c)
        d = ring.promote(d)
        x = a*d - b*c
        assert x == ring.one
        return cls(ring, a, b, c, d)

    def check(self):
        a, b, c, d = self.pos
        x = a*d - b*c
        assert x==self.ring.one

    @classmethod
    def rand(cls, ring, n=5):
        assert n>0
        while 1:
            a = ring.promote(randint(-n, n))
            b = ring.promote(randint(-n, n))
            if a != ring.zero or b != ring.zero: # easy
                break
        if a != ring.zero:
            c = ring.promote(randint(-n, n))
            d = ( ring.one + b*c ) / a
        else:
            d = ring.promote(randint(-n, n))
            c = (a*d - 1)/b
        m = Mat(ring, a, b, c, d)
        m.check()
        return m

    def __mul__(self, other):
        assert isinstance(other, Mat)
        assert self.ring is other.ring
        a, b, c, d = self.pos
        e, f, g, h = other.pos
        m = Mat(self.ring, a*e + b*g, a*f+b*h, c*e+d*g, c*f+d*h)
        return m

    def __eq__(self, other):
        assert isinstance(other, Mat)
        assert self.ring is other.ring
        return self.pos == other.pos or self.pos == other.neg

    def __hash__(self):
        return hash(self.pos)

    def __str__(self):
        return "Mat(%s, %s, %s, %s)"%(self.pos)
    __repr__ = __str__


def test():

    ring = element.Q
    construct = lambda *args : Mat.construct(ring, *args)
    I = construct(1, 0, 0, 1)
    nI = construct(-1, 0, 0, -1)
    A = construct(1, 1, 0, 1)
    Ai = construct(1, -1, 0, 1)
    #print(construct(2, 0, 0, 2))
    #assert I == construct(2, 0, 0, 2)

    assert I == nI

    assert A == A
    assert A != I
    assert A*A == construct(1, 2, 0, 1)
    assert A*Ai == I

    for i in range(1000):
        A = Mat.rand(ring, 3)
        B = Mat.rand(ring, 3)
        C = A*B
        lhs = (A==B) 
        rhs = (hash(A) == hash(B))
        assert not lhs or rhs


if __name__ == "__main__":
    test()







