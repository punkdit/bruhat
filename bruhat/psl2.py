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

    @property
    def a(self):
        return self.pos[0]

    @property
    def b(self):
        return self.pos[1]

    @property
    def c(self):
        return self.pos[2]

    @property
    def d(self):
        return self.pos[3]


    def check(self):
        a, b, c, d = self.pos
        x = a*d - b*c
        assert x==self.ring.one

# um... need complex numbers here...
#    def send(self, z):
#        "Apply mobius transform"
#        z = self.ring.promote(z)
#        a, b, c, d = self.pos
#        z = (a*z + b) / (c*z + d)
#        return z
#    __call__ = send

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

    def __pow__(self, i):
        assert i>=0
        assert int(i)==i
        if i==0:
            return self.construct(1, 0, 0, 1)
        if i==1:
            return self
        A = self
        while i>1:
            A = A*self
            i -= 1
        return A

    def __eq__(self, other):
        assert isinstance(other, Mat)
        assert self.ring is other.ring
        return self.pos == other.pos or self.pos == other.neg

    def __hash__(self):
        return hash(self.pos)

    def __str__(self):
        return "Mat(%s, %s, %s, %s)"%(self.pos)
    __repr__ = __str__

    def is_modular(self, n):
        a, b, c, d = self.pos
        result = (a.top%n==a.bot%n)
        result = result and (d.top%n==d.bot%n)
        result = result and (d.top%n==d.bot%n)

    is_modular = lambda m : m.a%n==1 and m.d%n==1 and m.b%n==0 and m.c%n==0

def mulclose_fast(gen, verbose=False, maxsize=None):
    els = set(gen)
    bdy = list(els)
    changed = True
    while bdy:
        if verbose:
            print(len(els), end=" ", flush=True)
        _bdy = []
        for A in gen:
            for B in bdy:
                C = A*B
                if C not in els:
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
        bdy = _bdy
    return els


def mulclose_subgroup(gen, test, verbose=False, maxsize=None):
    "test is a callback: is the element in the subgroup ?"
    els = set(g for g in gen if not test(g))
    bdy = list(els)
    changed = True
    while bdy:
        if verbose:
            print(len(els), end=" ", flush=True)
        _bdy = []
        for A in gen:
            for B in bdy:
                C = A*B
                if C not in els and not test(C):
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
        bdy = _bdy
    return els




def test():

    # --------------------------------------------------------------
    # PSL(2, Q)

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

    # --------------------------------------------------------------
    # PSL(2, Z)

    # https://en.wikipedia.org/wiki/Modular_group
    # G = <S, T | S*S==I, (S*T)**3 == I>

    ring = element.Z
    construct = lambda *args : Mat.construct(ring, *args)

    I = construct(1, 0, 0, 1)
    S = construct(0, -1, 1, 0)
    T = construct(1, 1, 0, 1)
    assert S*S == I
    assert (S*T)**3 == I

    G = mulclose_fast([S, T], maxsize=10000)
    print(len(G))
    assert len(G) >= 10000

    n = 3
    is_modular = lambda m : m.a%n==1 and m.d%n==1 and m.b%n==0 and m.c%n==0

#    # wah XXX
#    J = mulclose_subgroup([S, T], is_modular, verbose=True)
#    print(len(J))


if __name__ == "__main__":
    test()







