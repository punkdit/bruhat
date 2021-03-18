#!/usr/bin/env python3

import numpy

from bruhat.element import FiniteField, PolynomialRing, GaloisField
from bruhat.action import Perm, Group
from bruhat.util import factorize, cross, all_primes
from bruhat.argv import argv


def GF(q):
    ps = factorize(q)
    assert len(set(ps)) == 1
    p = ps[0]
    r = len(ps)
    assert q==p**r
    #print(p, r)
    
    field = FiniteField(p)
    if r==1:
        return field

    ring = PolynomialRing(field)

    zero = ring.zero
    one = ring.one
    x = ring.x

    itemss = [tuple(range(p)) for i in range(r)]
    for idxs in cross(itemss):
        poly = x**r
        for i, idx in enumerate(idxs):
            poly += idx*(x**i)
        #print(poly)
        for i in range(p):
            if poly(i) == zero:
                #print("i=", i)
                break
        else:
            break
    #print("poly:", poly)
    #print([str(poly(i)) for i in range(p)])

    F = GaloisField(ring, poly)
    def frobenius(a):
        return a**p
    def hermitian(a, b):
        return (a**p)*b + a*(b**p)
    F.frobenius = frobenius
    F.hermitian = hermitian
    return F


class Code(object):
    def __init__(self, field, gen):
        self.field = field


def test():

    # -------------------------

    # check we can build galois field GF(8)

    field = FiniteField(2)
    ring = PolynomialRing(field)

    one = ring.one
    x = ring.x

    f = x**3 - x - 1
    assert f == x**3+x+1
    assert f(0) != 0
    assert f(1) != 0

    b = x**5
    div, rem = f.reduce(b)
    assert f*div + rem == b

    group = []
    for i in range(2):
     for j in range(2):
      for k in range(2):
        a = i + j*x + k*x**2
        if a != 0:
            group.append(a)

    div = {}
    for a in group:
      for b in group:
        c = f.reduce(a*b)[1]
        div[c, a] = b
        div[c, b] = a

    # all non-zero pairs elements of GF(8) should be divisable:
    assert len(div) == 49

    # --------------------------------------------------

    # GF(4)
    x = ring.x
    F = GaloisField(ring, (x**2 + x + 1))

    omega = F.x
    one = F.one
    a, b = omega, omega+one
    assert a*a == b
    assert a*b == b*a
    assert a*b == one
    assert b*b == a

    assert one/omega == one+omega

    frobenius = lambda a : a**2
    assert frobenius(one) == one
    assert frobenius(one+one) == one+one
    assert frobenius(omega) == omega+1
    assert frobenius(omega+1) == omega

#    # --------------------------------------------------
#
#    for p in all_primes(20):
#        assert p>1
#        field = FiniteField(p)
#        ring = PolynomialRing(field)
#        zero = ring.zero
#        one = ring.one
#        x = ring.x
#        poly = 1+x+x**2
#        for i in range(1, p):
#            assert poly(i) != zero, "p=%d, poly=%s, i=%d"%(p, poly, i)

    # --------------------------------------------------

    p = argv.get("p", 2)
    r = argv.get("r", 2)

    print("q =", p**r)
    F = GF(p**r)
    print(F.mod)
    zero = 0

    els = F.elements
    assert len(els) == len(set(els)) == p**r

    for a in els:
      for b in els:
        if b!=0:
            c = a/b # XXX fails for p=3, r=4
            assert c*b==a

    # build the hexacode
    w = F.x
    w2 = w**2
    words = []
    for a in els:
     for b in els:
      for c in els:
        f = lambda x : a*x**2 + b*x + c
        v = [a, b, c, f(1), f(w), f(w2)]
        words.append(v)
    code = numpy.array(words)
    for word in code:
      print(' '.join("%4s"%(c) for c in word))

    print(code.shape)
    assert len(code)==4**3

    def inner(w0, w1):
        assert len(w0)==len(w1)
        n = len(w0)
        #r = sum([F.hermitian(a, b) for (a, b) in zip(w0, w1)])
        r = sum([a*F.frobenius(b) for (a, b) in zip(w0, w1)])
        return r

    for w0 in code:
      for w1 in code:
        print(inner(w0, w1), end=" ")
      print()


def test_galois():
    p = argv.get("p", 2)
    r = argv.get("r", 2)

    print("q =", p**r)
    F = GF(p**r)
    print(F.mod)
    zero = 0

    els = F.elements
    assert len(els) == len(set(els)) == p**r

    for a in els:
      for b in els:
        if b!=0:
            c = a/b # XXX fails for p=3, r=4
            assert c*b==a

    def orbit(a):
        items = set([a])
        while 1:
            a = F.frobenius(a)
            if a in items:
                break
            items.add(a)
        return items

    lookup = {e:i for (i,e) in enumerate(els)}
    g = {lookup[e] : lookup[F.frobenius(e)] for e in els}
    #print(["%s:%s"%(k,v) for (k,v) in g.items()])
    values = set(g.values())
    assert len(values)==len(g)
    items = list(range(len(els)))
    g = Perm(g, items)
    I = Perm.identity(items)
    G = [I]
    h = g
    print(g)
    while h != I:
        G.append(h)
        h = g*h
        #print(len(G))
        assert len(G) <= p**r
    print("|G| =", len(G))
    G = Group.generate([g])
    print("|G| =", len(G))

    for i in els:
        n = len(orbit(i))
        if n==1:
            print("*", end="")
        print(n, end=" ")
    print()

    print("OK")


if __name__ == "__main__":

    test()

