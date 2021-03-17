#!/usr/bin/env python3

from bruhat.element import FiniteField, PolynomialRing, GaloisField
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
            if poly(i) == 0:
                #print("i=", i)
                break
        else:
            break
        continue
    #print("poly:", poly)

    F = GaloisField(ring, poly)
    F.frobenius = lambda a : a**p
    return F


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

    p = argv.get("p", 3)
    r = argv.get("r", 2)

    print("q =", p**r)
    F = GF(p**r)
    zero = 0

    a = F.x
    def orbit(a):
        items = set([a])
        while 1:
            a = F.frobenius(a)
            if a in items:
                break
            items.add(a)
        return items

    els = F.elements
    assert len(els) == len(set(els)) == p**r

    for i in els:
        n = len(orbit(i))
        if n==1:
            print("*", end="")
        print(n, end=" ")
    print()

    print("OK")


if __name__ == "__main__":

    test()

