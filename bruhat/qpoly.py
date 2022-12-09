#!/usr/bin/env python
"""
q polynomials for grassmanians, etc.
"""

import string, os
from time import sleep, time
from functools import reduce
from operator import matmul

import numpy
from numpy import alltrue, zeros, dot

from bruhat.argv import argv
from bruhat.smap import SMap
from bruhat.element import Z
from bruhat.poly import Poly

# See:
# https://en.wikipedia.org/wiki/Q-analog
# https://math.ucr.edu/home/baez/week187.html


ring = Z
zero = Poly({}, ring)
one = Poly({():1}, ring)
q = Poly("q", ring)

def shortstr(p):
    assert isinstance(p, Poly)
    keys = p.keys
    items = {}
    m = 0
    for key in keys:
        if len(key)==0:
            items[0] = p.cs[key]
        else:
            c, n = key[0]
            items[n] = p.cs[key]
            m = max(n, m)
    items = [items.get(i, 0) for i in range(m+1)]
    return ''.join(str(c) for c in items)
    

def qbracket(n):
    p = zero
    for i in range(n):
        p = p + q**i
    return p

def qfactorial(n):
    p = one
    for i in range(n):
        p = p * qbracket(i+1)
    return p

def A(n):
    return qfactorial(n+1)

def B(n):
    p = one
    for i in range(n):
        p = p * qbracket(2*(i+1))
    return p
C = B

def D(n):
    p = one
    for i in range(n-1):
        p = p * qbracket(2*(i+1))
    p = p * qbracket(n)
    return p


def test():
    assert q**0 == one

    p = (1+q)**5

    assert qbracket(5) == 1 + q + q**2 + q**3 + q**4

    assert qfactorial(3) == 1 + 2*q + 2*q**2 + q**3

    # These count points of the maximal flag variety over the
    # field with q elements:
    assert shortstr(A(0)) == "1"
    assert shortstr(A(1)) == "11"
    assert shortstr(A(2)) == "1221"
    assert shortstr(A(3)) == "1356531"

    assert shortstr(B(1)) == "11"
    assert shortstr(B(2)) == "12221"
    assert shortstr(B(3)) == "1357887531"

    assert shortstr(D(1)) == "1"
    assert shortstr(D(2)) == "121"
    assert shortstr(D(3)) == "1356531"

    p = D(4)
    p = p / (A(1)**3)
    assert shortstr(p) == "1133443311"

    get = lambda p : p.substitute((('q',2),))
    print(get(p))

    assert get(B(1) / A(0)) == 3

    assert get(B(2) / A(1)) == 15

    assert get(B(3) / B(2)) == 63
    assert get(B(3) / (A(1)**2)) == 315
    assert get(B(3) / A(2)) == 135

    assert get(B(4) / B(3)) == 255
    assert get(B(4) / (A(1) * B(2))) == 5355
    assert get(B(4) / (A(2) * A(1))) == 11475
    assert get(B(4) / A(3)) == 2295


if __name__ == "__main__":

    start_time = time()


    profile = argv.profile
    name = argv.next()
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        test()

    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")


