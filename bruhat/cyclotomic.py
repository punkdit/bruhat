#!/usr/bin/env python


from time import time
start_time = time()
from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)

import operator
import itertools
from random import randint, seed

import numpy

from bruhat.argv import argv
#argv.fast = True

from bruhat.util import is_prime, all_primes, cross
from bruhat.elim import zeros, shortstr, pseudo_inverse, dot
from bruhat.element import Z, Q, cyclotomic, FiniteField, PolynomialRing
from bruhat.rational import Poly, Rational


def split_field():
    deg, N = 3, 3
    deg, N = 5, 5
    deg, N = 8, 8
    deg, N = 10, 5
    deg, N = 20, 63
    moduli = set()
    for p in all_primes(2000):
        if p < 10:
            continue
        print("p=%d:\t"%p, end="")
        ring = FiniteField(p)
        ring = PolynomialRing(ring)
        #f = cyclotomic(ring, deg)
        #print(f)
        x = ring.x
        #f = x**4 + 3*x**2 + 1
        f = x**2 -x - 1
        f = x**3 -7*x - 7
        splits = False
        for i in range(p):
            val = f(i)
            if val==0:
                #print("*", end=" ")
                splits = True
                break
            #print(f(i), end=" ")
        line = ['.' for i in range(N)]
        #print("*" if splits else "")
        line[p%N] = "X" if splits else "|"
        print(' '.join(line))
        if splits:
            moduli.add(p%N)
    moduli = list(moduli)
    moduli.sort()
    print(moduli)


def get_primes(ring, N=3, units=[]):
    # argh, this is too slow 
    deg = ring.mod.deg
    zero, one, x = ring.zero, ring.one, ring.x
    basis = [one, x]
    assert deg>1
    for i in range(2, deg):
        basis.append(x**i)

    items = set()
    a = tuple(range(-N, N+1))
    for cs in cross([a]*deg):
        b = zero
        for (c, v) in zip(cs, basis):
            b += c*v
        if b not in units:
            items.add(b)
    items.remove(zero)
    print("items:", len(items))

    units = set()
    for a in items:
     for b in items:
        if a*b == one:
            units.add(a)
            units.add(b)
    print("units:", [str(u) for u in units], len(units))
    for u in units:
        items.remove(u)

    found = set(items)
    for a in items:
     for b in items:
        ab = a*b
        if ab in found:
            found.remove(ab)
    return found


def main_artin():

    base = PolynomialRing(Z)
    R = base // cyclotomic(base, 12)
    G = base // cyclotomic(base, 4) # Gaussian integers
    E = base // cyclotomic(base, 3) # Eisenstein integers
    #print(R, type(R))
    #print(E.degree)

    ps = set(get_primes(E))
    assert (E.one-E.x) in ps

    #def hom(p):
    #    q = eval(str(p), {"x":

    for p in get_primes(R, 2):
        print(p)

    one = R.one
    x = R.x
    zero = R.zero

    w = x**4
    i = x**3
    assert i**2 == -one
    assert w**2 + w + one == zero
    


def main_generating():
    "some generating function ideas..."
    n = argv.get("n", 20)
    ring = Q

    divs = [i for i in range(1, n+1) if n%i == 0]
    print(divs)
    N = len(divs)
    A = zeros(ring, N, N)
    for ii, i in enumerate(divs):
     for jj, j in enumerate(divs):
        if i%j == 0:
            A[ii, jj] = ring.one
            
    print(shortstr(A))
    B = pseudo_inverse(ring, A)
    print(shortstr(B))
    print(B.sum())

    #print(shortstr(dot(ring, A, B)))

    n = 20
    M = numpy.zeros((n, n), dtype=int)
    for i in range(n):
        M[i, (i+1)%n] = 1
        M[i, (i-1)%n] = 1
    v = numpy.array([0 for i in range(n)])
    v[0] = 1
    v[2] = 1
    v[10] = -1
    v[12] = -1
    v0 = v.copy()
    for i in range(2*n):
        print(len([i for i in v if i!=0]), end=' ')
        #print(list(v))
        v = numpy.dot(M, v)

    #return

    class Ring: pass
    ring = Ring()
    ring.zero = Poly({}, Q)
    ring.one = Poly({():1}, Q)
    ring.x = Poly("x", Q)

    print()

    #for n in [2*3, 7, 3*3, 3*5, 2*2*5, 2*3*5, 2*3*7, 5*7, 3*5*7]:
    for n in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 30, 5*7, 3*5*11]:
        print("n =", n)
        p = cyclotomic(ring, n)
        print(p)
        f = Rational(Q, ring.one, p)
    
        pos, neg = [], []
        for i in range(2*n):
            val = f[i]
            if val == +1:
                #print("+%s "%i, end="")
                pos.append(i)
            elif val == -1:
                #print("-%s "%i, end="")
                neg.append(i)
        print()
        print(",".join(str(i) for i in pos))
        print(",".join(str(i) for i in neg))
        print()



if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds"%(time() - start_time))
    print()


