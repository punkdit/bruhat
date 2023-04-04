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

from bruhat.util import is_prime
from bruhat.elim import zeros, shortstr, pseudo_inverse, dot
from bruhat.element import Q, cyclotomic
from bruhat.rational import Poly, Rational

def main():
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


