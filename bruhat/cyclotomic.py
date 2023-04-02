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

    class Ring: pass
    ring = Ring()
    ring.zero = Poly({}, Q)
    ring.one = Poly({():1}, Q)
    ring.x = Poly("x", Q)

    for n in [2*3, 3*5, 2*3*5, 2*3*7, 5*7]:
        print("n =", n)
        p = cyclotomic(ring, n)
        print(p)
        f = Rational(Q, ring.one, p)
    
        for i in range(100):
            val = f[i]
            if val == +1:
                print("+%s "%i, end="")
            elif val == -1:
                print("-%s "%i, end="")
        print()
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


