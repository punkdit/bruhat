#!/usr/bin/env python

"""
central group extentions via 2-cocycles: unfinished
"""


from random import randint
from time import time
start_time = time()

import numpy

from bruhat.argv import argv
from bruhat.action import Perm, Group


def is_cocycle(f, G):

    n = len(G)
    for i, gi in enumerate(G):
     for j, gj in enumerate(G):
      for k, gk in enumerate(G):
        # finish me... XXX
        lhs = f[i, h*k] * f[h, k]
        rhs = f[g*h, k] * f[g, h]
        if lhs != rhs:
            return False
    return True



def main():
    p = argv.get("p", 2)

    G = Group.cyclic(5)
    n = len(G)

    f = numpy.zeros((n, n), dtype=int)
    for i in range(n):
     for j in range(n):
        f[i, j] = randint(0, p-1)

    print(f)
    assert is_cocycle(f, G)


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




