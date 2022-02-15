#!/usr/bin/env python3


"""
"""

from random import random, randint
from functools import reduce
from operator import mul

import numpy

from bruhat.solve import array2, zeros2, dot2, shortstr, rank, find_kernel, span
from bruhat.solve import linear_independent, parse, pseudo_inverse, eq2, rand2, rank
from bruhat.action import mulclose
from bruhat.comm import Poly
from bruhat.argv import argv
from bruhat.util import choose, cross, all_perms


def symplectic(n):
    F = zeros2(2*n, 2*n)
    for i in range(n):
        F[i, n+i] = 1
        F[n+i, i] = 1
    return F


def main():

    n = 8


#   ZZZZZZZZ|XXXXXXXX
#   12345678|12345678
    H = parse("""
    111..1..|........
    1..1....|1.1.....
    ........|11.11...
    .1..1...|.1...1..
    ..1...1.|...1..1.
    ...11.11|........
    .....1.1|....1..1
    ........|..1..111
    """.replace("|",""))

    print()
    print("H=")
    print(shortstr(H))

    F = symplectic(n)
    C = dot2(H, dot2(F, H.transpose()))

    for i in range(n):
      for j in range(i+1, n):
        if C[i, j]:
            print("fail:", i+1, j+1)

    print(rank(H))
    H = linear_independent(H)
    print("H=")
    print(shortstr(H))

    HF = dot2(H, F)
    K = array2(find_kernel(HF))
    print("K=")
    print(shortstr(K))

    HK = numpy.concatenate((H, K))
    L = linear_independent(HK)
    print()
    print(shortstr(L))



if __name__ == "__main__":

    main()


