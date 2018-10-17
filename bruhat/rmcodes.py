#!/usr/bin/env python3

"""
Build Reed-Muller codes.
"""


import numpy

from solve import array2, zeros2, dot2
from argv import argv
from util import choose


def main():

    r = argv.get("r", 1) # degree
    m = argv.get("m", 3)

    n = 2**m # length

    basis = [[1]*n]

    vs = [[] for i in range(m)]
    for i in range(2**m):
        for j in range(m):
            vs[j].append(i%2)
            i >>= 1
        assert i==0

    for v in vs:
        print(v)

    for k in range(r):
        for items in choose(vs, k+1):
            print(items)



if __name__ == "__main__":

    main()



