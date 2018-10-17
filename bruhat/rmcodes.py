#!/usr/bin/env python3

"""
Build Reed-Muller codes.
"""


import numpy

from solve import array2, zeros2, dot2, shortstr
from argv import argv
from util import choose


def main():

    r = argv.get("r", 1) # degree
    m = argv.get("m", 3)
    assert 0<r<=m

    n = 2**m # length

    one = array2([1]*n)
    basis = [one]

    vs = [[] for i in range(m)]
    for i in range(2**m):
        for j in range(m):
            vs[j].append(i%2)
            i >>= 1
        assert i==0

    vs = [array2(v) for v in vs]

    for k in range(r):
        for items in choose(vs, k+1):
            v = one
            #print(items)
            for u in items:
                v = v*u
            basis.append(v)
        
    G = numpy.array(basis)
    print(shortstr(G))



if __name__ == "__main__":

    main()



