#!/usr/bin/env python

import sys, os

import numpy

from bruhat.argv import argv 
from bruhat.util import cross, factorial
#from bruhat import action
#from bruhat.action import Perm, Group, mulclose
from bruhat.gelim import solve, array, identity, zeros, dot, shortstr, eq, dotx, kernel
from bruhat.gelim import Subspace


def mulclose(gen, maxsize=None, verbose=False):

    found = set(str(g) for g in gen)
    els = list(gen)
    changed = True 
    while changed:
        if verbose:
            print "mulclose:", len(els)
        changed = False
        _els = list(els)
        for A in _els:
            for B in _els:
                C = dot(A, B)
                if str(C) not in found:
                    els.append(C)
                    found.add(str(C))
                    if maxsize and len(els)>=maxsize:
                        return list(els)
                    changed = True 
    return els




def build(n):
    "Return list of generators for orthogonal reflection group B_n"
    assert n>=2

    gen = []

    # basis-swaps, "controlled bitflips"
    for i in range(n-1):
        X = zeros(n, n)
        for j in range(n):
            if j==i:
                X[j, j+1] = 1
            elif j==i+1:
                X[j, j-1] = 1
            else:
                X[j, j] = 1
        gen.append(X)

    # sign-swap, "controlled phase-flip"
    Z = identity(n)
    Z[n-1, n-1] = -1
    gen.append(Z)

    return gen




def main():

    for g in build(2):
        print shortstr(g)
        print

    B4 = build(4)
    for g in B4:
        print shortstr(g)
        print

    #assert len(mulclose(B4)) == 384 # slow...





if __name__ == "__main__":

    main()


