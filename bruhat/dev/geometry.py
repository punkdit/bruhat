#!/usr/bin/env python3

"""
higher genus weight enumerators.
Check _they satisfy various identities.
"""

from random import random, randint, shuffle
from functools import reduce
from operator import mul

import numpy

from bruhat.solve import parse, shortstr, shortstrx, span, int_scalar, find_kernel
from bruhat.solve import row_reduce, dot2
from bruhat.poly import Poly
from bruhat import element
from bruhat.gset import Group
from bruhat import algebraic

from bruhat.argv import argv
from bruhat.util import choose, cross, all_perms


def zeros(m, n):
    A = numpy.zeros((m, n), dtype=int_scalar)
    return A

def identity(m):
    A = numpy.identity(m, dtype=int_scalar)
    return A


#class Subspace(object):

def normal_form(G):
    #print("normal_form")
    #print(G)
    G = row_reduce(G) # makes copy
    #print(G)
    m, n = G.shape
    if m*n==0:
        return G
    for i in range(m):
        for j in range(n):
            if G[i, j]:
                break
        else:
            assert 0
        k = i-1
        while k>=0:
            if G[k, j]:
                G[k] += G[i]
            k -= 1
        G %= 2
    return G


def intersect(G1, G2):
    G = numpy.concatenate((G1, G2))
    #print("intersect")
    #print(G1, G2)
    #print(G)
    G = G.transpose()
    #print("find_kernel", G.shape)
    K = find_kernel(G)
    if not K:
        K = numpy.array(K)
        K.shape = 0, G.shape[1]
    else:
        K = numpy.array(K)
    #print("K:")
    #print(K, K.shape)
    G = dot2(K[:, :len(G1)], G1)
    #print("G:")
    #print(G, G.shape)
    #print()
    G = normal_form(G)
    return G




def get_cell(row, col, p=2):
    """
        return all matrices in bruhat cell at (row, col)
        These have shape (col, col+row).
    """

    if col == 0:
        yield zeros(0, row)
        return

    if row == 0:
        yield identity(col)
        return

    # recursive steps:
    m, n = col, col+row
    for left in get_cell(row, col-1, p):
        A = zeros(m, n)
        A[:m-1, :n-1] = left
        A[m-1, n-1] = 1
        yield A

    els = list(range(p))
    vecs = list(cross((els,)*m))
    for right in get_cell(row-1, col, p):
        for v in vecs:
            A = zeros(m, n)
            A[:, :n-1] = right
            A[:, n-1] = v
            yield A

#for row in range(3):
#  for col in range(4):
#    print(len(list(get_cell(row, col))), end=" ")
#  print()

def all_codes(m, n, q=2):
    """
        All full-rank generator matrices of shape (m, n)
    """
    assert m<=n
    col = m
    row = n-m
    return get_cell(row, col, q)


def get_codespace(G):
    space = list(span(G))
    return space


def all_auto_codes(m, n):
    for G in all_codes(2, 3):
        print(G)
        space = get_codespace(G)
        #shuffle(space)
        space.sort(key = lambda v : v.tostring())
        space = numpy.array(space)
        print(space)


def main():

    q = 2
    m = argv.get("m", 2)
    n = argv.get("n", 4)
    cells = list(all_codes(m, n))

#    for G in cells:
#        print(G)
    print(len(cells))

    GL = algebraic.GL(n, q)

    names = ".PLSH"
    for G1 in cells:
      for G2 in cells:
        H = intersect(G1, G2)
        i = H.shape[0]
        print(names[i], end="")
      print()


if __name__ == "__main__":

    main()



