#!/usr/bin/env python3

"""
higher genus weight enumerators.
Check they satisfy various identities.
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


ring = element.Z

def named_poly(cs, rank, names):
    items = {}
    for k, v in cs.items():
        tpl = tuple((names[i], exp) for (i, exp) in enumerate(k) if exp)
        items[tpl] = v
    return Poly(items, ring)


xpoly1 = lambda cs : named_poly(cs, 2, 
    "x_0 x_1".split())
ypoly1 = lambda cs : named_poly(cs, 2, 
    "y_0 y_1".split())
xpoly2 = lambda cs : named_poly(cs, 4, 
    "x_{00} x_{01} x_{10} x_{11}".split())
xpoly3 = lambda cs : named_poly(cs, 8,
    "x_{000} x_{001} x_{010} x_{011} x_{100} x_{101} x_{110} x_{111}".split())



def genus_enum1(G, verbose=False):
    m, n = G.shape

    cs = {}
    for v0 in span(G):
        exp = [0, 0]
        for i in range(n):
            exp[v0[i]] += 1
        exp = tuple(exp)
        cs[exp] = cs.get(exp, 0) + 1
    p = xpoly1(cs)

    return p


def genus_enum2(G, verbose=False):
    m, n = G.shape
    cs = {}
    for v0 in span(G):
        #print(".",end='',flush=True)
        for v1 in span(G):
            exp = [0, 0, 0, 0]
            #vv = numpy.array([2*v0, v1])
            for i in range(n):
                exp[2*v0[i] + v1[i]] += 1
            exp = tuple(exp)
            cs[exp] = cs.get(exp, 0) + 1
        #break
    #print()
    p = xpoly2(cs)

    return p


def genus_enum3(G, verbose=False):
    #print(shortstr(G))
    m, n = G.shape

    rank = 8
    cs = {}
    for v0 in span(G):
        #print(".",end='',flush=True)
        for v1 in span(G):
          for v2 in span(G):
            exp = [0, 0, 0, 0, 0, 0, 0, 0]
            #vv = numpy.array([2*v0, v1])
            for i in range(n):
                exp[4*v0[i] + 2*v1[i] + v2[i]] += 1
            exp = tuple(exp)
            cs[exp] = cs.get(exp, 0) + 1
        #break
    #print()
    p = xpoly3(cs)

    return p


def test_enum(G, verbose=False):

    x_0 = xpoly1({(1,0):1})
    x_1 = xpoly1({(0,1):1})
    y_0 = ypoly1({(1,0):1})
    y_1 = ypoly1({(0,1):1})
    x_00 = xpoly2({(1,0,0,0):1})
    x_01 = xpoly2({(0,1,0,0):1})
    x_10 = xpoly2({(0,0,1,0):1})
    x_11 = xpoly2({(0,0,0,1):1})

    w1 = genus_enum1(G)
    w2 = genus_enum2(G)
    w3 = genus_enum3(G)

    if verbose:
        print(w1)
        print(w2)
        print(w3)

    zero1 = xpoly1({})
    zero2 = xpoly2({})

    w2_10 = w2.substitute({ "x_00" : x_0, "x_10" : x_1, "x_01":zero1, "x_11":zero1 })
    assert w2_10 == w1

    w2_01 = w2.substitute({ "x_00" : x_0, "x_10" : zero1, "x_01":x_1, "x_11":zero1 })
    assert w2_01 == w1

    w2_11 = w2.substitute({ "x_00" : x_0, "x_10" : zero1, "x_01":zero1, "x_11":x_1 })
    assert w2_11 == w1

    lhs = w2.substitute({ 
        "x_00" : x_0*y_0, 
        "x_10" : x_1*y_0, 
        "x_01" : x_0*y_1, 
        "x_11" : x_1*y_1,
    })
    yw1 = w1.substitute({ "x_0" : y_0, "x_1" : y_1, })
    rhs = w1 * yw1
    assert lhs == rhs

    lhs = w3.substitute({ 
        "x_000" : x_00*y_0, 
        "x_010" : x_01*y_0, 
        "x_001" : x_00*y_1, 
        "x_011" : x_01*y_1,
        "x_100" : x_10*y_0, 
        "x_110" : x_11*y_0, 
        "x_101" : x_10*y_1, 
        "x_111" : x_11*y_1,
    })
    rhs = w2 * yw1
    assert lhs == rhs

    lhs = w3.substitute({ 
        "x_000" : y_0*x_00, 
        "x_010" : y_0*x_10, 
        "x_001" : y_0*x_01, 
        "x_011" : y_0*x_11,
        "x_100" : y_1*x_00, 
        "x_110" : y_1*x_10, 
        "x_101" : y_1*x_01, 
        "x_111" : y_1*x_11,
    })
    rhs = w2 * yw1
    assert lhs == rhs

    if argv.verbose:
        print("G =")
        print(shortstr(G))
        print("w1 =", w1)
        print("w2 =", w2)
        print("w3 =", w3)
        print()


def test():

    G = parse("11 ..")
    test_enum(G)

    G = parse("11. .11")
    test_enum(G, verbose=True)

    #G = parse("11.. .11. 1..1")
    #test_enum(G)

    for m in range(1, 5):
      for n in range(m, 6):
        for G in all_codes(m, n):
            test_enum(G)

    print("OK")


def hecke():

    n = 5 # cols

    codes = []
    lookup = {}
    for m in range(0, n): # rows
        for G in all_codes(m, n):
            codes.append(G)
            lookup[G.tostring()] = G
    N = len(codes)
    pairs = []
    for i in range(N):
      for j in range(i+1, N):
        GH = intersect(codes[i], codes[j])
        assert GH.tostring() in lookup, GH


if __name__ == "__main__":

    test()


