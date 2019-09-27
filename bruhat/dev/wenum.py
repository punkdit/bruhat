#!/usr/bin/env python3

"""
higher genus weight enumerators.
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


poly1 = lambda cs : Poly(cs, 2, "x_0 x_1".split())
poly2 = lambda cs : Poly(cs, 4, "x_{00} x_{01} x_{10} x_{11}".split())
poly3 = lambda cs : Poly(cs, 8,
    "x_{000} x_{001} x_{010} x_{011} x_{100} x_{101} x_{110} x_{111}".split())

def genus_enum1(G, verbose=False):
    m, n = G.shape

    #x_1 = poly1({(0, 1) : 1})
    #x_0 = poly1({(1, 0) : 1})
    #xs = [x0, x1]
    cs = {}
    for v0 in span(G):
        exp = [0, 0]
        #vv = numpy.array([2*v0, v1])
        for i in range(n):
            exp[v0[i]] += 1
        exp = tuple(exp)
        cs[exp] = cs.get(exp, 0) + 1
    p = poly1(cs)

    return p


def genus_enum2(G, verbose=False):
    m, n = G.shape

#    x_11 = poly2({(0, 0, 0, 1) : 1})
#    x_10 = poly2({(0, 0, 1, 0) : 1})
#    x_01 = poly2({(0, 1, 0, 0) : 1})
#    x_00 = poly2({(1, 0, 0, 0) : 1})
#    xs = [[x_00, x_01], [x_10, x_11]]
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
    p = poly2(cs)

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
    p = poly3(cs)

    return p


def main():

    G = parse("11 ..")
    test(G)

    G = parse("11. .11")
    test(G)


def test(G):

    print("G =")
    print(shortstr(G))

    w1 = genus_enum1(G)
    print("w1 =", w1)
    w2 = genus_enum2(G)
    print("w2 =", w2)
    w3 = genus_enum3(G)
    print("w3 =", w3)

    x_0 = poly1({(1,0):1})
    x_1 = poly1({(0,1):1})
    x_00 = poly2({(1,0,0,0):1})
    x_01 = poly2({(0,1,0,0):1})
    x_10 = poly2({(0,0,1,0):1})
    x_11 = poly2({(0,0,0,1):1})

    zero1 = poly1({})
    zero2 = poly2({})

    w2_10 = w2.substitute({ "x_00" : x_0, "x_10" : x_1, "x_01":zero1, "x_11":zero1 })
    assert w2_10 == w1

    w2_01 = w2.substitute({ "x_00" : x_0, "x_10" : zero1, "x_01":x_1, "x_11":zero1 })
    assert w2_01 == w1

    w2_11 = w2.substitute({ "x_00" : x_0, "x_10" : zero1, "x_01":zero1, "x_11":x_1 })
    assert w2_11 == w1

    q0 = w1.substitute({
        "x_0" : x_00,
        "x_1" : x_10,
    })
    q1 = w1.substitute({
        "x_0" : x_01,
        "x_1" : x_11,
    })
    #print(q0+q1)

    print()



if __name__ == "__main__":

    main()


