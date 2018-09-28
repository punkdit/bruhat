#!/usr/bin/env python3

"""
Exercise 13.4 from Serre "Linear Representations of Finite Groups ", 1977.
"""

import sys, os

import numpy

from element import Linear, Z, Q
from action import mulclose

M = Linear(4, Q)


def quaternion(a, b, c, d):
    # build matrix representation of quaternion
    A = M.get([
        [a, -b, -c, -d],
        [b, a, -d, c],
        [c, d, a, -b],
        [d, -c, b, a]])
    return A


e = quaternion(1, 0, 0, 0)
i = quaternion(0, 1, 0, 0)
j = quaternion(0, 0, 1, 0)
k = quaternion(0, 0, 0, 1)
basis = [e, i, j, k]


def dot(a, b):
    a = numpy.array(a.value)
    b = numpy.array(b.value)
    c = a*b
    return c.sum()


def get_rep(A, left=True):
    Arep = []
    for v in basis:
        row = []
        for u in basis:
            if left:
                B = A*u # left action
            else:
                B = u*A # right action
            r = dot(B, v)
            row.append(r/4)
        Arep.append(row)
    Arep = M.get(Arep)
    return Arep


def test():

    assert i*i == -e
    assert j*j == -e
    assert k*k == -e
    for a in [i, j, k]:
      for b in [i, j, k]:
        if a!=b:
            assert a*b == -b*a

    assert i*j*k == -e

    one = Q.one
    A = (one/2)*(i+j+k-e)
    assert A!=e
    assert A**2!=e
    assert A**3==e

    Q_8 = mulclose([i, j, k])
    assert len(Q_8)==8

    # Q_8 acts by right multiplication, C_3 by left multiplication

    Arep = get_rep(A)

    Qrep = [get_rep(V, False) for V in [i, j, k]]
    for V in Qrep:
        #print(V)
        assert V*Arep == Arep*V

    G = mulclose(Qrep + [Arep])
    assert len(G) == 24

    chi = []
    G = list(G)
    G.sort(key = get_order)
    for g in G:
        print(str(get_order(g)).rjust(3), end=" ")
        chi.append(g.trace())
    print()
    for x in chi:
        print(str(x).rjust(3), end=" ")
    print()
        

def get_order(g):
    n = 1
    a = g
    while a*a != a: # identity
        a = a*g
        n += 1

    return n


if __name__ == "__main__":

    test()
    
