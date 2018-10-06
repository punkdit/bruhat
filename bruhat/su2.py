#!/usr/bin/env python3

"""
Finite subgroups of SU(2).
"""

import sys, os

import numpy

from element import Linear, Z, Q
from action import mulclose

field = Q
M = Linear(4, field)


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


def get_order(g):
    n = 1
    a = g
    while a*a != a: # identity
        a = a*g
        n += 1

    return n


def test():

    assert i*i == -e
    assert j*j == -e
    assert k*k == -e
    for a in [i, j, k]:
      for b in [i, j, k]:
        if a!=b:
            assert a*b == -b*a

    assert i*j*k == -e

    one = field.one

    Q_8 = mulclose([i, j, k])
    assert len(Q_8)==8

    items = [i, -i, j, -j, k, -k]
    half = one/2
    for a in [-one, one]:
     for b in [-one, one]:
      for c in [-one, one]:
       for d in [-one, one]:
        items.append(half * quaternion(a, b, c, d))

    tetra = mulclose(items)
    assert len(tetra)==24
        

if __name__ == "__main__":

    test()
    
