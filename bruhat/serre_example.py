#!/usr/bin/env python3

"""
Exercise 13.4 from Serre "Linear Representations of Finite Groups ", 1977.
"""

import sys, os

from element import *
from action import mulclose

Z = IntegerRing()
M = Linear(4, Z)


def quaternion(a, b, c, d):
    # build matrix representation of quaternion
    A = M.get([
        [a, -b, -c, -d],
        [b, a, -d, c],
        [c, d, a, -b],
        [d, -c, b, a]])
    return A


def test():
    e = quaternion(1, 0, 0, 0)
    i = quaternion(0, 1, 0, 0)
    j = quaternion(0, 0, 1, 0)
    k = quaternion(0, 0, 0, 1)

    assert i*i == -e
    assert j*j == -e
    assert k*k == -e
    for a in [i, j, k]:
      for b in [i, j, k]:
        if a!=b:
            assert a*b == -b*a

    assert i*j*k == -e

    a = i+j+k-e
    assert a!=e
    assert a**2!=e
    assert a**3==8*e

    Q_8 = mulclose([i, j, k])
    assert len(Q_8)==8

    print(a**2)

    b = a**2
    print(b/3)

    if 0:
        search = [-1, 0, 1]
        for xe in search:
         for xi in search:
          for xj in search:
           for xk in search:
            if abs(xi)+abs(xj)+abs(xk)==0:
                continue
            a = xe*e + xi*i + xj*j + xk*k
            assert a!=0
            b = a*a # a**2
            if b==0 or b==e or b==8*e:
                continue
            b = b*a # a**3
            if b==e or b==8*e:
                print("found:", xe, xi, xj, xk)
      


if __name__ == "__main__":

    test()
    
