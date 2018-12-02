#!/usr/bin/env python3

"""
Algebraic groups: matrix groups over Z/pZ.

incomplete / abandoned for now.
"""


import sys, os
import random

import numpy

scalar = numpy.int64

from bruhat.action import mulclose
from bruhat.spec import isprime
from bruhat.argv import argv


class Op(object):
    def __init__(self, A, p):
        self.A = A.astype(scalar) # copy
        n = A.shape[0]
        assert A.shape == (n, n)
        assert int(p) == p
        assert p>=0
        self.p = p
        self.n = n
        if p>0:
            self.A %= p
        self.key = (self.p, self.A.tostring())
        self._hash = hash(self.key)

    def __str__(self):
        return str(self.A)

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        assert self.p == other.p
        return self.key == other.key

    def __ne__(self, other):
        assert self.p == other.p
        return self.key != other.key

    def __add__(self, other):
        A = self.A + other.A
        return Op(A, self.p)

    def __sub__(self, other):
        A = self.A - other.A
        return Op(A, self.p)

    def __neg__(self, other):
        A = -self.A
        return Op(A, self.p)

    def __mul__(self, other):
        A = numpy.dot(self.A, other.A)
        return Op(A, self.p)


# https://math.stackexchange.com/questions/34271/
# order-of-general-and-special-linear-groups-over-finite-fields

def order_gl(n, q):
    order = 1
    for i in range(n):
        order *= (q**n - q**i)
    return order

def order_sl(n, q):
    order = order_gl(n, q)
    assert order % (q-1) == 0
    return order//(q-1)


def order_sp(n, q):
    # n = 2*m
    assert n%2==0
    m = n//2
    N = q**(m**2)
    for i in range(m):
        N *= (q**(2*(i+1)) - 1)
    return N

assert order_sp(2, 2)==6     # 3!
assert order_sp(4, 2)==720   # 6!


def SL(n, p):
    "special linear group"
    assert int(n)==n
    assert int(p)==p
    assert n>0
    assert isprime(p)

    I = numpy.identity(n, scalar)
    gen = []
    for i in range(n):
        for j in range(n):
            if i==j:
                continue
            A = I.copy()
            A[i, j] = 1
            gen.append(Op(A, p))
#    for op in gen:
#        print(op)
#        print(numpy.linalg.det(op.A))
    order = order_sl(n, p)
#    G = mulclose(gen, maxsize=order)
    G = mulclose(gen)
    assert len(G)==order
    return G


def main():

    n = argv.get("n", 2)
    p = argv.get("p", 2)

    print("|GL(%d, %d)| = %d"%(p, n, order_gl(n, p)))
    print("|SL(%d, %d)| = %d"%(p, n, order_sl(n, p)))
    if n%2 == 0:
        print("|Sp(%d, %d)| = %d"%(p, n, order_sp(n, p)))

    if argv.SL:
        G = SL(n, p)
        print("|G| =", len(G))




if __name__ == "__main__":
    main()

    

    




