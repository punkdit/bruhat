#!/usr/bin/env python3

"""
The Heisenberg-Weyl group & clifford hierarchy for qudits.

"""

from operator import matmul
from functools import reduce

import numpy

from bruhat import element
from bruhat.action import mulclose
from bruhat.util import cross
from bruhat.element import CyclotomicRing, test, FiniteField, PolynomialRing
from bruhat.vec import Space, Hom, Map
from bruhat.argv import argv

def atpow(A, n):
    return reduce(matmul, [A]*n)


def astr(v):
    vs = numpy.empty(v.shape, dtype=object)
    for idx in numpy.ndindex(v.shape):
        vs[idx] = str(v[idx])
    s = str(vs)
    s = s.replace("'", "")
    return s


def main():

    d = argv.get("d", 3)
    assert d>=2

    if d==2:
        ring = CyclotomicRing(4)
        #ring = element.Z
        gamma = ring.x
        w = gamma**2
        assert w == -1
        argv.double = True
    elif 1:
        ring = CyclotomicRing(d**2)
        gamma = ring.x
        w = gamma**d
    else:
        ring = CyclotomicRing(d)
        w = ring.x

    I = numpy.zeros((d, d), dtype=object)
    wI = numpy.zeros((d, d), dtype=object)
    X = numpy.zeros((d, d), dtype=object)
    Z = numpy.zeros((d, d), dtype=object)
    Zdag = numpy.zeros((d, d), dtype=object)
    S = numpy.zeros((d, d), dtype=object)
    Sdag = numpy.zeros((d, d), dtype=object)
    T = numpy.zeros((d, d), dtype=object)
    Tdag = numpy.zeros((d, d), dtype=object)
    for j in range(d):
        I[j, j] = 1
        wI[j, j] = w
        X[j, (j+1)%d] = 1
        Z[j, j] = w**j
        Zdag[j, j] = w**(d-j)

        if d==2:
            val = gamma**(j**2)
            ival = gamma**(d**2 - (j**2)%(d**2))
        else:
            val = w**(j**2)
            ival = w**(d - (j**2)%d)
        assert val*ival == 1

        S[j, j] = val
        Sdag[j, j] = ival

        if d in [2, 3, 6]:
            val = gamma**(j**3)
            ival = gamma**(d**2 - (j**3)%(d**2))
        else:
            val = w**(j**3)
            ival = w**(d - (j**3)%d)
        assert val*ival == 1

        T[j, j] = val
        Tdag[j, j] = ival

    qu = Space(d, ring)
    hom = Hom(qu, qu)

    I = Map.from_array(I, hom)
    wI = Map.from_array(wI, hom)
    Xdag = Map.from_array(X.transpose(), hom)
    X = Xs = Map.from_array(X, hom)
    Z = Map.from_array(Z, hom)
    Zdag = Map.from_array(Zdag, hom)
    S = Map.from_array(S, hom)
    Sdag = Map.from_array(Sdag, hom)
    T = Map.from_array(T, hom)
    Tdag = Map.from_array(Tdag, hom)

    Y = w*X*Z # ?

    assert S*Sdag == I

    assert Z*Zdag == I
    assert X*Xdag == I

    if d>2:
        for j in range(1, d+1):
            assert (j==d) == (X**j==I)
            assert (j==d) == (Z**j==I)
    else:
        assert X != I
        assert Z != I
        assert X*X == I
        assert Z*Z == I

    assert Z*X == (w**(d-1))*X*Z

    if argv.double:
        pauli = mulclose([X, Z, gamma*I])
    else:
        pauli = mulclose([X, Z])
    pauli = set(pauli)
    print("pauli:", len(pauli))
    assert Zdag in pauli
    assert Xdag in pauli

    if d<6:
        # slow..
        lhs = atpow(X, d)
        rhs = atpow(Z, d)
        assert lhs * rhs == rhs * lhs
    
    lhs = X@X
    rhs = Z@Zdag
    assert lhs * rhs == rhs * lhs

    lhs = X@(X**(d-1))
    rhs = Z@Z
    assert lhs * rhs == rhs * lhs

    if 0:
        n = 5
        ops = [X, Xdag]
        lhs = []
        for op in cross([ops]*n):
            op = reduce(matmul, op)
            lhs.append(op) 
    
        ops = [Z, Zdag]
        rhs = []
        for op in cross([ops]*n):
            op = reduce(matmul, op)
            rhs.append(op) 
    
        for l in lhs:
          for r in rhs:
            if l*r == r*l:
                print("/", end=" ", flush=True)
            else:
                print(".", end=" ", flush=True)
        print()

    #print(Z@Z@Z)
    
    #print(Y)
    #print(S*X*Sdag)
    #print(Y == S*X*Sdag)
    if d==3:
        assert (Y == S*X*Sdag)
    
    inverse = {}
    for g in pauli:
      for h in pauli:
        if g*h == I:
            inverse[g] = h
            inverse[h] = g

    def is_cliff(A, Ai):
        assert A*Ai == I
        for g in pauli:
            h = A*g*Ai
            if h not in pauli:
                return False
        return True
    print(S)
    print(Sdag)
    assert is_cliff(S, Sdag)

    #print("is_cliff:", (S in pauli), is_cliff(S, Sdag))

    def is_third_level(A, Ai):
        assert A*Ai == I
        for g in pauli:
            gi = inverse[g]
            h = A*g*Ai*gi
            hi = g*A*gi*Ai
            if not is_cliff(h, hi):
                return False
        return True

    print("is_pauli(S)", S in pauli)
    print("is_cliff(S)", is_cliff(S, Sdag))
    #print("is_third_level(S)", is_third_level(S, Tdag))

    print("is_pauli(T)", T in pauli)
    print("is_cliff(T)", is_cliff(T, Tdag))
    print("is_third_level(T)", is_third_level(T, Tdag))
        

    print("OK")


def test():

    d = argv.get("d", 3)
    r = argv.get("r", 3)

    field = FiniteField(d)
    ring = PolynomialRing(field)

    funcs = []

    one = ring.one
    x = ring.x
    f = x**r # if r==d and is prime then this is the Frobenius: it's the identity function on the field.

    for i in range(d):
        print(f(i))




if __name__ == "__main__":

    main()
    #test()


