#!/usr/bin/env python3

from operator import matmul
from functools import reduce

import numpy

from bruhat import element
from bruhat.action import mulclose
from bruhat.element import CyclotomicRing, test
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
        #ring = CyclotomicRing(4)
        ring = element.Z
        w = -1
    else:
        ring = CyclotomicRing(d)
        w = ring.x

    I = numpy.zeros((d, d), dtype=object)
    wI = numpy.zeros((d, d), dtype=object)
    X = numpy.zeros((d, d), dtype=object)
    Z = numpy.zeros((d, d), dtype=object)
    Zs = numpy.zeros((d, d), dtype=object)
    S = numpy.zeros((d, d), dtype=object)
    Sdag = numpy.zeros((d, d), dtype=object)
    for j in range(d):
        I[j, j] = 1
        wI[j, j] = w
        X[j, (j+1)%d] = 1
        Z[j, j] = w**j
        Zs[j, j] = w**(d-j)
        S[j, j] = w**(j**2)
        Sdag[j, j] = w**(d - (j**2)%d)

    qu = Space(d, ring)
    hom = Hom(qu, qu)

    I = Map.from_array(I, hom)
    wI = Map.from_array(wI, hom)
    X = Xs = Map.from_array(X, hom)
    Z = Map.from_array(Z, hom)
    Zs = Map.from_array(Zs, hom)
    S = Map.from_array(S, hom)
    Sdag = Map.from_array(Sdag, hom)
    Y = w*X*Z # ?

    assert S*Sdag == I

    pauli = mulclose([X, Z])
    pauli = set(pauli)
    print(len(pauli))

    if d>2:
        for j in range(1, d+1):
            assert (j==d) == (X**j==I)
            assert (j==d) == (Z**j==I)
    else:
        assert X != I
        assert Z != I
        assert X*X == I
        assert Z*Z == I

    assert Z*X != X*Z

    if d<6:
        # slooooooow
        lhs = atpow(X, d)
        rhs = atpow(Z, d)
        assert lhs * rhs == rhs * lhs
    
    lhs = X@Xs
    rhs = Z@Zs
    assert lhs * rhs == rhs * lhs
    
    print(Y)
    print(S*X*Sdag)
    print(Y == S*X*Sdag)
    if d==3:
        assert (Y == S*X*Sdag)
    
    op = S*X*Sdag
    for g in pauli:
        if op == g:
            print("found")
            break
    else:
        print("not found")
    

    def is_cliff(A, Adag):
        for g in pauli:
            h = A*g*Adag
            if h not in pauli:
                return False
        return True

    print("is_cliff:", (S in pauli), is_cliff(S, Sdag))

    #def is_third_level(A, Adag):
        


if __name__ == "__main__":

    main()


