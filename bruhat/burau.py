#!/usr/bin/env python3

"""
Calculate the hyper-elliptic representation of
the braid group.
This is the reduced Burau representation at t=-1.
"""

import numpy


from argv import argv
import element
import vec
from action import mulclose


def zeros(m):
    return numpy.zeros((m, m), dtype=int)

from_array = vec.Map.from_array

def main():

    p = argv.p
    if p is None:
        ring = element.Z
    else:
        ring = element.FiniteField(p)

    n = argv.get("n", 3)
    assert n%2 == 1
    assert n>2
    m = n-1

    space = vec.Space(m, ring)
    hom = vec.Hom(space, space)

    J = zeros(m)
    for i in range(m):
        if i%2==0:
            J[i+1, i] = -1
        else:
            J[i-1, i] = 1
    J = from_array(J, hom)
    print("J =")
    print(J)

    make = lambda items : from_array(numpy.array(items), hom)
    is_symplectic = lambda A : (A.transpose() * J * A == J)


    if 0:
        a1 = make([
            [1, -1, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]])
    
        assert is_symplectic(a1)
    
        b1 = make([
            [1, 0, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]])
    
        assert is_symplectic(b1)
    
        a2 = make([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, -1],
            [0, 0, 0, 1]])
    
        assert is_symplectic(a2)
    
        b2 = make([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 1, 1]])
    
        assert is_symplectic(b2)
    
        c = make([
            [1, 0, 0, 0],
            [1, 1, -1, 0],
            [0, 0, 1, 0],
            [-1, 0, 1, 1]])
    
        assert is_symplectic(c)
    
        gen = [b1, a1, c, a2, b2]


    gen = []
    I = zeros(m)
    for i in range(m):
        I[i, i] = 1
    A = zeros(m) + I
    A[1, 0] = 1
    gen.append(A)

    for i in range(m//2):
        A = zeros(m) + I
        A[2*i, 2*i+1] = -1
        gen.append(A)

        if i+1 < m//2:
            A = zeros(m) + I
            A[2*i+1, 2*i] = 1
            A[2*i+3, 2*i] = -1
            A[2*i+1, 2*i+2] = -1
            A[2*i+3, 2*i+2] = 1
            gen.append(A)

    if m//2>1:
        A = zeros(m) + I
        A[m-1, m-2] = 1
        gen.append(A)

    gen = [make(A) for A in gen]

    print("gen:", len(gen))
    for A in gen:
        print(A)
        assert is_symplectic(A)

    for i in range(m):
      for j in range(i, m):
        a, b = gen[i], gen[j]
        if abs(i-j)==1:
            assert a*b*a == b*a*b # Yang-Baxter
        elif abs(i-j)>1:
            assert a*b == b*a

    A = zeros(m) + I
    A[3, 2] = 1
    A = make(A)
    assert is_symplectic(A)
    gen.append(A)

    if argv.p==2 and m==4 or argv.mulclose:
        assert p is not None, "not a finite group!"
        G = mulclose(gen)
        print("|G| =", len(G))
        if argv.p==2 and m==4:
            assert len(G)==720


if __name__ == "__main__":

    main()




