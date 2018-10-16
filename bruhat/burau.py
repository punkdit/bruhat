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

    assert a1*b1*a1 == b1*a1*b1
    assert a2*b2*a2 == b2*a2*b2

    assert b1*c == c*b1
    assert a1*c*a1 == c*a1*c

    gen = [b1, a1, c, a2, b2]

    assert c*b2 == b2*c

    return

    # getting desparate....
    Js = []
    for i in range(m-1):
        A = zeros(m)
        for j in range(m):
            if j==i:
                A[i, i+1] = 1
            elif j==i+1:
                A[i+1, i] = 1
            else:
                A[j, j] = 1
        A = from_array(A, hom)
        Js.append(A)
    Js = mulclose(Js)
    #print(len(Js))

    gen = []
    ngen = []
    for i in range(m):

        A = zeros(m)
#        for j in range(m):
#            A[j, j] = 1
#
#        if i>0:
#            A[i-1, i] = -1
#
#        if i+1 < m:
#            A[i+1, i] = 1

        for j in range(m):
            if j==i-1:
                A[i-1, j] = 1
                A[i, j] = 1
            elif j==i:
                A[i, j] = 1
            elif j==i+1:
                A[i, j] = -1
                A[i+1, j] = 1
            else:
                A[j, j] = 1

        #print("det:", numpy.linalg.det(A))
        ngen.append(A)
        A = from_array(A, hom)
        #print("A=")
        #print(A)
        #print()
        #print(A.transpose() * J * A)
        #print(A * J * A.transpose())
        #assert A.transpose() * J * A == J
        gen.append(A)

    print("gen:", len(gen))

    for i in range(m):
      for j in range(i, m):
        a, b = gen[i], gen[j]
        if abs(i-j)==1:
            assert a*b*a == b*a*b # Yang-Baxter
        elif abs(i-j)>1:
            assert a*b == b*a

    #G = mulclose(gen)

    if 0:
        from solve import Unknown, System, dot2
        U = Unknown(m, m)
        system = System(U)
        for A in ngen:
            system.append(dot2(A, dot2(U, A.transpose())), A)
        V = system.solve()
        print(V)

    return

    for J in Js:
        count = 0
        for A in gen:
            if A.transpose() * J * A == J:
                count += 1
        if count==len(gen):
            print("found:")
            print(J)
        else:
            print(count, end=" ")
    print()

    if argv.generate:
        assert p is not None, "not a finite group!"
        G = mulclose(gen)
        print("|G| =", len(G))


if __name__ == "__main__":

    main()




