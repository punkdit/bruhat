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

        print("det:", numpy.linalg.det(A))
        A = from_array(A, hom)
        print("A=")
        print(A)
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

    for a in gen:
        print(a)

    #G = mulclose(gen)

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




