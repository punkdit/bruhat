#!/usr/bin/env python3


import numpy

from bruhat.element import Z
from bruhat.vec import Space, Hom, Map


def main():

    #U = numpy.array([[0, 1], [-1, -1]])

    V = Space(2, Z)
    hom = Hom(V, V)
    U = Map.from_array([[0, 1], [-1, -1]], hom)

    A = Map.from_array([[1, 0], [0, 1]], hom)

    A0 = A
    for i in range(2):
        print(A)
        A = U * A * U.transpose()
        A0 = A0 + A
    print(A)
    print("sum:")
    print(A0)


if __name__=="__main__":
    main()


