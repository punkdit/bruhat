#!/usr/bin/env python3

"""
Kirchhoff's theorem:  You take the Laplacian matrix of
the graph (degree matrix - adjacency matrix), delete
an arbitrary row and its corresponding column, and then
find the determinant of the matrix. The value of the
determinant will be the number of spanning trees for the graph.
"""

import numpy

from bruhat import element
from bruhat.poly import Poly
from bruhat.util import allsignedperms
from bruhat.argv import argv


def get_submatrix(L, i):
    L1 = numpy.concatenate((L[:i], L[i+1:]))
    L2 = numpy.concatenate((L1[:, :i], L1[:, i+1:]), axis=1)
    return L2


def det(A, ring):
    n = len(A)
    assert A.shape == (n, n)

    items = list(range(n))
    result = ring.zero
    for sgn, perm in allsignedperms(items):
        x = sgn*ring.one
        for idx, i in enumerate(perm):
            x *= A[idx, i]

        #if x != 0:
        #    print(x)
        result += x

    return result


class Graph(object):
    "undirected graph"

    def __init__(self, n):
        self.A = numpy.zeros((n, n), dtype=int)
        self.n = n

    def add_edge(self, i, j, directed=False):
        A = self.A
        n = self.n
        A[i, j] = 1
        if not directed:
            A[j, i] = 1
    
    def get_laplacian(self):
        A = self.A
        n = self.n
    
        D = numpy.zeros((n, n), dtype=int) # degree matrix
        for i in range(n):
            D[i, i] = A[i].sum()
    
        return D - A

    def get_poly(self, ring=element.Z):
        A = self.A
        n = self.n

        B = numpy.zeros((n, n), dtype=object)
        ijs = numpy.transpose(numpy.nonzero(A))
        for (i, j) in ijs:
            assert i!=j
            ii, jj = (j, i) if i>j else (i,j)
            #ii, jj = i, j
            name = "a[%d,%d]"%(ii, jj)
            p = Poly(name, ring)
            B[i, j] = -p

        for i in range(n):
            B[i,i] = -B[i].sum()

        print("B =")
        print(B)
        #for row in B:
        #  for col in row:
        #    print("%6s"%col, end=" ")
        #  print()

        B = get_submatrix(B, 0)

        trees = det(B, ring)
        return trees


    def get_nspan(self):
        A = self.A
        n = self.n
    
        L = self.get_laplacian()
        #print(L)
    
        result = None
        for i in range(n):
            L2 = get_submatrix(L, i)
            r = numpy.linalg.det(L2)
            r = int(round(r))
            assert result is None or r==result, (r, result)
            result = r

        #print("same:", r, det(L2, element.Z))
    
        return r



def main():

    G = Graph(3)
    G.add_edge(0, 1)
    G.add_edge(1, 2)
    G.add_edge(2, 0)
    p = G.get_poly()
    assert len(p) == 3

    G = Graph(5)
    for (i, j) in [(0, 1), (1, 3), (3, 2), (2, 0), (3, 4)]:
        G.add_edge(i, j)

    n = G.get_nspan()
    assert n==4

    G = Graph(7)
    for (i, j) in [(0, 1), (1, 3), (3, 2), (2, 0), (3, 4),
        (4, 5), (5, 6), (6, 3)]:
        G.add_edge(i, j)

    n = G.get_nspan()
    assert n==4*4

    G = Graph(6)
    for (i, j) in [(0, 1), (1, 3), (3, 2), (2, 0),
        (2, 4), (4, 5), (5, 3)]:
        G.add_edge(i, j)

    n = G.get_nspan()
    assert n==15, n

    p = G.get_poly()
    print(p)

    print("OK")



if __name__ == "__main__":

    main()

    


