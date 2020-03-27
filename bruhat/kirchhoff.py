#!/usr/bin/env python3

"""
Kirchhoff's theorem:  You take the Laplacian matrix of
the graph (degree matrix - adjacency matrix), delete
an arbitrary row and its corresponding column, and then
find the determinant of the matrix. The value of the
determinant will be the number of spanning trees for the graph.
"""

import numpy

from bruhat.argv import argv


#class Graph(object):
    

def get_laplacian(A):
    shape = A.shape
    n = shape[0]
    assert shape == (n, n)

    D = numpy.zeros((n, n), dtype=int) # degree matrix
    for i in range(n):
        D[i, i] = A[i].sum()

    return D - A


def get_nspan(A):

    shape = A.shape
    n = shape[0]
    L = get_laplacian(A)
    #print(L)

    result = None
    for i in range(n):
        L1 = numpy.concatenate(
            (L[:i], L[i+1:]))
        #print(L1)
        L2 = numpy.concatenate(
            (L1[:, :i], L1[:, i+1:]), axis=1)
        #print(L2)
    
        r = numpy.linalg.det(L2)
        r = int(round(r))
        assert result is None or r==result, (r, result)
        result = r

    return r



def main():

    n = 5
    A = numpy.zeros((n, n), dtype=int)
    for (i, j) in [(0, 1), (1, 3), (3, 2), (2, 0), (3, 4)]:
        A[i, j] = 1.
        A[j, i] = 1.

    n = get_nspan(A)
    assert n==4

    n = 7
    A = numpy.zeros((n, n), dtype=int)
    for (i, j) in [(0, 1), (1, 3), (3, 2), (2, 0), (3, 4),
        (4, 5), (5, 6), (6, 3)]:
        A[i, j] = 1.
        A[j, i] = 1.

    n = get_nspan(A)
    assert n==4*4

    n = 6
    A = numpy.zeros((n, n), dtype=int)
    for (i, j) in [(0, 1), (1, 3), (3, 2), (2, 0),
        (2, 4), (4, 5), (5, 3)]:
        A[i, j] = 1.
        A[j, i] = 1.

    n = get_nspan(A)
    assert n==15, n



if __name__ == "__main__":

    main()

    


