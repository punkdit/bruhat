#!/usr/bin/env python

import sys, os

from solve import array2, zeros2, dot2, shortstrx, eq2

def main():

    I = array2([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    A = array2([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    B = array2([[1, 0, 0], [0, 1, 1], [0, 0, 1]])
    C = array2([[1, 0, 0], [1, 1, 0], [0, 0, 1]])
    D = array2([[1, 0, 0], [0, 1, 0], [0, 1, 1]])

    assert eq2(dot2(A, A), I)
    assert eq2(dot2(B, B), I)
    assert eq2(dot2(C, C), I)
    assert eq2(dot2(D, D), I)
    
    assert eq2(dot2(A, D), dot2(D, A))

    assert eq2(dot2(A, C, A), dot2(C, A, C))
    assert eq2(dot2(B, D, B), dot2(D, B, D))
    assert eq2(dot2(A, B, A, B), dot2(B, A, B, A))
    assert eq2(dot2(C, D, C, D), dot2(D, C, D, C))



if __name__ == "__main__":

    main()


