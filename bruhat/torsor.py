#!/usr/bin/env python3

"""
Find monoid, group, torsor using abstract tables.
"""

from bruhat.util import cross

import numpy


def table(*shape):
    return numpy.zeros(shape, dtype=int)

def all_mul(n):
    itemss = [tuple(range(n))]*(n-1)*(n-1)
    for idxs in cross(itemss):
        M = table(n, n)
        for i in range(n):
            M[0, i] = i # identity is always zero
            M[i, 0] = i # identity is always zero
        for i in range(n-1):
          for j in range(n-1):
            M[1+i, 1+j] = idxs[i + (n-1)*j]
        yield M


def is_assoc(M):
    n = M.shape[0]
    for a in range(n):
      for b in range(n):
        for c in range(n):
            if M[M[a, b], c] != M[a, M[b, c]]:
                return False
    return True


#def is_group(M):
    


def main():

    n = 3

    count = 0
    for M in all_mul(n):
        if not is_assoc(M):
            continue
        print(M)
        count += 1
    print("found:", count)



if __name__=="__main__":

    main()

