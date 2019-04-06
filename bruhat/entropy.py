#!/usr/bin/env python3

import sys
from math import log

import numpy

#from bruhat.gelim import row_reduce, shortstr, kernel
from qupy.dev import linalg
from bruhat.argv import argv

EPSILON = 1e-8

def W(*items):

    sitems = sum(items)
    r = 0.
    for n in items:
        r += n * log(n)
    return -1*(r - sitems*log(sitems))


if 0:
    def getname(*items):
        items = list(items)
        items.sort(reverse=True)
        return "W(%s)"%(','.join(str(i) for i in items))
    
    def namerank(key):
        assert key[0] == "W"
        i, j = eval(key[1:])
        return (i+j, -i, -j)

else:
    def getname(*items):
        items = list(items)
        items.sort(reverse=True)
        return tuple(items)
    
    def namerank(key):
        i, j = key
        return (i+j, -i, -j)


def fstr(x):
    s = "%.4f"%x 
    s = s.rjust(8)
    return s


def astr(A):
    m, n = A.shape
    lines = []
    for i in range(m):
        row = ' '.join(fstr(x) for x in A[i, :])
        lines.append(row)
    return '\n'.join(lines)


N0 = argv.get("N", 9)

for N in range(4, N0+1):
    
    # build a linear system of equations
    rows = []
    keys = set()
    
    #keys.add((1, 1))
    #rows.append({(1, 1) : 1})
    
    for n in range(4, N+1):
    
        # a + b + c == n
        triples = set()
        for a in range(1, n):
         for b in range(1, n):
          for c in range(1, n):
            if a+b+c != n:
                continue
            items = [a, b, c]
            items.sort(reverse=True)
            triples.add(tuple(items))
        triples = list(triples)
        triples.sort()
    #    print("n =", n)
    #    print(triples)
    
        for a, b, c in triples:
    
            row = {}
            lhs = [getname(a, b), getname(a+b, c)]
            lhs.sort()
            rhs = [getname(a, b+c), getname(b, c)]
            rhs.sort()
    
            if lhs == rhs:
                continue
    
            #print(lhs, "=", rhs)
    
            for key in lhs:
                row[key] = 1
                keys.add(key)
    
            for key in rhs:
                row[key] = -1
                keys.add(key)
        
    #        print(row)
            rows.append(row)
    
    for key in keys:
        i, j = key
        m = 2
        while 1:
            keym = (m*i, m*j)
            if keym not in keys:
                break
            row = {key : m, keym : -1}
            #print(key, keym, row)
            rows.append(row)
            m += 1
    
    keys = list(keys)
    keys.sort(key = namerank)
    #print("keys:", keys)
    
    m = len(rows)
    n = len(keys)
    A = numpy.zeros((m, n), dtype=object)
    for i, row in enumerate(rows):
        for j, key in enumerate(keys):
            value = row.get(key, 0)
            A[i, j] = value
    
    print("A.shape:", A.shape)
    print(A)

    v = numpy.array([W(i, j) for (i, j) in keys])

    #print(numpy.dot(A, v))

    A = A.astype(float)
    B = linalg.row_reduce(A, truncate=True)
    print("B.shape:", B.shape)
    print(B)

    assert numpy.abs(numpy.dot(B, v)).sum() < EPSILON

    K = linalg.kernel(B)
    BK = numpy.dot(B, K)
    assert numpy.abs(BK).sum() < EPSILON
    K = K.transpose()

    K = linalg.row_reduce(K)

    print(keys)
    print("kernel:")
    print(repr(K))
    print(astr(K))

    if N==6:
        v = [1, 1, 2, 2, 1, 2, 2,2, 3]
        print(numpy.dot(B, v))

    continue
    
    # _number of free variables follows
    # https://oeis.org/A056171
    # """ the _number of unitary prime divisors of n!. A prime
    # divisor of n is unitary iff its exponent is 1 in the
    # prime power factorization of n. In general, gcd(p, n/p) = 1 or p.
    # Here we count the cases when gcd(p, n/p) = 1. """
    # Starting at N=4:
    # 1, 2, 1, 2, 2, 2, 1, 2, 2, 3, 2, 2, 2, 3, 3, 4, 4, 4, 3, 4, 4, ...
    print("N = %d, free vars = %d" % (N, B.shape[1]-B.shape[0]))
    #print("N = %d" % N)

    s = ' '.join(str(k) for k in keys)
    s = s.replace(', ', ',')
    print(s)
    
    K = kernel(B)
    K = row_reduce(K)
    
    B = B.astype(int)
    K = K.astype(int)
    
    print(shortstr(K))

    if 0:
        
        
        leading = []
        m, n = K.shape
        j = 0
        for i in range(m):
            while K[i,j]==0 and j<n:
                #print("K[%d, %d] = %d" % (i, j, K[i, j]))
                j += 1
            #print(i, j)
            val = K[i, j]
            if val < 0:
                K[i, :] *= -1
            #print("K[%d, %d] = %d" % (i, j, K[i, j]))
            assert n>j>=i
            leading.append(j)
        
        #print("kernel:")
    
    if 0:
        
        # Rewrite kernel so that every entry is positive
        # while maintaining row reduced form. This means we 
        # can add lower rows to higher rows only.
        i = 0 # row
        j = 0 # col
        while i < m and j < n:
        
            if 0:
                print(i, j)
                print(shortstr(K))
        
            # loop invariant:
            for i0 in range(i):
              for j0 in range(j):
                assert K[i0, j0] >= 0, (i0, j0)
        
            # look for a positive entry in this column
            i0 = i
            while i0 >= 0 and K[i0, j] <= 0:
                if K[i0, j] < 0:
                    assert 0
                i0 -= 1
        
            if i0 >= 0:
                assert K[i0, j] > 0
            
            i1 = i0-1
            while i1 >= 0:
                while K[i1, j] < 0: # HACK
                    K[i1, :] += K[i0, :]
                i1 -= 1
        
        
            if i+1 < m:
                i += 1
                j += 1
            else:
                j += 1
                
            
        print("kernel:")
        s = shortstr(K)
        s = s.replace(" ", "")
        print(s)
        
        BK = numpy.dot(B, K.transpose())
        assert numpy.abs(BK).sum() == 0
        
        #B = B.astype(int)
    
    print(flush=True)


