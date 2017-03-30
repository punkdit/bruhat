#!/usr/bin/env python

"""
Gaussian elimination over the rationals.
"""

import sys, os
from random import randint, seed
from fractions import Fraction


import numpy
from numpy import dot

from smap import SMap
from argv import argv


def write(s):
    sys.stdout.write(str(s)+' ')
    sys.stdout.flush()


def fstr(x):
    x = Fraction(x)
    a, b = x.numerator, x.denominator
    if b==1:
        return str(a)
    if a==0:
        return "."
    return "%s/%s"%(a, b)


def shortstr(*items, **kw):

    if len(items)>1:
        return shortstrx(*items, **kw)

    A = items[0]
    if type(A) in [list, tuple]:
        A = array(A)

    shape = kw.get("shape")
    if shape:
        A = A.view()
        A.shape = shape

    if len(A.shape)==1:
        A = A.view()
        A.shape = (1, len(A))
    m, n = A.shape
    items = {}
    dw = 3
    for i in range(m):
      for j in range(n):
        x = A[i, j]
        s = fstr(x)
        dw = max(dw, len(s)+1)
        items[i, j] = s

    smap = SMap()
    for i in range(m):
      smap[i, 0] = "["
      smap[i, n*dw+1] = "]"
      for j in range(n):
        s = items[i, j]
        s = s.rjust(dw-1)
        smap[i, j*dw+1] = s

    return smap


def shortstrx(*items, **kw):
    smaps = [shortstr(item, **kw) for item in items]
    smap = SMap()
    col = 0
    for A in items:
        s = shortstr(A)
        smap[0, col] = s
        col += s.cols + 1

    return smap


def zeros(m, n):
    A = numpy.empty((m, n), dtype=object)
    A[:] = 0
    return A


def array(items):
    return numpy.array(items, dtype=object)


def identity(m):
    I = zeros(m, m)
    for i in range(m):
        I[i, i] = 1
    return I


def eq(A, B):
    r = numpy.abs(A-B).sum()
    return r==0


def dotx(*items):
    idx = 0
    A = items[idx]
    while idx+1 < len(items):
        B = items[idx+1]
        A = dot(A, B)
        idx += 1 
    return A


def compose(*items):
    items = list(reversed(items))
    A = dotx(*items)
    return A



#def shortstr(A):
#    return str(A)


#def shortstrx(*args):
#    return '\n'.join(str(A) for A in args)



def swap_row(A, j, k):
    row = A[j, :].copy()
    A[j, :] = A[k, :]
    A[k, :] = row

def swap_col(A, j, k):
    col = A[:, j].copy()
    A[:, j] = A[:, k]
    A[:, k] = col


def row_reduce(A, truncate=False, inplace=False, check=False, verbose=False):
    """ Remove zero rows if truncate==True
    """

    assert len(A.shape)==2, A.shape
    m, n = A.shape
    if not inplace:
        A = A.copy()

    if m*n==0:
        return A

    if verbose:
        print "row_reduce"
        #print "%d rows, %d cols" % (m, n)

    i = 0
    j = 0
    while i < m and j < n: 
        if verbose:
            print "i, j = %d, %d" % (i, j)
            print "A:"
            print shortstrx(A)

        assert i<=j 
        if i and check:
            assert (A[i:,:j]!=0).sum() == 0

        # first find a nonzero entry in this col
        for i1 in range(i, m):
            if A[i1, j]:
                break
        else:
            j += 1 # move to the next col
            continue # <----------- continue ------------

        if i != i1:
            if verbose:
                print "swap", i, i1
            swap_row(A, i, i1)

        assert A[i, j]
        for i1 in range(i+1, m):
            if A[i1, j]:
                if verbose: 
                    print "add row %s to %s" % (i, i1)
                r = -Fraction(A[i1, j], A[i, j])
                A[i1, :] += r*A[i, :]
                assert A[i1, j] == 0

        i += 1 
        j += 1 

    if truncate:
        m = A.shape[0]-1
        #print "sum:", m, A[m, :], A[m, :].sum()
        while m>=0 and (A[m, :]!=0).sum()==0:
            m -= 1
        A = A[:m+1, :]

    if verbose:
        print

    return A




def plu_reduce(A, truncate=False, check=False, verbose=False):
    """
    Solve PLU = A, st. P is permutation, L is lower tri, U is upper tri.
    Remove zero rows from U if truncate=True.
    """

    m, n = A.shape
    P = identity(m)
    L = identity(m)
    U = A.copy()

    assert m*n, (m, n)

    if verbose:
        print "plu_reduce:"
        print "%d rows, %d cols" % (m, n)

    i = 0
    j = 0
    while i < m and j < n:
        if verbose:
            print "i, j = %d, %d" % (i, j)
            print "P, L, U:"
            print shortstrx(P, L, U)

        assert i<=j
        if i and check:
            assert U[i:,:j].max() == 0 # XX rm

        # first find a nonzero entry in this col
        for i1 in range(i, m):
            if U[i1, j]:
                break
        else:
            j += 1 # move to the next col
            continue # <----------- continue ------------

        if i != i1:
            if verbose:
                print "swap", i, i1
            swap_row(U, i, i1)
            swap_col(P, i, i1)
            swap_col(L, i, i1)
            swap_row(L, i, i1)

        if check:
            A1 = dot(P, dot(L, U))
            assert eq(A1, A)

        r = U[i, j]
        assert r != 0
        for i1 in range(i+1, m):
            s = U[i1, j]
            if s != 0:
                t = Fraction(s, r)
                if verbose: 
                    print "add %s to %s" % (i, i1)
                L[i1, i] = t
                U[i1, :] -= t*U[i, :]
                assert U[i1, j] == 0

        if check:
            A1 = dot(P, dot(L, U))
            assert eq(A1, A)

        i += 1
        j += 1

    if check:
      for i in range(m):
        for j in range(i):
            assert U[i, j] == 0
      for i in range(m):
        for j in range(i+1, m):
            assert L[i, j] == 0

    if truncate:
        m = U.shape[0]-1
        #print "sum:", m, U[m, :], U[m, :].sum()
        while m>=0 and U[m, :].sum()==0:
            m -= 1
        U = U[:m+1, :]

    return P, L, U


def u_inverse(U, check=False, verbose=False):
    """invert a row reduced U
    """

    m, n = U.shape

    #items = []
    leading = []
    for row in range(m):
        cols = numpy.where(U[row, :])[0]
        if not len(cols):
            break
        col = cols[0]
        leading.append(col)

    assert sorted(leading) == leading
    assert len(set(leading)) == len(leading)

    U1 = zeros(n, m)

    #print shortstr(U)

    # Work backwards
    i = len(leading)-1 # <= m
    while i>=0:

        j = leading[i]
        #print "i=%d, j=%d"%(i, j)
        r = Fraction(1, U[i, j])
        U1[j, i] = r

        #print "U, U1, U*U1:"
        #print shortstrx(U, U1, dot(U, U1))

        k = i-1
        while k>=0:
            #print "dot"
            #print shortstr(U[k,:])
            #print shortstr(U1[:,i])
            r = dot(U[k, :], U1[:, i])
            #print "=", r
            if r != 0:
                j = leading[k]
                s = U[k, j]
                #print "set", j, i
                U1[j, i] = -Fraction(r, s)
            #print shortstr(U1[:,i])
            assert dot(U[k, :], U1[:, i]) == 0
            k -= 1
        i -= 1

    return U1


def l_inverse(L, check=False, verbose=False):
    """invert L (lower triangular, 1 on diagonal)
    """

    m, n = L.shape
    assert m==n
    L1 = identity(m)

    # Work forwards
    for i in range(m):
        #u = L1[:, i]
        assert L[i, i] == 1
        for j in range(i+1, m):
            #print "i=%d, j=%d"%(i, j)
            #print "L, L1, L*L1:"
            #print shortstrx(L, L1, dot(L, L1))
            r = dot(L[j, :], L1[:, i])
            #print "r =", r
            if r != 0:
                assert L1[j, i] == 0
                L1[j, i] = -r
            r = dot(L[j, :], L1[:, i])
            #print "r =", r
            #print shortstrx(L, L1, dot(L, L1))
            assert dot(L[j, :], L1[:, i]) == 0

    assert eq(dot(L, L1), identity(m))
    return L1


def pseudo_inverse(A, check=False):
    m, n = A.shape
    if m*n == 0:
        A1 = zeros(n, m)
        return A1
    P, L, U = plu_reduce(A, verbose=False, check=check)
    L1 = l_inverse(L, check=check)
    U1 = u_inverse(U, check=check)
    #print "P, L, U, PLU:"
    #print shortstr(P, L, U, dot(dot(P, L), U))
    
    A1 = dot(U1, dot(L1, P.transpose()))
    #print shortstr(dot(A1, A))
    return A1


def solve(H, u, force=False, verbose=False, check=False):
    "Solve Hv = u"
    assert len(H) == len(u)
    A = pseudo_inverse(H, check)
    #print "pseudo_inverse"
    #print shortstr(A)
    v = dot(A, u)
    if eq(dot(H, v), u) or force:
        return v


# https://en.wikipedia.org/wiki/Kernel_%28linear_algebra%29#Computation_by_Gaussian_elimination
def kernel(A, check=False, verbose=False):
    """
        return largest K such that dot(A, K) == 0.
    """

    if verbose:
        print "kernel"

    m, n = A.shape
    A, A0 = A.copy(), A

    K = identity(n)

    # Column reduce A, while duplicating operations onto K

    i = 0 # row
    for j in range(n): # col

        if verbose:
            print "A, K (i=%d, j=%d)" % (i, j)
            print shortstr(A)
            print "----------"
            print shortstr(K)
            print

        # look for a row
        while i<m and (A[i, j:]!=0).sum()==0:
            i += 1

        if i==m:
            break

        if A[i, j] == 0:
            k = j
            while A[i, k]==0:
                k += 1
            swap_col(A, j, k)
            swap_col(K, j, k)

        for k in range(j+1, n):
            r = -Fraction(A[i, k], A[i, j])
            A[:, k] += r * A[:, j]
            K[:, k] += r * K[:, j]

        i += 1
        if i==m:
            break

    if verbose:
        print "A, K (i=%d, j=%d)" % (i, j)
        print shortstr(A)
        print "----------"
        print shortstr(K)
        print

    j = K.shape[1] - 1
    while j>=0 and (A[:, j]!=0).sum() == 0:
        j -= 1
    j += 1

    #K = K[:, j+1:]
    K = K[:, j:]
    if check:
        B = dot(A0, K)
        assert numpy.abs(B).sum()==0

    return K.transpose()


def rank(A, **kw):
    A = row_reduce(A, truncate=True, **kw)
    return len(A)


def nullity(A, **kw):
    K = kernel(A, **kw)
    return len(K)


class Subspace(object):
    """ Subspace represented as the rowspace of a matrix.
    """
    def __init__(self, W):
        assert rank(W) == len(W)
        self.W = W
        self.m = W.shape[0]
        self.n = W.shape[1] # dimension

    def __len__(self):
        return self.m

    def __str__(self):
        s = shortstr(self.W)
        return str(s)

    def __eq__(self, other):
        W1 = self.W.transpose()
        W2 = other.W.transpose()
        if solve(W1, W2) is None:
            return False
        if solve(W2, W1) is None:
            return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def intersect(self, other):
        W1 = self.W
        W2 = other.W
        W = numpy.concatenate((W1, W2))
        #print "intersect"
        #print shortstr(W)
        #print
        K = kernel(W.transpose())#.transpose()
        #print "K:"
        #print shortstr(K)
        W = dot(K[:, :len(W1)], W1)
        return Subspace(W)


def rand(m, n, p=3, q=3):
    A = zeros(m, n)
    for i in range(m):
      for j in range(n):
        a = randint(-p, p)
        A[i, j] = Fraction(a, randint(1, q))
    return A



def test():

    _seed = argv.get("seed")
    if _seed is not None:
        seed(_seed)

    m, n = 4, 5
    A = zeros(m, n)
    U = zeros(m, n)

    for i in range(m):
      for j in range(i, n):
        a = randint(-3, 3)
        if i==j and a==0:
            a = 1
        U[i, j] = Fraction(a, randint(1, 3))

    V = u_inverse(U)

    L = identity(m)
    for i in range(m):
        for j in range(i):
            L[i, j] = Fraction(randint(-3, 3), randint(1, 3))
    V = l_inverse(L)

    A = rand(m, n)
    P, L, U = plu_reduce(A, check=True, verbose=False)

    #m, n = 2, 3
    #while 1:
    for i in range(100):
        A = rand(m, n)
        #A = row_reduce(A, check=True)
        if rank(A, check=True) + nullity(A, check=True) == n:
            #write('.')
            continue

        print "FAIL"
        print "A:"
        print shortstr(A)
        print "row_reduce:"
        B = row_reduce(A, verbose=True)
        print shortstr(B)
        print "kernel:"
        print shortstr(kernel(A))
        return


    A = zeros(3, 4)
    A[0, 0] = 1
    A[1, 1] = 1
    A[2, 2] = 1
    A = numpy.concatenate((A, A))
    A = A.transpose()
    #print "A:"
    #print shortstr(A)
    K = kernel(A)
    #print "K:"
    #print shortstr(K)
    assert len(K) == 3

    while 1:
        m, n = 3, 4
        A = zeros(m, n)
        for i in range(m):
          for j in range(n):
            a = randint(-2, 2)
            A[i, j] = Fraction(a, randint(1, 3))
    
        K = kernel(A, check=True)
        #print "kernel: A, K, A*K"
        #print shortstr(A, K, dot(A, K))
        B = dot(A, K.transpose())
        assert numpy.abs(B).sum()==0

        if K.shape[1]>1:
            break

    #while 1:
    for i in range(100):
        m, n = 3, 4
        W = rand(m, n, 2, 3)
    
        W = row_reduce(W, truncate=True)
        s = Subspace(W)
        #print "-"*79
        #print s
        #print
        assert s == s
        assert Subspace(2*W) == s
        assert Subspace(W[:-1]) != s
        ss = s.intersect(s)
        #print ss
        assert ss == s

    #print P
    #print L
    #print U

    print "OK"


if __name__ == "__main__":

    test()



