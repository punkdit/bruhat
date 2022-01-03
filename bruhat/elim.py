#!/usr/bin/env python3

"""
Gaussian elimination over a field or ring.
Modified from gelim.py
"""

import sys, os
from random import randint, seed


import numpy

from bruhat.smap import SMap
from bruhat.argv import argv




def write(s):
    sys.stdout.write(str(s)+' ')
    sys.stdout.flush()



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

    smap = SMap()
    smap[0, 0] = "["

    items = {}
    dw = 3
    for i in range(m):
      for j in range(n):
        x = A[i, j]
        s = str(x) if x is not None else "N"
        dw = max(dw, len(s)+1)
        items[i, j] = s

    for i in range(m):
      smap[i, 1] = "["
      smap[i, n*dw+2] = "]"
      smap[i, n*dw+3] = ","
      for j in range(n):
        s = items[i, j]
        s = s.rjust(dw-1)
        smap[i, j*dw+2] = s
    smap[max(0, m-1), n*dw+3] = "]"

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


def zeros(ring, m, n):
    A = numpy.empty((m, n), dtype=object)
    A[:] = ring.zero
    return A


def rand(ring, m, n, p=3, q=3):
    A = zeros(ring, m, n)
    for i in range(m):
      for j in range(n):
        a = randint(-p, p)
        a = ring.promote(a)
        if q > 1:
            a /= randint(1, q)
        A[i, j] = a
    return A


def array(items):
    return numpy.array(items, dtype=object)


def identity(ring, m):
    I = zeros(ring, m, m)
    for i in range(m):
        I[i, i] = ring.one
    return I


def eq(A, B):
    return numpy.alltrue(A==B)


def dot(ring, A, B):
    C = numpy.dot(A, B)
    if len(A.shape)==2 and A.shape[1] == 0:
        C[:] = ring.zero
    return C


def dotx(ring, *items):
    idx = 0
    A = items[idx]
    while idx+1 < len(items):
        B = items[idx+1]
        A = dot(ring, A, B)
        idx += 1 
    return A


def compose(ring, *items):
    #assert isinstance(ring, element.Ring)
    items = list(reversed(items))
    A = dotx(ring, *items)
    return A


def swap_row(A, j, k):
    row = A[j, :].copy()
    A[j, :] = A[k, :]
    A[k, :] = row

def swap_col(A, j, k):
    col = A[:, j].copy()
    A[:, j] = A[:, k]
    A[:, k] = col


def row_reduce(ring, A, truncate=False, inplace=False, check=False, verbose=False):
    """ Remove zero rows if truncate==True
    """

    zero = ring.zero
    one = ring.one

    assert len(A.shape)==2, A.shape
    m, n = A.shape
    if not inplace:
        A = A.copy()

    if m*n==0:
        if truncate and m:
            A = A[:0, :]
        return A

    if verbose:
        print("row_reduce")
        #print("%d rows, %d cols" % (m, n))

    i = 0
    j = 0
    while i < m and j < n:
        if verbose:
            print("i, j = %d, %d" % (i, j))
            print("A:")
            print(shortstrx(A))

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
                print("swap", i, i1)
            swap_row(A, i, i1)

        assert A[i, j] != zero, "FAIL: %s != %s"%(A[i,j], zero)
        for i1 in range(i+1, m):
            if A[i1, j]:
                if verbose:
                    print("add row %s to %s" % (i, i1))
                r = -A[i1, j] / A[i, j]
                A[i1, :] += r*A[i, :]
                assert A[i1, j] == zero

        i += 1
        j += 1

    if truncate:
        m = A.shape[0]-1
        #print("sum:", m, A[m, :], A[m, :].sum())
        while m>=0 and (A[m, :]!=0).sum()==0:
            m -= 1
        A = A[:m+1, :]

    if verbose:
        print()

    return A



def complement(ring, A, verbose=False):
    A = row_reduce(ring, A, truncate=True)
    zero = ring.zero
    one = ring.one
    m, n = A.shape
    idxs = []
    j = 0
    for i in range(m):
        while A[i, j] == zero:
            idxs.append(j)
            j += 1
        j += 1
    while j<n:
        idxs.append(j)
        j += 1
    B = zeros(ring, len(idxs), n)
    for i, j in enumerate(idxs):
        B[i, j] = one
    return B


def plu_reduce(ring, A, truncate=False, check=False, verbose=False):
    """
    Solve PLU = A, st. P is permutation, L is lower tri, U is upper tri.
    Remove zero rows from U if truncate=True.
    """

    zero = ring.zero
    one = ring.one

    m, n = A.shape
    P = identity(ring, m)
    L = identity(ring, m)
    U = A.copy()

    assert m*n, (m, n)

    if verbose:
        print("plu_reduce:")
        print("%d rows, %d cols" % (m, n))

    i = 0
    j = 0
    while i < m and j < n:
        if verbose:
            print("i, j = %d, %d" % (i, j))
            print("P, L, U:")
            print(shortstrx(P, L, U))

        assert i<=j
        if i and check:
            assert numpy.alltrue(U[i:,:j]==zero) # XX rm

        # first find a nonzero entry in this col
        for i1 in range(i, m):
            if U[i1, j]:
                break
        else:
            j += 1 # move to the next col
            continue # <----------- continue ------------

        if i != i1:
            if verbose:
                print("swap", i, i1)
            swap_row(U, i, i1)
            swap_col(P, i, i1)
            swap_col(L, i, i1)
            swap_row(L, i, i1)

        if check:
            A1 = dot(ring, P, dot(ring, L, U))
            assert eq(A1, A)

        r = U[i, j]
        assert r != zero
        for i1 in range(i+1, m):
            s = U[i1, j]
            if s != zero:
                t = s/r
                if verbose: 
                    print("add %s to %s" % (i, i1))
                L[i1, i] = t
                U[i1, :] -= t*U[i, :]
                assert U[i1, j] == zero

        if check:
            A1 = dot(ring, P, dot(ring, L, U))
            assert eq(A1, A)

        i += 1
        j += 1

    if check:
      for i in range(m):
        for j in range(i):
            assert U[i, j] == zero
      for i in range(m):
        for j in range(i+1, m):
            assert L[i, j] == zero

    if truncate:
        m = U.shape[0]-1
        #print("sum:", m, U[m, :], U[m, :].sum())
        while m>=0 and (U[m, :]!=0).sum()==0:
            m -= 1
        U = U[:m+1, :]

    return P, L, U


def u_inverse(ring, U, check=False, verbose=False):
    """invert a row reduced U
    """

    zero = ring.zero
    one = ring.one

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

    U1 = zeros(ring, n, m)

    # Work backwards
    i = len(leading)-1 # <= m
    while i>=0:

        j = leading[i]
        #print("i=%d, j=%d"%(i, j))
        r = 1 / U[i, j]
        U1[j, i] = r

        #print("U, U1, U*U1:")
        #print(shortstrx(U, U1, dot(ring, U, U1)))

        k = i-1
        while k>=0:
            #print("dot")
            #print(shortstr(U[k,:]))
            #print(shortstr(U1[:,i]))
            r = dot(ring, U[k, :], U1[:, i])
            #print("=", r)
            if r != 0:
                j = leading[k]
                s = U[k, j]
                #print("set", j, i)
                U1[j, i] = -r/s
            #print(shortstr(U1[:,i]))
            assert dot(ring, U[k, :], U1[:, i]) == 0
            k -= 1
        i -= 1

    return U1


def l_inverse(ring, L, check=False, verbose=False):
    """invert L (lower triangular, 1 on diagonal)
    """

    zero = ring.zero
    one = ring.one

    m, n = L.shape
    assert m==n
    L1 = identity(ring, m)

    # Work forwards
    for i in range(m):
        #u = L1[:, i]
        assert L[i, i] == 1
        for j in range(i+1, m):
            #print("i=%d, j=%d"%(i, j))
            #print("L, L1, L*L1:")
            #print(shortstrx(L, L1, dot(ring, L, L1)))
            r = dot(ring, L[j, :], L1[:, i])
            #print("r =", r)
            if r != 0:
                assert L1[j, i] == 0
                L1[j, i] = -r
            r = dot(ring, L[j, :], L1[:, i])
            #print("r =", r)
            #print(shortstrx(L, L1, dot(ring, L, L1)))
            assert dot(ring, L[j, :], L1[:, i]) == 0

    assert eq(dot(ring, L, L1), identity(ring, m))
    return L1


def pseudo_inverse(ring, A, check=False):
    m, n = A.shape
    if m*n == 0:
        A1 = zeros(ring, n, m)
        return A1
    P, L, U = plu_reduce(ring, A, verbose=False, check=check)
    L1 = l_inverse(ring, L, check=check)
    U1 = u_inverse(ring, U, check=check)
    #print("P, L, U, PLU:")
    #print(shortstr(P, L, U, dot(ring, dot(ring, P, L), U)))
    
    A1 = dot(ring, U1, dot(ring, L1, P.transpose()))
    #print(shortstr(dot(ring, A1, A)))
    return A1


def solve(ring, H, u, force=False, verbose=False, check=False):
    "Solve Hv = u"
    assert len(H) == len(u)
    A = pseudo_inverse(ring, H, check)
    #print("pseudo_inverse")
    #print(shortstr(A))
    v = dot(ring, A, u)
    if eq(dot(ring, H, v), u) or force:
        return v


# https://en.wikipedia.org/wiki/Kernel_%28linear_algebra%29#Computation_by_Gaussian_elimination
def kernel(ring, A, check=False, verbose=False):
    """
        return largest K such that dot(A, K) == 0.
    """

    if verbose:
        print("kernel")

    zero = ring.zero
    one = ring.one

    m, n = A.shape
    A, A0 = A.copy(), A

    K = identity(ring, n)

    # Column reduce A, while duplicating operations onto K

    i = 0 # row
    for j in range(n): # col

        if verbose:
            print("A, K (i=%d, j=%d)" % (i, j))
            print(shortstr(A))
            print("----------")
            print(shortstr(K))
            print()

        # look for a row
        while i<m and (A[i, j:]!=zero).sum()==0:
            i += 1

        if i==m:
            break

        if A[i, j] == zero:
            k = j
            while A[i, k]==zero:
                k += 1
            swap_col(A, j, k)
            swap_col(K, j, k)

        for k in range(j+1, n):
            r = -ring.promote(A[i, k]) / A[i, j]
            A[:, k] += r * A[:, j]
            K[:, k] += r * K[:, j]

        i += 1
        if i==m:
            break

    if verbose:
        print("A, K (i=%d, j=%d)" % (i, j))
        print(shortstr(A))
        print("----------")
        print(shortstr(K))
        print()

    j = K.shape[1] - 1
    while j>=0 and (A[:, j]!=zero).sum() == 0:
        j -= 1
    j += 1

    #K = K[:, j+1:]
    K = K[:, j:]
    if check:
        B = dot(ring, A0, K)
        assert numpy.alltrue(B==zero)

    return K


def projector(ring, A, check=False):

    """
        Find universal projector that kills the columns of A,
        ie. PP=P and PA = 0, st. given any other Q with
        QQ=Q and QA=0, then there exists R st. Q=RP.
    """

    """
        Alternatively
        Find universal projector that preserves the columns of A,
        ie. PP=P and PA=A, st. given any other Q with
        QQ=Q and QA=A, then there exists R st. P=RQ.
    """

    m, n = A.shape

    P = identity(ring, m) - dot(ring, A, pseudo_inverse(ring, A))

    return P


def pushout(ring, J, K, J1=None, K1=None, check=False):
    """  
    Return JJ,KK given J and K in the following diagram:

       J
    A ---> B
    |      |
    | K    | JJ
    v      v
    C ---> B+C/~
       KK

    if also given J1:B->T and K1:C->T (st. J1*J=K1*K)
    return unique arrow F : B+C/~ --> T (st. F*JJ=J1 and F*KK=K1).
    """
    assert J.shape[1] == K.shape[1]

    b, c = J.shape[0], K.shape[0]
    JJ = zeros(ring, b+c, b)
    JJ[:b] = identity(ring, b)

    KK = zeros(ring, b+c, c)
    KK[b:] = identity(ring, c)

    #print(K.shape)
    #print(KK.shape)
    kern = compose(ring, J, JJ) - compose(ring, K, KK)
    # We need to kill the columns of kern
    R = projector(ring, kern, check=check)
    R = row_reduce(ring, R, truncate=True, check=check)

    assert eq(compose(ring, J, JJ, R), compose(ring, K, KK, R))

    JJ = compose(ring, JJ, R)
    KK = compose(ring, KK, R)

    if J1 is not None:
        assert K1 is not None
        assert J1.shape[0] == K1.shape[0]
        assert eq(compose(ring, J, J1), compose(ring, K, K1))
        m = J1.shape[0]
        n = R.shape[0]
        F = zeros(ring, m, n)

#        print "F=", F.shape
#        print "R=", R.shape
        #print shortstr(R)

        Rinv = pseudo_inverse(ring, R, check=check)

#        print "Rinv=", Rinv.shape
        #print shortstr(Rinv)

        for i in range(n):
            r = Rinv[:, i]
            u, v = r[:b], r[b:]
            u = dot(ring, J1, u)
            v = dot(ring, K1, v)
            #assert eq(u, v)
            uv = (u+v)
            F[:, i] = uv

        assert eq(compose(ring, JJ, F), J1)
        assert eq(compose(ring, KK, F), K1)

        return JJ, KK, F

    else:
        return JJ, KK


def pullback(ring, J, K, J1=None, K1=None, check=False):
    assert J.shape[0] == K.shape[0]
    Jt = J.transpose()
    Kt = K.transpose()
    if J1 is not None:
        assert K1 is not None
        assert J1.shape[1] == K1.shape[1]
        J1t = J1.transpose()
        K1t = K1.transpose()
        JJt, KKt, Ft = pushout(ring, Jt, Kt, J1t, K1t, check)
        return JJt.transpose(), KKt.transpose(), Ft.transpose()

    JJt, KKt = pushout(ring, Jt, Kt, None, None, check)
    return JJt.transpose(), KKt.transpose()


def old_cokernel(ring, J, check=False):
    """  
    find f as a pushout of the following diagram:

      J
    k ---> n
    |      |
    |      | f
    v      v
    0 ---> n'

    """

    n, k = J.shape
    K = zeros(ring, 0, k)

    f, g = pushout(ring, J, K, check=check)
    return f


def cokernel(ring, J, P1=None, check=False):
    """
           J
        k ---> n
        |      |
        |      | P
        v      v
        0 ---> m
    """

    n, k = J.shape
    L = zeros(ring, 0, k)
    Q1 = None
    if P1 is not None:
        Q1 = zeros(ring, P1.shape[0], 0)
        assert eq(compose(ring, J, P1), 0)

    if P1 is None:
        f, g = pushout(ring, J, L, check=check)
        return f

    else:
        P, Q, R = pushout(ring, J, L, P1, Q1, check=check)
        assert Q.shape[1] == 0
        return P, R


def coequalizer(ring, J, K, JK1=None):
    JK = J-K
    result = cokernel(ring, JK, JK1)
    return result


def rank(ring, A, **kw):
    A = row_reduce(ring, A, truncate=True, **kw)
    return len(A)


def nullity(ring, A, **kw):
    K = kernel(ring, A, **kw)
    return len(K[0])


class Subspace(object):
    """ Subspace represented as the rowspace of a matrix.
    """
    def __init__(self, ring, W):
        assert rank(ring, W) == len(W)
        self.ring = ring
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
        if solve(self.ring, W1, W2) is None:
            return False
        if solve(self.ring, W2, W1) is None:
            return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def intersect(self, other):
        W1 = self.W
        W2 = other.W
        W = numpy.concatenate((W1, W2))
        #print("intersect")
        #print(shortstr(self.W))
        #print(shortstr(other.W))
        #print(shortstr(W))
        #print()
        K = kernel(self.ring, W.transpose()).transpose()
        #print("K:")
        #print(shortstr(K))
        W = dot(self.ring, K[:, :len(W1)], W1)
        return Subspace(self.ring, W)



def test():

    _seed = argv.get("seed")
    if _seed is not None:
        seed(_seed)

    if argv.fast:
        from bruhat import _element as element
    else:
        from bruhat import element

    ring = element.Q # field of rationals
    zero = ring.zero
    one = ring.one

    m, n = 4, 5
    A = zeros(ring, m, n)
    U = zeros(ring, m, n)

    for i in range(m):
      for j in range(i, n):
        a = randint(-3, 3)
        if i==j and a==0:
            a = 1
        U[i, j] = ring.promote(a) / randint(1, 3)

    V = u_inverse(ring, U)

    L = identity(ring, m)
    for i in range(m):
        for j in range(i):
            L[i, j] = ring.promote(randint(-3, 3)) / randint(1, 3)
    V = l_inverse(ring, L)

    A = rand(ring, m, n)
    P, L, U = plu_reduce(ring, A, check=True, verbose=False)

    for i in range(100):
        A = rand(ring, m, n)
        if rank(ring, A, check=True) + nullity(ring, A, check=True) == n:
            #write('.')
            continue

        print("FAIL")
        print("A:")
        print(shortstr(A))
        print("rank=%s, nullity=%s, n=%d"%(
            rank(ring, A, check=True) ,
            nullity(ring, A, check=True) , n))
        print("row_reduce:")
        B = row_reduce(ring, A, verbose=True)
        print(shortstr(B))
        print("kernel:")
        print(shortstr(kernel(ring, A)))
        return

    A = zeros(ring, 3, 4)
    A[0, 0] = 1
    A[1, 1] = 1
    A[2, 2] = 1
    A = numpy.concatenate((A, A))
    A = A.transpose()
    #print("A:")
    #print(shortstr(A))
    K = kernel(ring, A)
    #print("K:")
    #print(shortstr(K))
    assert len(K[0]) == 3, len(K)

    while 1:
        m, n = 3, 4
        A = zeros(ring, m, n)
        for i in range(m):
          for j in range(n):
            a = randint(-2, 2)
            A[i, j] = ring.promote(a) / randint(1, 3)
    
        K = kernel(ring, A, check=True)
        #print("kernel: A, K, A*K")
        #print(shortstr(A, K, dot(ring, A, K)))
        B = dot(ring, A, K)
        assert numpy.alltrue(B==zero)

        if K.shape[1]>1:
            break

    #while 1:
    for i in range(100):
        m, n = 3, 4
        W = rand(ring, m, n, 2, 3)
    
        W = row_reduce(ring, W, truncate=True)
        s = Subspace(ring, W)
        #print("-"*79)
        #print(s)
        #print()
        assert s == s
        assert Subspace(ring, 2*W) == s
        assert Subspace(ring, W[:-1]) != s
        ss = s.intersect(s)
        #print(ss)
        assert ss == s

    #print(P)
    #print(L)
    #print(U)


    for trial in range(10):

        m, n = 5, 10
    
        A = rand(ring, m, n, 0, 1)
        #print(shortstr(A))
    
        B = row_reduce(ring, A)
        #print(shortstr(B))
     
        B = complement(ring, A)
        #print(shortstr(B))
    
        C = numpy.concatenate((A, B)) 
        assert rank(ring, C) == n

        m, n = 10, 5
        f = rand(ring, m, n, 1, 1)
        #print(shortstr(f))
        g = cokernel(ring, f)
        #print(shortstr(g))
        assert rank(ring, g)==m-n
        gf = dot(ring, g, f)
        #print(shortstr(gf))
        assert numpy.alltrue(gf==ring.zero)


    print("OK")


if __name__ == "__main__":

    test()



