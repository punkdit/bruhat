#!/usr/bin/env python

"""
Gaussian elimination over the rationals. Sparse version.
"""

from __future__ import print_function

import sys, os
from random import randint, seed
from fractions import Fraction

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


class Vec(object):
    def __init__(self, n, data=None):
        self.n = n
        if data is None:
            data = {}
        self.data = data

    def __str__(self):
        return "Vec(%d, %s)"%(self.n, self.data)
    __repr__ = __str__

    def copy(self):
        return Vec(self.n, dict(self.data))

    def dot(A, B):
        assert A.n==B.n
        r = 0
        for idx, value in A.data.items():
            r += value*B.data.get(idx, 0)
        return r

    def __len__(self):
        return self.n

    def __getitem__(self, i):
        return self.data.get(i, 0)

    def __setitem__(self, i, v):
        self.data[i] = v

    def __eq__(self, other):
        assert self.n == other.n
        fail

    def __ne__(self, other):
        assert self.n == other.n
        fail

    def __rmul__(self, r):
        data = dict((i, r*value) for (i, value) in self.data.items())
        return Vec(self.n, data)

    def __add__(A, B):
        A = A.copy()
        for (i, v) in B.data.items():
            A.data[i] = A.data.get(i, 0) + v
        return A

    def __sub__(A, B):
        A = A.copy()
        for (i, v) in B.data.items():
            A.data[i] = A.data.get(i, 0) - v
        return A

    def is_zero(A):
        for v in A.data.values():
            if v != 0:
                return False
        return True

assert Fraction(0, 1) == 0

def mkslice(slc, n):
    assert slc.step == None, "not implemented"
    stop = slc.stop
    if stop is None:
        stop = n
    elif stop < 0:
        stop = stop%n
    start = slc.start
    if start is None:
        start = 0
    elif start < 0:
        start = start%n
    slc = slice(start, stop)
    return slc


class Sparse(object):
    def __init__(self, m, n, data=None):
        self.shape = m, n # rows, cols
        if data is None:
            data = {}
        self.data = data
        self.rows = [[] for i in range(m)] # nonzero col idx for each row
        self.cols = [[] for i in range(n)] # nonzero row idx for each col
        self.check()

    def __str__(self):
        return "Sparse(%d, %d, %s)"%(self.shape[0], self.shape[1], self.data)
    __repr__ = __str__

    def copy(self):
        A = Sparse(*self.shape)
        A.data.update(self.data)
        A.rows = [list(idxs) for idxs in self.rows]
        A.cols = [list(idxs) for idxs in self.cols]
        A.check() # XX
        return A

    def check(self):
        for key, value in self.data.items():
            row, col = key
            assert col in self.rows[row]
            assert row in self.cols[col]
        data = self.data
        for row, cols in enumerate(self.rows):
            for col in cols:
                assert data.get((row, col), 0) != 0
        for col, rows in enumerate(self.cols):
            for row in rows:
                assert data.get((row, col), 0) != 0

    def __len__(self):
        return self.shape[0]

    def __setitem__(self, idx, value):
        row, col = idx
        if isinstance(row, slice) and isinstance(col, slice):
            assert 0, "not implemented"
        if isinstance(row, slice):
            assert row == slice(None)
            assert isinstance(value, Vec)
            for i in list(self.cols[col]):
                self[i, col] = 0 # recurse
            assert not self.cols[col]
            for i, x in value.data.items():
                self[i, col] = x # recurse
            return
        if isinstance(col, slice):
            assert col == slice(None)
            assert isinstance(value, Vec)
            #print("__setitem__", row, col, value)
            #print(shortstr(self))
            #print("rows:", self.rows[row])
            for j in list(self.rows[row]):
                #print("clear", row, j)
                self[row, j] = 0 # recurse
            #print(" =>")
            #print(shortstr(self))
            assert not self.rows[row], self.rows[row]
            for j, x in value.data.items():
                self[row, j] = x # recurse
            return
        n, m = self.shape
        assert 0<=row<n
        assert 0<=col<m
        if value != 0:
            self.data[idx] = value
            if col not in self.rows[row]:
                self.rows[row].append(col)
            if row not in self.cols[col]:
                self.cols[col].append(row)
        elif idx in self.data:
            del self.data[idx]
            if col in self.rows[row]:
                self.rows[row].remove(col)
            if row in self.cols[col]:
                self.cols[col].remove(row)
        #self.check()

    def __getitem__(self, idx):
        if type(idx) != tuple:
            idx = (idx, slice(None))
        row, col = idx
        if isinstance(row, slice) and isinstance(col, slice):
            row = mkslice(row, self.shape[0])
            col = mkslice(col, self.shape[1])
            A = Sparse(row.stop - row.start, col.stop-col.start)
            for key, value in self.data.items():
                i, j = key
                if row.start <= i < row.stop and col.start <= j < col.stop:
                    A[i-row.start, j-col.start] = value
            return A
        elif isinstance(row, slice):
            # return col A[:, col]
            row = mkslice(row, self.shape[0])
            idxs = self.cols[col]
            data = dict((i, self.data[i, col]) for i in idxs if row.start<=i<row.stop)
            A = Vec(self.shape[0], data)
        elif isinstance(col, slice):
            # return row A[row, :]
            col = mkslice(col, self.shape[1])
            idxs = self.rows[row]
            data = dict((j, self.data[row, j]) for j in idxs if col.start<=j<col.stop)
            A = Vec(self.shape[1], data)
        else:
            A = self.data.get(idx, 0)
        return A

    @classmethod
    def identity(cls, m):
        I = cls(m, m)
        for i in range(m):
            I[i, i] = 1
        return I

    def __eq__(self, other):
        assert self.shape==other.shape
        return self.data==other.data

    def __ne__(self, other):
        assert self.shape==other.shape
        return self.data!=other.data

    def __add__(self, other):
        A = self.copy()
        for key, value in other.data.items():
            A[key] = A[key] + value
        return A

    def __sub__(self, other):
        A = self.copy()
        for key, value in other.data.items():
            A[key] = A[key] - value
        return A

    def __rmul__(self, r):
        A = Sparse(*self.shape)
        for key, value in self.data.items():
            A[key] = r*value
        return A

    def __neg__(self):
        return (-1)*self

    def dot(self, other):
        #data = {}
        shape = self.shape[0], other.shape[1]
        A = Sparse(*shape)
        for key, value in self.data.items():
            i, j = key
            for k in other.rows[j]:
                A[i, k] = A[i, k] + value * other.data[j, k]
        return A
    __mul__ = dot

    def is_zero(self):
        return not self.data

    def get_rows(self, col):
        rows = self.cols[col]
        rows.sort()
        return list(rows)

    def get_cols(self, row):
        cols = self.rows[row]
        cols.sort()
        return list(cols)

    def concatenate(A, B):
        assert A.shape[1]==B.shape[1]
        C = Sparse(A.shape[0]+B.shape[0], A.shape[1])
        for idx, value in A.data.items():
            C[idx] = value
        offset = A.shape[0]
        for idx, value in B.data.items():
            row, col = idx
            C[row+offset, col] = value
        return C

    def transpose(A):
        B = Sparse(A.shape[1], A.shape[0])
        for idx, value in A.data.items():
            i, j = idx
            B[j, i] = value
        return B


def array(items):
    assert items
    m = len(items)
    n = len(items[0])
    A = Sparse(m, n)
    for i, row in enumerate(items):
        for j, value in enumerate(row):
            A[i, j] = value
    return A


def test_sparse():

    m, n = 5, 5
    A = Sparse(m, n)
    A[1, 2] = 3

    B = Sparse.identity(m)

    assert A==A
    assert A!=B
    assert (A*B) == A

    I = Sparse.identity(2)
    A = array([[1, 1], [1, 0]])
    B = array([[0, 1], [1, 1]])
    C = array([[1, 2], [0, 1]])
    assert A*B == C

    Z = array([[1, 0], [0, -1]])
    X = array([[0, 1], [1, 0]])

    assert (Z*Z)==I
    assert (A-B)==Z
    assert B+Z==A


test_sparse()

#zeros = lambda n, m : Sparse((n, m))
zeros = Sparse
eq = lambda A, B : A==B
identity = Sparse.identity


def dotx(*items):
    idx = 0
    A = items[idx]
    while idx+1 < len(items):
        B = items[idx+1]
        A = A.dot(B)
        idx += 1 
    return A

dot = lambda A, B : A.dot(B)


def compose(*items):
    items = list(reversed(items))
    A = dotx(*items)
    return A



#def shortstr(A):
#    return str(A)


#def shortstrx(*args):
#    return '\n'.join(str(A) for A in args)



def swap_row(A, j, k):
    #print("swap_row", j, k)
    #print(shortstr(A))
    #print(" =>")
    row = A[j, :].copy()
    A[j, :] = A[k, :]
    A[k, :] = row
    #print(shortstr(A))
    #print()

def swap_col(A, j, k):
    #print("swap_col", j, k)
    #print(shortstr(A))
    #print(" =>")
    col = A[:, j].copy()
    A[:, j] = A[:, k]
    A[:, k] = col
    #print(shortstr(A))
    #print()


def row_reduce(A, truncate=False, inplace=False, check=False, verbose=False):
    """ Remove zero rows if truncate==True
    """

    assert len(A.shape)==2, A.shape
    m, n = A.shape
    if not inplace:
        A = A.copy()

    if m*n==0:
        if truncate and m:
            A = A[:0, :]
        return A

    if verbose:
        print( "row_reduce")
        #print( "%d rows, %d cols" % (m, n))

    i = 0
    j = 0
    while i < m and j < n: 
        if verbose:
            print( "i, j = %d, %d" % (i, j))
            print( "A:")
            print( shortstrx(A))

        assert i<=j 
        if i and check:
#            assert (A[i:,:j]!=0).sum() == 0
            assert A[i:,:j].is_zero()

        # first find a nonzero entry in this col
        for i1 in range(i, m):
            if A[i1, j]:
                break
        else:
            j += 1 # move to the next col
            continue # <----------- continue ------------

        if i != i1:
            if verbose:
                print( "swap", i, i1)
            swap_row(A, i, i1)

        assert A[i, j]
        for i1 in range(i+1, m):
            if A[i1, j]:
                if verbose: 
                    print( "add row %s to %s" % (i, i1))
                r = -Fraction(A[i1, j], A[i, j])
                A[i1, :] += r*A[i, :]
                assert A[i1, j] == 0

        i += 1 
        j += 1 

    if truncate:
        m = A.shape[0]-1
        #print( "sum:", m, A[m, :], A[m, :].sum())
#        while m>=0 and (A[m, :]!=0).sum()==0:
        while m>=0 and A[m, :].is_zero():
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
        print( "plu_reduce:")
        print( "%d rows, %d cols" % (m, n))

    i = 0
    j = 0
    while i < m and j < n:
        if verbose:
            print( "i, j = %d, %d" % (i, j))
            print( "P, L, U:")
            print( shortstrx(P, L, U))

        assert i<=j
        if i and check:
#            assert U[i:,:j].max() == 0 # XX rm
            assert U[i:,:j].is_zero()

        # first find a nonzero entry in this col
        for i1 in range(i, m):
            if U[i1, j]:
                break
        else:
            j += 1 # move to the next col
            continue # <----------- continue ------------

        if i != i1:
            if verbose:
                print( "swap", i, i1)
            swap_row(U, i, i1)
            swap_col(P, i, i1)
            swap_col(L, i, i1)
            swap_row(L, i, i1)

        if check:
            A1 = dot(P, dot(L, U))
            if not eq(A1, A):
                print("FAIL:")
                print("A")
                print(shortstr(A))
                print("A1")
                print(shortstr(A1))
                assert 0

        r = U[i, j]
        assert r != 0
        for i1 in range(i+1, m):
            s = U[i1, j]
            if s != 0:
                t = Fraction(s, r)
                if verbose: 
                    print( "add %s to %s" % (i, i1))
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
        #print( "sum:", m, U[m, :], U[m, :].sum())
#        while m>=0 and U[m, :].sum()==0:
        while m>=0 and U[m, :].is_zero():
            m -= 1
        U = U[:m+1, :]

    return P, L, U


def u_inverse(U, check=False, verbose=False):
    """invert a row reduced U
    """

    m, n = U.shape

    if verbose:
        print("u_inverse")
        print(shortstr(U))

    #items = []
    leading = []
    for row in range(m):
        #cols = numpy.where(U[row, :])[0]
        cols = U.get_cols(row)
        #print("row %d, cols %s"%(row, cols))
        if not len(cols):
            break
        col = cols[0]
        assert U[row, col]
        leading.append(col)

    #print("leading:", leading)
    assert sorted(leading) == leading
    assert len(set(leading)) == len(leading)

    U1 = zeros(n, m)

    #print( shortstr(U))

    # Work backwards
    i = len(leading)-1 # <= m
    while i>=0:

        j = leading[i]
        #print( "i=%d, j=%d"%(i, j))
        r = Fraction(1, U[i, j])
        U1[j, i] = r

        #print( "U, U1, U*U1:")
        #print( shortstrx(U, U1, dot(U, U1)))

        k = i-1
        while k>=0:
            #print( "dot")
            #print( shortstr(U[k,:]))
            #print( shortstr(U1[:,i]))
            r = dot(U[k, :], U1[:, i])
            #print( "=", r)
            if r != 0:
                j = leading[k]
                s = U[k, j]
                #print( "set", j, i)
                U1[j, i] = -Fraction(r, s)
            #print( shortstr(U1[:,i]))
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
            #print( "i=%d, j=%d"%(i, j))
            #print( "L, L1, L*L1:")
            #print( shortstrx(L, L1, dot(L, L1)))
            r = dot(L[j, :], L1[:, i])
            #print( "r =", r)
            if r != 0:
                assert L1[j, i] == 0
                L1[j, i] = -r
            r = dot(L[j, :], L1[:, i])
            #print( "r =", r)
            #print( shortstrx(L, L1, dot(L, L1)))
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
    #print( "P, L, U, PLU:")
    #print( shortstr(P, L, U, dot(dot(P, L), U)))
    
    A1 = dot(U1, dot(L1, P.transpose()))
    #print( shortstr(dot(A1, A)))
    return A1


def solve(H, u, force=False, verbose=False, check=False):
    "Solve Hv = u"
    assert len(H) == len(u)
    A = pseudo_inverse(H, check)
    #print( "pseudo_inverse")
    #print( shortstr(A))
    v = dot(A, u)
    if eq(dot(H, v), u) or force:
        return v


# https://en.wikipedia.org/wiki/Kernel_%28linear_algebra%29#Computation_by_Gaussian_elimination
def kernel(A, check=False, verbose=False):
    """
        return largest K such that dot(A, K) == 0.
    """

    if verbose:
        print( "kernel")

    m, n = A.shape
    A, A0 = A.copy(), A

    K = identity(n)

    # Column reduce A, while duplicating operations onto K

    i = 0 # row
    for j in range(n): # col
        #if j%10==0:
        #    write("(%d/%d)"%(j, n))

        if verbose:
            print( "A, K (i=%d, j=%d)" % (i, j))
            print( shortstr(A))
            print( "----------")
            print( shortstr(K))
            print()

        # look for a row
        #while i<m and A[i, j:].is_zero():
        while i<m and max(A.get_cols(i) or [-1]) < j:
            i += 1

        if i==m:
            break

        if A[i, j] == 0:
            k = j
            #while A[i, k]==0:
            #    k += 1
            for k in A.get_cols(i):
                if k >= j:
                    break
            else:
                assert 0
            swap_col(A, j, k)
            swap_col(K, j, k)

        #for k in range(j+1, n):
        #    if A[i, k]==0:
        #        continue
        for k in A.get_cols(i):
            if k<=j:
                continue
            r = -Fraction(A[i, k], A[i, j])
            A[:, k] += r * A[:, j]
            K[:, k] += r * K[:, j]

        i += 1
        if i==m:
            break

    if verbose:
        print( "A, K (i=%d, j=%d)" % (i, j))
        print( shortstr(A))
        print( "----------")
        print( shortstr(K))
        print()

    j = K.shape[1] - 1
    while j>=0 and A[:, j].is_zero():
        j -= 1
    j += 1

    #K = K[:, j+1:]
    K = K[:, j:]
    if check:
        B = dot(A0, K)
        #assert numpy.abs(B).sum()==0
        assert B.is_zero(), repr(B)

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
    def __init__(self, W, check=False):
        if check:
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
        #W = numpy.concatenate((W1, W2))
        W = W1.concatenate(W2)
        #print( "intersect")
        #print( shortstr(W))
        #print()
        K = kernel(W.transpose())#.transpose()
        #print( "K:")
        #print( shortstr(K))
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

    for trial in range(100):
        A = rand(m, n)
        row_reduce(A, check=True)

    for trial in range(100):
        for i in range(m):
          for j in range(i, n):
            a = randint(-3, 3)
            if i==j and a==0:
                a = 1
            U[i, j] = Fraction(a, randint(1, 3))
    
        V = u_inverse(U, check=True)

        
    for trial in range(100):
        L = identity(m)
        for i in range(m):
            for j in range(i):
                L[i, j] = Fraction(randint(-3, 3), randint(1, 3))
        V = l_inverse(L, check=True)
    
    for trial in range(100):
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

        print( "FAIL")
        print( "A:")
        print( shortstr(A))
        print( "row_reduce:")
        B = row_reduce(A, verbose=True)
        print( shortstr(B))
        print( "kernel:")
        print( shortstr(kernel(A)))
        return


    A = zeros(3, 4)
    A[0, 0] = 1
    A[1, 1] = 1
    A[2, 2] = 1
    #A = numpy.concatenate((A, A))
    A = A.concatenate(A)
    A = A.transpose()
    #print( "A:")
    #print( shortstr(A))
    K = kernel(A)
    #print( "K:")
    #print( shortstr(K))
    assert len(K) == 3

    while 1:
        m, n = 3, 4
        A = zeros(m, n)
        for i in range(m):
          for j in range(n):
            a = randint(-2, 2)
            A[i, j] = Fraction(a, randint(1, 3))
    
        K = kernel(A, check=True)
        #print( "kernel: A, K, A*K")
        #print( shortstr(A, K, dot(A, K)))
        B = dot(A, K.transpose())
        #assert numpy.abs(B).sum()==0
        assert B.is_zero()

        if K.shape[1]>1:
            break

    #while 1:
    for i in range(100):
        m, n = 3, 4
        W = rand(m, n, 2, 3)
    
        W = row_reduce(W, truncate=True)
        s = Subspace(W)
        #print( "-"*79)
        #print( s)
        #print()
        assert s == s
        assert Subspace(2*W) == s
        assert Subspace(W[:-1]) != s
        ss = s.intersect(s)
        #print( ss)
        assert ss == s

    #print( P)
    #print( L)
    #print( U)

    print( "OK")


if __name__ == "__main__":

    test()



