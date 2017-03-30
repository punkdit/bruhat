#!/usr/bin/env python

import sys
from random import random, randint, shuffle, seed

import numpy
import numpy.random as ra
from numpy import dot, concatenate


def array2(items):
    return numpy.array(items, dtype=numpy.int32)


def zeros2(*shape):
    if len(shape)==1 and type(shape[0]) is tuple:
        shape = shape[0]
    return numpy.zeros(shape, dtype=numpy.int32)


def identity2(n):
    return numpy.identity(n, dtype=numpy.int32)


def dot2(*items):
    idx = 0
    A = items[idx]
    while idx+1 < len(items):
        B = items[idx+1]
        A = dot(A, B)
        idx += 1
    A = A%2
    return A


def compose2(*items):
    items = list(reversed(items))
    A = dot2(*items)
    return A


def eq2(A, B):
    A = (A-B)%2
    return numpy.abs(A).sum() == 0


def pop2(A, i):
    m, n = A.shape
    B = zeros2(m-1, n)
    B[:i] = A[:i]
    B[i:] = A[i+1:]
    return B


def insert2(A, i, u):
    m, n = A.shape
    B = zeros2(m+1, n)
    B[:i] = A[:i]
    B[i] = u
    B[i+1:] = A[i:]
    return B


def append2(A, u):
    A = insert2(A, A.shape[0], u)
    return A


def binomial2(p, *shape):
    error = ra.binomial(1, p, shape)
    error = error.astype(numpy.int32)
    return error



def rand2(m, n, p=0.5, weight=None):
    if weight is not None:
        r = zeros2(m, n)
        for i in range(m):
            idxs = range(n)
            for j in range(weight):
                idx = idxs.pop(randint(0, len(idxs)-1))
                assert r[i, idx] == 0
                r[i, idx] = 1
    else:
        r = binomial2(p, m, n)
    return r


def randexpo(n, C=10.):
    "use exponential distribution to choose 0,...,n-1"
    assert n>0, n
    x = n
    while x>n-1:
        x = ra.exponential(n/float(C))
    i = int(x)
    assert 0<=i<n
    return i


def parse(s):
    s = s.replace('.', '0')
    lines = s.split()
    lines = [l.strip() for l in lines if l.strip()]
    rows = [list(int(c) for c in l) for l in lines]
    if rows:
        n = len(rows[0])
        for row in rows:
            assert len(row)==n, "rows have varying lengths"
    a = numpy.array(rows, dtype=numpy.int32)
    return a


def shortstr(A, deco=False, zero='.'):
    A = array2(A)
    A = A.view()
    assert len(A.shape)<=2
    if 1 in A.shape:
        A.shape = (1, A.shape[0]*A.shape[1])
    if len(A.shape)==1:
        A.shape = (1, A.shape[0])
    m, n = A.shape
    rows = []
    for row in range(m):
        idx = row
        row = list(A[row, :])
        row = ''.join(str(int(i)) for i in row)
        row = row.replace('0', zero)
        if deco:
            row = "%3d "%idx + row
        rows.append(row)
    if deco:
        row = ''.join(str(i%10) for i in range(n))
        rows.append('    '+row)
    return '\n'.join(rows)


def shortstr2(H, u):
    m, n = H.shape
    rows = []
    for idx in range(m):
        row = list(H[idx, :])
        row = ''.join(str(i) for i in row)
        row += ' ' + str(u[idx])
        row = row.replace('0', '.')
        rows.append("%3d "%idx + row)
    return '\n'.join(rows)


def hbox2(a, b, space=''):
    a = a.split('\n')
    b = b.split('\n')
    if len(a)<len(b):
        #a += ['']*(len(b)-len(a))
        a = ['']*(len(b)-len(a)) + a
    else:
        #b += ['']*(len(a)-len(b))
        b = ['']*(len(a)-len(b)) + b
    assert len(a)==len(b), (len(a), len(b), repr(a), repr(b))
    c = []
    w = max(len(aa) for aa in a)
    c = [a[i].ljust(w) + space + b[i] for i in range(len(a))]
    return '\n'.join(c)


def hbox(items, space=''):
    a = items[0]
    for i in range(1, len(items)):
        a = hbox2(a, items[i], space)
    return a



def enum2(n):
    "enumerate through each vector of dimension n"
    #assert n < 20, "too big"
    N = 2**n
    i = 0
    rn = range(n)
    while i < N:
        a = [(i>>j)%2 for j in range(n)]
        a = array2(a)
        yield a
        i += 1


def swap_row(A, j, k):
    row = A[j, :].copy()
    A[j, :] = A[k, :]
    A[k, :] = row

def swap_col(A, j, k):
    col = A[:, j].copy()
    A[:, j] = A[:, k]
    A[:, k] = col


def shortstrx(*As, **kw):
    from smap import SMap
    smap = SMap()

    zero = kw.get('zero', '.')

    col = 0
    for A in As:
        if A is None:
            smap[0, col] = 'None'
            col += 5
            continue
        if len(A.shape)==1:
            A = A.view()
            A.shape = A.shape[0], 1
        m, n = A.shape
        for i in range(m):
          for j in range(n):
            smap[i, col + j] = zero if A[i, j]==0 else str(A[i, j])
        col += n + 1

    return str(smap)


def lu_decompose(A, verbose=False):
    """
        Solve LU = A, st. L is lower, U is upper.
        Reorder A rows *inplace* as needed.
    """

    m, n = A.shape
    if verbose:
        print "lu_decompose: shape=", m, n
    assert m<=n, "umm.."

    for k in range(m-1):
        if A[k, k]:
            continue
        for j in range(k+1, m):
            if A[j, k]:
                swap_row(A, k, j)
                break

    L = identity2(m)
    U = A.copy()

    for k in range(m-1):
        
        if verbose:
            print "_"*79
            print "k =", k
            print shortstrx(L, U, A)
            print "LU == A", eq2(dot(L, U), A)
            print

        if U[k, k] == 0:
            continue

        for j in range(k+1, m):
            L[j, k] = U[j, k]
            U[j, k:] = (U[j, k:] - L[j, k] * U[k, k:]) % 2

    if verbose:
        print shortstrx(L, U, A)
        print "LU == A", eq2(dot(L, U), A)
        print

    return L, U
        

def row_reduce(H, truncate=True, inplace=False, check=False, debug=False):
    """Remove zero rows if truncate=True
    """

    assert len(H.shape)==2, H.shape
    m, n = H.shape
    if not inplace:
        H = H.copy()

    if m*n==0:
        return H

    if debug:
        print "solve:"
        print "%d rows, %d cols" % (m, n)

    i = 0
    j = 0
    while i < m and j < n:
        if debug:
            print "i, j = %d, %d" % (i, j)

        assert i<=j
        if i and check:
            assert H[i:,:j].max() == 0 # XX rm

        # first find a nonzero entry in this col
        for i1 in range(i, m):
            if H[i1, j]:
                break
        else:
            j += 1 # move to the next col
            continue # <----------- continue ------------

        if i != i1:
            if debug:
                print "swap", i, i1
            swap_row(H, i, i1)

        assert H[i, j]
        for i1 in range(i+1, m):
            if H[i1, j]:
                if debug: 
                    print "add %s to %s" % (i, i1)
                H[i1, :] += H[i, :]
                H[i1, :] %= 2

        assert 0<=H.max()<=1

        i += 1
        j += 1

    if truncate:
        m = H.shape[0]-1
        #print "sum:", m, H[m, :], H[m, :].sum()
        while m>=0 and H[m, :].sum()==0:
            m -= 1
        H = H[:m+1, :]

    return H


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

    U1 = zeros2(n, m)

    #print shortstrx(U, U1)

    # Work backwards
    i = len(leading)-1 # <= m
    while i>=0:

        j = leading[i]
        #print "i=", i, "j=", j
        U1[j, i] = 1

        #print shortstrx(U, U1)

        k = i-1
        while k>=0:
            #print "dot2"
            #print shortstr(U[k,:])
            #print shortstr(U1[:,i])
            if dot2(U[k, :], U1[:, i]):
                j = leading[k]
                #print "set", j, i
                U1[j, i] = 1
            #print shortstr(U1[:,i])
            assert dot2(U[k, :], U1[:, i]) == 0
            k -= 1
        i -= 1

    return U1


def l_inverse(L, check=False, verbose=False):
    """invert L (lower triangular, 1 on diagonal)
    """

    m, n = L.shape
    assert m==n
    L1 = identity2(m)

    # Work forwards
    for i in range(m):
        #u = L1[:, i]
        for j in range(i+1, m):
            r = dot2(L[j, :], L1[:, i])
            if r:
                L1[j, i] = 1
            assert dot2(L[j, :], L1[:, i]) == 0

    assert eq2(dot2(L, L1), identity2(m))
    return L1


def plu_reduce(A, truncate=False, check=False, verbose=False):
    """
    Solve PLU = A, st. P is permutation, L is lower tri, U is upper tri.
    Remove zero rows from U if truncate=True.
    """

    m, n = A.shape
    P = identity2(m)
    L = identity2(m)
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
            A1 = dot2(P, dot2(L, U))
            assert eq2(A1, A)

        assert U[i, j]
        for i1 in range(i+1, m):
            if U[i1, j]:
                if verbose: 
                    print "add %s to %s" % (i, i1)
                L[i1, i] = 1
                U[i1, :] += U[i, :]
                U[i1, :] %= 2

        if check:
            A1 = dot2(P, dot2(L, U))
            assert eq2(A1, A)

        i += 1
        j += 1

    if truncate:
        m = U.shape[0]-1
        #print "sum:", m, U[m, :], U[m, :].sum()
        while m>=0 and U[m, :].sum()==0:
            m -= 1
        U = U[:m+1, :]

    return P, L, U


def linear_independent(A, check=False, verbose=False):
    """ Remove linear dependant rows of A
    """

    if A.shape[0] == 0:
        return A

    P, L, U = plu_reduce(A)

    if verbose:
        print "linear_independent"
        print "P L U = A"
        print P.shape, L.shape, U.shape
        print shortstrx(P, L, U, A)
        print

    m, n = U.shape

    # find the zero rows of U
    k = m-1
    while k>=0 and U[k, :].sum() == 0:
        k -= 1

    if k==m-1:
        return A

    L1 = l_inverse(L)
    W = dot2(L1, P.transpose())

    if check:
        assert eq2(U, dot2(W, A))

    if verbose:
        #print m, n
        #print k
        print "U = W A"
        print shortstrx(U, W, A)
        print

    # Key idea: rows of W[k+1:] denote linearly dependant rows of A.
    # So we use the (column) span of W[:k+1]:

    W1 = W[:k+1, :]
    W1 = row_reduce(W1)
    #print shortstr(W1)
    idxs = []
    for i in range(k+1):
        nz = numpy.where(W1[i])[0]
        assert len(nz), "doh!"
        idxs.append(nz[0])

    #print "u:", u
    #idxs = [i for i in range(m) if u[i]]
    #print "idxs:", idxs
    A = A[idxs]

    return A



def find_kernel(A, inplace=False, check=False, verbose=False):
    """return a list of vectors that span the nullspace of A
    """

    if check:
        A0 = A.copy() # save

#    L, U = lu_decompose(A)
#    assert eq2(dot(L, U), A)

    U = row_reduce(A, inplace=inplace)

    # We are looking for a basis for the nullspace of A

    m, n = U.shape

    if verbose:
        print "find_kernel: shape", m, n
        #print shortstr(U, deco=True)
        print

    items = []
    for row in range(m):
        cols = numpy.where(U[row, :])[0]
        if not len(cols):
            break
        col = cols[0]
        items.append((row, col))

    #items.sort(key = lambda item : item[1])
    #print items
    #rows = [row for (row, col) in items]
    #U = U[rows]
    leading = [col for (row, col) in items]
    degeneracy = m - len(leading)

    if verbose:
        print "leading:", leading
        print "degeneracy:", degeneracy

    # Look for the free variables
    vars = []
    row = 0
    col = 0
    while row < m and col < n:
        #print row, col
        if U[row:, col].max() == 0: # XXX optimize this XXX
            #print "*"
            assert U[row:, col].max() == 0, U[row:, col]
            vars.append(col)
        else:
            #print U[row:, col]
            while row<m and U[row:, col].max():
                row += 1
                #print "row", row
                #if row<m:
                #    print U[row:, col]
        col += 1
    for k in range(col, n):
        vars.append(k)

    if verbose:
        print "found %d free vars:" % len(vars), vars

    basis = []
    for var in vars:

        #print "var:", var
        v = numpy.zeros((n,), dtype=numpy.int32)
        v[var] = 1
        row = min(var-1, m-1)
        while row>=0:
            u = dot(U[row], v)
            if u.sum()%2:
                col = leading[row]
                #print "\trow", row, "leading:", col
                v[col] = 1
                #print '\t', shortstr(v)
            assert dot(U[row], v).sum()%2==0, row
            row -= 1
        #print '\t', shortstr(v)
        if check:
            assert dot(A0, v).sum()%2 == 0, shortstr(v)
        basis.append(v)

    return basis


class RowReduction(object):
    "deprecated: use get_reductor"
    def __init__(self, H):
        self.H = H
        self.H1 = H1 = row_reduce(H, truncate=True)

        m, n = H1.shape
        leading = []
        for i in range(m):
            j = i
            while j<n and H1[i, j] == 0:
                j += 1
            if j==n:
                break
            leading.append(j)
        self.leading = leading

    def reduce(self, v):
        """Generate a canonical form for v modulo image of H"""
        H1 = self.H1
        leading = self.leading
        v = v.copy()
        for i, j in enumerate(leading):
            if v[j]:
                v += H1[i]
                v %= 2
        return v


def get_reductor(H):
    """
        Construct a projector P onto the complement of the rowspace of H.
        Multiply with P on the left: u = dot2(P, v).
    """
    H = row_reduce(H, truncate=True)

    A = zeros2(H.shape)
    m, n = H.shape
    for i in range(m):
        j = i
        while j<n and H[i, j] == 0:
            j += 1
        if j==n:
            break
        A[i, j] = 1
        assert H[i, j] == 1
        for i1 in range(i):
            if H[i1, j]:
                H[i1] = (H[i1] + H[i])%2
    #P = identity2(n) + dot2(A.transpose(), H)
    P = identity2(n) + dot2(H.transpose(), A)
    P %= 2
    return P


def kernel(H):
    m, n = H.shape
    for v in enum2(n):
        if dot2(H, v).sum()==0:
            yield v


def image(H):
    m, n = H.shape
    image = set()
    for v in enum2(n):
        u = dot2(H, v)
        su = u.tostring()
        if su in image:
            continue
        image.add(su)
        yield u


def log2(x):
    i = 0
    while 2**i < x:
        i += 1
    return i


def span(vs):
    H = array2(vs)
    return image(H.transpose())


def contains(vs, v):
    for v1 in vs:
        if eq2(v1, v):
            return True
    return False


def find_logops(Hx, Hz, check=False, verbose=False):
    """
        Find kernel of Hx modulo image of Hz^t.
        These are the z type logical operators.
    """

    basis = find_kernel(Hx.copy(), check=check, verbose=verbose)
    if verbose:
        print "find_logops: basis len", len(basis)

    logops = basis

    if not logops:
        shape = 0, Hx.shape[1]
        Lz = zeros2(*shape)
        return Lz

    Lz = array2(logops)
    Lz = row_reduce(Lz)

    redHz = RowReduction(Hz)

    Hz = row_reduce(Hz)

    U = []
    for i in range(Lz.shape[0]):
        #u = reduce(Hz, Lz[i, :])
        #assert eq2(redHz.reduce(Lz[i, :]), u)
        u = redHz.reduce(Lz[i, :])
        assert dot2(Hx, u).sum()==0
        U.append(u)
    Lz = array2(U)
    Lz = linear_independent(Lz, check=check)

    if verbose:
        print "Lz="
        print shortstr(Lz)

    return Lz


def find_stabilizers(Gx, Gz):
    n = Gx.shape[1]
    A = dot2(Gx, Gz.transpose())
    vs = find_kernel(A)
    vs = list(vs)
    #print "kernel GxGz^T:", len(vs)
    Hz = zeros2(len(vs), n)
    for i, v in enumerate(vs):
        Hz[i] = dot2(v.transpose(), Gz)  
    Hz = linear_independent(Hz)
    return Hz


def find_errors(Hx, Lx, Rx=None):
    "find inverse of Hx commuting with Lx"

    if Rx is not None:
        Lx = concatenate((Lx, Rx))

    # find Tz
    n = Hx.shape[1]

    Lx = row_reduce(Lx)
    k = len(Lx)
    mx = len(Hx)

    HL = row_reduce(concatenate((Lx, Hx)))
    assert len(HL) == mx+k
    assert k+mx <= n, (k, mx, n)

    U = zeros2(mx+k, n)
    U[:mx] = Hx
    U[mx:mx+k] = Lx

    B = zeros2(mx+k, mx)
    B[:mx] = identity2(mx)

    Tz_t = solve(U, B)
    assert Tz_t is not None, "no solution"
    Tz = Tz_t.transpose()
    assert len(Tz) == mx

    check_conjugate(Hx, Tz)
    check_commute(Lx, Tz)

    return Tz





def check_conjugate(A, B):
    if A is None or B is None:
        return
    assert A.shape == B.shape
    I = numpy.identity(A.shape[0], dtype=numpy.int32)
    assert eq2(dot2(A, B.transpose()), I)


def check_commute(A, B):
    if A is None or B is None:
        return
    C = dot2(A, B.transpose())
    assert C.sum() == 0, "\n%s"%shortstr(C)





"""
def orthogonal(A):
    A = row_reduce(A, truncate=True)
    m, n = A.shape
    B = zeros2(n-m, n)
    i, j = 0, 0
    for row in A:
        while row[i] == 0:
            B[j, i] = 1
            i += 1
            j += 1
        else:
            i += 1
    while j<n-m:
        B[j, i] = 1
        i += 1
        j += 1
#    print "orthogonal"
#    print shortstrx(A)
#    print
#    print shortstrx(B)
    AB = numpy.concatenate((A, B))
    assert len(row_reduce(AB, truncate=True))==n
    return B
"""


def pseudo_inverse(A, check=False):
    m, n = A.shape
    if m*n == 0:
        A1 = zeros2(n, m)
        return A1
    P, L, U = plu_reduce(A, verbose=False, check=check)
    L1 = l_inverse(L, check=check)
    U1 = u_inverse(U, check=check)
    A1 = dot2(U1, dot2(L1, P.transpose()))
    return A1


def solve(H, u, force=False, debug=False):
    "Solve Hv = u"
    A = pseudo_inverse(H)
    v = dot2(A, u)
    if eq2(dot2(H, v), u) or force:
        return v


def randsolve(A, b, C=10, **kw):
    "Solve Ax = b"
    "http://rjlipton.wordpress.com/2012/08/09/a-new-way-to-solve-linear-equations/"

    m, n = A.shape
    #assert b.shape == (m,)

    N = C*n

    S = rand2(N, n)
    #print S.shape

    for i in range(m):

        #for j in range(i):
        #    assert eq2(dot2(A[i], S[0]), b[i])

        # Selection
        x1 = dot2(A[i], S.transpose())

        T = S[(x1 == b[i])]

        #print T.shape
        #print (dot2(A[i], T.transpose())), b[i]

        if T.shape[0] == 0:
            return None

        # Recombination
        for j in range(N):
            rows = [randint(0, T.shape[0]-1) for k in range(3)]
            S[j] = T[rows[0]] + T[rows[1]] + T[rows[2]]

        S %= 2

        #print dot2(A, S.transpose())
        #print

    if S.shape[0] == 0:
        return None

    x = S[0]

    return x

        

def test_randsolve():

    m, n = 8, 12

    for i in range(100):
        A = rand2(m, n)
        b = rand2(m, 1)
    
        x = randsolve(A, b)

        #print "b:", b
        #print "x:", x
        #print "Ax:", dot2(A, x)

        if x is not None:
            Ax = dot2(A, x)
            Ax.shape = b.shape
            assert eq2(Ax, b)
        else:
            assert solve(A, b) is None # may not always succeed.


def in_rowspace(A, u):
    #if len(u.shape)==1:
    #n = u.shape[0]*u.shape[1]

    At = A.transpose()
    ut = u.transpose()
    B = pseudo_inverse(At)
    v = dot2(B, ut)
    #v = solve(At, ut)

    return eq2(dot2(At, v), ut)


class LinearCombination(object):
    def __init__(self, data={}):
        self.data = dict(data) # map ob to coefficient

    @classmethod
    def promote(cls, item):
        if type(item) is cls:
            lin = item
        elif type(item) in (int, long, float): #, complex, numpy.complex128): # wah...
            if item==0:
                lin = cls()
            else:
                lin = cls({1 : item})
        elif type(item) is tuple:
            lin = cls({item : 1})
        else:
            raise TypeError, (item, type(item))
        return lin

    def __str__(self):
        if not self.data:
            return "0"
        #s = ["%d*(%d,%d)" % (value, key[0], key[1]) 
        #    for (key, value) in self.data.items()]
        s = ["%s*%s" % (value, key)
            for (key, value) in self.data.items()]
        s = "+".join(s)
        s = s.replace("1*(", "(")
        return s
        #return "lin(%s)"%(self.data)
    __repr__ = __str__

    def __add__(self, other):
        other = self.promote(other)
        keys = set(self.data.keys() + other.data.keys())
        data = {}
        for key in keys:
            value = self.data.get(key, 0) + other.data.get(key, 0)
            if value != 0:
                data[key] = value
        return self.__class__(data)
    __radd__ = __add__

    def __neg__(self):
        return 0-self

    def __pos__(self):
        return self

    def __sub__(self, other):
        other = self.promote(other)
        keys = set(self.data.keys() + other.data.keys())
        data = {}
        for key in keys:
            value = self.data.get(key, 0) - other.data.get(key, 0)
            if value != 0:
                data[key] = value
        return self.__class__(data)

    def __rsub__(self, other):
        other = self.promote(other)
        keys = set(self.data.keys() + other.data.keys())
        data = {}
        for key in keys:
            value = other.data.get(key, 0) - self.data.get(key, 0)
            if value != 0:
                data[key] = value
        return self.__class__(data)

    def __mul__(self, other):
        #other = self.promote(other)
        #if type(other) in (int, long):
        if other==0:
            return self.__class__()
        if other==1:
            return self
        data = {}
        for key, value in self.data.items():
            data[key] = other*value
        return self.__class__(data)
    __rmul__ = __mul__

    def __mod__(self, i):
        data = dict((key, value%i) for (key, value) in self.data.items())
        return self.__class__(data)

    def __eq__(self, other):
        other = self.promote(other)
        return self.data == other.data

    def __ne__(self, other):
        return not self==other



def MatrixEntry(i, j):
    return LinearCombination({(i, j) : 1})


def Unknown(*args):
    if len(args)==2:
        m, n = args
    elif len(args)==1:
        m, n = args[0]
    else:
        raise TypeError
    data = [[MatrixEntry(i, j) for j in range(n)] for i in range(m)]
    A = numpy.array(data)
    return A


class System(object):
    "system of linear equations"
    def __init__(self, *Ts):
        assert Ts
        m = sum(T.shape[0] for T in Ts)
        n = sum(T.shape[1] for T in Ts)
        T = Unknown(m, n)
        self.lhss = []
        self.rhss = []
        i = 0
        blocks = []
        for T0 in Ts:
            assert T0.dtype == object # the Unknown
            j = 0
            row = []
            for T1 in Ts:
                row.append( [[i, i+T1.shape[0]], [j, j+T1.shape[1]]] )
                j += T1.shape[1]
            blocks.append(row)
            i += T0.shape[0]
        self.blocks = blocks = numpy.array(blocks, dtype=numpy.int32)
        n = self.n = len(Ts)
        for i in range(n):
          for j in range(n):
            blk = blocks[i, j]
            T1 = T[blk[0, 0]:blk[0, 1], blk[1, 0]:blk[1, 1]]
            if i==j:
                Ts[i][:] = T1
            else:
                self.append(T1, 0)
        self.Ts = Ts
        self.T = T

    def append(self, lhs, rhs):
        if type(rhs) in (int, float, long) and rhs==0:
            rhs = zeros2(*lhs.shape)
        if type(rhs) in (int, float, long) and rhs==1:
            n = min(*lhs.shape)
            rhs = identity2(n)
        if lhs.dtype==numpy.int32:
            lhs, rhs = rhs, lhs # swap
        elif lhs.dtype==object and rhs.dtype==object:
            lhs = lhs-rhs
            rhs = zeros2(*rhs.shape)
        assert lhs.dtype == object, lhs.dtype
        assert rhs.dtype == numpy.int32, rhs.dtype
        if len(lhs.shape)==1:
            lhs = lhs.view()
            lhs.shape = (lhs.shape[0], 1)
        if len(rhs.shape)==1:
            rhs = rhs.view()
            rhs.shape = (rhs.shape[0], 1)
        assert lhs.shape == rhs.shape, (lhs.shape, rhs.shape)
        self.lhss.append(lhs)
        self.rhss.append(rhs)

    def build(self):
        T, lhss, rhss = self.T, self.lhss, self.rhss
        m = sum(lhs.shape[0]*lhs.shape[1] for lhs in lhss)
        n = T.shape[0]*T.shape[1]
        H = zeros2(
            m, # this many constraints
            n) # this many unknowns
        v = zeros2(m, 1)
        row = 0
        for lhs, rhs in zip(lhss, rhss):
          for i0 in range(lhs.shape[0]):
            for j0 in range(lhs.shape[1]):
              lin = lhs[i0, j0]
              assert isinstance(lin, LinearCombination), repr(lhs)
              for (i1, j1), w in lin.data.items():
                  H[row, T.shape[1]*i1+j1] = w
              v[row] = rhs[i0, j0]
              row += 1
        H %= 2
        return H, v

    def kernel(self):
        H, v = self.build()
        basis = find_kernel(H)
        return basis

    def solve(self, unpack=True, check=False):
        H, v = self.build()
        Hinv = pseudo_inverse(H, check=check)
        u = dot2(Hinv, v)
        #print shortstrx(H, u, v)
        if not eq2(dot2(H, u), v):
            # No solution
            return None
        u.shape = self.T.shape
        #print shortstrx(T)
        if unpack:
            u = self.unpack(u)
        return u

    def unpack(self, u):
        if len(self.Ts)>1:
            us = []
            i, j = 0, 0
            for T in self.Ts:
                u0 = u[i:i+T.shape[0], j:j+T.shape[1]]
                i += T.shape[0]
                j += T.shape[1]
                us.append(u0)
            return us
        else:
            return u

    def all_solutions(self):
        u = self.solve(unpack=False)
        kern = self.kernel()
        if not kern:
            return
        for v in span(kern):
            v.shape = u.shape
            yield self.unpack((u+v)%2)

    def random_solution(self):
        u = self.solve(unpack=False)
        kern = self.kernel()
        if not kern:
            return
        for v in kern:
            v.shape = u.shape
            if random()<=0.5:
                u += v
        return self.unpack(u%2)


def rank(H):
    if len(H)==0:
        return 0
    H = row_reduce(H, truncate=True)
    return H.shape[0]


def projector(A, check=False):

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

    P = identity2(m) - dot2(A, pseudo_inverse(A))
    P %= 2

    return P


def fromkernel(J, check=False):
    """
    find f as a pushout of the following diagram:

      J^T
    k ---> n
    |      |
    |      | f
    v      v
    0 ---> n'

    """

    k, n = J.shape
    K = zeros2(0, k)

    f, g = pushout(J.transpose(), K, check=check)

    return f


def pushout(J, K, J1=None, K1=None, check=False):
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
    JJ = zeros2(b+c, b)
    JJ[:b] = identity2(b)

    KK = zeros2(b+c, c)
    KK[b:] = identity2(c)

    kern = (compose2(J, JJ) - compose2(K, KK))%2
    # We need to kill the columns of kern
    R = projector(kern, check=check)
    R = row_reduce(R, truncate=True, check=check)

    assert eq2(compose2(J, JJ, R), compose2(K, KK, R))

    JJ = compose2(JJ, R)
    KK = compose2(KK, R)

    if J1 is not None:
        assert K1 is not None
        assert J1.shape[0] == K1.shape[0]
        assert eq2(compose2(J, J1), compose2(K, K1))
        m = J1.shape[0]
        n = R.shape[0]
        F = zeros2(m, n)

#        print "F=", F.shape
#        print "R=", R.shape
        #print shortstr(R)

        Rinv = pseudo_inverse(R, check=check)

#        print "Rinv=", Rinv.shape
        #print shortstr(Rinv)

        for i in range(n):
            r = Rinv[:, i]
            u, v = r[:b], r[b:]
            u = dot2(J1, u)
            v = dot2(K1, v)
            #assert eq2(u, v)
            uv = (u+v)%2
            F[:, i] = uv

        assert eq2(compose2(JJ, F), J1)
        assert eq2(compose2(KK, F), K1)

        return JJ, KK, F

    else:
        return JJ, KK


# _________________________________________________________ #
#                                                           #
#                       test code                           #
#                                                           #

def test_solve():

    H = numpy.array([
       [0, 1, 0, 0, 0, 0, 1, 0],
       [0, 0, 1, 0, 0, 1, 0, 0],
       [0, 0, 1, 0, 0, 1, 0, 0],
       [1, 0, 1, 1, 0, 1, 0, 0],
       [0, 0, 0, 0, 1, 0, 0, 1],
       [1, 1, 0, 1, 0, 0, 1, 0],
       [0, 0, 0, 0, 1, 0, 0, 1],
       [0, 1, 0, 0, 0, 0, 1, 0],
       [1, 0, 0, 1, 1, 0, 0, 1]], dtype=numpy.int32)

    #u = numpy.array([0, 0, 0, 1, 0, 1, 0, 0, 1], dtype=numpy.int32)

    v = numpy.array([0, 1, 1, 0, 1, 0, 0, 1], dtype=numpy.int32)

    v = numpy.array([
        [0, 1, 1, 0, 1, 0, 0, 1],
        [0, 1, 1, 0, 0, 0, 0, 0]],
        dtype=numpy.int32)

    #print H.shape
    #print v.shape
    v = v.transpose()

    u = dot(H, v) % 2

    #print u

    #v = search(H, u)
    v = solve(H, u)

    #print "v =", v
    #print v.shape

    u1 = dot(H, v) % 2
    #print "u =", u1

    #assert numpy.abs(u-u1)


def test_solve_rand():

    numpy.random.seed(0)
    #m, n = 7, 20
    m, n = 4, 6

    p = 0.2
    for i in range(100):
        H = numpy.random.binomial(1, p, (m, n))
        #print shortstr(H), '\n'
        v = numpy.random.binomial(1, p, (n,))
        u = dot(H, v) % 2

#        H = numpy.array([[0, 0, 1, 0, 0, 0],
#               #[0, 0, 0, 0, 0, 0],
#               #[0, 0, 0, 0, 0, 0],
#               [0, 0, 1, 0, 0, 1]])
#        v = numpy.array([1, 0, 0, 1, 0, 1])
#        #u = numpy.array([0, 0, 0, 1])
#        u = dot(H, v) % 2
#
#        print "H =", repr(H)
#        print "v =", repr(v)
#        print "u =", repr(u)
#
        v1 = solve(H, u, debug=False)
#        print "v1 =", repr(v)
#        return

        assert v1 is not None
        #if v1 is None:
        #    assert 0
        #    continue
        delta = dot(H, v+v1) % 2
        assert numpy.abs(delta).max() == 0
        u1 = dot(H, v1) % 2
        assert ((u1+u)%2).max() == 0
        #print ".",
    #print
    

    p = 0.1
    for i in range(1000):
        H = numpy.random.binomial(1, p, (m, n))
        #v = numpy.random.binomial(1, p, (n,))
        #u = dot(H, v) % 2
        u = numpy.random.binomial(1, p, (m,))
        v = solve(H, u, debug=False)
        if v is None:
            continue
        u1 = dot(H, v) % 2
        delta = (u+u1)%2
        assert numpy.abs(delta).max() == 0
        #print ".",
    #print
    

def test_logops_toric():

    # Toric code Hz
    Hz = parse("""
    111......1......................
    ..111......1....................
    ....111......1..................
    1.....11.......1................
    ........111......1..............
    ..........111......1............
    ............111......1..........
    ........1.....11.......1........
    ................111......1......
    ..................111......1....
    ....................111......1..
    ................1.....11.......1
    .1......................111.....
    ...1......................111...
    .....1......................111.
    """)

    Hx = parse("""
    11.....1................1.......
    .111......................1.....
    ...111......................1...
    .....111......................1.
    1.......11.....1................
    ..1......111....................
    ....1......111..................
    ......1......111................
    ........1.......11.....1........
    ..........1......111............
    ............1......111..........
    ..............1......111........
    ................1.......11.....1
    ..................1......111....
    ....................1......111..
    """)

    z_ops = find_logops(Hx, Hz, check=True)
    assert len(z_ops) == 2

    #z_ops = [shortstr(op) for op in z_ops]
    #z_ops.sort()
    #print z_ops
    #assert z_ops == [
    #    '.1.1.1.1........................',
    #    '1.......1.......1.......1.......']


def check_logops(Hx):

    Hz = Hx.copy()
    m, n = Hx.shape # 6, 16

    assert (dot(Hx, Hz.transpose())%2).sum() == 0

    # logops is kern(Hx)/Im(Hz^t)
    z_ops = find_logops(Hx, Hz, check=True)
    assert len(z_ops) == n-2*m # 4

    # kernel Hx has dimension n-m
    kernHx = list(kernel(Hx))
    d = log2(len(kernHx))
    assert d == n-m

    imHzt = list(image(Hz.transpose()))
    d = log2(len(imHzt))
    assert d == m

    for v in imHzt:
        assert dot2(Hx, v).sum()==0

    for v in span(z_ops):
        if v.sum()==0:
            continue
        assert dot2(Hx, v).sum()==0 # in the kernel
        assert not contains(imHzt, v), shortstr(v) # not in the image

    # linear independent:
    assert len(list(image(z_ops))) == 2**len(z_ops)


def test_logops_bicycle():

    Hx = parse("""
    .1..111.1.11.1..
    ...1.1111...111.
    1.1...11.1...111
    .11..1.111.1..1.
    11.11.....111..1
    1.111....11.1..1
    """)

    #Hx = row_reduce(Hx)
    #print shortstr(Hx)

    Hx = parse("""
    1..1...1.111.11.
    .11.1.1...111..1
    ..11.11.1.11.1..
    ...1...1..111111
    .....1..1.11....
    """)

    check_logops(Hx)


def test_plu():

    m, n = 9, 9
    rows = [[int(i==j) for i in range(n)] for j in range(n)]
    shuffle(rows)
    A = numpy.array(rows, dtype=numpy.int32)

    A = parse("""
    1........
    .1.......
    ..1......
    ...1.....
    .....1...
    ....1....
    ......1..
    .......1.
    ........1
    """)

    A = parse("""
    1...1....
    .1...1...
    ..1...1..
    ...1...1.
    ....1...1
    .....1...
    ...1..1..
    ....1..1.
    .....1..1
    """)

    A = parse("""
    1...1.....
    .1...1....
    ..1...1...
    ...1...1..
    ....1...1.
    .....1....
    ...1..1...
    ....1..1..
    .....1..1.
    """)

    # Toric code Hz
    A = parse("""
    111......1......................
    ..111......1....................
    ....111......1..................
    1.....11.......1................
    ........111......1..............
    ..........111......1............
    ............111......1..........
    ........1.....11.......1........
    ................111......1......
    ..................111......1....
    ....................111......1..
    ................1.....11.......1
    .1......................111.....
    ...1......................111...
    .....1......................111.
    """)

    _A = parse("""
    111....
    ..111..
    ....111
    1.....1
    """)

    _A = parse("""
    ..11
    1...
    .111
    """)

    _A = parse("""
    111.
    ..11
    ...1
    1...
    """)

    #print "A="
    #print shortstr(A)

    #P, L, U = plu_decompose(A)
    P, L, U = plu_reduce(A, verbose=False, check=True)

    L1 = l_inverse(L)
    #print shortstrx(L, L1, dot2(L, L1))

    U1 = u_inverse(U)
    #print shortstrx(U, U1, dot2(U, U1))
    assert eq2(dot2(U, U1), identity2(U.shape[0]))

    is_lower = True
    m, n = L.shape
    for i in range(m):
      for j in range(i+1, m):
        if L[i, j]:
            is_lower = False
    #print "L is_lower", is_lower
    assert is_lower

    is_upper = True
    m, n = L.shape
    for i in range(m):
      for j in range(i):
        if U[i, j]:
            is_upper = False
    #print "U is_upper", is_upper
    assert is_upper

#test_plu()
#sys.exit(0)


def test_lu():

    # Toric code Hz
    A = parse("""
    111......1......................
    ..111......1....................
    ....111......1..................
    1.....11.......1................
    ........111......1..............
    ..........111......1............
    ............111......1..........
    ........1.....11.......1........
    ................111......1......
    ..................111......1....
    ....................111......1..
    ................1.....11.......1
    .1......................111.....
    ...1......................111...
    .....1......................111.
    """)

    L, U = lu_decompose(A)

    #print shortstrx(L, U, A)
    #print "LU == A", 

    assert eq2(dot(L, U), A)


def test_conjugate_ops():

    # Toric code Lz and Hz
    Uz = parse("""
    .1.1.1.1........................
    1.......1.......1.......1.......
    111......1......................
    ..111......1....................
    ....111......1..................
    1.....11.......1................
    ........111......1..............
    ..........111......1............
    ............111......1..........
    ........1.....11.......1........
    ................111......1......
    ..................111......1....
    ....................111......1..
    ................1.....11.......1
    .1......................111.....
    ...1......................111...
    .....1......................111.
    """)

    _Uz = parse("""
    .1...........1111...
    ......1.....1..1.1..
    .........11.11111.11
    ..........1.1.11..1.
    ............11..1111
    .............1..1.11
    111.1.1.1....111...1
    ..11.1..1.1...1.1..1
    ...1.11..11.1..1111.
    ....1.11.1.1.1.1..1.
    .....111.1.1...11.1.
    .......111.11....111
    ........11..1....111
    """)

    Tx = pseudo_inverse(Uz, check=True)
    Tx = Tx.transpose()
    assert Tx.shape==Uz.shape
    I = dot(Uz, Tx.transpose())%2
    assert eq2(identity2(Tx.shape[0]), I)


def test_linear_independent():

    A = parse("""
    ...1111.......
    .11..11.......
    .1111.........
    1.1.1.1.......
    1.11.1........
    11..11........
    11.1..1.......
    ..........1111
    ........11..11
    ........1111..
    .......1.1.1.1
    .......1.11.1.
    .......11..11.
    """)

    spA = list(span(A))
    AA = linear_independent(A, check=True, verbose=False)
    #print shortstr(AA)
    assert AA.shape == row_reduce(A, truncate=True).shape
    #print shortstr(row_reduce(A))
    spAA = list(span(AA))

    #print "LI:"
    #print shortstr(AA)

    #print "reduced:"
    #print shortstr(row_reduce(AA))

    assert len(spA) == len(spAA)
    for x in spA:
        assert contains(spAA, x)
    for x in spAA:
        assert contains(spA, x)

def test_entry():

#    m, n = 5, 12
#
#    print A
#
#    B = zeros2(n, m)
#    B[1,1] = 1
#
#    print numpy.dot(A, B)

    H = parse("""
    11.11...
    .111..1.
    1...11.1
    """)

    m, n = H.shape

    T = Unknown(n, m)

    #print numpy.dot(H, T)

    system = System(T)
    system.append(dot2(H, T), identity2(m))

    T = system.solve()

    assert eq2(dot2(H, T), identity2(m))


def test_pushout():

    J = zeros2(2, 1)
    J[0, 0] = 1
    K = zeros2(2, 1)
    K[1, 0] = 1

    JJ, KK = pushout(J, K)

    #print shortstrx(JJ, KK)

    assert eq2(compose2(J, JJ), compose2(K, KK))

    J1, K1 = JJ, KK

    JJ, KK, F = pushout(J, K, J1, K1)

    #print
    #print shortstr(F)
    assert eq2(F, identity2(3))


def test_fromkernel():

    n = 54
    J = zeros2(1, n)
    for i in [0, 2, 8, 10, 14, 16, 22, 28, 32, 34, 37, 40, 42, 44, 49, 51]:
        J[0, i] = 1

    f = fromkernel(J)

    assert f.shape == (53, 54)

    assert compose2(J.transpose(), f).sum() == 0


def test_reductor():
    #seed(0)
    #numpy.random.seed(0)
    m = 10
    n = 20
    for trial in range(10):
        H = rand2(m, n)
        rr = RowReduction(H)
        P = get_reductor(H)
        v = rand2(1, n)
        v.shape = (n,)
        vP = dot2(P, v)
        rrv = rr.reduce(v)
        vP.shape = rrv.shape
        #print shortstr(vP)
        #print shortstr(rrv)
        assert eq2(vP, rrv)



if __name__=="__main__":

    test_solve()
    test_solve_rand()
    test_randsolve()
    test_lu()
    test_plu()
    test_logops_bicycle()
    test_logops_toric()
    test_conjugate_ops()
    test_linear_independent()
    test_entry()
    test_pushout()
    test_fromkernel()
    test_reductor()

    print "OK"


