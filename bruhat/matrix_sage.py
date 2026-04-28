#!/usr/bin/env python

"""

"""

from random import choice
from operator import mul, matmul, add
from functools import reduce
#from functools import cache
from functools import lru_cache
cache = lru_cache(maxsize=None)

import numpy

from sage.all_cmdline import (FiniteField, CyclotomicField, latex, block_diagonal_matrix,
    PolynomialRing)
from sage import all_cmdline as sage

from bruhat.action import mulclose, mulclose_names, mulclose_find
from bruhat.argv import argv
from bruhat.smap import SMap


def unify(R, S):
    if R.coerce_map_from(S) is not None:
        return R
    if S.coerce_map_from(R) is not None:
        return S
    assert 0


K = CyclotomicField(8)
w8 = K.gen()
w4 = w8**2
r2 = w8 + w8**7
assert r2**2 == 2


def simplify_latex(self):
    M = self.M
    m, n = self.shape
    idxs = [(i,j) for i in range(m) for j in range(n)]
    for idx in idxs:
        if M[idx] != 0:
            break
    else:
        assert 0
    scale = M[idx]
    if scale != 1 and scale != -1:
        M = (1/scale) * M
        s = {
            r2 : r"\sqrt{2}",
            1/r2 : r"\frac{1}{\sqrt{2}}",
            #2/r2 : r"\frac{2}{\sqrt{2}}",
            2/r2 : r"\sqrt{2}",
            #r2/2 : r"\sqrt{2}/2",
        }.get(scale, latex(scale))
        if "+" in s:
            s = "("+s+")"
        s = "%s %s"%(s, latex(M))
    else:
        s = latex(M)
    s = s.replace(r"\zeta_{8}^{2}", "i")
    return s



class Matrix(object):
    def __init__(self, ring, rows, name=(), shape=None):
        if shape is not None:
            m, n = shape
            M = sage.Matrix(ring, m, n, rows)
            assert M.nrows() == m
            assert M.ncols() == n
        else:
            M = sage.Matrix(ring, rows)
        M.set_immutable()
        self.M = M
        self.ring = ring
        self.shape = (M.nrows(), M.ncols())
        if type(name) is str:
            name = (name,)
        self.name = name

    @classmethod
    def promote(cls, ring, rows, name=()):
        if isinstance(rows, cls):
            return rows
        return cls(ring, rows, name)

    def __eq__(self, other):
        assert isinstance(other, Matrix)
        #assert self.ring == other.ring
        #assert self.shape == other.shape
        return self.M == other.M

    def __hash__(self):
        return hash(self.M)

    def __str__(self):
        lines = str(self.M).split("\n")
        lines[0] = "[" + lines[0]
        lines[-1] = lines[-1] + "]"
        lines[1:] = [" "+l for l in lines[1:]]
        s = '\n'.join(lines)
        s = s.replace(" 0", " .")
        s = s.replace("[0 ", "[. ")
        s = s.replace(" 0]", " .]")
        return s
    __repr__ = __str__

    def __len__(self):
        return self.shape[0]

    def __mul__(self, other):
        assert isinstance(other, Matrix), type(other)
        #assert self.ring == other.ring
        ring = unify(self.ring, other.ring)
        assert self.shape[1] == other.shape[0], (
            "cant multiply %sx%s by %sx%s"%(self.shape + other.shape))
        M = self.M * other.M
        name = self.name + other.name
        return Matrix(ring, M, name)

    def __add__(self, other):
        assert isinstance(other, Matrix)
        #assert self.ring == other.ring
        ring = unify(self.ring, other.ring)
        M = self.M + other.M
        return Matrix(ring, M)

    def __sub__(self, other):
        assert isinstance(other, Matrix)
        #assert self.ring == other.ring
        ring = unify(self.ring, other.ring)
        M = self.M - other.M
        return Matrix(ring, M)

    def __neg__(self):
        M = -self.M
        return Matrix(self.ring, M)

    def __pow__(self, n):
       assert n>=0
       if n==0:
           return Matrix.identity(self.ring, self.shape[0])
       return reduce(mul, [self]*n)

    def __rmul__(self, r):
        M = r*self.M
        return Matrix(self.ring, M)

    def __matmul__(self, other):
        assert isinstance(other, Matrix)
        #assert self.ring == other.ring
        ring = unify(self.ring, other.ring)
        M = self.M.tensor_product(other.M)
        return Matrix(ring, M)
    tensor_product = __matmul__

    def direct_sum(self, other):
        assert isinstance(other, Matrix)
        #assert self.ring == other.ring
        ring = unify(self.ring, other.ring)
        #M = self.M.direct_sum(other.M)
        M = block_diagonal_matrix(self.M, other.M)
        return Matrix(ring, M)
    __lshift__ = direct_sum

    def stack(self, other):
        "concatenate rows of self & other"
        #assert self.ring == other.ring
        ring = unify(self.ring, other.ring)
        return Matrix(ring, self.M.stack(other.M))

    def augment(self, other):
        "concatenate cols of self & other"
        #assert self.ring == other.ring
        ring = unify(self.ring, other.ring)
        return Matrix(ring, self.M.augment(other.M))

    def __getitem__(self, idx):
        #if type(idx) is int:
        #    return self.M[idx]
        M = self.M[idx]
        #print(type(M), type(self.M))
        return Matrix(self.ring, M)

    def reshape(self, m, n):
        shape = self.shape
        assert m*n == shape[0]*shape[1]
        M = self.M
        els = [M[i,j] for i in range(shape[0]) for j in range(shape[1])]
        rows = [[els[n*i+j] for j in range(n)] for i in range(m)]
        M = Matrix(self.ring, rows)
        assert M.shape == (m,n)
        return M

    def _latex_(self):
        M = self.M
        s = M._latex_()
        if "zeta_" not in s:
            return s
        return simplify_latex(self)

    @classmethod
    def identity(cls, ring, n):
        rows = []
        for i in range(n):
            row = [0]*n
            row[i] = 1
            rows.append(row)
        return Matrix(ring, rows)

    @classmethod
    def zero(cls, ring, m, n):
        rows = [[0]*n for i in range(n)]
        return Matrix(ring, rows)

    @classmethod
    def get_perm(cls, ring, perm):
        n = len(perm)
        cols = []
        for i in perm:
            col = [0]*n
            col[i] = 1
            cols.append(col)
        M = Matrix(ring, cols)
        return M.t

    def order(self):
        I = self.identity(self.ring, len(self))
        count = 1
        A = self
        while A != I:
            count += 1
            A = self*A
        return count

    def inverse(self):
        M = self.M.inverse()
        return Matrix(self.ring, M)
    __invert__ = inverse

    def pseudoinverse(self, algorithm=None):
        """
        Return the Moore-Penrose pseudoinverse of this matrix.

        INPUT:

        - ``algorithm`` -- (default: guess) one of the following:

          - ``'numpy'`` -- use numpy's ``linalg.pinv()`` which is
            suitable over real or complex fields

          - ``'exact'`` -- use a simple algorithm which is not
            numerically stable but useful over exact fields. Assume that
            no conjugation is needed, that the conjugate transpose is
            just the transpose.

          - ``'exactconj'`` -- like ``exact`` but use the conjugate
            transpose
        """

        M = self.M.pseudoinverse(algorithm=algorithm)
        return Matrix(self.ring, M)

    def to_numpy(self):
        u = numpy.empty(self.shape, dtype=object)
        m, n = self.shape
        for i in range(m):
          for j in range(n):
            u[i,j] = self.M[i,j]
        return u

    def transpose(self):
        M = self.M.transpose()
        n, m = self.shape
        return Matrix(self.ring, M, shape=(m,n))

    def trace(self):
        return self.M.trace()

    @property
    def t(self):
        return self.transpose()

    def dagger(self):
        M = self.M.conjugate_transpose()
        name = tuple(name+".d" for name in reversed(self.name))
        return Matrix(self.ring, M, name)

    def real(self):
        ring = self.ring
        M = self.M + self.M.conjugate()
        M = (ring.one()/2)*M
        return Matrix(self.ring, M)

    def conjugate(self):
        M = self.M.conjugate()
        return Matrix(self.ring, M)

    def imag(self):
        ring = self.ring
        M = self.M - self.M.conjugate()
        M = (ring.one()/2)*M
        return Matrix(self.ring, M)

    def coerce(self, ring):
        return Matrix(ring, self.M)

    @property
    def d(self):
        return self.dagger()

    def is_diagonal(self):
        M = self.M
        return M.is_diagonal()

    def is_zero(self):
        return self == -self

    def is_identity(self):
        return self == Matrix.identity(self.ring, len(self))

    def rank(self):
        return self.M.rank()

    def determinant(self):
        return self.M.determinant()

    def eigenvectors(self):
        evs = self.M.eigenvectors_right()
        spaces = []
        for val,vecs,dim in evs:
            vecs = sage.Matrix(self.ring, vecs)
            vecs = vecs.transpose() 
            vecs = Matrix(self.ring, vecs)
            assert vecs.shape[1] == dim
            spaces.append((val, vecs, dim))
        #print(type(self.M))
        return spaces

    def solve(self, other):
        assert len(self) == len(other)
        try:
            A = self.M.solve_right(other.M)
        except ValueError:
            return None
        return Matrix(self.ring, A)

    def row_reduce(self):
        M = self.M.echelon_form()
        #print(self.M.echelon_form.__doc__)
        return Matrix(self.ring, M)

    def cokernel(self):
        #print("cokernel")
        #print(self, self.shape)
        K = sage.kernel(self.M)
        B = K.basis()
        #print(K)
        #print(K.degree(), type(K))
        n = K.degree()
        m = len(B)
        M = Matrix(self.ring, B, shape=(m,n))
        assert M.shape == (m,n), M.shape
        return M

    def kernel(self):
        K = self.t.cokernel().t
        return K

    def sum(self):
        m, n = self.shape
        u = sum(self.M[i,j] for i in range(m) for j in range(n))
        return u

    def intersect_rowspace(self, other):
        assert self.shape[1] == other.shape[1]
        V1, V2 = self, other
        V = V1.stack(-V2)
        K = V.cokernel()
        if len(K)==0:
            return Matrix.zero(self.ring, 0, self.shape[1])
        #print()
        #print(V1, V1.shape)
        #print(V2, V2.shape)
        #print(K, K.shape)
        W = K[:, :len(V1)]*V1
        #u = V1.t.solve(W.t)
        #assert u is not None
        ##print("solve:", u.shape)
        #u = V2.t.solve(W.t)
        #assert u is not None
        ##print("solve:", u.shape)
        return W

    def puncture(A, i):
        m, n = A.shape
        assert 0<=i<n
        A0 = A[:, :i]
        A1 = A[:, i+1:]
        return A0.augment(A1)
    delete = puncture

    def contract(A, i):
        m, n = A.shape
        assert 0<=i<n
        B = A.t.cokernel()
        Bi = B.puncture(i)
        C = Bi.kernel().t
        return C

    def is_loop(self, i):
        m, n = self.shape
        assert 0<=i<n
        #print(self, "is_loop", i, self.shape)
        #print(self[:,i])
        #result = self[:, i].sum() == 0
        #print(result)
        for j in range(m):
            if self.M[j,i] != 0:
                return False
        return True

    def is_coloop(self, i):
        m, n = self.shape
        assert 0<=i<n
        #print(self, "is_coloop", i)
        H = self.t.cokernel()
        #print("H =")
        #print(H, H.shape)
        return H.is_loop(i)

    def _tutte(self, x, y, depth=0):
        m, n = self.shape
        if n == 0:
            return 1
        i = 0
        #print(" "*depth+ "_tutte", i)
        #smap = SMap()
        #smap[0,depth] = str(self)
        #print(smap)
        if self.is_loop(i) and self.is_coloop(i):
            print("FAIL")
            print(self)
            assert 0
        if self.is_loop(i):
            assert not self.is_coloop(i)
            #print(" "*depth+ "is_loop")
            lhs = self.delete(i)
            p = y*lhs._tutte(x, y, depth+1)
        elif self.is_coloop(i):
            assert not self.is_loop(i)
            #print(" "*depth+ "is_coloop")
            rhs = self.contract(i)
            p = x*rhs._tutte(x, y, depth+1)
        else:
            #print(" "*depth+ "else")
            lhs = self.delete(i)
            rhs = self.contract(i)
            lhs = lhs._tutte(x, y, depth+1)
            rhs = rhs._tutte(x, y, depth+1)
            p = lhs+rhs
        return p

    def get_tutte(self):
        R = sage.PolynomialRing(sage.ZZ, list("xy"))
        x, y, = R.gens()
        p = self._tutte(x, y)
        return p


def test_tutte():
    R = sage.PolynomialRing(sage.ZZ, list("xy"))
    x, y, = R.gens()

    from bruhat.dev.geometry import all_codes

    F = sage.FiniteField(2)
    u = Matrix(F, [[1]])
    assert u.is_coloop(0)
    assert not u.is_loop(0)

    M = Matrix(F, [[1,0,1],[0,1,1]])
    assert not M.is_coloop(0)
    assert not M.is_loop(0)

    M1 = M.contract(0)
    assert M1 == Matrix(F, [[1,1]])

    p = M.get_tutte()
    assert p == x**2 + x + y

    # --------------------------

    q = argv.get("q", 3)
    #m = argv.get("m", 2)
    n = argv.get("n", 6)

    F = sage.FiniteField(q)

    row = 0
    for m in range((n+2)//2):
        count = 0
        found = set()
        for Gt in all_codes(m, n, q):
            Gt = Matrix(F, Gt)
            #print(Gt, type(Gt))
            p = Gt.get_tutte()
            found.add(p)
            count += 1
        print(m, count, len(found))
        #row += len(found)
    #print("row =", row)


def get_orbits(F, n, found):
    gen = []
    for i in range(n-1):
        idxs = list(range(n))
        idxs[i:i+2] = [i+1, i]
        rows = [[0]*n for i in range(n)]
        for i,j in enumerate(idxs):
            rows[i][j] = 1
        #g = lambda M,idxs=idxs : M[:, idxs]
        g = Matrix(F, rows)
        gen.append(g)

    for i in range(n):
        rows = [[0]*n for i in range(n)]
        for j in range(n):
            if i==j:
                rows[j][j] = 2
            else:
                rows[j][j] = 1
        #print(rows)
        g = Matrix(F, rows)
        #print(g)
        gen.append(g)


    orbits = []
    remain = set(found)
    while remain:
        M = remain.pop()
        orbit = {M}
        bdy = list(orbit)
        while bdy:
            _bdy = []
            for M1 in bdy:
              #print("="*n)
              #print(M1)
              for g in gen:
                #print(idxs)
                #M2 = M1[:, idxs]
                #M2 = g(M1)
                M2 = M1 * g
                #print(M2)
                M2 = M2.row_reduce()
                #print(M2)
                if M2 in orbit:
                    #print()
                    continue
                #print("new!")
                orbit.add(M2)
                _bdy.append(M2)
                remain.remove(M2)
            bdy = _bdy
        orbits.append(orbit)
    return orbits


def test_enum():
    from bruhat.dev.geometry import all_codes
    from bruhat.util import allperms

    q = argv.get("q", 3)
    #m = argv.get("m", 2)
    n = argv.get("n", 6)

    #idxs = list(range(n))
    #perms = list(allperms(idxs))
    #print(len(perms))

    F = sage.FiniteField(q)

    row = 0
    for m in range((n+2)//2):
        count = 0
        found = []
        for Gt in all_codes(m, n, q):
            Gt = Matrix(F, Gt)
            found.append(Gt)
        print(m, len(found), end=" ", flush=True)

        orbits = get_orbits(F, n, found)
        print(len(orbits))
        #for o in orbits:
        #    print(list(o)[0], len(o))







    


def test():
    "Orthogonal matrixes over GF(2)[x] / x**2 "
    from util import cross
    from random import shuffle
    from sage.all_cmdline import GF, PolynomialRing
    K = GF(2)
    R = PolynomialRing(K, "x")
    x = R.gens()[0]
    S = R.quo(x**2)
    print(S)

    x = S.gens()[0]
    zero = x**2
    one = (x+1)**2
    print("one:", one)
    print("x:", x)
    #return

    els = [zero, one, one+x, x]
    found = []
    n = argv.get("n", 2)
    I = Matrix.identity(S, n)
    rows = list(cross([els]*n))
    shuffle(rows)
    #for cols in cross( [rows]*n ):
    print("search:")
    while 1:
        cols = [choice(rows) for i in range(n)]
        M = Matrix(S, cols)
        if M*M.t == I:
            found.append(M)
            print('/', end='', flush=True)
            if len(found) > 1:
                G = mulclose(found, verbose=True)
                print("|G|=", len(G))
                #del G # free mem
                #break
    print()
    #for g in G:
    #    print(g)


def test_linear():

    R = sage.QQ
    M = Matrix(R, [
        [1,1,1,1,0],
        [1,2,3,4,0],
    ])

    K = M.kernel()
    print(K.t)

    print(M*K)


if __name__ == "__main__":

    from time import time
    start_time = time()

    profile = argv.profile
    name = argv.next() or "test"
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        test()


    t = time() - start_time
    print("OK! finished in %.3f seconds\n"%t)




