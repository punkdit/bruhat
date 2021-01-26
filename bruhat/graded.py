#!/usr/bin/env python3

"""
Earlier version: qupy.ldpc.cell
Used by: bruhat.morse
See also: bruhat.vec 
"""

from random import shuffle, seed
from math import sin, cos, pi

import numpy

from bruhat.argv import argv
from bruhat import element
from bruhat import elim
from bruhat.elim import shortstr
from bruhat.solve import parse


class Cell(object):
    "The basis elements for the vector space at a grade in a Chain complex."

    # XXX enforce key as a tuple ???

    def __init__(self, grade, key=None, name=None, **attrs):
        self.grade = grade
        #if not isinstance(key, tuple):
        #    key = (key,)
        self.key = key
        self.__dict__.update(attrs)
        self._name = name

    pos = None # for layout code
    infty = False # for layout code

    @property
    def name(self):
        if self._name is not None:
            return self._name
        c = 'vef'[self.grade] # vertex, edge, face, ...
        key = self.key
        if type(key) is int:
            key = str(key)
        elif type(key) is tuple:
            key = ''.join(str(k) for k in key)
        else:
            key = str(key)
        return "%s_{%s}"%(c, key)

    def __str__(self):
        return self.name

    def __repr__(self):
        #return "Cell(%d, %s)"%(self.grade, self.key)
        return "<%d:%r>"%(self.grade, self.key)

    def __lt__(self, other):
        return self.key < other.key

    def __le__(self, other):
        return self.key <= other.key

    def tensor(self, other):
        "graded tensor product"
        assert isinstance(other, Cell)
        left = self.key
        right = other.key
        left = left if isinstance(left, tuple) else (left,) 
        right = right if isinstance(right, tuple) else (right,)
        key = left + right
        grade = self.grade + other.grade
        return Cell(grade, key)
    __matmul__ = tensor

    def dual(self):
        key = self.key 
        key = key if type(key) is tuple else (key,)
        #if not key or key[0] != "*":
        #    key = ("*",)+key
        #elif key[0]=="*":
        #    key == key[1:]
        return Cell(-self.grade, key)

    def __eq__(self, other):
        return self.grade==other.grade and self.key==other.key

    def __ne__(self, other):
        return self.grade!=other.grade or self.key!=other.key

    def __hash__(self):
        return hash((self.grade, self.key))


the_star = Cell(0, ("*",)) # tensor unit

class Matrix(object):
    """ 
        Matrix over arbitrary indexing Cell objects (not just numbers).
        rows and cols : list of hashable items
        elements : map (row, col) -> value, such that value != 0
    """

    def __init__(self, rows, cols, elements, ring, src_grade=None, tgt_grade=None):
        for item in rows:
            assert isinstance(item, Cell)
            assert tgt_grade==None or tgt_grade==item.grade
            tgt_grade = item.grade
        for item in cols:
            assert isinstance(item, Cell)
            assert src_grade==None or src_grade==item.grade
            src_grade = item.grade
        rows = list(rows)
        cols = list(cols)
        self.rows = rows
        self.cols = cols
        self.set_rows = set(rows)
        self.set_cols = set(cols)
        self.elements = dict(elements)
        self.src_grade = src_grade
        self.tgt_grade = tgt_grade
        self.grade = None
        if tgt_grade is not None and src_grade is not None:
            self.grade = tgt_grade - src_grade
        self.ring = ring
        self.one = ring.one
        self.zero = ring.zero
        self.check() # REMOVE ME TODO <-----------------

    @property
    def shape(self):
        return (len(self.rows), len(self.cols))

    @classmethod
    def unit(cls, ring):
        one = ring.one
        return Matrix([the_star], [the_star], {(the_star,the_star):one}, ring, 0, 0)

    @classmethod
    def zero(cls, ring):
        return Matrix([], [], {}, ring, 0, 0)

    def check(self):
        elements = self.elements
        set_rows = self.set_rows
        set_cols = self.set_cols
        for ((i,j),v) in elements.items():
            assert i in set_rows
            assert j in set_cols
            assert v != 0

    @classmethod
    def from_array(cls, M, ring, rows=None, cols=None, src_grade=0, tgt_grade=0):
        if isinstance(M, list):
            M = numpy.array(M)
        m, n = M.shape
        # WARNING:
        # by default we reuse labels for cols & rows
        if cols is None:
            cols = [i+1 for i in range(n)]
            cols = [Cell(src_grade, col) for col in cols]
        if rows is None:
            rows = [i+1 for i in range(m)]
            rows = [Cell(tgt_grade, row) for row in rows]
        A = cls(rows, cols, {}, ring)
        for i in range(m):
          for j in range(n):
            if M[i, j]:
                A[rows[i], cols[j]] = ring.promote(M[i, j])
        return A

    def to_array(self):
        zero = self.ring.zero
        elements = self.elements
        rows, cols = self.rows, self.cols
        items = [[elements.get((row, col), zero) 
            for row in rows] for col in cols]
        A = numpy.array(items)
        return A

    @property
    def shape(self):
        return len(self.rows), len(self.cols)

    def copy(self):
        return Matrix(self.rows, self.cols, self.elements, self.ring, 
            self.src_grade, self.tgt_grade)

    def __len__(self):
        return len(self.rows)

    def __str__(self):
        #return "Matrix(%d, %d)"%(len(self.rows), len(self.cols))
        s = str(shortstr(self.todense()))
        s = s.replace(" 0 ", " . ")
        return s

    def longstr(self):
        elements = self.elements
        keys = list(elements.keys())
        keys.sort()
        items = ',\n'.join("%s:%s*%s"%(col, elements[row,col], row) for (row,col) in keys)
        return "{%s}"%items

    def __eq__(self, other):
        return self.elements == other.elements

    def __ne__(self, other):
        return self.elements != other.elements

    def dual(self):
        elements = dict()
        rows = [row.dual() for row in self.rows]
        cols = [col.dual() for col in self.cols]
        for ((i, j), v) in self.elements.items():
            elements[j.dual(), i.dual()] = v
        M = Matrix(cols, rows, elements, self.ring, -self.tgt_grade, -self.src_grade)
        return M

    @classmethod
    def identity(cls, items, ring, grade=None):
        one = ring.one
        elements = dict(((i, i), one) for i in items)
        A = cls(items, items, elements, ring, grade, grade)
        return A

    @classmethod
    def inclusion(cls, rows, cols, ring, grade=None):
        A = cls(rows, cols, {}, ring, grade, grade)
        one = ring.one
        for i in cols:
            A[i, i] = one
        return A

    @classmethod
    def retraction(cls, rows, cols, ring, grade=None):
        A = cls(rows, cols, {}, ring, grade, grade)
        one = ring.one
        for i in rows:
            A[i, i] = one
        return A

    def keys(self):
        return self.elements.keys()

    def __setitem__(self, key, value):
        row, col = key
        #idx = self.row_lookup[row]
        #jdx = self.col_lookup[col]
        assert row in self.set_rows
        assert col in self.set_cols
        elements = self.elements
        value = self.ring.promote(value)
        if value != self.zero:
            elements[row, col] = value
        elif (row, col) in elements:
            del elements[row, col]
        self.check() # REMOVE ME TODO <-----------------

    def __getitem__(self, item):
        row, col = item

        if isinstance(row, slice) and isinstance(col, slice):
            assert row == slice(None)
            assert col == slice(None)
            value = self

        elif isinstance(col, slice):
            assert col == slice(None)
            value = []
            assert row in self.set_rows
            for (r,c) in self.elements.keys():
                if r == row:
                    value.append(c)

        elif isinstance(row, slice):
            assert row == slice(None)
            value = []
            assert col in self.set_cols
            for (r,c) in self.elements.keys():
                if c == col:
                    value.append(r)

        else:
            assert row in self.set_rows
            assert col in self.set_cols
            value = self.elements.get((row, col), 0)
        return value

    def __add__(self, other):
        assert self.ring == other.ring
        assert self.set_cols == other.set_cols
        assert self.set_rows == other.set_rows
        A = self.copy()
        for i,j in other.elements.keys():
            A[i, j] += other[i, j]
        return A

    def __sub__(self, other):
        assert self.ring == other.ring
        assert self.set_cols == other.set_cols
        assert self.set_rows == other.set_rows
        A = self.copy()
        for i,j in other.elements.keys():
            A[i, j] -= other[i, j]
        return A

    def __neg__(self):
        els = dict()
        elements = self.elements
        for ((i,j),v) in elements.items():
            els[i, j] = -v
        return Matrix(self.rows, self.cols, els, self.ring)

    def __rmul__(self, r):
        r = self.ring.promote(r)
        els = dict()
        elements = self.elements
        for ((i,j),v) in elements.items():
            els[i, j] = r*v
        return Matrix(self.rows, self.cols, els, self.ring)

    def __mul__(self, other):
        assert self.ring == other.ring
        assert self.set_cols == other.set_rows, (
            "%s != %s" % (self.set_cols, other.set_rows))
        A = Matrix(self.rows, other.cols, {}, self.ring, other.src_grade, self.tgt_grade)
        for (i, j) in self.elements.keys():
            for k in other.cols:
                A[i, k] = A[i, k] + self[i, j] * other[j, k]
        return A

    def tensor(left, right):
        assert left.ring == right.ring
        ring = left.ring
        elements = dict()
        rows = [li@ri for li in left.rows for ri in right.rows]
        cols = [lj@rj for lj in left.cols for rj in right.cols]
        lels = left.elements
        rels = right.elements
        for ((li,lj),lv) in lels.items():
          for ((ri,rj),rv) in rels.items():
            v = lv*rv
            if v != 0:
                row = li@ri
                col = lj@rj
                elements[row, col] = v
        return Matrix(rows, cols, elements, ring, 
            left.src_grade+right.src_grade, left.tgt_grade+right.tgt_grade)
    __matmul__ = tensor

    def cokernel(a): # TODO
        A = a.to_array()
        B = elim.cokernel(a.ring, A)
        n = B.shape[0]
        #X = Space(n, a.ring)
        #hom = Hom(a.tgt, X)
        b = Matrix.from_array(B, a.ring)
        return b

    def kernel(a): # TODO
        A = a.to_array()
        B = elim.kernel(a.ring, A)
        n = B.shape[1] # src
        #X = Space(n, a.ring)
        #hom = Hom(X, a.src)
        b = Matrix.from_array(B, a.ring)
        return b

    def image(a): # TODO
        A = a.to_array()
        At = A.transpose()
        At = elim.row_reduce(a.ring, At, truncate=True)
        A = At.transpose()
        #X = Space(A.shape[1], a.ring)
        #hom = Hom(X, a.tgt)
        b = Matrix.from_array(A, a.ring)
        return b

    def rank(a):
        A = a.to_array()
        d = elim.rank(a.ring, A)
        return d

    def sum(self):
        return sum(self.elements.values(), self.zero)

    def nonzero(self):
        zero = self.zero
        for v in self.elements.values():
            if v != zero:
                return True
        assert not self.elements
        return False

    def is_zero(self):
        return not self.nonzero()

    def todense(self, rows=None, cols=None):
        if rows is None:
            rows = self.rows
        if cols is None:
            cols = self.cols
        assert set(rows) == self.set_rows, "%s\n%s"%(set(rows), self.set_rows)
        assert set(cols) == self.set_cols, "%s\n%s"%(set(cols), self.set_cols)
        row_lookup = dict((row, idx) for (idx, row) in enumerate(rows))
        col_lookup = dict((col, idx) for (idx, col) in enumerate(cols))
        A = numpy.zeros((len(rows), len(cols)), dtype=object)
        elements = self.elements
        for (row, col) in elements.keys():
            value = elements[row, col]
            A[row_lookup[row], col_lookup[col]] = value
        return A



class Chain(object):
    "chain complex"
    def __init__(self, cells, bdys, ring, check=False):
        self.ring = ring
        self.cells = dict(cells) # map grade -> list of Cell's
        self.bdys = dict(bdys) # map grade -> Matrix(_, grade)
        lookup = {}
        for grade, cells in cells.items():
            for cell in cells:
                lookup[cell.grade, cell.key] = cell
        self.lookup = lookup
        if check:
            self.check()

    def get_grades(self):
        "keys of bdys, reverse sorted"
        grades = list(self.bdys.keys())
        grades.sort(reverse=True)
        return grades

    def get_mingrade(self):
        g = min(self.bdys.keys(), default=0)
        return g-1

    def get_maxgrade(self):
        g = max(self.bdys.keys(), default=0)
        return g

    @classmethod
    def from_array_i(cls, M, ring, rows=None, cols=None): # arghhh...
        src_grade = 1
        tgt_grade = 0
        A = Matrix.from_array(M, ring, rows, cols, src_grade, tgt_grade)
        chain = cls({tgt_grade:A.rows, src_grade:A.cols}, {src_grade:A}, ring)
        return chain

    @classmethod
    def from_array(cls, M, ring):
        src_grade = 1
        tgt_grade = 0
        if isinstance(M, list):
            M = numpy.array(M)
        m, n = M.shape
        # We need distinct names for the row cells & col cells...
        cols = [i+1 for i in range(n)]
        cols = [Cell(src_grade, "e%s"%col, "e_{%s}"%col) for col in cols]
        rows = [i+1 for i in range(m)]
        rows = [Cell(tgt_grade, "v%s"%row, "v_{%s}"%row) for row in rows]
        A = Matrix.from_array(M, ring, rows, cols, src_grade, tgt_grade)
        chain = cls({tgt_grade:A.rows, src_grade:A.cols}, {src_grade:A}, ring)
        return chain

    def dump(self):
        print("=========== Chain ==============")
        for grade in self.get_grades():
            M = self.get_bdymap(grade)
            print("grade:", M.src_grade, "-->", M.tgt_grade)
            print(M.cols, "-->")
            print("  ", M.rows)
            print(M)
        print("================================")

    def check(self):
        try:
            self._check()
        except AssertionError:
            print("================================")
            print("******* AssertionError *********")
            self.dump()
            raise

    def _check(self):
        cells = self.cells
        grades = [grade for grade in cells if cells[grade]]
        grades.sort()
        if not grades:
            return
        g0 = grades[0]
        g1 = grades[-1]
        assert len(grades) == g1-g0+1, "must have contiguous grades"
        bdys = self.bdys
        if g0 in bdys:
            A = bdys[g0]
            assert A.rows == [], A.rows
        #assert not grades or min(grades)>=0
        for src_grade in grades:
            if src_grade > g0:
                assert src_grade in bdys
            A = self.get_bdymap(src_grade)
            A.check()
            tgt_grade = A.tgt_grade
            assert tgt_grade == src_grade - 1
            assert A.cols == cells[src_grade]
            assert A.rows == cells.get(tgt_grade, [])
            B = self.get_bdymap(tgt_grade)
            assert (B*A).is_zero(), (B*A)

    def get_cells(self, grade=None):
        if grade is None:
            cells = []
            for grade in self.cells.keys():
                cells += self.cells[grade]
        else:
            cells = self.cells.setdefault(grade, [])
        return cells

    def get_bdymap(self, grade):
        A = self.bdys.get(grade)
        if A is None: # default to zero
            A = Matrix(self.get_cells(grade-1), self.get_cells(grade), {}, self.ring, 
                grade, grade-1)
            if A.rows or A.cols:
                self.bdys[grade] = A
        return A.copy()

    def dual(self, check=False):
        bdys = dict()
        cells = {}
        for (grade, M) in self.bdys.items():
            Mt = M.dual()
            bdys[Mt.src_grade] = Mt
            cells[Mt.src_grade] = Mt.cols
            cells[Mt.tgt_grade] = Mt.rows
        chain = Chain(cells, bdys, self.ring, check=check)
        return chain

    def tensor(left, right):
        assert left.ring == right.ring
        ring = left.ring
        lgrades = left.get_grades()
        assert lgrades
        rgrades = right.get_grades()
        assert rgrades
        lids = {}
        for lg in lgrades:
            M = left.bdys[lg]
            I = Matrix.identity(M.cols, ring, lg)
            lids[lg] = I
        lg -= 1
        cols = M.rows
        I = Matrix.identity(cols, ring, lg)
        lids[lg] = I

        rids = {}
        for rg in rgrades:
            M = right.bdys[rg]
            I = Matrix.identity(M.cols, ring, rg)
            rids[rg] = I
        rg -= 1
        cols = M.rows
        I = Matrix.identity(cols, ring, rg)
        rids[rg] = I

        lgrades.append(lgrades[-1]-1)
        rgrades.append(rgrades[-1]-1)
        pairs = {}
        spaces = {}
        #print()
        for lg in lgrades:
          for rg in rgrades:
            assert lids[lg].cols == left.get_bdymap(lg).cols
            assert rids[rg].cols == right.get_bdymap(rg).cols
            H = lids[lg] @ right.get_bdymap(rg)
            V = left.get_bdymap(lg) @ rids[rg]
            key = lg, rg
            pairs[key] = (H, V)
            #print(H.cols, len(H.cols))
            #print(V.cols, len(V.cols))
            assert H.cols == V.cols

            src = H.cols
            if spaces.get(key) is not None:
                assert spaces[key] == src
            else:
                spaces[key] = src

            tgt = H.rows
            key = lg, rg-1
            if spaces.get(key) is not None:
                assert spaces[key] == tgt
            else:
                spaces[key] = tgt

            tgt = V.rows
            key = lg-1, rg
            if spaces.get(key) is not None:
                assert spaces[key] == tgt
            else:
                spaces[key] = tgt

        keys = list(spaces.keys())
        keys.sort(reverse=True)
        #print("keys:", keys)
        grades = list(set([sum(key) for key in keys]))
        grades.sort(reverse=True)
        #print("grades:", grades)
        cells = dict((g, []) for g in grades)
        for lg, rg in keys:
            g = lg + rg
            space = spaces[lg, rg]
            #print(lg, rg, space)
            assert set(space).intersection(cells[g]) == set(), \
                "non disjoint bases"
            cells[g] += space
        bdys = {}
        assert set(grades) == set(cells.keys())
        #grades = list(cells.keys())
        #grades.sort(reverse=True)
        for g in grades:
            if g-1 not in cells:
                continue
            rows = cells[g-1]
            cols = cells[g]
            M = Matrix(rows, cols, {}, ring, g, g-1)
            bdys[g] = M

        for lg,rg in pairs.keys():
            H, V = pairs[lg,rg]
            #print("lg, rg = ", lg, rg)
            #if not H.is_zero() or 1:
            #    print("H:", H.shape)
            #    print(H)
            #if not V.is_zero() or 1:
            #    print("V:", V.shape)
            #    print(V)
            #print()

            g = H.src_grade
            assert g is not None
            M = bdys[g]
            #print(H.cols)
            #print("------>")
            #print(M.rows)
            sign = -1 if (lg%2) else 1
            H = Matrix.inclusion(M.rows, H.rows, ring) \
                * H * Matrix.retraction(H.cols, M.cols, ring)
            M = M + sign * H

            V = Matrix.inclusion(M.rows, V.rows, ring) \
                * V * Matrix.retraction(V.cols, M.cols, ring)
            M = M + V

            bdys[g] = M

        chain = Chain(cells, bdys, ring, True)
        return chain
                
    __matmul__ = tensor

    def transvect(self, i, j, grade=1, check=True):
        """
            apply _transvection i-->j (a.k.a. CNOT gate) at grade
        """
        ring = self.ring
        one = ring.one
        C0 = self.get_cells(grade-1)
        C1 = self.get_cells(grade)
        C2 = self.get_cells(grade+1)
        
        #Sx = Matrix(C1, C2, {}, ring)
        #Sz = Matrix(C0, C1, {}, ring)

        if type(i) is int:
            i = C1[i]
        if type(j) is int:
            j = C1[j]
    
        f0 = Matrix.identity(C0, ring)
        f1 = Matrix.identity(C1, ring)
        f1[j, i] = one
        f1_inv = Matrix.identity(C1, ring)
        f1_inv[j, i] = -one
        assert f1*f1_inv == Matrix.identity(C1, ring)
        f2 = Matrix.identity(C2, ring)

        Sx = f1*self.get_bdymap(grade+1)
        Sz = self.get_bdymap(grade)*f1_inv
        tgt = Chain({0:C0,1:C1,2:C2}, {1:Sz,2:Sx}, ring, check=check)
        hom = Hom(self, tgt, {0:f0,1:f1,2:f2}, ring, check=check)

        return tgt, hom



class Hom(object):
    "Chain complex map"
    def __init__(self, src, tgt, fs, ring, check=True):
        assert isinstance(src, Chain)
        assert isinstance(tgt, Chain)
        self.ring = ring
        self.fs = dict(fs) # map grade --> Matrix
        self.src = src
        self.tgt = tgt
        self.g0 = min(src.get_mingrade(), tgt.get_mingrade())
        self.g1 = max(src.get_maxgrade(), tgt.get_maxgrade())
        if check:
            self.check()

    def get(self, grade):
        f = self.fs.get(grade)
        if f is None:
            f = Matrix(
                self.tgt.get_cells(grade), 
                self.src.get_cells(grade), 
                {}, self.ring, grade, grade)
            self.fs[grade] = f
        return f

    def check(self):
        fs = self.fs
        src, tgt = self.src, self.tgt
        for grade in range(self.g0, self.g1+1):
            f = self.get(grade)
            g = self.get(grade-1)

            lhs = g*src.get_bdymap(grade)
            rhs = tgt.get_bdymap(grade)*f
            assert lhs == rhs

    #def homology(self):




def test_matrix(ring):
    one = ring.one
    zero = ring.zero

    rows = cols = [Cell(0, i) for i in [0, 1]]
    I = Matrix.from_array([[1,0],[0,1]], ring, rows, cols)
    X = Matrix.from_array([[0,1],[1,0]], ring, rows, cols)
    Z = Matrix.from_array([[1,0],[0,-1]], ring, rows, cols)

    assert I*X == X
    assert X*X == I
    assert Z*Z == I
    assert ring.one==-ring.one or X*Z != Z*X
    assert X*Z == -Z*X
    assert (-1)*I == -I
    assert (I + -I).is_zero()

    II = I@I
    XI = X@I
    IZ = I@Z
    assert XI*XI == II
    assert XI*IZ == X@Z

    M = parse("""
    1.1..1
    11.1..
    .11.1.
    ...111
    """)
    A = Matrix.from_array(M, ring)
    #print(A)
    A.check()

    if 0:
        print("kernel:")
        print(A.kernel())
        print("image:")
        print(A.image())
        print("cokernel:")
        print(A.cokernel())
    
    At = A.dual()
    #print(A*At)
    At.check()

    M = [[1,0,1],[1,1,0]]
    C = Chain.from_array(M, ring)
    Ct = C.dual()

    #Ct.dump()
    #print(Ct.cells)
    #print(Ct.get_cells(-1))
    #M = Ct.get_bdymap(-1)
    #print(M.cols, M.rows)

    Ct.check()

    #C.dump()
    #Ct.dump()
    A = C @ Ct
    #test_chain(A)

    A.check()


def test_main():

    for ring in [element.Q, element.FiniteField(2), element.FiniteField(7)]:
    #for ring in [element.FiniteField(2)]:

        test_matrix(ring)

        #break


if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        seed(_seed)

    test_main()

    print("OK")


