#!/usr/bin/env python3

"""
compute homology using algebraic morse theory.

Earlier version: qupy.ldpc.cell
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
        if cols is None:
            cols = [i+1 for i in range(n)]
            cols = [Cell(src_grade, "c%s"%col) for col in cols]
        if rows is None:
            rows = [i+1 for i in range(m)]
            rows = [Cell(tgt_grade, "r%s"%row) for row in rows]
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
        X = Space(n, a.ring)
        hom = Hom(a.tgt, X)
        b = Map.from_array(B, hom)
        return b

    def kernel(a): # TODO
        A = a.to_array()
        B = elim.kernel(a.ring, A)
        n = B.shape[1] # src
        X = Space(n, a.ring)
        hom = Hom(X, a.src)
        b = Map.from_array(B, hom)
        return b

    def image(a): # TODO
        A = a.to_array()
        At = A.transpose()
        At = elim.row_reduce(a.ring, At, truncate=True)
        A = At.transpose()
        X = Space(A.shape[1], a.ring)
        hom = Hom(X, a.tgt)
        b = Map.from_array(A, hom)
        return b

    def rank(a): # TODO
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

    def get_degree(self): # XXX REMOVE THIS XXX
        cells = self.cells
        grades = [grade for grade in cells if cells[grade]]
        grades.sort()
        return max(grades)+1 if grades else 0

    @classmethod
    def from_array(cls, M, ring, rows=None, cols=None):
        src_grade = 1
        tgt_grade = 0
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
        # bit of a shit-show here... 
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
        if check:
            self.check()

    def check(self):
        fs = self.fs
        src, tgt = self.src, self.tgt
        for grade in fs.keys(): # um, what grades to check ?? XXX
            f = fs[grade]
            g = fs.get(grade-1)
            if g is None:
                continue

            lhs = g*src.get_bdymap(grade)
            rhs = tgt.get_bdymap(grade)*f
            assert lhs == rhs




class Assembly(object):
    "An assembly of cells, used for constructing Chain complexes."

    def __init__(self, lookup, ring):
        # lookup: map grade to map of key->cell
        self.lookup = dict(lookup)
        self.bdys = {} # map cell to bdy map
        self.ring = ring
        self.zero = ring.zero
        self.one = ring.one

    def set(self, key, cell, bdy={}):
        assert isinstance(cell, Cell)
        lookup = self.lookup
        cells = lookup.setdefault(cell.grade, {})
        assert cells.get(key) is None
        cells[key] = cell
        assert cell.key is None
        cell.key = key
        bdys = self.bdys
        assert bdys.get(cell) is None
        bdys[cell] = dict(bdy)

    def mk_vert(self, key):
        cell = Cell(0)
        self.set(key, cell)
        return cell

    def mk_edge(self, key, src, tgt):
        one = self.one
        cell = Cell(1)
        bdy = {self[0, src]:-one, self[0, tgt]:one}
        self.set(key, cell, bdy)
        return cell

    def __getitem__(self, key):
        grade, key = key
        return self.lookup[grade][key]

    def get_cells(self, grade):
        cells = list(self.lookup.setdefault(grade, {}).values())
        cells.sort()
        return cells

    def get_bdymap(self, grade): # warning: needs to return new object
        # grade -> grade-1
        cols = self.get_cells(grade)
        rows = self.get_cells(grade-1)
        bdys = self.bdys
        A = Matrix(rows, cols, {}, self.ring, grade, grade-1)
        zero = self.zero
        one = self.one
        for col in cols:
            assert col.grade == grade
            bdy = bdys[col]
            #for row in col:
            for row, value in bdy.items():
                assert row.grade == grade-1
                assert A[row, col] == zero
                A[row, col] = value
        return A

    def get_degree(self):
        lookup = self.lookup
        grades = [grade for grade in lookup if lookup[grade]]
        grades.sort()
        return max(grades)+1 if grades else 0

    def get_chain(self):
        ring = self.ring

        n = self.get_degree()
        cells = {}
        bdys = {}
        for grade in range(n+1):
            cells[grade] = self.get_cells(grade)
            bdys[grade] = self.get_bdymap(grade)
        chain = Chain(cells, bdys, ring)
        return chain

    @classmethod
    def build_tetrahedron(cls, ring):
        cx = Assembly({}, ring)
        one = cx.one
    
        radius = 1.5 # for layout
        for idx in range(4):
            cell = cx.mk_vert(idx+1)
            if idx == 0:
                x, y = 0, 0
            else:
                theta = 2*pi*idx/3.
                x, y = radius*sin(theta), radius*cos(theta)
            cell.pos = x, y
    
        cx.mk_edge(1, 3, 1) # inner edge
        cx.mk_edge(2, 1, 2) # inner edge
        cx.mk_edge(3, 1, 4) # inner edge
        cx.mk_edge(4, 3, 2) # outer edge
        cx.mk_edge(5, 4, 2) # outer edge
        cx.mk_edge(6, 4, 3) # outer edge
    
        cx.set(1, Cell(2), {cx[1, 1]:-one, cx[1, 6]:-one, cx[1, 3]:-one})
        cx.set(2, Cell(2), {cx[1, 1]: one, cx[1, 2]: one, cx[1, 4]:-one})
        cx.set(3, Cell(2), {cx[1, 2]:-one, cx[1, 3]: one, cx[1, 5]: one})
        cell = Cell(2)
        cell.infty = True # note for layout code
        theta = 2*pi/6.
        cell.pos = -radius*sin(theta), radius*cos(theta)
        cx.set(4, cell, {cx[1, 4]: one, cx[1, 5]:-one, cx[1, 6]: one})
    
        return cx

    @classmethod
    def build_surface(cls, ring, top_left, bot_right,
            open_top=False, open_bot=False, 
            open_left=False, open_right=False):
        # open boundary is "rough". default is smooth

        cx = Assembly({}, ring)
        one = cx.one

        i0, j0 = top_left
        i1, j1 = bot_right

        # build verts
        for i in range(i0, i1):
          for j in range(j0, j1):
            if i==i0 and open_top:
                continue
            if i==i1-1 and open_bot:
                continue
            if j==j0 and open_left:
                continue
            if j==j1-1 and open_right:
                continue
            cell = Cell(0)
            cell.pos = j, -i
            key = (i, j)
            cx.set(key, cell)

        verts = cx.lookup[0]

        # build edges
        for i in range(i0, i1):
          for j in range(j0, j1):
            ks = "hv" # horizontal, _vertical edge
            if i==i1-1 and j==j1-1:
                continue
            elif i==i1-1:
                ks = "h"
            elif j==j1-1:
                ks = "v"
            for k in ks:
                bdy = {}
                vert = verts.get((i, j))
                if vert is not None:
                    bdy[vert] = one
                if k == "h":
                    vert = verts.get((i, j+1))
                    pos = j+0.5, -i
                else:
                    vert = verts.get((i+1, j))
                    pos = j, -(i+0.5)
                if vert is not None:
                    bdy[vert] = -one
                if len(bdy)==0:
                    continue
                cell = Cell(1, pos=pos)
                key = (i, j, k)
                cx.set(key, cell, bdy)

        edges = cx.lookup[1]

        # build faces
        for i in range(i0, i1-1):
          for j in range(j0, j1-1):
            top = edges.get((i, j, "h"))
            left = edges.get((i, j, "v"))
            bot = edges.get((i+1, j, "h"))
            right = edges.get((i, j+1, "v"))
            bdy = {top:one, left:-one, bot:-one, right:one}
            bdy = dict((cell, value) 
                for (cell, value) in bdy.items() if cell)
            if len(bdy)==0:
                continue
            cell = Cell(2)
            cx.set((i, j), cell, bdy)

        return cx

    @classmethod
    def build_torus(cls, ring, rows, cols):
        cx = Assembly({}, ring)
        assert not cx.lookup
        one = cx.one

        # build verts
        for i in range(rows):
          for j in range(cols):
            cell = Cell(0)
            key = (i, j)
            cell.pos = i, j
            cx.set(key, cell)

        verts = cx.lookup[0]

        # build edges
        for i in range(rows):
          for j in range(cols):
            ks = "hv" # horizontal, _vertical edge
            for k in ks:
                a_vert = verts.get((i, j))
                assert a_vert is not None
                if k == "h":
                    b_vert = verts.get((i, (j+1)%cols))
                else:
                    b_vert = verts.get(((i+1)%rows, j))
                assert b_vert is not None
                key = (i, j, k)
                cx.set(key, Cell(1), {a_vert:-one, b_vert:one})

        edges = cx.lookup[1]

        # build faces
        for i in range(rows):
          for j in range(cols):
            top = edges[(i, j, "h")]
            left = edges[(i, j, "v")]
            bot = edges[((i+1)%rows, j, "h")]
            right = edges[(i, (j+1)%cols, "v")]
            cell = Cell(2)
            bdy = {top:one, left:-one, bot:-one, right:one}
            cx.set((i, j), cell, bdy)

        return cx


class Field(object):
    "height field on the cells of a chain complex"

    def __init__(self, chain):
        assert isinstance(chain, Chain)

        nbd = {} # map cell --> list of neighbour cells
        cells = chain.get_cells()
        for cell in cells:
            nbd[cell] = set()

        for grade in chain.get_grades():
            A = chain.get_bdymap(grade)
            for row, col in A.elements.keys():
                assert A[row, col] != 0
                nbd[row].add(col)
                nbd[col].add(row)

        self.nbd = nbd
        self._clamp = {}
        self.cells = cells
        self.chain = chain

    def clamp(self, item, value):
        if isinstance(item, Cell):
            cell = item
        else:
            for cell in self.cells:
                #print(cell.name)
                if cell.name == item:
                    break
                if (cell.grade, cell.key) == item:
                    break
            else:
                assert 0, "cell %r not found"%item
        #print("clamp", cell, value)
        self._clamp[cell] = value

    def solve(self):
        nbd = self.nbd
        cells = self.cells
        clamp  = self._clamp
        field = dict((cell, 0.) for cell in cells)
        field.update(clamp)

        interior = [cell for cell in cells if cell not in clamp]
        for _ in range(1000):
            shuffle(interior)
            for cell in interior:
                bdy = nbd[cell]
                assert bdy
                value = sum(field[c] for c in bdy) / len(bdy)
                field[cell] = value
        return field

    def get_flow(self):
        nbd = self.nbd
        chain = self.chain
        field = self.solve()
        flow = Flow(chain)
        #pairs = flow.get_all_pairs()
        #shuffle(pairs)
        #for (row, col) in pairs:
        #    assert row.grade == col.grade-1
        #    if field[row] > field[col] and flow.accept(row, col):
        #        flow.add(row, col)

        done = False
        while not done:
            done = True
            for grade in range(2):
                cells = chain.get_cells(grade)
                for c0 in cells:
                    cs = [c for c in nbd[c0] if c.grade==grade+1]
                    cs = [c for c in cs if field[c] < field[c0]]
                    cs = [c1 for c1 in cs if flow.accept(c0, c1)]
                    if len(cs)!=1:
                        continue
                    c1 = cs[0]
                    flow.add(c0, c1)
                    done = False

        for grade in range(2):
            cells = chain.get_cells(grade)
            for c0 in cells:
                cs = [c for c in nbd[c0] if c.grade==grade+1]
                cs = [c for c in cs if field[c] < field[c0]]
                cs = [c1 for c1 in cs if flow.accept(c0, c1)]
                if len(cs)<2:
                    continue
                cs.sort(key = lambda c1:field[c1] - field[c0])
                c1 = cs[0]
                flow.add(c0, c1)

        return flow

    def show(self):
        "use polyscope : https://polyscope.run/py/ "

        assert 0, "TODO"
        import polyscope as ps

        ps.init()

        ### Register a point cloud
        # `my_points` is a Nx3 numpy array
        #ps.register_point_cloud("my points", my_points)

        ### Register a mesh
        # `verts` is a Nx3 numpy array of vertex positions
        # `faces` is a Fx3 array of indices, or a nested list
        # ...

        ps.register_surface_mesh("my mesh", verts, faces, smooth_shade=True)

        ps.show()



class Flow(object):
    "morse matching of a chain complex"

    def __init__(self, chain):
        assert isinstance(chain, Chain)
        self.chain = chain
        self.ring = chain.ring
        self.zero = chain.ring.zero
        self.one = chain.ring.one

        # These are the pairs of the Morse matching we are building.
        # Map grade -> list of pairs of cells (a,b) 
        # with a.grade==grade and b.grade==grade+1 .
        self.pairs = {} 

    def add(self, src, tgt, check=False):
        "Add a match "
        if isinstance(src, tuple):
            src = self.chain.lookup[src]
        if isinstance(tgt, tuple):
            tgt = self.chain.lookup[tgt]
        assert isinstance(src, Cell)
        assert isinstance(tgt, Cell)
        assert src.grade == tgt.grade-1
        #assert src in tgt
        #print("add", src.grade, tgt.grade, list(self.pairs.keys()))
        if check:
            assert self.accept(src, tgt), "%s --> %s DISSALLOWED" % (src, tgt)
        pairs = self.pairs.setdefault(src.grade, [])
        #print("add", list(self.pairs.keys()))
        #print()
        pairs.append((src, tgt))

    def add_match(self, grade, src, tgt, check=True):
        src = self.chain.lookup[grade, src]
        tgt = self.chain.lookup[grade+1, tgt]
        self.add(src, tgt, check)

    def remove(self, src, tgt):
        self.pairs[src.grade].remove((src, tgt))

    def get_pairs(self, grade=None):
        if grade is None:
            pairs = []
            for grade, _pairs in self.pairs.items():
                pairs += _pairs
            return pairs
        pairs = self.pairs.get(grade, [])
        return pairs

    def all_grades(self):
        grades = set()
        for grade in self.pairs.keys():
            grades.add(grade)
            grades.add(grade+1)
        grades = list(grades)
        grades.sort()
        return grades

    def get_critical(self, grade=None):
        " return list of critical cells at grade "
        if grade is None:
            critical = []
            for grade in self.all_grades():
                critical += self.get_critical(grade) # recurse
            return critical

        chain = self.chain
        remove = set()
        for _grade, pairs in self.pairs.items():
          for src, tgt in pairs:
            remove.add(src)
            remove.add(tgt)
        critical = []
        cells = chain.get_cells(grade)
        for cell in cells:
            if cell not in remove:
                critical.append(cell)
        return critical

    def get_match(self, grade):
        " return linear op: C_grade --> C_{grade+1} "
        chain = self.chain
        zero = self.zero
        one = self.one
        bdy = chain.get_bdymap(grade+1) # grade+1 --> grade
        src = chain.get_cells(grade)
        tgt = chain.get_cells(grade+1)
        A = Matrix(tgt, src, {}, self.ring) # C_grade --> C_{grade+1}
        pairs = self.get_pairs(grade)
        for src, tgt in pairs: # grade --> grade+1
            assert src.grade == grade
            assert A[tgt, src] == zero
            value = bdy[src, tgt]
            assert value != zero, "cannot match non-bdy"
            A[tgt, src] = -one/value
        return A

    def get_down_match(self, grade):
        # C_grade --> C_grade-1 --> C_grade
        chain = self.chain
        zero = self.zero
        one = self.one
        bdy = chain.get_bdymap(grade) # grade --> grade-1
        pairs = self.get_pairs(grade-1)
        A = self.get_match(grade-1) # C_{grade-1} --> C_grade
        for src, tgt in pairs: # grade-1 --> grade
            assert bdy[src, tgt] != zero
            bdy[src, tgt] = zero
        return A*bdy

    def get_up_match(self, grade):
        # C_grade --> C_grade+1 --> C_grade
        chain = self.chain
        zero = self.zero
        one = self.one
        bdy = chain.get_bdymap(grade+1) # grade+1 --> grade
        pairs = self.get_pairs(grade)
        A = self.get_match(grade)
        for src, tgt in pairs: # grade --> grade+1
            assert bdy[src, tgt] != zero
            bdy[src, tgt] = zero
        return bdy*A

    def accept(self, src, tgt):
        assert isinstance(src, Cell)
        assert isinstance(tgt, Cell)
        assert src.grade == tgt.grade-1
        pairs = self.get_pairs(src.grade)
        for pair in pairs:
            if pair[0] == src or pair[1] == tgt:
                return False
        pairs = self.get_pairs(src.grade-1)
        for pair in pairs:
            if pair[1] == src:
                return False
        pairs = self.get_pairs(src.grade+1)
        for pair in pairs:
            if pair[0] == tgt:
                return False

        chain = self.chain
        zero = self.zero
        self.add(src, tgt)
        for grade in self.all_grades():
            A = self.get_down_match(grade)
            B = A.copy()
            cells = chain.get_cells(grade)
            while B.nonzero():
                B = A*B
                for j in cells:
                    if B[j, j] != zero: # found a loop
                        self.remove(src, tgt)
                        return False 
        self.remove(src, tgt)
        return True

    def __str__(self):
        grades = list(self.pairs.keys())
        grades.sort()
        lines = []
        for grade in grades:
            pairs = self.pairs[grade]
            if not pairs:
                continue
            line = ', '.join("%s->%s" % (src.key, tgt.key) for (src, tgt) in pairs)
            lines.append("grade=%d: %s"%(grade, line))
        return "Flow(%s)"%(', '.join(lines),)

    def get_all_pairs(self):
        pairs = []
        chain = self.chain
        n = chain.get_degree()
        for grade in range(1, n):
            bdy = chain.get_bdymap(grade)
            for row, col in bdy.keys(): # random order in keys
                assert row.key
                assert col.key
                pairs.append((row, col))
        pairs.sort()
        return pairs

    def build(self):
        "randomly build a Morse matching "
        pairs = self.get_all_pairs()
        shuffle(pairs)

        idx = 0
        while idx < len(pairs):
            src, tgt = pairs[idx]
            if self.accept(src, tgt):
                self.add(src, tgt)
            idx += 1

    def get_bdymap(self, grade): # grade --> grade-1
        chain = self.chain
        src = self.get_critical(grade) # CM_grade
        tgt = self.get_critical(grade-1) # CM_{grade-1}
    
        bdy = Matrix(tgt, src, {}, self.ring) # CM_grade --> CM_{grade-1}
        A = self.get_down_match(grade)  # C_grade --> C_grade-1 --> C_grade
        B = chain.get_bdymap(grade) # C_grade --> C_grade-1
        cells = chain.get_cells(grade)
        cells_1 = chain.get_cells(grade-1)
        L = Matrix.retraction(tgt, cells_1, self.ring) # C_grade-1 --> CM_grade-1
        B = L*B
        R = Matrix.inclusion(cells, src, self.ring) # CM_grade --> C_grade
        while R.nonzero():
            bdy += B*R
            R = A*R
        return bdy

    def get_f(self, grade):
        # chain map, from cells --> critical cells
        chain = self.chain
        src = chain.get_cells(grade)
        tgt = self.get_critical(grade)
    
        f = Matrix(tgt, src, {}, self.ring)
        A = self.get_up_match(grade)  # grade --> grade+1 --> grade
        C = Matrix.retraction(tgt, src, self.ring)
        while C.nonzero():
            f += C
            C = C*A
        return f

    def get_g(self, grade):
        # chain map, from critical cells --> cells
        chain = self.chain
        src = self.get_critical(grade)
        tgt = chain.get_cells(grade)
    
        g = Matrix(tgt, src, {}, self.ring)
        A = self.get_down_match(grade)  # grade --> grade-1 --> grade
        C = Matrix.inclusion(tgt, src, self.ring)
        while C.nonzero():
            g += C
            C = A*C
        return g

    def get_homotopy(self, grade):
        # grade --> grade+1
        ring = self.ring
        chain = self.chain
        tgt = chain.get_cells(grade+1)
        src = chain.get_cells(grade)
        A = self.get_match(grade)
        B = self.get_down_match(grade+1)
        chi = Matrix(tgt, src, {}, ring)
        while A.nonzero():
            chi += A
            A = B*A
        return chi




def test_chain(chain):

    chain.check()
    ring = chain.ring

    #print(chain.get_bdymap(1))

#    return
#
#    # FIX TODO TODO 
#    # FIX TODO TODO 
#    # FIX TODO TODO 
#    # FIX TODO TODO 
#    # FIX TODO TODO 
#    # FIX TODO TODO 
#    # FIX TODO TODO 

    flow = Flow(chain)
    flow.build()

    #print(flow)

    A = flow.get_bdymap(2) # 2 --> 1
    B = flow.get_bdymap(1) # 1 --> 0

    #print(A)
    #print(B)
    C = B*A
    assert C.is_zero()

    #for grade in [0, 1, 2]:
    #  for cell in flow.get_critical(grade):
    #    print(cell)

    f0 = flow.get_f(0)
    f1 = flow.get_f(1)
    f2 = flow.get_f(2)

    assert f0*chain.get_bdymap(1) == B*f1
    assert f1*chain.get_bdymap(2) == A*f2

    g0 = flow.get_g(0)
    g1 = flow.get_g(1)
    g2 = flow.get_g(2)

    assert chain.get_bdymap(1)*g1 == g0*B
    assert chain.get_bdymap(2)*g2 == g1*A

    fs = [f0, f1, f2]
    gs = [g0, g1, g2]

    #print()
    #for grade in range(3):
    #    print("grade =", grade)
    #    f = flow.get_f(grade)
    #    print(f)

    chi0 = flow.get_homotopy(0)
    chi1 = flow.get_homotopy(1)
    chi2 = flow.get_homotopy(2)

    chis = [chi0, chi1, chi2]

    for i in [1, 2]:
        cells = chain.get_cells(i)
        I = Matrix.identity(cells, ring)
        lhs = gs[i] * fs[i] - I
        rhs = chain.get_bdymap(i+1)*chis[i] + chis[i-1]*chain.get_bdymap(i)
        assert lhs == rhs
        #print(lhs)

        cells = flow.get_critical(i)
        I = Matrix.identity(cells, ring)
        lhs = fs[i] * gs[i]
        #print(lhs)
        #print(I)
        assert lhs == I
        

    #return

    #for idx in [0, 1, 2]:
        #print(fs[idx]*gs[idx])
        #print(gs[idx]*fs[idx])
        #print()
    

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

    test_chain(Ct)

    #C.dump()
    #Ct.dump()
    A = C @ Ct
    #test_chain(A)



def test_main():

    for ring in [element.Q, element.FiniteField(2), element.FiniteField(7)]:
    #for ring in [element.FiniteField(2)]:

        test_matrix(ring)

        cx = Assembly.build_tetrahedron(ring)
        chain = cx.get_chain()
        #chain.dump()
        test_chain(chain)
    
        M = parse("""
        1.1..1
        11.1..
        .11.1.
        ...111
        """)
        chain = Chain.from_array(M, ring)
        test_chain(chain)
    
        #cx = Assembly.build_torus(ring, 2, 2)
        cx = Assembly.build_surface(ring, (0, 0), (2, 2))
        chain = cx.get_chain()
        test_chain(chain)
    



def main():
    ring = element.Q
    #ring = element.FiniteField(2)

    # tetra-hedron face map:
    M = parse("""
    1.1..1
    11.1..
    .11.1.
    ...111
    """)
    chain = Chain.from_array(M, ring)
    chain.check()

    M = chain.get_bdymap(1)
    print(M)
    rows, cols = M.rows, M.cols

    flow = Flow(chain)
    #flow.build()
    
    # match goes from src:grade-1 --> tgt:grade
    # rows are faces, cols are edges
    assert flow.accept(rows[1], cols[0])
    flow.add(rows[1], cols[0])

    assert flow.accept(rows[2], cols[1])
    flow.add(rows[2], cols[1])

    assert flow.accept(rows[3], cols[4])
    flow.add(rows[3], cols[4])

    print(flow)

    A = flow.get_bdymap(2) # 2 --> 1
    B = flow.get_bdymap(1) # 1 --> 0

    print(A, A.shape)
    print(B, B.shape)
    print(B.cols, B.rows)
    C = B*A
    assert C.is_zero()


    f0 = flow.get_f(0)
    f1 = flow.get_f(1)
    f2 = flow.get_f(2)

    assert f0*chain.get_bdymap(1) == B*f1
    assert f1*chain.get_bdymap(2) == A*f2

    g0 = flow.get_g(0)
    g1 = flow.get_g(1)
    g2 = flow.get_g(2)

    assert chain.get_bdymap(1)*g1 == g0*B
    assert chain.get_bdymap(2)*g2 == g1*A

    fs = [f0, f1, f2]
    gs = [g0, g1, g2]

    chi0 = flow.get_homotopy(0)
    chi1 = flow.get_homotopy(1)
    chi2 = flow.get_homotopy(2)

    chis = [chi0, chi1, chi2]

    for i in [1, 2]:
        cells = chain.get_cells(i)
        I = Matrix.identity(cells, ring)
        lhs = gs[i] * fs[i] - I
        rhs = chain.get_bdymap(i+1)*chis[i] + chis[i-1]*chain.get_bdymap(i)
        assert lhs == rhs
        #print(lhs)

        cells = flow.get_critical(i)
        I = Matrix.identity(cells, ring)
        lhs = fs[i] * gs[i]
        #print(lhs)
        #print(I)
        assert lhs == I
        
    print()
    for grade in range(2):
        print("grade =", grade)
        print("f =")
        print(fs[grade])
        print("g =")
        print(gs[grade])
        print("chi =")
        print(chis[grade])
        #print("rows:", chis[grade].rows)
        #print("cols:", chis[grade].cols)
        print("fg =")
        print(fs[grade]*gs[grade])
        print("gf =")
        print(gs[grade]*fs[grade])
        print()

    print(flow.get_bdymap(1))
    #print(chain.get_bdymap(1) * chis[0])


def test_surface():

    ring = element.FiniteField(2)
    one = ring.one

    C0 = [Cell(0, i) for i in [0]]
    C1 = [Cell(1, i) for i in [0,1]]
    C2 = [Cell(2, i) for i in [0]]

    Sx = Matrix(C1, C2, {}, ring)
    Sx[C1[0], C2[0]] = one
    Sz = Matrix(C0, C1, {}, ring)
    Sz[C0[0], C1[1]] = one
    src = Chain({0:C0,1:C1,2:C2}, {1:Sz,2:Sx}, ring, check=True)

    Sx = Matrix(C1, C2, {}, ring)
    Sx[C1[0], C2[0]] = one
    Sx[C1[1], C2[0]] = one
    Sz = Matrix(C0, C1, {}, ring)
    Sz[C0[0], C1[0]] = one
    Sz[C0[0], C1[1]] = one
    tgt = Chain({0:C0,1:C1,2:C2}, {1:Sz,2:Sx}, ring, check=True)

    f0 = Matrix.identity(C0, ring)
    f1 = Matrix.identity(C1, ring)
    f1[C1[1], C1[0]] = one
    f2 = Matrix.identity(C2, ring)
    hom = Hom(src, tgt, {0:f0,1:f1,2:f2}, ring, check=True)

    #tgt, cnot = src.cnot(0, 1)

    # ---- Higgott encoder ----

    C0 = [Cell(0, i) for i in range(2)]
    C1 = [Cell(1, i) for i in range(5)]
    C2 = [Cell(2, i) for i in range(2)]

    Sx = Matrix(C1, C2, {}, ring)
    Sx[C1[1], C2[0]] = one
    Sx[C1[3], C2[1]] = one
    Sz = Matrix(C0, C1, {}, ring)
    Sz[C0[0], C1[0]] = one
    Sz[C0[1], C1[4]] = one
    src = Chain({0:C0,1:C1,2:C2}, {1:Sz,2:Sx}, ring, check=True)
    #src.dump()

    chain = src
    for i, j in [(2, 0), (1, 4), (2, 4), (3, 0), (3, 2), (1, 2)]:
        chain, hom = chain.transvect(i, j)
    #chain.dump()
    


    # -------------------------


    m, n = 3, 2
    
    ambly = Assembly.build_surface(
        ring, (0, 0), (m, n),
        open_top=True, open_bot=True)
    
    chain = ambly.get_chain()
    #chain.dump()

    tgt, hom = chain.transvect(1,2)
    #tgt.dump()

    tgt, hom = tgt.transvect(1,2)
    #tgt.dump()




if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        seed(_seed)

    test_main()
    #main()
    test_surface()

    print("OK")


