#!/usr/bin/env python3

"""
compute homology using algebraic morse theory.

Earlier version: qupy.ldpc.cell
"""

from random import shuffle, seed

import numpy

from bruhat.argv import argv
from bruhat import element
from bruhat.elim import shortstr



class Matrix(object):
    """ 
        Matrix over arbitrary indexing objects (not just numbers).
        rows and cols : list of hashable items
        elements : map (row, col) -> value
    """

    def __init__(self, rows, cols, elements, ring):
        rows = list(rows)
        cols = list(cols)
        self.rows = rows
        self.cols = cols
        self.set_rows = set(rows)
        self.set_cols = set(cols)
        self.elements = dict(elements)
        self.ring = ring
        self.one = ring.one
        self.zero = ring.zero

    def copy(self):
        return Matrix(self.rows, self.cols, self.elements, self.ring)

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

    @classmethod
    def identity(cls, items, ring):
        one = ring.one
        elements = dict(((i, i), one) for i in items)
        A = cls(items, items, elements, ring)
        return A

    @classmethod
    def inclusion(cls, rows, cols, ring):
        A = cls(rows, cols, {}, ring)
        one = ring.one
        for i in cols:
            A[i, i] = one
        return A

    @classmethod
    def retraction(cls, rows, cols, ring):
        A = cls(rows, cols, {}, ring)
        one = ring.one
        for i in rows:
            A[i, i] = one
        return A

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
            #idx = self.row_lookup[row]
            #jdx = self.col_lookup[col]
            assert row in self.set_rows
            assert col in self.set_cols
            value = self.elements.get((row, col), 0)
        return value

    def __add__(self, other):
        assert self.ring == other.ring
        assert self.set_cols == other.set_rows
        A = self.copy()
        for i,j in other.elements.keys():
            A[i, j] += other[i, j]
        return A

    def __mul__(self, other):
        assert self.ring == other.ring
        assert self.set_cols == other.set_rows
        A = Matrix(self.rows, other.cols, {}, self.ring)
        for (i, j) in self.elements.keys():
            for k in other.cols:
                A[i, k] = A[i, k] + self[i, j] * other[j, k]
        return A

    def sum(self):
        return sum(self.elements.values(), self.zero)

    def nonzero(self):
        zero = self.zero
        for v in self.elements.values():
            if v != zero:
                return True
        assert not self.elements
        return False

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



class Cell(object):

    def __init__(self, grade, bdy={}, key=None):
        self.grade = grade
        self.bdy = dict(bdy) # map Cell -> scalar
        self.key = key

    def __repr__(self):
        return "Cell(%d, %s, %s)"%(self.grade, self.bdy, self.key)

    def __str__(self):
        return "Cell(%d, %s)"%(self.grade, self.key)

    def check(self):
        grade = self.grade
        for cell in self.bdy:
            assert cell.grade == grade-1

    def append(self, cell, w):
        assert isinstance(cell, Cell)
        assert cell.grade == self.grade-1
        self.bdy[cell] = w

    def __len__(self):
        return len(self.bdy)

    def __iter__(self):
        return iter(self.bdy)

    def __getitem__(self, cell):
        return self.bdy[cell]

    def __lt__(self, other):
        return self.key < other.key

    def __le__(self, other):
        return self.key <= other.key

    #def get_dual(self):


class Complex(object):
    "chain complex"

    def __init__(self, lookup, ring):
        # lookup: map grade to map of key->cell
        #self.lookup = {0:{}, 1:{}, 2:{}}
        self.lookup = dict(lookup)
        self.ring = ring
        self.zero = ring.zero
        self.one = ring.one

    def set(self, key, cell):
        assert isinstance(cell, Cell)
        lookup = self.lookup
        cells = lookup.setdefault(cell.grade, {})
        assert cells.get(key) is None
        cells[key] = cell
        assert cell.key is None
        cell.key = key

    def mk_vert(self, key):
        cell = Cell(0, {})
        self.set(key, cell)
        return cell

    def mk_edge(self, key, src, tgt):
        one = self.one
        cell = Cell(1, {self[0, src]:-one, self[0, tgt]:one})
        self.set(key, cell)
        return cell

    def __getitem__(self, key):
        grade, key = key
        return self.lookup[grade][key]

    def get_cells(self, grade):
        cells = list(self.lookup[grade].values())
        cells.sort()
        return cells

    def get_bdymap(self, grade): # warning: needs to return unique object
        # grade -> grade-1
        cols = self.get_cells(grade)
        rows = self.get_cells(grade-1)
        A = Matrix(rows, cols, {}, self.ring)
        zero = self.zero
        one = self.one
        for col in cols:
            assert col.grade == grade
            for row in col:
                assert row.grade == grade-1
                assert A[row, col] == zero
                A[row, col] = col[row] # confused?
        return A

    def build_tetrahedron(cx):
        one = cx.one
    
        for idx in range(4):
            cx.mk_vert(idx+1)
    
        cx.mk_edge(1, 3, 1)
        cx.mk_edge(2, 1, 2)
        cx.mk_edge(3, 1, 4)
        cx.mk_edge(4, 3, 2)
        cx.mk_edge(5, 4, 2)
        cx.mk_edge(6, 4, 3)
    
        cx.set(1, Cell(2, {cx[1, 1]:-one, cx[1, 6]:-one, cx[1, 3]:-one}))
        cx.set(2, Cell(2, {cx[1, 1]: one, cx[1, 2]: one, cx[1, 4]:-one}))
        cx.set(3, Cell(2, {cx[1, 2]:-one, cx[1, 3]: one, cx[1, 5]: one}))
        cx.set(4, Cell(2, {cx[1, 4]: one, cx[1, 5]:-one, cx[1, 6]: one}))
    
        return cx

    def build_torus(self, rows, cols):
        one = self.one

        # build verts
        for i in range(rows):
          for j in range(cols):
            cell = Cell(0)
            key = (i, j)
            self.set(key, cell)

        verts = self.lookup[0]

        # build edges
        for i in range(rows):
          for j in range(cols):
            ks = "hv" # horizontal, _vertical edge
            for k in ks:
                bdy = []
                a_vert = verts.get((i, j))
                assert a_vert is not None
                if k == "h":
                    b_vert = verts.get((i, (j+1)%cols))
                else:
                    b_vert = verts.get(((i+1)%rows, j))
                assert b_vert is not None
                cell = Cell(1, {a_vert:-one, b_vert:one})
                key = (i, j, k)
                self.set(key, cell)

        edges = self.lookup[1]

        # build faces
        for i in range(rows):
          for j in range(cols):
            top = edges[(i, j, "h")]
            left = edges[(i, j, "v")]
            bot = edges[((i+1)%rows, j, "h")]
            right = edges[(i, (j+1)%cols, "v")]
            cell = Cell(2, {top:one, left:-one, bot:-one, right:one})
            self.set((i, j), cell)

        return self



class Flow(object):
    "morse matching of a chain complex"

    def __init__(self, cx):
        self.cx = cx
        self.ring = cx.ring
        self.zero = cx.ring.zero
        self.one = cx.ring.one
        # map grade -> list of pairs of cells (a,b) with a.grade==grade and b.grade==grade+1
        self.pairs = {} 

    def add(self, src, tgt):
        assert isinstance(src, Cell)
        assert isinstance(tgt, Cell)
        assert src.grade == tgt.grade-1
        assert src in tgt
        pairs = self.pairs.setdefault(src.grade, [])
        pairs.append((src, tgt))

    def remove(self, src, tgt):
        self.pairs[src.grade].remove((src, tgt))

    def get_pairs(self):
        pairs = []
        cx = self.cx
        edges = cx.get_cells(1)
        faces = cx.get_cells(2)
        for edge in edges:
            for vert in edge:
                pairs.append((vert, edge))
        for face in faces:
            for edge in face:
                pairs.append((edge, face))
        return pairs

    def get_critical(self, grade):
        cx = self.cx
        remove = set()
        for _grade, pairs in self.pairs.items():
          for src, tgt in pairs:
            remove.add(src)
            remove.add(tgt)
        critical = []
        cells = cx.get_cells(grade)
        for cell in cells:
            if cell not in remove:
                critical.append(cell)
        return critical

#    def get_flowmap(self, grade):
#        # grade --> grade+1
#        cx = self.cx
#        zero = self.zero
#        one = self.one
#        bdy = cx.get_bdymap(grade+1) # grade+1 --> grade
#        edges = cx.get_cells(grade)
#        faces = cx.get_cells(grade+1)
#        A = Matrix(faces, edges, {}, self.ring)
#        pairs = self.pairs.setdefault(grade, [])
#        for src, tgt in pairs:
#            assert src.grade == grade
#            assert A[tgt, src] == zero
#            value = bdy[src, tgt]
#            assert value != zero
#            A[tgt, src] = -one/value
#        return A

#    def get_adj(self, grade):
#        "return flow (adjacency) matrix for this grade"
#        cx = self.cx
#        A = cx.get_bdymap(grade) # grade --> grade-1
#        B = self.get_flowmap(grade-1) # grade-1 --> grade
#        C = B*A # grade --> grade
#        for cell in cx.get_cells(grade):
#            C[cell, cell] = 0
#        return C

    def get_adj(self, grade):
        # grade --> grade
        cx = self.cx
        bdy = cx.get_bdymap(grade) # grade --> grade-1
        zero = self.zero
        one = self.one
        edges = cx.get_cells(grade)
        verts = cx.get_cells(grade-1)
        A = Matrix(edges, verts, {}, self.ring) # grade-1 --> grade
        pairs = self.pairs.setdefault(grade-1, [])
        for src, tgt in pairs: # grade-1 --> grade
            assert src.grade == grade-1
            assert A[tgt, src] == zero
            value = bdy[src, tgt]
            assert value != zero
            A[tgt, src] = -one/value
        for src, tgt in pairs: # grade-1 --> grade
            assert bdy[src, tgt] != 0
            bdy[src, tgt] = 0
        return A*bdy

    def accept(self, src, tgt):
        assert isinstance(src, Cell)
        assert isinstance(tgt, Cell)
        assert src.grade == tgt.grade-1
        pairs = self.pairs.setdefault(src.grade, [])
        for pair in pairs:
            if pair[0] == src or pair[1] == tgt:
                return False
        pairs = self.pairs.setdefault(src.grade-1, [])
        for pair in pairs:
            if pair[1] == src:
                return False
        pairs = self.pairs.setdefault(src.grade+1, [])
        for pair in pairs:
            if pair[0] == tgt:
                return False

        cx = self.cx
        zero = cx.zero
        self.add(src, tgt)
        for grade in (1, 2):
            A = self.get_adj(grade)
            #print(A)
            #N = len(A)
            B = A.copy()
            cells = cx.get_cells(grade)
            while B.nonzero():
                #B = numpy.dot(A, B)
                B = A*B
                #for j in range(N):
                for j in cells:
                    if B[j, j] != zero:
                        self.remove(src, tgt)
                        return False
        self.remove(src, tgt)
        return True

#    def __str__(self):
#        items = ['%s->%s' % (src, tgt) for (src, tgt) in self.pairs[0]]
#        items += ['%s->%s' % (src, tgt) for (src, tgt) in self.pairs[1]]
#        s = ", ".join(items)
#        return "Flow(%s)"%(s,)

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

    def build(self):
        pairs = self.get_pairs()
        shuffle(pairs)

        idx = 0
        while idx < len(pairs):
            src, tgt = pairs[idx]
            if self.accept(src, tgt):
                self.add(src, tgt)
            idx += 1

#    def XX_morse_homology(self, grade=1):
#        cx = self.cx
#        cells = {}
#        crit = {}
#        for deg in [grade-1, grade, grade+1]:
#            cells[deg] = cx.get_cells(deg)
#            crit[deg] = self.get_critical(deg)
#
#        chain = []
#        for deg in [grade+1, grade]:
#            bdy = Matrix(crit[deg-1], crit[deg], {}, self.ring)
#            A = self.get_adj(deg)  # deg --> deg-1 --> deg
#            B = cx.get_bdymap(deg) # deg --> deg-1
#            C = Matrix.identity(cells[deg], self.ring)
#            while C.sum():
#                #print()
#                #print(shortstr(C.todense()))
#                D = B*C # deg --> deg-1
#                for a in crit[deg-1]:
#                  for b in crit[deg]:
#                    if D[a, b]:
#                        bdy[a, b] += D[a, b]
#
#                C = A*C
#
#            #print(bdy.todense())
#            chain.append(bdy)
#        return chain

    def get_bdymap(self, grade=1): # grade --> grade-1
        cx = self.cx

        tgt = self.get_critical(grade-1)
        src = self.get_critical(grade)
    
        bdy = Matrix(tgt, src, {}, self.ring)
        A = self.get_adj(grade)  # grade --> grade-1 --> grade
        B = cx.get_bdymap(grade) # grade --> grade-1
        cells = cx.get_cells(grade)
        C = Matrix.identity(cells, self.ring)
        while C.sum():
            #print()
            #print(shortstr(C.todense()))
            D = B*C # grade --> grade-1
            for a in tgt:
              for b in src:
                if D[a, b]:
                    bdy[a, b] += D[a, b]

            C = A*C

        #print(bdy.todense())
        return bdy





def main():
    ring = element.Q

    cx = Complex({}, ring).build_tetrahedron()
    #cx = Complex({}, ring).build_torus(2, 2)
    
    A = cx.get_bdymap(1)
    #print(shortstr(A.todense()))

    B = cx.get_bdymap(2)
    #print(shortstr(B.todense()))

    C = A*B
    assert C.sum() == 0
    #print(shortstr(C.todense()))

    flow = Flow(cx)
    flow.build()

    print(flow)

#    
#    A = flow.get_adj(1)
#    print(A.longstr())
#
#    B = A
#    for i in range(100):
#        #print(B)
#        B = A*B
#        assert B.sum()
#        
#
#    return

    A = flow.get_bdymap(2)
    B = flow.get_bdymap(1)

    print(shortstr(A.todense()))
    print(shortstr(B.todense()))
    C = B*A
    assert C.sum() == 0



if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        seed(_seed)

    main()



