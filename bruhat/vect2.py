#!/usr/bin/env python3

"""
Kapranov-Voevodsky 2-vector spaces.

This is a bicategory with 0,1,2-cells, 
or objects, morphisms, 2-morphisms.
Here we call these types Cell0, Cell1, Cell2.

"""

from functools import reduce
import operator
import itertools
from random import randint, seed

import numpy

from bruhat import element
from bruhat import elim
from bruhat.chain import Space, Lin, AddSpace, MulSpace
from bruhat.argv import argv

def array(A, cols=None):
    #A = numpy.array(A) # not good enough...
    for row in A:
        assert type(row) is list, row
        assert cols is None or cols==len(row)
        cols = len(row)
    rows = len(A)
    A1 = numpy.empty((rows, cols or 0), dtype=object)
    if len(A):
        A1[:] = A
    return A1


class Matrix(object):
    """
    An object numpy array with type'd entries.
    """
    def __init__(self, tp, A, cols=None):
        orig = A
        if type(A) is list:
            A = array(A, cols)
        assert type(A) is numpy.ndarray, type(A)
        assert cols is None or A.shape[1] == cols, (cols, A.shape, orig)
        self.tp = tp
        self.shape = A.shape
        rows, cols = A.shape
        for row in range(rows):
          for col in range(cols):
            assert isinstance(A[row, col], tp)
        self.A = A.copy()
        #self.rows = rows
        #self.cols = cols

    def __getitem__(self, idx):
        return self.A[idx]

    def __eq__(self, other):
        assert self.tp == other.tp
        return numpy.alltrue(self.A == other.A)

    def indexes(self):
        rows, cols = self.shape
        for i in range(rows):
            for j in range(cols):
                yield (i, j)

    def send(self, f):
        rows, cols = self.shape
        A = self.A
        A = [[f(A[i,j]) for j in range(cols)] for i in range(rows)]
        A = array(A, cols)
        return A

#    # not used:
#    def __mul__(self, other):
#        assert self.tp == other.tp
#        assert self.src == other.tgt
#        A = numpy.dot(self.A, other.A)
#        return Matrix(self.tp, A)


class Rig(object):
    def __init__(self, ring):
        self.ring = ring
        self.zero = Space(ring, 0, name="N") # null
        self.zero._dual = self.zero
        self.one = Space(ring, 1, name="K")  # the underlying field is the tensor unit
        self.one._dual = self.one

    def dot(self, A, B):
        assert type(A) is numpy.ndarray
        assert type(B) is numpy.ndarray
        rows = A.shape[0]
        cols = B.shape[1]
        l = A.shape[1]
        assert l==B.shape[0]
        C = numpy.empty((rows, cols), dtype=object)
        for i in range(rows):
          for j in range(cols):
            items = [A[i,k] @ B[k, j] for k in range(l)] or [self.zero]
            C[i,j] = reduce(operator.add, items)
        return C


# ----------------------------------------------------------------
# Here we build the 0,1,2-cells of a bicategory.
# Note use of classmethod/staticmethod: these could also be
# method's of the first argument, but we like seeing the type when calling.
# It helps to keep track of what level (0,1 or 2) we are on.


class Cell0(object):
    """
    The 0-cell's (object's) are determinued by a natural _number dimension.
    """
    def __init__(self, rig, n=0, name="?"):
        # should we use a list (a basis) instead of n ??
        assert n>=0
        self.rig = rig
        self.ring = rig.ring
        self.n = n
        self.name = name
        self.dual = self # apparently we are self-dual ...

    def __str__(self):
        return "%s(%d, %r)"%(self.__class__.__name__, self.n, self.name)
    __repr__ = __str__

    def __hash__(self):
        return id(self)

    def __getitem__(self, idx):
        if idx<0 or idx >= self.n:
            raise IndexError
        return idx

#    # __eq__ is object identity.
#    # __add__, __matmul__ for (bi-)monoidal structure
#    # See: chain.Space 
#    # XXX i think this is silly... should at least make __add__ strict
#
#    def __add__(self, other):
#        return AddCell0(self.rig, (self, other))
#
#    def __matmul__(self, other):
#        return MulCell0(self.rig, (self, other))

    def __add__(self, other):
        assert isinstance(other, Cell0)
        return Cell0(self.rig, self.n+other.n, "(%s+%s)"%(self.name, other.name))

    def __matmul__(self, other):
        assert isinstance(other, Cell0)
        return Cell0(self.rig, self.n*other.n, "%s@%s"%(self.name, other.name))

    def __eq__(self, other):
        assert isinstance(other, Cell0)
        return self.n == other.n


#
#class AddCell0(Cell0):
#
#    # XXX use https://docs.python.org/3/library/weakref.html
#    cache = {}
#    def __new__(cls, rig, _items):
#        assert _items
#        items = []
#        for item in _items:
#            assert item.rig == rig
#            if type(item) is AddCell0:
#                items += item.items # _associative on the nose
#            else:
#                items.append(item)
#        key = tuple(items)
#        if key in cls.cache:
#            return cls.cache[key]
#        cell0 = object.__new__(cls)
#        n = sum(item.n for item in items)
#        name = "("+"+".join(item.name for item in items)+")"
#        Cell0.__init__(cell0, rig, n, name)
#        cell0.items = items
#        cls.cache[key] = cell0
#        return cell0
#
#    def __init__(self, *_items):
#        pass
#
#    def get_swap(self, perm):
#        rig = self.rig
#        perm = tuple(perm)
#        items = self.items
#        assert len(perm) == len(items), ("len(%s)!=%d"%(perm, len(items)))
#        assert set(perm) == set(range(len(items)))
#        tgt = reduce(operator.add, [items[i] for i in perm])
#        N = len(items)
#        rows = []
#        for i in perm:
#          row = []
#          for j in range(N):
#            if i==j:
#                lin = rig.one
#            else:
#                lin = rig.zero
#            row.append(lin)
#          rows.append(row)
#        #rows = [numpy.concatenate(row, axis=1) for row in rows]
#        #A = numpy.concatenate(rows)
#        return Cell1(tgt, self, rows)
#
#
#
#
#class MulCell0(Cell0):
#
#    # XXX use https://docs.python.org/3/library/weakref.html
#    cache = {}
#    def __new__(cls, rig, _items):
#        assert _items
#        items = []
#        for item in _items:
#            assert item.rig == rig
#            if type(item) is AddCell0:
#                items += item.items # _associative on the nose
#            else:
#                items.append(item)
#        key = tuple(items)
#        if key in cls.cache:
#            return cls.cache[key]
#        cell0 = object.__new__(cls)
#        n = reduce(operator.mul, [item.n for item in items], 1)
#        name = "("+"@".join(item.name for item in items)+")"
#        Cell0.__init__(cell0, rig, n, name)
#        cell0.items = items
#        cls.cache[key] = cell0
#        return cell0
#
#    def __init__(self, *_items):
#        pass



class Cell1(Matrix):
    """
    The 1-cell's (morphism's) : a Matrix of Space's
    """
    def __init__(self, tgt, src, A):
        assert isinstance(tgt, Cell0)
        assert isinstance(src, Cell0)
        assert tgt.rig == src.rig
        self.rig = tgt.rig
        self.tgt = tgt
        self.src = src
        self.hom = (tgt, src) # yes it's backwards, just like shape is.
        self._dual = None # cache this
        Matrix.__init__(self, Space, A, src.n)
        assert self.shape == (tgt.n, src.n), ("%s != %s" % (self.shape, (tgt.n, src.n)))

    def __str__(self):
        A = self.A
        rows, cols = A.shape
        s = "["+''.join(
            ["["+' '.join(A[i,j].name for j in range(cols))+"]" for i in range(rows)])+"]"
        return s

    def dimstr(self):
        A = self.A
        rows, cols = A.shape
        s = "["+''.join(
            ["["+' '.join(str(A[i,j].n) for j in range(cols))+"]" for i in range(rows)])+"]"
        return s

    def __repr__(self):
        s = str(self.A)
        s = s.replace("\n", " ")
        return "Cell1(%s, %s, %s)"%(self.tgt, self.src, s)

    @classmethod
    def promote(cls, cell):
        if isinstance(cell, Cell0):
            I = Cell1.identity(cell)
            return I
        if isinstance(cell, Cell1):
            return cell
        assert 0, "what's this: %s"%(type(cell),)

    def __eq__(self, other):
        if isinstance(other, Cell2):
            self = Cell2.promote(self)
            return self.__eq__(other) # <-------- return
        other = Cell1.promote(other)
        assert self.hom == other.hom
        return numpy.alltrue(self.A == other.A)

    def __mul__(self, other):
        "composition of 1-morphism's"
        if isinstance(other, Cell2):
            self = Cell2.promote(self)
            return self.__mul__(other) # <-------- return
        other = Cell1.promote(other)
        assert self.rig == other.rig
        assert self.src == other.tgt
        rig = self.rig
        A = rig.dot(self.A, other.A)
        return Cell1(self.tgt, other.src, A)

    def __add__(self, other):
        "direct sum of 1-morphism's"
        if isinstance(other, Cell2):
            self = Cell2.promote(self)
            return self.__add__(other) # <-------- return
        other = Cell1.promote(other)
        assert self.rig == other.rig
        rig = self.rig
        src = self.src + other.src
        tgt = self.tgt + other.tgt
        A = numpy.empty((tgt.n, src.n), dtype=object)
        A[:, :] = [[rig.zero]]
        A[:self.tgt.n, :self.src.n] = self.A
        A[self.tgt.n:, self.src.n:] = other.A
        return Cell1(tgt, src, A)

    def __matmul__(self, other):
        "tensor (monoidal) product of 1-morphism's"
        if isinstance(other, Cell2):
            self = Cell2.promote(self)
            return self.__matmul__(other) # <-------- return
        other = Cell1.promote(other)
        assert self.rig == other.rig
        src = self.src @ other.src
        tgt = self.tgt @ other.tgt
        A = numpy.empty((tgt.n, src.n), dtype=object)
        m, n = other.shape
        for i, j in numpy.ndindex(self.shape): # big
          for k, l in numpy.ndindex(other.shape): # little
            A[i*m+k, j*n+l ] = self[i,j] @ other[k,l]
        return Cell1(tgt, src, A)

    def __lshift__(self, other):
        if isinstance(other, Cell1):
            return self.__mul__(other) # <-------- return
        self = Cell2.promote(self)
        return self << other

    @property
    def transpose(self):
        A = self.A.transpose()
        tgt, src = self.src, self.tgt
        cell1 = Cell1(tgt, src, A)
        return cell1

    @property
    def dual(self):
        if self._dual is not None:
            return self._dual
        src, tgt = self.hom
        A = [[self[j,i].dual for j in src] for i in tgt]
        self._dual = Cell1(tgt, src, A)
        return self._dual

    @staticmethod
    def identity(cell0):
        rig = cell0.rig
        A = numpy.empty((cell0.n, cell0.n), dtype=object)
        for row in cell0:
          for col in cell0:
            A[row, col] = rig.one if row==col else rig.zero
        return Cell1(cell0, cell0, A)

    @staticmethod
    def zero(tgt, src): # tgt <---- src
        assert tgt.rig == src.rig
        rig = tgt.rig
        A = numpy.empty((tgt.n, src.n), dtype=object)
        for row in tgt:
          for col in src:
            A[row, col] = rig.zero
        return Cell1(tgt, src, A)

    @staticmethod
    def all_ones(tgt, src): # tgt <---- src
        assert tgt.rig == src.rig
        rig = tgt.rig
        A = numpy.empty((tgt.n, src.n), dtype=object)
        for row in tgt:
          for col in src:
            A[row, col] = rig.one
        return Cell1(tgt, src, A)

#    @staticmethod
#    def get_addswap(cell0, perm):
#        assert len(perm) == cell0.n
        
    @property
    def normalized(self):
        rig = self.rig
        N, K = rig.zero, rig.one
        n = Cell2.send(self, lambda space : space.get_normal(N, K, inverse=False, force=False))
        return n.tgt

    def to_normal(self, force=False):
        rig = self.rig
        N, K = rig.zero, rig.one
        n = Cell2.send(self, lambda space : space.get_normal(N, K, False, force))
        return n

    def from_normal(self, force=False):
        rig = self.rig
        N, K = rig.zero, rig.one
        n = Cell2.send(self, lambda space : space.get_normal(N, K, True, force))
        return n

    def get_normal(self, force=False):
        rig = self.rig
        N, K = rig.zero, rig.one
        n = Cell2.send(self, lambda space : space.get_normal(N, K, False, force))
        ni = n.transpose2() # inverse is transpose 
        return n*ni

    @staticmethod
    def rand(tgt, src, mindim=0, maxdim=4, name="A"): # tgt <---- src
        assert tgt.rig == src.rig
        rig = tgt.rig
        ring = rig.ring
        A = numpy.empty((tgt.n, src.n), dtype=object)
        for row in tgt:
          for col in src:
            n = randint(mindim, maxdim)
            A[row, col] = Space(ring, n, name="%s_%d%d"%(name, row, col))
        return Cell1(tgt, src, A)

    @staticmethod
    def fold(i, cell0): # the unit of a higher adjunction
        "cell0.dual @ cell0 <---- i"
        assert isinstance(i, Cell0)
        assert isinstance(cell0, Cell0)
        rig = i.rig
        assert i.n == 1, "not tensor identity"
        tgt = cell0.dual @ cell0
        src = i
        A = numpy.empty((cell0.n, cell0.n), dtype=object)
        for row in cell0:
          for col in cell0:
            A[row, col] = rig.one if row==col else rig.zero
        A.shape = tgt.n, src.n
        return Cell1(tgt, src, A)

    @staticmethod
    def unfold(i, cell0): # the counit of a higher adjunction
        "i <---- cell0 @ cell0.dual"
        assert isinstance(i, Cell0)
        assert isinstance(cell0, Cell0)
        rig = i.rig
        assert i.n == 1, "not tensor identity"
        tgt = i
        src = cell0 @ cell0.dual
        A = numpy.empty((cell0.n, cell0.n), dtype=object)
        for row in cell0:
          for col in cell0:
            A[row, col] = rig.one if row==col else rig.zero
        A.shape = tgt.n, src.n
        return Cell1(tgt, src, A)


class Cell2(Matrix):
    """
    The 2-cell's (2-morphism's) : a Matrix of Lin's
    """
    def __init__(self, tgt, src, linss):
        assert isinstance(tgt, Cell1)
        assert isinstance(src, Cell1)
        assert tgt.rig == src.rig
        assert tgt.hom == src.hom
        self.rig = tgt.rig
        self.tgt = tgt
        self.src = src
        self.hom = (tgt, src) # yes it's backwards, just like shape is.
        Matrix.__init__(self, Lin, linss, tgt.hom[1].n)
        assert self.shape[1] == tgt.hom[1].n
        self.check()

    def check(self):
        rows, cols = self.shape
        tgt, src = self.hom
        #print( self.shape , (tgt.hom[0].n, tgt.hom[1].n) )
        assert self.shape == (tgt.hom[0].n, tgt.hom[1].n)
        for i in range(rows):
          for j in range(cols):
            f = self[i,j]
            assert f.tgt == tgt[i,j], ("%s != %s = tgt[%d,%d]"%(f.tgt.name, tgt[i,j].name, i, j))
            assert f.src == src[i,j], ("%s != %s"%(f.src.name, src[i,j].name))

    def homstr(self):
        A = self.A
        rows, cols = A.shape
        s = "["+''.join(
            ["["+' '.join(A[i,j].homstr()
            for j in range(cols))+"]" 
            for i in range(rows)])+"]"
        return s
    __str__ = homstr

    @classmethod
    def promote(cls, cell):
        if isinstance(cell, Cell0):
            I = Cell1.identity(cell)
            i = Cell2.identity(I)
            return i
        if isinstance(cell, Cell1):
            i = Cell2.identity(cell)
            return i
        if isinstance(cell, Cell2):
            return cell
        assert 0, "what's this: %s"%(cell,)

    @classmethod
    def zero(cls, tgt, src): # tgt <---- src
        assert isinstance(tgt, Cell1)
        assert isinstance(src, Cell1)
        assert tgt.rig == src.rig
        assert tgt.hom == src.hom
        rig = tgt.rig
        rows, cols = tgt.shape
        linss = [[Lin.zero(tgt[i,j], src[i,j]) for j in range(cols)] for i in range(rows)]
        return Cell2(tgt, src, linss)

#    def transpose(self): # what am i trying to do here....
#        A = self.A.transpose()
#        tgt, src = self.src, self.tgt
#        cell2 = Cell2(tgt, src, A)
#        return cell2

    def transpose2(self):
        rows, cols = self.shape
        tgt, src = self.hom
        #A = self.A.transpose()
        A = self.A
        A = [[lin.transpose() for lin in row] for row in A]
        return Cell2(src, tgt, A)

    def __repr__(self):
        s = str(self.A)
        s = s.replace("\n", " ")
        return "Cell2(%s, %s, %s)"%(self.tgt, self.src, s)

    def __eq__(self, other):
        other = Cell2.promote(other)
        assert self.hom == other.hom
        return numpy.alltrue(self.A == other.A)

    def __mul__(self, other):
        "(vertical) composition of 2-morphism's"
        other = Cell2.promote(other)
        assert self.rig == other.rig
        if self.src != other.tgt:
            self = self * self.src.from_normal() # recurse
            other = other.tgt.to_normal() * other # recurse
        assert self.src == other.tgt
        rig = self.rig
        A = self.A * other.A # compose Lin's elementwise
        return Cell2(self.tgt, other.src, A)

    def __lshift__(left, right):
        "horizontal composition of 2-morphism's"
        right = Cell2.promote(right)
        tgt = left.tgt * right.tgt
        src = left.src * right.src
        lins = []
        for i in range(left.shape[0]):
          row = []
          for k in range(right.shape[1]):
            f = reduce(Lin.direct_sum,
                [left[i,j]@right[j,k] for j in range(left.shape[1])])
            row.append(f)
          lins.append(row)
        return Cell2(tgt, src, lins)

    def __add__(left, right):
        "direct_sum of 2-morphisms's"
        right = Cell2.promote(right)
        tgt = left.tgt + right.tgt
        src = left.src + right.src
        rig = left.rig
        zero = Lin(rig.zero, rig.zero)
        A = numpy.empty(tgt.shape, dtype=object)
        A[:] = [[zero]]
        A[:left.shape[0], :left.shape[1]] = left.A
        A[left.shape[0]:, left.shape[1]:] = right.A
        return Cell2(tgt, src, A)

    def __matmul__(self, other):
        "tensor (monoidal) product of 2-morphism's"
        other = Cell2.promote(other)
        src = self.src @ other.src
        tgt = self.tgt @ other.tgt
        A = numpy.empty(tgt.shape, dtype=object)
        m, n = other.shape
        for i, j in numpy.ndindex(self.shape): # big
          for k, l in numpy.ndindex(other.shape): # little
            A[i*m+k, j*n+l ] = self[i,j] @ other[k,l]
        return Cell2(tgt, src, A)

    @property
    def normalized(self):
        self = self.tgt.to_normal() * self * self.src.from_normal()
        return self

    @classmethod
    def send(cls, cell1, f):
        "apply f component-wise to construct a Cell2 tgt or src"
        #print("send", cell1.shape)
        A = Matrix.send(cell1, f)
        rows, cols = cell1.shape
        tgt = [[A[i,j].tgt for j in range(cols)] for i in range(rows)]
        tgt = Cell1(cell1.tgt, cell1.src, tgt)
        src = [[A[i,j].src for j in range(cols)] for i in range(rows)]
        src = Cell1(cell1.tgt, cell1.src, src)
        lin2 = Cell2(tgt, src, A)
        return lin2

    @classmethod
    def identity(cls, cell1):
        tgt, src = cell1, cell1
        rows, cols = cell1.shape
        A = [[cell1[i,j].identity() for j in range(cols)] for i in range(rows)]
        return Cell2(tgt, src, A)

    @classmethod
    def left_unitor(cls, cell1, inverse=False):
        m, n = cell1.hom
        I_m = Cell1.identity(m)
        tgt, src = cell1, I_m * cell1
        if inverse:
            tgt, src = src, tgt
        # Bit of a hack just using .iso here!
        # Should use MulSpace.unitor, etc. etc. XXX
        rows, cols = cell1.shape
        A = [[Lin.iso(tgt[i,j], src[i,j])
            for j in range(cols)] for i in range(rows)]
        return Cell2(tgt, src, A)

    @classmethod
    def right_unitor(cls, cell1, inverse=False):
        m, n = cell1.hom
        I_n = Cell1.identity(n)
        tgt, src = cell1, cell1 * I_n
        if inverse:
            tgt, src = src, tgt
        # Bit of a hack just using .iso here!
        # Should use MulSpace.unitor, etc. etc. XXX
        rows, cols = cell1.shape
        A = [[Lin.iso(tgt[i,j], src[i,j])
            for j in range(cols)] for i in range(rows)]
        return Cell2(tgt, src, A)

    @classmethod
    def rand(cls, tgt, src):
        assert tgt.hom == src.hom
        shape = tgt.shape
        lins = [[Lin.rand(tgt[i,j], src[i,j]) 
            for j in range(shape[1])] for i in range(shape[0])]
        return Cell2(tgt, src, lins)

    @staticmethod
    def reassociate(A, B, C, inverse=False):
        "A*(B*C) <--- (A*B)*C"

        assert C.tgt == B.src
        assert B.tgt == A.src
        m, n = A.shape
        p, q = C.shape

        if n*p == 0:
            tgt, src = A*(B*C), (A*B)*C
            if inverse:
                tgt, src = src, tgt
            return Cell2.zero(tgt, src)
            
        rig = A.rig
        ring = rig.ring
        def add(items):
            items = list(items)
            assert len(items)
            if len(items) == 1:
                return items[0]
            return AddSpace(ring, *items)

        #ABC = numpy.empty((m, n, p, q), dtype=object)
        #for i in range(m):
        # for j in range(n):
        #  for k in range(p):
        #   for l in range(q):
        #    ABC[i,j,k,l] = A[i,j]@B[j,k]@C[k,l]

        # src
        AB_C = numpy.empty((m, q), dtype=object)
        ab_c = numpy.empty((m, q), dtype=object)
        for i in range(m):
          for l in range(q):
            AB_C[i, l] = add(add(A[i,j]@B[j,k] for j in range(n))@C[k,l] for k in range(p))
            items = [Lin.right_distributor(
                add(A[i,j]@B[j,k] for j in range(n)), C[k,l]) for k in range(p)]
            #rd = reduce(Lin.direct_sum, items or [rig.zero.identity()]) # will need nullitor's ?
            rd = reduce(Lin.direct_sum, items)
            assert rd.src == AB_C[i, l], ("%s != %s"%(rd.src, AB_C[i,l]))
            ab_c[i, l] = rd
        src = Cell1(A.tgt, C.src, AB_C)
        cell1 = Cell1(A.tgt, C.src, [[ab_c[i,l].tgt for l in range(q)] for i in range(m)])
        ab_c = Cell2(cell1, src, ab_c)

        # tgt
        A_BC = numpy.empty((m, q), dtype=object)
        a_bc = numpy.empty((m, q), dtype=object)
        for i in range(m):
          for l in range(q):
            A_BC[i, l] = add(A[i,j]@add(B[j,k]@C[k,l] for k in range(p)) for j in range(n))
            items = [Lin.left_distributor(
                A[i,j], add(B[j,k]@C[k,l] for k in range(p)), inverse=True) for j in range(n)]
            #ld = reduce(Lin.direct_sum, items or [rig.zero.identity()]) # will need nullitor's ?
            ld = reduce(Lin.direct_sum, items)
            assert ld.tgt == A_BC[i, l]
            a_bc[i, l] = ld
        tgt = Cell1(A.tgt, C.src, A_BC)
        cell1 = Cell1(A.tgt, C.src, [[a_bc[i,l].src for l in range(q)] for i in range(m)])
        a_bc = Cell2(tgt, cell1, a_bc)

        lookup = {}
        for k in range(p):
            for j in range(n):
                lookup[(j,k)] = len(lookup)
        perm = tuple(lookup[j,k] for j in range(n) for k in range(p))
        def get_swap(u):
            assert isinstance(u, AddSpace), (u, perm)
            f = u.get_swap(perm)
            return f

        if len(perm) > 1:
            s = Cell2.send(ab_c.tgt, get_swap)
            assert s.tgt == a_bc.src
        else:
            s = Cell2.identity(ab_c.tgt)
            assert s.tgt == a_bc.src
    
        f = a_bc*s*ab_c
        assert f.src == src
        assert f.tgt == tgt
        if inverse:
            f = f.transpose2()
        return f

    @staticmethod
    def interchanger(A, B, C, D, inverse=False):
        """ 
            (A @ C) << (B @ D) <------- (A << B) @ (C << D) 
        """
        assert isinstance(A, Cell1)
        assert isinstance(B, Cell1)
        assert isinstance(C, Cell1)
        assert isinstance(D, Cell1)
        rig = A.rig
        N, K = rig.zero, rig.one
        tgt = (A @ C) << (B @ D)
        src = (A << B) @ (C << D)
        linss = []
        for i in src.hom[0]:
          lins = []
          for j in src.hom[1]:
            s, t = src[i,j], tgt[i,j]
            #print(s.name)
            #print(t.name)
            #print()
            f = s.get_normal() # XXX too much in general
            s = f.tgt
            assert isinstance(f.tgt, AddSpace)
            gs = [space.get_swap((0,2,1,3)) for space in f.tgt.items]
            g = reduce(Lin.direct_sum, gs)
            lin = g*f
            #lin = Lin.zero(t, s)
            lins.append(lin)
          linss.append(lins)
        mu = Cell2(tgt, src, linss)
        if inverse:
            mu = mu.transpose2()
        return mu

    @staticmethod
    def unit(A):
        assert isinstance(A, Cell1)
        src = Cell1.identity(A.src)
        tgt = A.dual * A
        assert tgt.hom == src.hom
        shape = tgt.shape
        n = shape[0]
        assert n == shape[1]
        rig = A.rig
        ring = rig.ring
        linss = []
        for i in range(n):
          lins = [] # row
          for k in range(n):
            t, s = tgt[i,k], src[i,k]
            if i!=k:
                lin = Lin.zero(t, s)
            else:
                # argh... why can't we direct_sum the Lin.unit's ?
                a = elim.zeros(ring, t.n, s.n)
                idx = 0
                for j in range(A.shape[0]):
                    unit = Lin.unit(rig.one, A[j,i])
                    a[idx:idx+unit.shape[0], :] = unit.A
                    idx += unit.shape[0]
                assert idx == t.n
                lin = Lin(t, s, a)
            lins.append(lin)
          linss.append(lins)
        return Cell2(tgt, src, linss)
        
    @staticmethod
    def counit(A):
        assert isinstance(A, Cell1)
        src = A * A.dual
        tgt = Cell1.identity(A.tgt)
        assert tgt.hom == src.hom
        shape = tgt.shape
        n = shape[0]
        assert n == shape[1]
        rig = A.rig
        ring = rig.ring
        linss = []
        for i in range(n):
          lins = [] # row
          for k in range(n):
            t, s = tgt[i,k], src[i,k]
            if i!=k:
                lin = Lin.zero(t, s)
            else:
                # argh... why can't we direct_sum the Lin.counit's ?
                a = elim.zeros(ring, t.n, s.n)
                idx = 0
                for j in range(A.shape[1]):
                    counit = Lin.counit(rig.one, A[i,j])
                    a[:, idx:idx+counit.shape[1]] = counit.A
                    idx += counit.shape[1]
                assert idx == s.n
                lin = Lin(t, s, a)
            lins.append(lin)
          linss.append(lins)
        return Cell2(tgt, src, linss)

    @staticmethod
    def fold(one, A):
        assert isinstance(A, Cell1)
        m, n = A.hom
        Im = Cell1.identity(m)
        In = Cell1.identity(n)
        Fm = Cell1.fold(one, m)
        Fn = Cell1.fold(one, n)
        tgt = (Im @ A.transpose) * Fm
        src = (A @ In) * Fn
        tgt = tgt.get_normal().tgt
        src = src.get_normal().tgt
        assert tgt == src
        i = Cell2.identity(src)
        return i

    @staticmethod
    def unfold(one, A):
        assert isinstance(A, Cell1)
        m, n = A.hom
        Im = Cell1.identity(m)
        In = Cell1.identity(n)
        Gm = Cell1.unfold(one, m)
        Gn = Cell1.unfold(one, n)
        tgt = Gn * (In @ A.transpose)
        src = Gm * (A @ Im)
        tgt = tgt.get_normal().tgt
        src = src.get_normal().tgt
        assert tgt == src
        i = Cell2.identity(src)
        return i

    @staticmethod
    def fold_counit(one, A): # fold the cap
        assert isinstance(A, Cell1)
        m, n = A.hom
        Im = Cell1.identity(m)
        In = Cell1.identity(n) # not needed
        At = A.transpose
        Atd = At.dual
        Fm = Cell1.fold(one, m)
        Fn = Cell1.fold(one, n) # not needed
        fA = Cell2.fold(one, A)
        assert fA.src == ((A @ In) << Fn).normalized
        cap = Cell2.counit(Atd)
        assert cap.src == Atd * At
        assert cap.tgt == Im
        bot = (Im @ Atd) << fA
        lhs, rhs = bot.src.normalized , ((A @ Atd)<<Fn).normalized
        #assert ((Im @ Atd) << (A @ In)).normalized == ((A @ Im) << (In @ Atd)).normalized # no
        #assert (A @ Atd).normalized == ((A @ Im) << (In @ Atd)).normalized # yes
            #(A @ C) << (B @ D) <------- (A << B) @ (C << D) 
        rhs = (Im @ Atd) << (A @ In) << Fn
        rhs = rhs.normalized
        assert lhs == rhs
        top = (Im @ cap) << Fm
        assert (bot.tgt.normalized) == (top.src.normalized)
        #mu = Cell2.interchanger(Im, A, Atd, In, inverse=True)
        mu = Cell2.interchanger(Im, A, Atd, In) << Fn
        #print(mu.tgt)
        #print(bot.src)
        left = (top.normalized * bot.normalized * mu.normalized).normalized
        return left

    @staticmethod
    def unfold_counit(one, A): # unfold the cap
        assert isinstance(A, Cell1)
        m, n = A.hom
        In = Cell1.identity(n)
        At = A.transpose
        Atd = At.dual
        Gn = Cell1.unfold(one, n)
        gA = Cell2.unfold(one, A)
        cap = Cell2.counit(At)
        assert cap.src == At * Atd
        assert cap.tgt == In
        bot = gA << (In @ Atd)
        top = Gn << (In @ cap)
        assert (bot.tgt.normalized) == (top.src.normalized)
        left = (top * bot).normalized
        return left

        

class Frobenius(object):
    def __init__(self, **kw):
        self.__dict__.update(kw)


def make_frobenius(A):
    assert isinstance(A, Cell1)

    m, n = A.hom
    I_m = Cell1.identity(m)
    I_n = Cell1.identity(n)
    l_cup = Cell2.unit(A)
    l_cap = Cell2.counit(A)
    r_cup = Cell2.unit(A.dual)
    r_cap = Cell2.counit(A.dual)
    i_A = Cell2.identity(A)
    l_A = Cell2.left_unitor(A)
    r_A = Cell2.right_unitor(A)
    i_dA = Cell2.identity(A.dual)
    l_dA = Cell2.left_unitor(A.dual)
    r_dA = Cell2.right_unitor(A.dual)

    # getting down to business
    X = A * A.dual # the 'object' supporting the algebra
    i_X = Cell2.identity(X)
    unit = r_cup
    counit = l_cap
    assert unit.tgt == X
    assert counit.src == X

    left_unitor = l_A << i_dA
    left_unitor = left_unitor * Cell2.reassociate(I_m, A, A.dual, inverse=True)
    assert left_unitor.tgt == X
    assert left_unitor.src == I_m * X

    right_unitor = i_A << r_dA
    right_unitor = right_unitor * Cell2.reassociate(A, A.dual, I_m)
    assert right_unitor.tgt == X
    assert right_unitor.src == X * I_m

    mul = i_A << r_cap << i_dA
    assert mul.tgt == (A * I_n) * A.dual
    mul = (r_A << i_dA)*mul
    assert mul.tgt == X
    assert mul.src == A * (A.dual * A) * A.dual
    a = Cell2.reassociate(A, A.dual, A) << i_dA
    mul = mul*a
    assert mul.src == ((A * A.dual) * A) * A.dual
    a = Cell2.reassociate(A*A.dual, A, A.dual, inverse=True)
    mul = mul*a
    assert mul.src == X*X

    comul = i_A << l_cup << i_dA
    comul = comul * (Cell2.right_unitor(A, inverse=True) << i_dA)
    a = Cell2.reassociate(A, A.dual, A, inverse=True) << i_dA
    comul = a * comul
    a = Cell2.reassociate(A*A.dual, A, A.dual)
    comul = a * comul
    assert comul.src == X
    assert comul.tgt == X*X

    a = Cell2.reassociate(X, X, X)                # X*(X*X) <--- (X*X)*X 
    ai = Cell2.reassociate(X, X, X, inverse=True) # (X*X)*X <--- X*(X*X)

    # Verify Frobenius algebra axioms
    # (1) monoid axioms
    assert mul * (i_X << unit) == right_unitor
    assert mul * (unit << i_X) == left_unitor
    assert mul * (i_X << mul) * a == mul * (mul << i_X)

    # (2) comonoid axioms
    assert right_unitor * (i_X << counit) * comul == i_X
    assert left_unitor * (counit << i_X) * comul == i_X
    assert (i_X << comul) * comul == a * (comul << i_X) * comul

    # (3) Frobenius law
    lhs = (mul << i_X) * ai * (i_X << comul)
    rhs = (i_X << mul) * a * (comul << i_X)
    assert lhs == rhs

    # a consequence of the Frobenius law
    mid = comul * mul
    assert lhs == mid

    s = mul * comul
    special = s == i_X

    frobenius = Frobenius(**locals())

    return frobenius


def test_monoidal():

    class Ring(element.Type):
        one = 1
        zero = 0
        @classmethod
        def promote(cls, a):
            return a
    ring = Ring()
    rig = Rig(ring)

    zero = Cell0(rig, 0, "z")
    one = Cell0(rig, 1, "i")
    m = Cell0(rig, 3, "m")
    n = Cell0(rig, 2, "n")
    Im = Cell1.identity(m)
    In = Cell1.identity(n)

    A = Cell1.rand(m, n, mindim=2, maxdim=2, name="A")
    B = Cell1.rand(m, n, mindim=2, maxdim=2, name="B")
    C = Cell1.rand(m, n, mindim=2, maxdim=2, name="C")
    D = Cell1.rand(m, n, mindim=2, maxdim=2, name="D")

    a = Cell2.rand(B, A)
    b = Cell2.rand(B, A)

    AB = (A+B)
    AB = (A @ B)
    ab = (a+b)
    ab = (a @ b)

    F = Cell1.unfold(one, m)
    G = Cell1.fold(one, m)

    assert F.dual == G

    frobenius = make_frobenius(F)
    assert frobenius.special

    lhs = (A @ Im) << (In @ B)
    mid = (A @ B)
    rhs = (Im @ B) << (A @ In)
    lhs = lhs.normalized
    rhs = rhs.normalized
    mid = mid.normalized
    assert(lhs == mid) # just happens to work
    assert(rhs != mid)
    assert(lhs != rhs)

    l, m, n, o, p, q = [Cell0(rig, d, name) for (d, name) in zip([2,4,3,2,2,3], 'lmnopq')]
    A = Cell1.rand(l, m, mindim=2, maxdim=2, name="A")
    B = Cell1.rand(m, n, mindim=2, maxdim=2, name="B")
    C = Cell1.rand(o, p, mindim=2, maxdim=2, name="C")
    D = Cell1.rand(p, q, mindim=2, maxdim=2, name="D")

    # interchanger
    tgt = (A @ C) << (B @ D)
    src = (A << B) @ (C << D)

    mu = Cell2.interchanger(A, B, C, D)
    assert mu.src == src
    assert mu.tgt == tgt


def test_hopf():

    class Ring(element.Type):
        one = 1
        zero = 0
        @classmethod
        def promote(cls, a):
            return a
    ring = Ring()

    #ring = element.Q # slooooow
    rig = Rig(ring)
    zero = Cell0(rig, 0, "z")
    one = Cell0(rig, 1, "i")

    if argv.trivial:
        m = Cell0(rig, 3, "m")
        n = m
        A = Cell1.identity(m)
        B = Cell1.identity(m)
    else:
        m = Cell0(rig, 2, "m")
        n = Cell0(rig, 2, "n")
        A = Cell1.rand(m, n, mindim=0, maxdim=1, name="A")

    Im = Cell1.identity(m)
    In = Cell1.identity(n)
    Ad = A.dual
    At = A.transpose
    Atd = At.dual

    Gm = Cell1.unfold(one, m)
    Fm = Cell1.fold(one, m)
    assert Gm.dual == Fm

    Gn = Cell1.unfold(one, n)
    Fn = Cell1.fold(one, n)
    assert Gn.dual == Fn

    f = Cell2.fold(one, A)
    g = Cell2.unfold(one, A)

    left = Cell2.unfold_counit(one, A)
    assert left.src == (Gm << (A @ Atd)).normalized

    right = Cell2.fold_counit(one, Ad)
    lhs, rhs = right.src , ((Ad @ At) << Fm).normalized
    assert lhs == rhs

    right = Cell2.fold_counit(one, Ad)
    lhs, rhs = right.src , ((Ad @ At) << Fm).normalized
    assert lhs == rhs

    stem = left << right
    lhs = stem.src
    left, right = (Gm << (A @ Atd)), ((Ad @ At) << Fm)
    rhs = left << right

    assert lhs.normalized == rhs.normalized

    cap = Cell2.counit(Gn).normalized

    cap = cap * stem


    if 0:

        # only works when the saddle is invertible

        frobenius = make_frobenius(Gm << (A @ Atd))
    
        lhs = cap * frobenius.mul 
        rhs = cap << cap
        print( lhs.normalized[0,0] )
        print( rhs.normalized[0,0] )
        assert lhs.normalized == rhs.normalized
    

def test_frobenius():

    # ----------------------------------------------------
    # Here we construct a Frobenius algebra object 
    # See: https://math.ucr.edu/home/baez/week174.html

    #ring = element.Z # too slow
    class Ring(element.Type):
        one = 1
        zero = 0
    ring = Ring()

    rig = Rig(ring)
    m = Cell0(rig, 2, "m")
    n = Cell0(rig, 2, "n")
    #A = Cell1.rand(m, n, name="A") # too slow
    S = lambda n, name : Space(ring, n, 0, name)
    A = Cell1(m, n, [[S(2, "A0"), S(1, "A1")], [S(0, "A2"), S(1, "A3")]])
    #print(A.dimstr())
    make_frobenius(A)


def test_reassociate():
    #ring = element.Z # too slow
    class Ring(element.Type):
        one = 1
        zero = 0
    ring = Ring()

    rig = Rig(ring)
    zero = Cell0(rig, 0, "z")
    one = Cell0(rig, 1, "i")
    m = Cell0(rig, 2, "m")
    n = Cell0(rig, 3, "n")

    N, K = rig.zero, rig.one
    cell0s = [zero, one, one+one]
    spaces = itertools.cycle([N, K, N+K, K@K, K@K@K, N@N@K@(K+K), K+K])
    for q in cell0s:
     for p in cell0s:
      for n in cell0s:
       for m in cell0s:
        for trial in range(10):
            A = Cell1(m, n, [[spaces.__next__() for j in n] for i in m])
            B = Cell1(n, p, [[spaces.__next__() for k in p] for j in n])
            C = Cell1(p, q, [[spaces.__next__() for l in q] for k in p])
            a = Cell2.reassociate(A, B, C)
            ai = Cell2.reassociate(A, B, C, inverse=True)
            left = a*ai
            assert left == Cell2.identity(A*(B*C))
            right = ai*a
            assert right == Cell2.identity((A*B)*C)

            A = Cell1.rand(m, n, name="A")
            B = Cell1.rand(n, p, name="B")
            C = Cell1.rand(p, q, name="C")
    
            a = Cell2.reassociate(A, B, C)


def test_bialgebra():
    class Ring(element.Type):
        one = 1
        zero = 0
    ring = Ring()

    rig = Rig(ring)
    zero = Cell0(rig, 0, "0")
    m = Cell0(rig, 1, "1")

    I = Cell1.all_ones(m, m)
    Null = Cell1.zero(zero, zero)
    NAdd = Cell1.zero(zero, zero+zero)
    NCopy = Cell1.zero(zero+zero, zero)
    Add = Cell1.all_ones(m, m+m)
    Zero = Cell1.zero(m, zero)
    Copy = Cell1.all_ones(m+m, m)
    Delete = Cell1.zero(zero, m)

    # assoc
    lhs = Add * (Add + I)
    rhs = Add * (I + Add)
    assert lhs != rhs # not strict

    lhs = (Copy + I) * Copy
    rhs = (I + Copy) * Copy
    assert lhs != rhs # not strict

    print( Add * (Zero + I), I)
    #assert Add * (Zero + I) == I

    lhs = (Delete * Add)
    rhs = NAdd*(Delete + Delete)
    assert lhs == rhs

    return

    # --------------------
    # bialgebrator

    lhs = Copy * Add
    print(lhs)

    CC = Copy + Copy
    S = CC.tgt.get_swap((0, 2, 1, 3))

    AA =Add + Add
    rhs = AA * S * CC
    print(rhs)

    #print((Delete*Zero).hom)

    


def test():

    #ring = element.Z # too slow
    class Ring(element.Type):
        one = 1
        zero = 0
        @classmethod
        def promote(cls, a):
            return a
    ring = Ring()

    rig = Rig(ring)

    l, m, n, o, p = 2, 3, 2, 2, 2

    l = Cell0(rig, l, "l")
    m = Cell0(rig, m, "m")
    n = Cell0(rig, n, "n")
    o = Cell0(rig, o, "o")
    p = Cell0(rig, p, "p")

    I_l = Cell1.identity(l)
    I_m = Cell1.identity(m)
    I_n = Cell1.identity(n)
    I_o = Cell1.identity(o)
    I_p = Cell1.identity(p)

    #     A      B      C      D
    # l <--- m <--- n <--- o <--- p
    A = Cell1.rand(l, m, name="A")
    B = Cell1.rand(m, n, name="B")
    C = Cell1.rand(n, o, name="C")
    D = Cell1.rand(o, p, name="D")

    assert A == A
    assert A*I_m == A*I_m
    assert A != A*I_m

    # Does not hold strictly, only up to 2-cell 
    #assert (A*I_m)*B == A*(I_m*B) # nope

    i_A = Cell2.identity(A)
    i_B = Cell2.identity(B)
    i_C = Cell2.identity(C)
    i_D = Cell2.identity(D)
    assert i_A * i_A == i_A

    tgt, src = A*(B*C), (A*B)*C
    a = Cell2.reassociate(A, B, C)  #  A*(B*C) <--- (A*B)*C
    assert a.tgt == tgt
    assert a.src == src

    if argv.reassociate:
        b = Cell2.reassociate(A, B, C, inverse=True)
        assert a*b == Cell2.identity(A*(B*C))
        assert b*a == Cell2.identity((A*B)*C)

    ru_A = Cell2.right_unitor(A)
    assert ru_A.src == A * I_m
    assert ru_A.tgt == A

    lu_B = Cell2.left_unitor(B)
    assert lu_B.src == I_m * B
    assert lu_B.tgt == B

    a = Cell2.reassociate(A, I_m, B)

    # check the triangle equation for identity 1-cell's
    lhs = ru_A << i_B
    rhs = (i_A << lu_B) * a
    assert lhs == rhs

    a = Cell2.reassociate(A, B, C)
    b = Cell2.reassociate(A, B, C, inverse=True)
    lhs = b*a
    assert lhs == Cell2.identity((A*B)*C)

    if argv.pentagon:
        # check the pentagon equation (slow !)
        a = Cell2.reassociate(A*B, C, D)
        b = Cell2.reassociate(A, B, C*D)
        lhs = b*a
    
        a = (i_A << Cell2.reassociate(B, C, D))
        b = Cell2.reassociate(A, B*C, D)
        c = (Cell2.reassociate(A, B, C) << i_D)
        rhs = a*b*c
        assert lhs==rhs


    # ----------------------------------------------------
    # For object's (Cell0's), m & n
    # the 1- and 2- cell's form a category on the nose

    m, n = Cell0(rig, 2, "m"), Cell0(rig, 3, "n")

    A = Cell1.rand(m, n, name="A")
    B = Cell1.rand(m, n, name="B")
    C = Cell1.rand(m, n, name="C")
    D = Cell1.rand(m, n, name="D")

    f = Cell2.rand(B, A)
    g = Cell2.rand(C, B)
    h = Cell2.rand(D, C)

    # identity
    assert (Cell2.identity(B) * f) == f
    assert (f * Cell2.identity(A)) == f

    # _assoc
    assert (h*g)*f == h*(g*f)

    # ----------------------------------------------------
    # snake equations (zig-zag) for adjoint 1-cell's

    cup = Cell2.unit(A)
    cap = Cell2.counit(A)

    i_A = Cell2.identity(A)
    l_A = Cell2.left_unitor(A)
    r_A = Cell2.right_unitor(A)

    a, b, c = l_A , (cap << i_A) , (i_A << cup)
    assoc = Cell2.reassociate(A, A.dual, A, inverse=True)
    lhs = l_A * (cap << i_A) * assoc * (i_A << cup)
    rhs = r_A
    assert lhs == rhs

    i_dA = Cell2.identity(A.dual)
    r_dA = Cell2.right_unitor(A.dual)
    l_dA = Cell2.left_unitor(A.dual)

    assoc = Cell2.reassociate(A.dual, A, A.dual)
    
    lhs = r_dA * (i_dA << cap) * assoc * (cup << i_dA)
    rhs = l_dA
    assert lhs == rhs

    # ----------------------------------------------------
    # _interchange law for 2-morphisms

    #     A      B    
    # l <--- m <--- n 
    A = Cell1.rand(l, m, mindim=1, maxdim=2, name="A")
    B = Cell1.rand(m, n, mindim=1, maxdim=2, name="B")

    a = Cell2.rand(A, A)
    c = Cell2.rand(A, A)
    b = Cell2.rand(B, B)
    d = Cell2.rand(B, B)

    lhs = (c << d) * (a << b)
    rhs = (c*a) << (d*b)

    assert lhs == rhs
    


def test_all():
    test_monoidal()
    test_frobenius()
    test_reassociate()
    test_bialgebra()
    test()



if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)
    
    profile = argv.profile
    fn = argv.next() or "test_all"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK")
    print()



