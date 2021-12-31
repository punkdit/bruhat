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
        self.one = Space(ring, 1, name="K")  # the underlying field is the tensor unit

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

    # __eq__ is object identity.

    def __getitem__(self, idx):
        if idx<0 or idx >= self.n:
            raise IndexError
        return idx

    # todo: __add__, __mul__ for (bi-)monoidal structure
    # See: chain.Space 

    def __add__(self, other):
        return AddCell0(self.rig, (self, other))

    def __matmul__(self, other):
        return MulCell0(self.rig, (self, other))



class AddCell0(Cell0):

    # XXX use https://docs.python.org/3/library/weakref.html
    cache = {}
    def __new__(cls, rig, _items):
        assert _items
        items = []
        for item in _items:
            assert item.rig == rig
            if type(item) is AddCell0:
                items += item.items # _associative on the nose
            else:
                items.append(item)
        key = tuple(items)
        if key in cls.cache:
            return cls.cache[key]
        cell0 = object.__new__(cls)
        n = sum(item.n for item in items)
        name = "("+"+".join(item.name for item in items)+")"
        Cell0.__init__(cell0, rig, n, name)
        cell0.items = items
        cls.cache[key] = cell0
        return cell0

    def __init__(self, *_items):
        pass



class MulCell0(Cell0):

    # XXX use https://docs.python.org/3/library/weakref.html
    cache = {}
    def __new__(cls, rig, _items):
        assert _items
        items = []
        for item in _items:
            assert item.rig == rig
            if type(item) is AddCell0:
                items += item.items # _associative on the nose
            else:
                items.append(item)
        key = tuple(items)
        if key in cls.cache:
            return cls.cache[key]
        cell0 = object.__new__(cls)
        n = reduce(operator.mul, [item.n for item in items], 1)
        name = "@".join(item.name for item in items)
        Cell0.__init__(cell0, rig, n, name)
        cell0.items = items
        cls.cache[key] = cell0
        return cell0

    def __init__(self, *_items):
        pass



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

    def __eq__(self, other):
        assert self.hom == other.hom
        return numpy.alltrue(self.A == other.A)

    def __mul__(self, other):
        "composition of 1-morphism's"
        assert self.rig == other.rig
        assert self.src == other.tgt
        rig = self.rig
        A = rig.dot(self.A, other.A)
        return Cell1(self.tgt, other.src, A)

    @property
    def dual(self):
        if self._dual is not None:
            return self._dual
        src, tgt = self.hom
        A = [[self[j,i].dual for j in src] for i in tgt]
        self._dual = Cell1(tgt, src, A)
        return self._dual

    @classmethod
    def identity(cls, cell0):
        rig = cell0.rig
        A = numpy.empty((cell0.n, cell0.n), dtype=object)
        for row in cell0:
          for col in cell0:
            A[row, col] = rig.one if row==col else rig.zero
        return Cell1(cell0, cell0, A)

    @classmethod
    def zero(cls, tgt, src): # tgt <---- src
        assert tgt.rig == src.rig
        rig = tgt.rig
        A = numpy.empty((tgt.n, src.n), dtype=object)
        for row in tgt:
          for col in src:
            A[row, col] = rig.zero
        return Cell1(tgt, src, A)

    @classmethod
    def rand(cls, tgt, src, maxdims=4, name="A"): # tgt <---- src
        assert tgt.rig == src.rig
        rig = tgt.rig
        ring = rig.ring
        A = numpy.empty((tgt.n, src.n), dtype=object)
        for row in tgt:
          for col in src:
            n = randint(0, maxdims)
            A[row, col] = Space(ring, n, name="%s_%d%d"%(name, row, col))
        return Cell1(tgt, src, A)

    @staticmethod
    def unfold(i, cell0): # the unit of a higher adjunction
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
    def fold(i, cell0): # the counit of a higher adjunction
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
            assert f.tgt == tgt[i,j], ("%s != %s"%(f.tgt.name, tgt[i,j].name))
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
    def zero(cls, tgt, src): # tgt <---- src
        assert isinstance(tgt, Cell1)
        assert isinstance(src, Cell1)
        assert tgt.rig == src.rig
        assert tgt.hom == src.hom
        rig = tgt.rig
        rows, cols = tgt.shape
        linss = [[Lin.zero(tgt[i,j], src[i,j]) for j in range(cols)] for i in range(rows)]
        return Cell2(tgt, src, linss)

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
        assert self.hom == other.hom
        return numpy.alltrue(self.A == other.A)

    def __mul__(self, other):
        "(vertical) composition of 2-morphism's"
        assert self.rig == other.rig
        assert self.src == other.tgt
        rig = self.rig
        A = self.A * other.A # compose Lin's elementwise
        return Cell2(self.tgt, other.src, A)

    def __lshift__(left, right):
        "horizontal composition of 2-morphism's"
        src = left.src * right.src
        tgt = left.tgt * right.tgt
        lins = []
        for i in range(left.shape[0]):
          row = []
          for k in range(right.shape[1]):
            f = reduce(Lin.direct_sum,
                [left[i,j]@right[j,k] for j in range(left.shape[1])])
            row.append(f)
          lins.append(row)
        return Cell2(tgt, src, lins)

    @classmethod
    def send(cls, cell1, f):
        "apply f component-wise to construct a Cell2"
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

#        print("reassociate")
#        print("A =", A, A.shape)
#        print("B =", B, B.shape)
#        print("C =", C, C.shape)
        assert C.tgt == B.src
        assert B.tgt == A.src
        m, n = A.shape
        p, q = C.shape

        if n*p == 0:
            tgt, src = A*(B*C), (A*B)*C
            #print(tgt)
            #print(src)
            #assert tgt == src
            return Cell2.zero(tgt, src)
            
        rig = A.rig
        ring = rig.ring
        def add(items):
            items = list(items)
            assert len(items)
            if len(items) == 1:
                return items[0]
            return AddSpace(ring, *items)
        #add = lambda items: AddSpace(ring, *list(items))
        #add = lambda items: reduce(operator.add, list(items) or [rig.zero]) # nullitor ?

        ABC = numpy.empty((m, n, p, q), dtype=object)
        for i in range(m):
         for j in range(n):
          for k in range(p):
           for l in range(q):
            ABC[i,j,k,l] = A[i,j]@B[j,k]@C[k,l]

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
            #print(rd.tgt)
            #print(AB_C[i,l])
            assert rd.src == AB_C[i, l], ("%s != %s"%(rd.src, AB_C[i,l]))
            ab_c[i, l] = rd
        src = Cell1(A.tgt, C.src, AB_C)
        cell1 = Cell1(A.tgt, C.src, [[ab_c[i,l].tgt for l in range(q)] for i in range(m)])
        ab_c = Cell2(cell1, src, ab_c)

#        print("ab_c")
#        print(ab_c)

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
#        print("a_bc")
#        print(a_bc)

        lookup = {}
        for k in range(p):
            for j in range(n):
                lookup[(j,k)] = len(lookup)
        perm = tuple(lookup[j,k] for j in range(n) for k in range(p))
#        print("perm:", perm)
        def get_swap(u):
            assert isinstance(u, AddSpace), (u, perm)
#            print("get_swap")
#            print(u.name)
            #print(perm)
            f = u.get_swap(perm)
            return f

        if len(perm) > 1:
            #s = rd.tgt.send(get_swap)
            s = Cell2.send(ab_c.tgt, get_swap)
#            if s.tgt != a_bc.src:
#                print("get_swap:")
#                print(s.homstr())
#                print(s.tgt)
#                print(a_bc.src)
            assert s.tgt == a_bc.src
        else:
            s = Cell2.identity(ab_c.tgt)
            assert s.tgt == a_bc.src
    
        f = a_bc*s*ab_c
        #print(f.homstr())
        assert f.src == src
        assert f.tgt == tgt
        if inverse:
            f = f.transpose2()
        return f


    @staticmethod
    def _XXX_reassociate(V, W, U, inverse=False):
        "V*(W*U) <--- (V*W)*U"
        def right_distributor(u):
            u = AddSpace.promote(u)
            #print("right_distributor:", u)
            fs = [v.right_distributor() for v in u.items]
            f = reduce(Lin.direct_sum, fs)
            return f
        def left_distributor(u):
            u = AddSpace.promote(u)
            fs = [v.left_distributor(inverse=True) for v in u.items]
            f = reduce(Lin.direct_sum, fs)
            return f
        lhs, rhs = (V*W)*U, V*(W*U)
        print("reassociate")
        print("V", V)
        print("W", W)
        print("U", U)
        print("lhs:", lhs)
        print("rhs:", rhs)
        rd = Cell2.send(lhs, right_distributor)
        #print("right_distributor")
        #print(rd.homstr())
        ld = Cell2.send(rhs, left_distributor)
        #print("left_distributor")
        #print(ld.homstr())
    
        lookup = {}
        for k in range(U.shape[0]):
            for j in range(V.shape[1]):
                lookup[(j,k)] = len(lookup)
        perm = tuple(lookup[j,k] for j in range(V.shape[1]) for k in range(U.shape[0]))
        def get_swap(u):
            assert isinstance(u, AddSpace)
            #print("get_swap")
            #print(u.name)
            #print(perm)
            f = u.get_swap(perm)
            return f
        #s = rd.tgt.send(get_swap)
        s = Cell2.send(rd.tgt, get_swap)
        #print("get_swap:")
        #print(s.homstr())
        #print(s.tgt)
        #print(ld.src)
        assert s.tgt == ld.src
    
        f = ld*s*rd
        #print(f.homstr())
        assert f.src == lhs
        assert f.tgt == rhs
        if inverse:
            f = f.transpose2()
        return f

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
        

def test_reassociate():
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

            A = Cell1.rand(m, n, name="A")
            B = Cell1.rand(n, p, name="B")
            C = Cell1.rand(p, q, name="C")
    
            a = Cell2.reassociate(A, B, C)


def test_monoidal():
    class Ring(element.Type):
        one = 1
        zero = 0
    ring = Ring()

    rig = Rig(ring)
    zero = Cell0(rig, 0, "z")
    one = Cell0(rig, 1, "i")
    m = Cell0(rig, 2, "m")
    n = Cell0(rig, 3, "n")

    print(m+n)
    print(m@n)

    F = Cell1.fold(one, m)
    G = Cell1.unfold(one, m)

    print(F)
    print(F*G)
    print(G*F)

    make_frobenius(F)


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



def test():

    ring = element.Z
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





if __name__ == "__main__":

    seed(1)
    
    profile = argv.profile
    fn = argv.next() or "test"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK")
    print()



