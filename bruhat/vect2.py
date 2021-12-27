#!/usr/bin/env python3

"""
Kapranov-Voevodsky 2-vector spaces.

This is a bicategory with 0,1,2-cells, 
or objects, morphisms, 2-morphisms.
Here we call these types Cell0, Cell1, Cell2.

"""

from functools import reduce
import operator
from random import randint, seed

import numpy

from bruhat import element
from bruhat.chain import Space, Lin, AddSpace, MulSpace
from bruhat.argv import argv

def array(A):
    #A = numpy.array(A) # not good enough...
    cols = None
    for row in A:
        assert type(row) is list, row
        assert cols is None or cols==len(row)
        cols = len(row)
    rows = len(A)
    A1 = numpy.empty((rows, cols or 0), dtype=object)
    A1[:] = A
    return A1


class Matrix(object):
    """
    An object numpy array with type'd entries.
    """
    def __init__(self, tp, A):
        if type(A) is list:
            A = array(A)
        assert type(A) is numpy.ndarray, type(A)
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
        A = array(A)
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


class Cell0(object):
    """
    The 0-cell's (object's) are determinued by a natural _number dimension.
    """
    def __init__(self, rig, n=0, name="?"):
        assert n>=0
        self.rig = rig
        self.ring = rig.ring
        self.n = n
        self.name = name

    def __str__(self):
        return "%s(%d, %r)"%(self.__class__.__name__, self.n, self.name)
    __repr__ = __str__

    def __hash__(self):
        return id(self)

    # __eq__ is object identity.

    # todo: __add__, __mul__ for (bi-)monoidal structure
    # See: chain.Space 

    def __getitem__(self, idx):
        if idx<0 or idx >= self.n:
            raise IndexError
        return idx

    def __add__(self, other):
        return AddCell0(self, other)

    def __mul__(self, other):
        return MulCell0(self, other)

    # should we use a list (a basis) instead of n ??

    def identity(self):
        rig = self.rig
        A = numpy.empty((self.n, self.n), dtype=object)
        for row in self:
          for col in self:
            A[row, col] = rig.one if row==col else rig.zero
        return Cell1(self, self, A)

    def zero(tgt, src): # tgt <---- src
        assert tgt.rig == src.rig
        rig = tgt.rig
        A = numpy.empty((tgt.n, src.n), dtype=object)
        for row in tgt:
          for col in src:
            A[row, col] = rig.zero
        return Cell1(tgt, src, A)

    def rand(tgt, src, maxdims=4, name="A"): # tgt <---- src
        assert tgt.rig == src.rig
        rig = tgt.rig
        ring = rig.ring
        A = numpy.empty((tgt.n, src.n), dtype=object)
        for row in tgt:
          for col in src:
            n = randint(0, maxdims)
            A[row, col] = Space(ring, n, name="%s_%d%d"%(name, row, col))
        return Cell1(tgt, src, A)


class AddCell0(Cell0):

    cache = {} # XXX use https://docs.python.org/3/library/weakref.html
    def __new__(cls, *_items):
        assert _items
        #if len(_items)==1:
        #    return _items[0] # ???
        items = []
        for item in _items:
            if type(item) is AddCell0:
                items += item.items # _associative on the nose
            else:
                items.append(item)
        key = tuple(items)
        if key in cls.cache:
            #print("cache hit", key)
            return cls.cache[key]
        #print("cache miss", key)
        lin0 = object.__new__(cls)
        rig = items[0].rig
        n = sum(item.n for item in items)
        name = "("+"+".join(item.name for item in items)+")"
        Cell0.__init__(lin0, rig, n, name)
        lin0.items = items
        cls.cache[key] = lin0
        return lin0

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
        Matrix.__init__(self, Space, A)

    def __str__(self):
        A = self.A
        rows, cols = A.shape
        s = "["+''.join(
            ["["+' '.join(A[i,j].name for j in range(cols))+"]" for i in range(rows)])+"]"
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

    # are we doing classmethod's of Cell1, or method's of Cell0 ???
    @classmethod
    def rand(cls, tgt, src, *args, **kw):
        return tgt.rand(src, *args, **kw)

    def send(self, f):
        "apply f component-wise to construct a Cell2"
        A = Matrix.send(self, f)
        rows, cols = self.shape
        tgt = [[A[i,j].tgt for j in range(cols)] for i in range(rows)]
        tgt = Cell1(self.tgt, self.src, tgt)
        src = [[A[i,j].src for j in range(cols)] for i in range(rows)]
        src = Cell1(self.tgt, self.src, src)
        lin2 = Cell2(tgt, src, A)
        return lin2

    def identity(self):
        tgt, src = self, self
        rows, cols = self.shape
        A = [[self[i,j].identity() for j in range(cols)] for i in range(rows)]
        return Cell2(tgt, src, A)

    def left_unitor(self, inverse=False):
        m, n = self.hom
        i_m = m.identity()
        tgt, src = self, i_m * self
        if inverse:
            tgt, src = src, tgt
        # Bit of a hack just using .iso here!
        # Should use MulSpace.unitor, etc. etc. XXX
        rows, cols = self.shape
        A = [[Lin.iso(tgt[i,j], src[i,j])
            for j in range(cols)] for i in range(rows)]
        return Cell2(tgt, src, A)

    def right_unitor(self, inverse=False):
        m, n = self.hom
        i_n = n.identity()
        tgt, src = self, self * i_n
        if inverse:
            tgt, src = src, tgt
        # Bit of a hack just using .iso here!
        # Should use MulSpace.unitor, etc. etc. XXX
        rows, cols = self.shape
        A = [[Lin.iso(tgt[i,j], src[i,j])
            for j in range(cols)] for i in range(rows)]
        return Cell2(tgt, src, A)


class Cell2(Matrix):
    """
    The 2-cell's (2-morphism's) : a Matrix of Lin's
    """
    def __init__(self, tgt, src, A):
        assert isinstance(tgt, Cell1)
        assert isinstance(src, Cell1)
        assert tgt.rig == src.rig
        assert tgt.hom == src.hom
        self.rig = tgt.rig
        self.tgt = tgt
        self.src = src
        self.hom = (tgt, src) # yes it's backwards, just like shape is.
        Matrix.__init__(self, Lin, A)
        self.check()

    def check(self):
        rows, cols = self.shape
        tgt, src = self.hom
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
    def rand(cls, tgt, src):
        assert tgt.hom == src.hom
        shape = tgt.shape
        lins = [[Lin.rand(tgt[i,j], src[i,j]) 
            for j in range(shape[1])] for i in range(shape[0])]
        return Cell2(tgt, src, lins)

    @staticmethod
    def reassociate(V, W, U):
        "V*(W*U) <--- (V*W)*U"
        def right_distributor(u):
            assert isinstance(u, AddSpace)
            fs = [v.right_distributor() for v in u.items]
            f = reduce(Lin.direct_sum, fs)
            return f
        def left_distributor(u):
            assert isinstance(u, AddSpace)
            fs = [v.left_distributor(inverse=True) for v in u.items]
            f = reduce(Lin.direct_sum, fs)
            return f
        lhs, rhs = (V*W)*U, V*(W*U)
        rd = lhs.send(right_distributor)
        #print("right_distributor")
        #print(rd.homstr())
        ld = rhs.send(left_distributor)
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
        s = rd.tgt.send(get_swap)
        #print("get_swap:")
        #print(s.homstr())
        #print(s.tgt)
        #print(ld.src)
        assert s.tgt == ld.src
    
        f = ld*s*rd
        #print(f.homstr())
        assert f.src == lhs
        assert f.tgt == rhs
        return f


def test():

    ring = element.Z
    rig = Rig(ring)

    l, m, n, o, p = 1, 3, 2, 2, 2

    l = Cell0(rig, l, "l")
    m = Cell0(rig, m, "m")
    n = Cell0(rig, n, "n")
    o = Cell0(rig, o, "o")
    p = Cell0(rig, p, "p")

    I_l = l.identity()
    I_m = m.identity()
    I_n = n.identity()
    I_o = o.identity()
    I_p = p.identity()

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

    i_A = A.identity() # Cell2
    i_B = B.identity() # Cell2
    i_C = C.identity() # Cell2
    i_D = D.identity() # Cell2
    assert i_A * i_A == i_A

    tgt, src = A*(B*C), (A*B)*C
    a = Cell2.reassociate(A, B, C)
    assert a.tgt == tgt
    assert a.src == src

    ru_A = A.right_unitor()
    assert ru_A.src == A * I_m
    assert ru_A.tgt == A

    lu_B = B.left_unitor()
    assert lu_B.src == I_m * B
    assert lu_B.tgt == B

    a = Cell2.reassociate(A, I_m, B)

    # check the triangle equation for identity 1-cell's
    lhs = ru_A << i_B
    rhs = (i_A << lu_B) * a
    assert lhs == rhs

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
    assert (B.identity() * f) == f
    assert (f * A.identity()) == f

    # _assoc
    assert (h*g)*f == h*(g*f)


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



