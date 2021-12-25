#!/usr/bin/env python3

"""
Kapranov-Voevodsky 2-vector spaces.

"""

from functools import reduce
import operator

import numpy

from bruhat import element
from bruhat.chain import Space, Lin


class Matrix(object):
    """
    An object numpy array with type'd entries.
    """
    def __init__(self, tp, A):
        assert type(A) is numpy.ndarray
        self.tp = tp
        self.shape = A.shape
        rows, cols = A.shape
        for row in range(rows):
          for col in range(cols):
            assert isinstance(A[row, col], tp)
        self.A = A.copy()

    def __getitem__(self, idx):
        return self.A[idx]

    def __mul__(self, other):
        assert self.tp == other.tp
        assert self.src == other.tgt
        A = numpy.dot(self.A, other.A)
        return Matrix(self.tp, A)


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


class Lin0(object):
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
        return "Lin0(%d, %r)"%(self.n, self.name)
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
        return AddLin0(self, other)

    def __mul__(self, other):
        return MulLin0(self, other)

    # should we use a list (a basis) instead of n ??

    def identity(self):
        rig = self.rig
        A = numpy.empty((self.n, self.n), dtype=object)
        for row in self:
          for col in self:
            A[row, col] = rig.one if row==col else rig.zero
        return Lin1(self, self, A)

    def zero(tgt, src): # tgt <---- src
        assert tgt.rig == src.rig
        rig = tgt.rig
        A = numpy.empty((tgt.n, src.n), dtype=object)
        for row in tgt:
          for col in src:
            A[row, col] = rig.zero
        return Lin1(tgt, src, A)


class AddLin0(Lin0):

    cache = {} # XXX use https://docs.python.org/3/library/weakref.html
    def __new__(cls, *_items):
        assert _items
        #if len(_items)==1:
        #    return _items[0] # ???
        items = []
        for item in _items:
            if type(item) is AddLin0:
                items += item.items # associative on the nose
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
        Lin0.__init__(lin0, rig, n, name)
        lin0.items = items
        cls.cache[key] = lin0
        return lin0

    def __init__(self, *_items):
        pass




class Lin1(Matrix):
    """
    The 1-cell's (morphism's) : a Matrix of Space's
    """
    def __init__(self, tgt, src, A):
        assert isinstance(tgt, Lin0)
        assert isinstance(src, Lin0)
        assert tgt.rig == src.rig
        self.rig = tgt.rig
        self.tgt = tgt
        self.src = src
        self.hom = (tgt, src) # yes it's backwards, just like shape is.
        Matrix.__init__(self, Space, A)

    def __str__(self):
        s = str(self.A)
        s = s.replace("\n", " ")
        return "Lin1(%s, %s, %s)"%(self.tgt, self.src, s)
    __repr__ = __str__

    def __mul__(self, other):
        assert self.rig == other.rig
        assert self.src == other.tgt
        rig = self.rig
        A = rig.dot(self.A, other.A)
        return Lin1(self.tgt, other.src, A)


        



def main():

    ring = element.Q
    rig = Rig(ring)

    V0 = Lin0(rig, 0, "V0")
    V1 = Lin0(rig, 1, "V1")
    V2 = V1 + V1
    print(V2)

    iV0 = V0.identity()
    iV1 = V1.identity()

    f = Lin0.zero(V2, V1)

    print(f)
    print(f * iV1)


if __name__ == "__main__":

    main()

    print()



