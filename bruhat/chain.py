#!/usr/bin/env python3

"""



"""

from random import shuffle, seed
from functools import reduce
from operator import mul

import numpy

from bruhat.argv import argv
from bruhat import element
from bruhat import elim
from bruhat.elim import eq
#from bruhat.elim import (
#    zeros, rand, dot, identity, coequalizer, compose,
#    rank, pseudo_inverse)
#from bruhat.solve import parse
from bruhat import solve


def shortstr(A):
    s = elim.shortstr(A)
    s = str(s)
    s = s.replace(" 0 ", " . ")
    return s


class Space(object):
    def __init__(self, ring, n=0):
        assert isinstance(ring, element.Ring)
        assert type(n) is int
        assert 0<=n
        self.ring = ring
        self.n = n

    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__, self.n)
    __repr__ = __str__

#    def __eq__(self, other):
#        return self.n == other.n
#
#    def __ne__(self, other):
#        return self.n != other.n

    def __hash__(self):
        return id(self)

    def __len__(self):
        return self.n

    def __getitem__(self, idx):
        if idx<0 or idx >= self.n:
            raise IndexError
        return idx

    # we should cache these ... ?

    def __add__(self, other):
        return AddSpace(self, other)

    def __matmul__(self, other):
        return MulSpace(self, other)

    def __inv__(self):
        return DualSpace(self)

    def parse(self, decl):
        A = solve.parse(decl)
        decl = decl.replace('.', '0')
        lines = decl.split()
        lines = [l.strip() for l in lines if l.strip()]
        promote = self.ring.promote
        rows = [list(promote(int(c)) for c in l) for l in lines]
        assert rows
        n = len(rows[0])
        for row in rows:
            assert len(row)==n, "rows have varying lengths"
        A = numpy.array(rows, dtype=object)
        m, n = A.shape
        src = Space(self.ring, n)
        tgt = Space(self.ring, m)
        return Lin(tgt, src, A)

    def identity(self):
        A = elim.identity(self.ring, self.n)
        return Lin(self, self, A)

    def sym2(U):
        UU = U@U
        I = UU.identity()
        s = UU.get_swap()
        f = I.coequalizer(s)
        return f


class AddSpace(Space):
    "direct sum of vector spaces"

    cache = {}
    def __new__(cls, *_items):
        assert _items
        items = []
        for item in _items:
            if type(item) is AddSpace:
                items += item.items
            else:
                items.append(item)
        key = tuple(items)
        if key in cls.cache:
            #print("cache hit", key)
            return cls.cache[key]
        #print("cache miss", key)
        space = object.__new__(cls)
        ring = items[0].ring
        n = sum(item.n for item in items)
        space.n = n
        space.ring = ring
        space.items = items
        cls.cache[key] = space
        return space

    def __init__(self, *_items):
        pass

    def _old__init__(self, *_items):
        items = []
        for item in _items:
            if type(item) is AddSpace:
                items += item.items
            else:
                items.append(item)
        ring = items[0].ring
        n = sum(item.n for item in items)
        Space.__init__(self, ring, n)
        self.items = items

    def get_swap(self):
        assert len(self.items) == 2, "not implemented.."
        U, V = self.items
        other = V+U
        f = Lin(other, self)
        one = self.ring.one
        a, b = len(U), len(V)
        for i in range(b):
            f.A[i, i+a] = one
        for i in range(a):
            f.A[i+b, i] = one
        return f


class MulSpace(Space):
    "tensor product of vector spaces"

    cache = {}
    def __new__(cls, *_items):
        assert _items
        items = []
        for item in _items:
            if type(item) is MulSpace:
                items += item.items
            else:
                items.append(item)
        key = tuple(items)
        if key in cls.cache:
            #print("cache hit", key)
            return cls.cache[key]
        #print("cache miss", key)
        space = object.__new__(cls)
        ring = items[0].ring
        n = reduce(mul, [item.n for item in items])
        space.n = n
        space.ring = ring
        space.items = items
        cls.cache[key] = space
        return space

    def __init__(self, *_items):
        pass

    def _old__init__(self, *_items):
        items = []
        for item in _items:
            if type(item) is MulSpace:
                items += item.items
            else:
                items.append(item)
        ring = items[0].ring
        n = reduce(mul, [item.n for item in items])
        Space.__init__(self, ring, n)
        self.items = items

    def get_swap(self):
        assert len(self.items) == 2, "not implemented.."
        U, V = self.items
        other = V@U
        f = Lin(other, self)
        one = self.ring.one
        a, b = len(U), len(V)
        for i in range(a):
          for j in range(b):
            f.A[j*a + i, i*b + j] = one
        return f


class DualSpace(Space):
    def __init__(self, item):
        Space.__init__(self, item.ring, item.n)
        self.item = item

    def __inv__(self):
        return self.item


class Lin(object):
    def __init__(self, tgt, src, A=None):
        assert tgt.ring is src.ring
        self.ring = tgt.ring
        self.tgt = tgt
        self.src = src
        if A is None:
            A = elim.zeros(self.ring, tgt.n, src.n)
        assert A.shape == (tgt.n, src.n), "%s != %s" % ( A.shape , (tgt.n, src.n) )
        self.shape = (tgt, src)
        self.A = A.copy()

    @classmethod
    def zeros(cls, tgt, src):
        assert tgt.ring is src.ring
        ring = tgt.ring
        A = elim.zeros(ring, tgt.n, src.n)
        return Lin(tgt, src, A)

    @classmethod
    def rand(cls, tgt, src, a=1, b=1):
        assert tgt.ring is src.ring
        ring = tgt.ring
        A = elim.rand(ring, tgt.n, src.n, a, b)
        return Lin(tgt, src, A)

    def is_zero(self):
        B = Lin.zeros(self.tgt, self.src)
        return eq(self.A, B)

    def is_identity(self):
        assert self.tgt == self.src, "wah?"
        B = self.tgt.identity()
        return eq(self.A, B)

    def __str__(self):
        return shortstr(self.A)

    def __repr__(self):
        return "Lin( %s <--- %s )"%(self.tgt, self.src)

    def weak_eq(self, other):
        assert self.A.shape == other.A.shape
        return eq(self.A, other.A)

    def __eq__(self, other):
        assert self.shape == other.shape
        return eq(self.A, other.A)

    def __ne__(self, other):
        assert self.shape == other.shape
        return not eq(self.A, other.A)

    def __getitem__(self, idx):
        return self.A[idx]

    def __setitem__(self, idx, val):
        self.A[idx] = val

    def __add__(self, other):
        assert self.shape == other.shape
        A = self.A + other.A
        return Lin(*self.shape, A)

    def __sub__(self, other):
        assert self.shape == other.shape
        A = self.A - other.A
        return Lin(*self.shape, A)

    def __mul__(self, other):
        assert other.tgt == self.src
        A = numpy.dot(self.A, other.A)
        return Lin(self.tgt, other.src, A)

    def __rmul__(self, r):
        r = self.ring.promote(r)
        A = r*self.A
        return Lin(self.tgt, self.src, A)

    def __neg__(self):
        A = -self.A
        return Lin(self.tgt, self.src, A)

    def __matmul__(self, other):
        src = self.src @ other.src
        tgt = self.tgt @ other.tgt
        ring = self.ring
        A, B = self.A, other.A
        if 0 in A.shape or 0 in B.shape:
            C = self.zeros(A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])
        else:
            #print("kron", A.shape, B.shape)
            C = numpy.kron(A, B)
            #print("\t", C.shape)
        return Lin(tgt, src, C)

    def rank(self):
        return elim.rank(self.ring, self.A)

    def direct_sum(self, other):
        tgt = self.tgt + other.tgt
        src = self.src + other.src
        f = Lin.zeros(tgt, src)
        m, n = self.A.shape
        f.A[:m, :n] = self.A
        f.A[m:, n:] = other.A
        return f

    def coequalizer(self, other):
        assert self.shape == other.shape
        ring = self.ring
        A = elim.coequalizer(ring, self.A, other.A)
        src, _ = self.shape
        tgt = Space(ring, len(A))
        return Lin(tgt, src, A)



def test():

    p = argv.get("p", 3)
    ring = element.FiniteField(p)

    space = Space(ring)
    f = space.parse("11. .11")
    assert -f == -1*f

    assert f.coequalizer(f).weak_eq(f.tgt.identity())

    U, V, W = Space(ring, 3), Space(ring, 4), Space(ring, 7)

    assert hash(U) is not None
    assert U+V == U+V
    assert U+V is U+V
    assert U+V+W is U+V+W
    assert len((U+V+W).items) == 3

    sUV = (U+V).get_swap()
    sVU = (V+U).get_swap()
    assert sVU*sUV == (U+V).identity()

    sUV = (U@V).get_swap()
    sVU = (V@U).get_swap()
    assert sVU*sUV == (U@V).identity()

    A = Lin.rand(U, V, 1, 1)
    B = Lin.rand(U, V, 1, 1)
    assert A@B != B@A

    sUU = (U@U).get_swap()
    sVV = (V@V).get_swap()

    assert sUU*(A@B)*sVV == B@A

    for m in range(1, 5):
        U = Space(ring, m)
        UU = U@U
        I = UU.identity()
        s = UU.get_swap()
        f = I.coequalizer(s)
        assert eq(f, f*s)
        assert f.rank() == [1, 3, 6, 10][m-1]


if __name__ == "__main__":


    fn = argv.next() or "test"
    fn = eval(fn)
    fn()

    print("OK")


