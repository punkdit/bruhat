#!/usr/bin/env python3

"""
Category of finite dimensional vector spaces as a compact closed category.
We avoid using associators by treating tensor as a multi-arity operation.
(Even though python treats it as a binary operation.)
Based on rel.py 

See:
Physics, Topology, Logic and Computation: A Rosetta Stone
John C. Baez, Mike Stay
https://arxiv.org/abs/0903.0340

Also:
2-Hilbert spaces:
https://arxiv.org/pdf/q-alg/9609018.pdf
"""

import sys, os
from functools import lru_cache

from bruhat.element import Z


class Space(object):
    "Vector space with an *ordered* basis, over a ring."

    def __init__(self, basis, ring):
        basis = [(x if type(x) is tuple else tuple(x)) for x in basis]
        self.basis = tuple(basis)
        self.set_basis = set(basis)
        self.ring = ring

    def __contains__(self, x):
        return x in self.set_basis

    def __str__(self):
        return "Space(%s)"%(str(list(self.items)))
    __repr__ = __str__

    def __eq__(A, B):
        assert A.ring == B.ring
        return A.basis == B.basis # ordered basis

    def __ne__(A, B):
        assert A.ring == B.ring
        return A.basis != B.basis # ordered basis

    def __hash__(self):
        return hash(self.basis) # ordered basis

    def __mul__(A, B):
        "tensor"
        assert A.ring == B.ring
        basis = [x+y for x in A.basis for y in B.basis] # tuple addition
        return Space(basis, A.ring)
        
    def __matmul__(A, B):
        assert A.ring == B.ring
        if A==B:
            return A
        else:
            assert 0

    @classmethod
    @lru_cache()
    def zero(cls, ring):
        "zero object"
        return cls([], ring)
    
    the_star = ("*",)

    @classmethod
    @lru_cache()
    def unit(cls, ring):
        "Tensor unit object"
        return cls([cls.the_star], ring)
    
    #@lru_cache
    @property
    def ident(self):
        ring = self.ring
        one, zero = ring.one, ring.zero
        items = [((x, x), one) for x in self.basis]
        ident = Map(items, self, self)
        return ident

    @property
    def left_unitor(self):
        "IxA --> A"
        src = Space.unit(self.ring) * self
        tgt = self
        one = self.ring.one
        basis = [((Space.the_star+x, x), one) for x in self.basis] # tuple addition
        return Map(basis, src, tgt)

    @property
    def right_unitor(self):
        "AxI --> A"
        src = self * Space.unit(self.ring)
        tgt = self
        one = self.ring.one
        basis = [((x+Space.the_star, x), one) for x in self.basis] # tuple addition
        return Map(basis, src, tgt)

    @property
    def cap(self):
        "the unit, eta: I --> AxA"
        src = Space.unit(self.ring)
        the_star = Space.the_star
        tgt = self * self
        one = self.ring.one
        basis = [((the_star, x+x), one) for x in self.basis] # tuple addition
        return Map(basis, src, tgt)

    @property
    def cup(self):
        "the co-unit, epsilon: AxA --> I"
        tgt = Space.unit(self.ring)
        the_star = Space.the_star
        src = self * self
        one = self.ring.one
        basis = [((x+x, the_star), one) for x in self.basis] # tuple addition
        return Map(basis, src, tgt)

    #get_perm


class Map(object):

    def __init__(self, _items, src, tgt):
        assert isinstance(src, Space)
        assert isinstance(tgt, Space)

        ring = src.ring
        zero = ring.zero
        assert tgt.ring == ring
        items = []
        keys = []
        for item in _items:
            (i, j), val = item
            i = (i if type(i) is tuple else (i,))
            j = (j if type(j) is tuple else (j,))
            if val != zero: # make canonical. sparse
                items.append(((i, j), val))
                keys.append((i, j))
            assert val.tp == ring
            assert i in src, "%s not found in %s" % (i, src) # col
            assert j in tgt, "%s not found in %s" % (j, tgt) # row
        assert len(keys)==len(items), "duplicate key: %s" % (keys)
        items.sort() # Make canonical. careful with this...
        self.items = tuple(items)
        self.map_items = dict(items)
        self.src = src
        self.tgt = tgt
        self.shape = (src, tgt)
        self.ring = ring

    def __str__(self):
        return "Map(%s)"%(str(list(self.items)))
    __repr__ = __str__

    def __eq__(a, b):
        return a.items == b.items

    def __ne__(a, b):
        return a.items != b.items

    def __hash__(self):
        return hash(self.items)

    def __contains__(self, item):
        return item in self.set_items

    def __mul__(a, b):
        "tensor"
        src = a.src * b.src
        tgt = a.tgt * b.tgt
        items = [
            ((x+y, u+v), val*wal)
            for ((x, u), val) in a.items 
            for ((y, v), wal) in b.items]
        return Map(items, src, tgt)

    def __matmul__(a, b):
        assert b.tgt == a.src, "%s != %s" % (b.tgt, a.src)

        zero = a.ring.zero
        map_items = dict()
        for ((j, i), u) in a.items:
          for ((k, j1), v) in b.items:
            if j != j1:
                continue
            val = map_items.get((k, i), zero) + u*v
            map_items[k, i] = val
        items = list(map_items.items())
        return Map(items, b.src, a.tgt)

    def transpose(self):
        items = [((j, i), v) for ((i, j), v) in self.items]
        return Map(items, self.tgt, self.src)



def dot(*maps):
    idx = 0
    A = maps[idx]
    while idx+1 < len(maps):
        B = maps[idx+1]
        A = A@B
        idx += 1 
    A = A%2
    return A


def compose(*maps):
    maps = list(reversed(maps))
    A = dot(*maps)
    return A



def test():

    ring = Z

    zero = Space.zero(ring)
    I = Space.unit(ring)

    A = Space("abcd", Z)
    B = Space("uv", Z)
    C = Space("1234", Z)
    D = Space("678", Z)

    assert (A*B)*C == A*(B*C)
    assert zero*A == A*zero == zero

    one = Z.one
    f = Map([(('a', 'u'), one), (('a', 'v'), one)], A, B) # A--f-->B

    assert f@A.ident == f
    assert B.ident@f == f

    assert ((f * A.ident) @ A.cap).shape == (I, B*A)
    assert ((A.ident * f) @ A.cap).shape == (I, A*B)
    
    g = Map([(('2', '6'), one), (('3', '8'), one)], C, D) # C--g-->D

    h = Map([(('a', 'd'), one), (('b', 'c'), one), (('c', 'c'), one)], A, A) # A--g-->A

    assert f*(g*h) == (f*g)*h

    #print(h)

    # right zig-zag equation
    lhs = A.left_unitor @ ( A.cup * A.ident ) @ ( A.ident * A.cap )
    rhs = A.right_unitor
    assert lhs == rhs

    # left zig-zag equation
    lhs = A.right_unitor @ ( A.ident * A.cup ) @ ( A.cap * A.ident )
    rhs = A.left_unitor
    assert lhs == rhs

    # check dual is transpose
    for f in [h, g, f]:
        A, B = f.src, f.tgt
        lhs = A.left_unitor @ ( B.cup * A.ident ) @ (B.ident * f * A.ident ) @ ( B.ident * A.cap )
        rhs = f.transpose() @ B.right_unitor
        assert lhs == rhs

    #print(lhs)
    #print(rhs)



if __name__ == "__main__":

    test()



