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

#from functools import lru_cache # fail
from random import randint, seed

import numpy

from bruhat import element
from bruhat.element import Keyed, Type, Element, GenericElement


class Space(Keyed, Type):
    "Vector space with a basis, over a ring."
    "Vector space _generated by a set."

    def __init__(self, gen, ring):
        gen = [(x if type(x) is tuple else (x,)) for x in gen]
        gen.sort() # <------- canonical <-----
        self.gen = tuple(gen)
        self.set_gen = set(gen)
        self.dim = len(gen)
        self.ring = ring
        Type.__init__(self)
        #Keyed.__init__(self, (self.gen, self.ring))

    def __contains__(self, x):
        return x in self.set_gen

    def __str__(self):
        return "Space(%s, %s)"%(str(list(self.gen)), self.ring)
    __repr__ = __str__

    def __eq__(A, B):
        assert A.ring == B.ring
        return A.gen == B.gen

    def __ne__(A, B):
        assert A.ring == B.ring
        return A.gen != B.gen

    def __hash__(self):
        return hash(self.gen)

    the_star = ("*",)

    @classmethod
    def unit(cls, ring):
        "Tensor unit object"
        return cls([cls.the_star], ring)
    
    def basis(self, i):
        assert 0<=i<self.dim
        x = self.gen[i]
        ring = self.ring
        v = Map([((x, self.the_star), self.ring.one)], Hom(self.unit(self.ring), self))
        return v

    def __add__(A, B):
        assert A.ring == B.ring
        assert not A.set_gen.intersection(B.set_gen)
        gen = A.gen + B.gen
        return Space(gen, A.ring)

    def tensor(A, B):
        assert A.ring == B.ring
        gen = [x+y for x in A.gen for y in B.gen] # tuple addition
        return Space(gen, A.ring)
    __mul__ = tensor
        
    def __matmul__(A, B):
        assert A.ring == B.ring
        if A==B:
            return A
        else:
            assert 0

    @classmethod
    def zero(cls, ring):
        "zero object"
        return cls([], ring)
    
    @property
    def dual(self):
        # XX um...
        return self

    @property
    def ident(self): # XX cache me
        "ident: A --> A"
        ring = self.ring
        one, zero = ring.one, ring.zero
        items = [((x, x), one) for x in self.gen]
        ident = Map(items, Hom(self, self))
        return ident

    @property
    def left_unitor(self): # XX cache me
        "IxA --> A"
        src = Space.unit(self.ring) * self
        tgt = self
        one = self.ring.one
        items = [((x, Space.the_star+x), one) for x in self.gen] # tuple addition
        return Map(items, Hom(src, tgt))

    @property
    def right_unitor(self): # XX cache me
        "AxI --> A"
        src = self * Space.unit(self.ring)
        tgt = self
        one = self.ring.one
        items = [((x, x+Space.the_star), one) for x in self.gen] # tuple addition
        return Map(items, Hom(src, tgt))

    @property
    def cap(self): # XX cache me
        "the unit, eta: I --> A*xA"
        src = Space.unit(self.ring)
        the_star = Space.the_star
        tgt = self.dual * self
        one = self.ring.one
        items = [((x+x, the_star), one) for x in self.gen] # tuple addition
        return Map(items, Hom(src, tgt))

    @property
    def cup(self): # XX cache me
        "the co-unit, epsilon: AxA* --> I"
        tgt = Space.unit(self.ring)
        the_star = Space.the_star
        src = self * self.dual
        one = self.ring.one
        items = [((the_star, x+x), one) for x in self.gen] # tuple addition
        return Map(items, Hom(src, tgt))


class Hom(Keyed, Type):
    def __init__(self, src, tgt):
        assert isinstance(src, Space)
        assert isinstance(tgt, Space)
        assert src.ring == tgt.ring
        Keyed.__init__(self, (src, tgt))
        Type.__init__(self)
        self.ring = src.ring
        self.src = src
        self.tgt = tgt
        self.shape = (tgt.dim, src.dim)

    @property
    def zero(self):
        return Map([], self) # the zero map

    def __getitem__(self, i):
        assert 0<=i<1
        return [self.src, self.tgt][i]

    def tensor(a, b):
        # not sure if this makes sense mathematically..
        src = a.src * b.src
        tgt = a.tgt * b.tgt
        return Hom(src, tgt)
    __mul__ = tensor

    def __matmul__(a, b):
        assert b.tgt == a.src, "%s != %s" % (b.tgt, a.src)
        return Hom(b.src, a.tgt)

    def transpose(a):
        return Hom(a.tgt, a.src)


class Map(Element):

    def __init__(self, _items, hom):
        assert isinstance(hom, Hom)
        Element.__init__(self, hom)
        ring = hom.src.ring
        zero = ring.zero
        assert hom.tgt.ring == ring
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
            assert i in hom.tgt, "%s not found in %s" % (i, hom.tgt) # row
            assert j in hom.src, "%s not found in %s" % (j, hom.src) # col
        assert len(keys)==len(items), "duplicate key: %s" % (keys)
        items.sort() # Make canonical. careful with this...
        self.items = tuple(items)
        self.map_items = dict(items)
        self.src = hom.src
        self.tgt = hom.tgt
        self.hom = hom
        self.ring = ring

    def __repr__(self):
        return "Map(%s)"%(str(list(self.items)))

    def __str__(self):
        zero = self.ring.zero
        map_items = self.map_items
        rows = [[str(map_items.get((i, j), zero)) 
            for j in self.src.gen] for i in self.tgt.gen]
        w = 1
        for row in rows:
            for col in row: 
                w = max(w, len(col))
        rows = ['[%s]'%' '.join(s.rjust(w) for s in row) for row in rows]
        lines = [] 
        for i, row in enumerate(rows):
            if i==0:
                row = "["+row
            else:
                row = " "+row
            if i==len(rows)-1:
                row = row+"]"
            else:
                row = row+","
            lines.append(row)
        return '\n'.join(lines)

    def to_array(self):
        zero = self.ring.zero
        map_items = self.map_items
        rows = [[map_items.get((i, j), zero)
            for j in self.src.gen] for i in self.tgt.gen]
        A = numpy.array(rows)
        return A

    @classmethod
    def from_array(cls, A, hom):
        ring = hom.ring
        zero = ring.zero
        items = []
        assert A.shape == hom.shape
        for idx, i in enumerate(hom.tgt.gen):
          for jdx, j in enumerate(hom.src.gen):
            u = A[idx, jdx]
            if u != zero:
                items.append(((i, j), u))
        return cls(items, hom)

    def __eq__(a, b):
        return a.items == b.items

    def __ne__(a, b):
        return a.items != b.items

    def __hash__(self):
        return hash(self.items)

    def __contains__(self, item):
        return item in self.set_items

    # these operations could also live in the Hom:

    def __add__(a, b):
        assert a.hom == b.hom
        zero = a.ring.zero
        map_items = dict(a.items)
        for ((i, j), u) in b.items:
            val = map_items.get((i, j), zero) + u
            map_items[i, j] = val
        items = list(map_items.items())
        return Map(items, a.hom)

    def __sub__(a, b):
        assert a.hom == b.hom
        zero = a.ring.zero
        map_items = dict(a.items)
        for ((i, j), u) in b.items:
            val = map_items.get((i, j), zero) - u
            map_items[i, j] = val
        items = list(map_items.items())
        return Map(items, a.hom)

    def __neg__(a):
        items = [((i, j), -u) for ((i, j), u) in a.items]
        return Map(items, self.hom)

    def tensor(a, b):
        hom = a.hom*b.hom
        items = [
            ((x+y, u+v), val*wal)
            for ((x, u), val) in a.items 
            for ((y, v), wal) in b.items]
        return Map(items, hom)
    __mul__ = tensor # ?

#    def __rmul__(a, r):
#        return NotImplemented

    def __matmul__(a, b):
        hom = a.hom@b.hom
        zero = a.ring.zero
        map_items = dict()
        for ((i, j), u) in a.items:
          for ((j1, k), v) in b.items:
            if j != j1:
                continue
            val = map_items.get((i, k), zero) + u*v
            map_items[i, k] = val
        items = list(map_items.items())
        return Map(items, hom)

    def __rmatmul__(a, r):
        #assert isinstance(r, Element)
        #assert r.tp == self.ring
        r = a.ring.promote(r)
        items = [((i, j), u*r) for ((i, j), u) in a.items]
        return Map(items, a.hom)
    __rmul__ = __rmatmul__ # yes?

    def transpose(a):
        items = [((j, i), v) for ((i, j), v) in a.items] # conjugate v?
        return Map(items, a.hom.transpose())



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




def swap_row(A, j, k):
    row = A[j, :].copy()
    A[j, :] = A[k, :]
    A[k, :] = row

def swap_col(A, j, k):
    col = A[:, j].copy()
    A[:, j] = A[:, k]
    A[:, k] = col


def row_reduce(f, truncate=False, inplace=False, check=False, verbose=False):
    """ Remove zero rows if truncate==True
    """

    A = f.to_array()
    ring = f.ring
    zero = ring.zero

    assert len(A.shape)==2, A.shape
    m, n = A.shape
    if not inplace:
        A = A.copy()

    if m*n==0:
        if truncate and m:
            A = A[:0, :]
        return A

    if verbose:
        print("row_reduce")
        #print("%d rows, %d cols" % (m, n))

    i = 0
    j = 0
    while i < m and j < n: 
        if verbose:
            print("i, j = %d, %d" % (i, j))
            print("A:")
            print(shortstrx(A))

        assert i<=j 
        if i and check:
            assert (A[i:,:j]!=0).sum() == 0

        # first find a nonzero entry in this col
        for i1 in range(i, m):
            if A[i1, j]:
                break
        else:
            j += 1 # move to the next col
            continue # <----------- continue ------------

        if i != i1:
            if verbose:
                print("swap", i, i1)
            swap_row(A, i, i1)

        assert A[i, j] != zero
        for i1 in range(i+1, m):
            if A[i1, j]:
                if verbose: 
                    print("add row %s to %s" % (i, i1))
                r = -A[i1, j] / A[i, j]
                A[i1, :] += r*A[i, :]
                assert A[i1, j] == zero

        i += 1 
        j += 1 

    if truncate:
        m = A.shape[0]-1
        #print("sum:", m, A[m, :], A[m, :].sum())
        while m>=0 and (A[m, :]!=0).sum()==0:
            m -= 1
        A = A[:m+1, :]

    if verbose:
        print()

    f1 = Map.from_array(A, f.hom)

    return f1





def test_over_ring(ring):

    zero = Space.zero(ring)
    I = Space.unit(ring)

    A = Space("abcd", ring)
    B = Space("uv", ring)
    C = Space("1234", ring)
    D = Space("678", ring)

    assert str(A)
    assert (A*B)*C == A*(B*C)
    assert zero*A == A*zero == zero
    assert A*(B+C) == A*B + A*C

    a = A.basis(0)
    u = B.basis(0)
    v = B.basis(0)

    f = (a.transpose()*u) # kind of a mess...
    #print(f.hom)

    one = ring.one
    f = Map([(('u', 'a'), one), (('v', 'a'), one)], Hom(A, B)) # A--f-->B

    U = f.to_array()
    f1 = Map.from_array(U, f.hom)
    assert f==f1

    assert str(f)
    assert f-f == f.hom.zero
    assert f+f.hom.zero == f
    assert f@A.ident == f
    assert B.ident@f == f

    assert ((f * A.ident) @ A.cap).hom == Hom(I, B*A)
    assert ((A.ident * f) @ A.cap).hom == Hom(I, A*B)
    
    g = Map([(('6', '2'), one), (('8', '3'), one)], Hom(C, D)) # C--g-->D

    h = Map([(('d', 'a'), one), (('c', 'b'), one), (('c', 'c'), one)], Hom(A, A)) # A--g-->A

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

    #assert (2@f) == f+f # do we swap @ and * ? this looks ugly...
    assert 2*f == f+f



def zeros(ring, m, n):
    A = numpy.empty((m, n), dtype=object)
    A[:] = ring.zero
    return A


def rand(ring, m, n, p=3, q=3):
    A = zeros(ring, m, n)
    for i in range(m):
      for j in range(n):
        a = randint(-p, q)
        A[i, j] = ring.promote(a) / randint(1, q)
    return A


def test():

    Q = element.Q

    n = 5
    gen = list(range(n))
    X = Space(gen, Q)

    A = rand(Q, n, n, 0, 1)

    hom = Hom(X, X)
    f = Map.from_array(A, hom)

    print(f)

    f = row_reduce(f)
    
    print(f)
    print()




if __name__ == "__main__":


    test_over_ring(element.Z)
    test_over_ring(element.Q)
    test()
    print("OK")


