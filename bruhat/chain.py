#!/usr/bin/env python3

"""
previous version: schur.py



"""

from random import shuffle, seed, choice
from functools import reduce
from operator import mul, add, matmul

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
from bruhat.frobenius import GF
from bruhat.action import Perm, Group, mulclose, mulclose_hom
from bruhat.rep import Young
from bruhat.util import partitions, cross


def shortstr(A):
    s = elim.shortstr(A)
    s = str(s)
    s = s.replace(" 0 ", " . ")
    return s


def none_uniq(grades):
    if len(set(grades))==1:
        return grades[0]


none_add = lambda g0, g1 : g0+g1 if g0 is not None and g1 is not None else None
none_sub = lambda g0, g1 : g0-g1 if g0 is not None and g1 is not None else None
    

class Space(object):
    def __init__(self, ring, n=0, grade=0, name="?"):
        assert isinstance(ring, element.Ring)
        assert type(n) is int
        assert 0<=n
        self.ring = ring
        self.n = n
        self.grade = grade
        self.name = name

    def __str__(self):
        return "%s(%s, grade=%s, name=%r)"%(
            self.__class__.__name__, self.n, self.grade, self.name)
    __repr__ = __str__

    # Note:
    # __eq__ is object identity ! Yes.
    # This is a tricky business, see __new__'s below

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

    def __add__(self, other):
        #assert self.grade == other.grade
        return AddSpace(self, other)

    def __matmul__(self, other):
        return MulSpace(self, other)

    def __inv__(self):
        return DualSpace(self)

    def parse(self, decl):
        A = solve.parse(decl)
        decl = decl.replace('.', '0')
        rows = decl.split()
        rows = [l.strip() for l in rows if l.strip()]
        promote = self.ring.promote
        rows = [list(promote(int(c)) for c in l) for l in rows]
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

    cache = {} # XXX use https://docs.python.org/3/library/weakref.html
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
        #grade = items[0].grade
        grade = none_uniq([item.grade for item in items])
        n = sum(item.n for item in items)
        name = "("+"+".join(item.name for item in items)+")"
        Space.__init__(space, ring, n, grade, name)
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

    def get_swap(self, perm=(1, 0)):
        perm = tuple(perm)
        items = self.items
        assert len(perm) == len(items)
        assert set(perm) == set(range(len(items)))
        tgt = reduce(add, [items[i] for i in perm])
        N = len(items)
        rows = []
        for i in perm:
          row = []
          for j in range(N):
            if i==j:
                A = elim.identity(self.ring, items[i].n)
            else:
                A = elim.zeros(self.ring, items[i].n, items[j].n)
            row.append(A)
          rows.append(row)
        rows = [numpy.concatenate(row, axis=1) for row in rows]
        A = numpy.concatenate(rows)
        return Lin(tgt, self, A)

    def send_outof(self, lins):
        assert lins
        tgt = None
        As = []
        for idx, lin in enumerate(lins):
            assert lin.src == self.items[idx]
            assert tgt is None or tgt == lin.tgt
            tgt = lin.tgt
            As.append(lin.A)
        A = numpy.concatenate(As, axis=1)
        lin = Lin(tgt, self, A)
        return lin

    def send_into(self, lins):
        assert lins
        src = None
        As = []
        for idx, lin in enumerate(lins):
            assert lin.tgt == self.items[idx]
            assert src is None or src == lin.src
            src = lin.src
            As.append(lin.A)
        A = numpy.concatenate(As, axis=0)
        lin = Lin(self, src, A)
        return lin

    def get_slice(self, space):
        i = 0
        for idx, item in enumerate(self.items):
            if item == space:
                return slice(i, i+item.n)
            i += item.n
        assert 0, "space %s not found in %s"%(space, self)


class MulSpace(Space):
    "tensor product of vector spaces"

    cache = {} # XXX use https://docs.python.org/3/library/weakref.html
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
        #grade = sum(item.grade for item in items)
        grade = reduce(none_add, [item.grade for item in items])
        n = reduce(mul, [item.n for item in items])
        name = "".join(item.name for item in items)
        Space.__init__(space, ring, n, grade, name=name)
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

    def get_swap(self, perm=(1, 0)):
        perm = tuple(perm)
        items = self.items
        assert self.grade is not None
        assert len(perm) == len(items)
        assert set(perm) == set(range(len(items)))
        sign = self.ring.one
        N = len(items)
        for i in range(N):
          for j in range(i+1, N):
            if perm[i] > perm[j] and (items[i].grade * items[j].grade) % 2:
                sign *= -1
        A = sign*elim.identity(self.ring, self.n)
        shape = tuple(item.n for item in items)
        A.shape = shape + shape
        axes = list(perm)
        for i in range(N):
            axes.append(i+len(perm))
        A = A.transpose(axes)
        tgt = reduce(matmul, [self.items[i] for i in perm])
        A = A.reshape((tgt.n, self.n)) # its a square matrix
        return Lin(tgt, self, A)


class DualSpace(Space):
    def __init__(self, item):
        Space.__init__(self, item.ring, item.n, -item.grade)
        self.item = item

    def __inv__(self):
        return self.item


class Lin(object):
    def __init__(self, tgt, src, A=None):
        assert tgt.ring is src.ring
        self.ring = tgt.ring
        self.tgt = tgt
        self.src = src
        self.grade = none_sub(tgt.grade, src.grade)
        if A is None:
            A = elim.zeros(self.ring, tgt.n, src.n)
        assert A.shape == (tgt.n, src.n), "%s != %s" % ( A.shape , (tgt.n, src.n) )
        self.hom = (tgt, src) # yes it's backwards, just like shape is.
        self.shape = A.shape
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
        return eq(self.A, B.A)

    def is_identity(self):
        assert self.tgt == self.src, "wah?"
        B = self.tgt.identity()
        return eq(self.A, B.A)

    def __str__(self):
        return shortstr(self.A)

    def __repr__(self):
        return "Lin( %s <--- %s )"%(self.tgt, self.src)

    def weak_eq(self, other):
        assert self.A.shape == other.A.shape
        return eq(self.A, other.A)

    def __eq__(self, other):
        assert self.hom == other.hom
        return eq(self.A, other.A)

    def __ne__(self, other):
        assert self.hom == other.hom
        return not eq(self.A, other.A)

    def __getitem__(self, idx):
        return self.A[idx]

    def __setitem__(self, idx, val):
        self.A[idx] = val

    def __add__(self, other):
        assert self.hom == other.hom
        A = self.A + other.A
        return Lin(*self.hom, A)

    def __sub__(self, other):
        assert self.hom == other.hom
        A = self.A - other.A
        return Lin(*self.hom, A)

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

    def _coequalizer(self, other):
        assert self.src == other.src
        assert self.tgt == other.tgt
        ring = self.ring
        A = elim.coequalizer(ring, self.A, other.A)
        src = self.tgt
        tgt = Space(ring, len(A))
        return Lin(tgt, src, A)

    def coequalizer(self, *others):
        if len(others)==1:
            return self._coequalizer(*others)
        colim = self.tgt.identity()
        prev = colim
#        result = []
        #print("coequalizer")
        #print(self.hom)
        for other in others:
            #print(other.hom)
            assert prev.src == self.tgt
            assert self.hom == other.hom
            other = prev * other
            colim = (prev*self)._coequalizer(other)
#            if result:
#                assert colim.src == result[-1].tgt 
#            result.append(colim)
            prev = colim * prev
#        assert result[0].src == self.tgt
        #colim = reduce(mul, reversed(result))
        colim = prev
        assert colim.src == self.tgt
        return colim


# ------------------------------------------------------------

# ------------------------ graded stuff ----------------------


class Seq(object):
    def __init__(self, lins):
        assert lins
        self.ring = lins[0].ring
        self.lins = list(lins)



class Chain(Seq):
    def __init__(self, lins):
        Seq.__init__(self, lins)
        prev = None
        for lin in lins:
            tgt, src = lin.hom
            if prev is not None:
                assert prev.tgt == src
                assert (lin*prev).is_zero()
            prev = lin

    def __str__(self):
        spaces = [lin.src for lin in self.lins] + [self.lins[-1].tgt]
        return "%s(%s)"%(self.__class__.__name__,
            "-->".join(str(space) for space in spaces))

    def __matmul__(self, other):
        assert isinstance(other, Chain)
        return MulChain(self, other)


class MulChain(Chain):
    "tensor product of Chain complexes"
    def _init(self):
        self.spaces = set()  # found spaces
        self.srcs = {}       # Space -> [Lin]
        self.tgts = {}       # Space -> [Lin]
        self.grades = {}     # int -> [Space]

    def _addspace(self, space):
        spaces = self.spaces
        grades = self.grades
        if space not in spaces:
            spaces.add(space)
            grades.setdefault(space.grade, []).append(space)

    def _addlin(self, lin):
        tgt, src = lin.hom
        self._addspace(src)
        self.srcs.setdefault(lin.src, []).append(lin)
        self._addspace(tgt)
        self.tgts.setdefault(lin.tgt, []).append(lin)

    def __init__(self, lhs, rhs):
        assert lhs.lins
        assert rhs.lins
        assert lhs.ring is rhs.ring
        ring = lhs.ring
        self._init()
        for g in rhs.lins: # cols
            for f in lhs.lins: # rows
                self._addlin( f @ g.src.identity() ) # vertical arrow
                sign = -1 if f.src.grade % 2 else 1
                self._addlin( sign * f.src.identity() @ g ) # _horizontal arrow
            sign = -1 if f.tgt.grade % 2 else 1
            self._addlin( sign * f.tgt.identity() @ g ) # _horizontal arrow
  
        for f in lhs.lins: # rows
            self._addlin( f @ g.tgt.identity() ) # vertical arrow

        keys = list(self.grades.keys())
        keys.sort(reverse=True)
        #print(keys)

        N = len(keys)
        lins = []
        for idx in range(N-1):
            i = keys[idx]
            assert keys[idx+1] == i-1, keys
            tgt = AddSpace(*self.grades[i-1])
            src = AddSpace(*self.grades[i])
            #print(tgt)
            #print(src)
            A = elim.zeros(ring, tgt.n, src.n)
            #print(shortstr(A))
            for s in src.items:
                for lin in self.srcs[s]:
                    assert lin.src is s
                    cols = src.get_slice(lin.src)
                    rows = tgt.get_slice(lin.tgt)
                    A[rows, cols] = lin.A
            #print(shortstr(A))
            lin = Lin(tgt, src, A)
            #print(repr(lin))
            lins.append(lin)
            #print()
        Chain.__init__(self, lins)


# ------------------------------------------------------------

# ------------------------ Young symmetrizers ----------------

# code ripped from qu.py

def test_young():

    d = argv.get("d", 2)
    n = argv.get("n", 3)
    p = argv.get("p")

    if p is None:
        ring = element.Q
    else:
        ring = element.FiniteField(p)
    print("ring:", ring)

    space = Space(ring, d)

    # tensor power of the space
    tspace = reduce(matmul, [space]*n)

    # build action of symmetric group on the tensor power
    items = list(range(n))
    gen1 = []
    gen2 = []
    for i in range(n-1):
        perm = dict((item, item) for item in items)
        perm[items[i]] = items[i+1]
        perm[items[i+1]] = items[i]
        g = Perm(perm, items)
        gen1.append(g)
        lin = tspace.get_swap([g[i] for i in items])
        #print(lin.hom)
        gen2.append(lin)

    perms = mulclose(gen1)
    G = Group(perms, items)

    #print(G)

    action = mulclose_hom(gen1, gen2)
    for g in G:
      for h in G:
        assert action[g*h] == action[g]*action[h] # check it's a group hom

    # Build the young symmetrizers
    projs = []
    part = argv.get("part")
    parts = list(partitions(n)) if part is None else [part]
    for part in parts:
        assert sum(part) == n

        t = Young(G, part)

        rowG = t.get_rowperms()
        colG = t.get_colperms()
        horiz = None
        for g in rowG:
            P = action[g]
            horiz = P if horiz is None else (horiz + P)

        vert = None
        for g in colG:
            P = action[g]
            s = g.sign()
            P = ring.promote(s)*P
            vert = P if vert is None else (vert + P)
        A = horiz * vert

        assert vert*vert == len(colG) * vert
        assert horiz*horiz == len(rowG) * horiz
        #A = A.transpose()
        projs.append(A)

        print("part:", part)
        print(t)
        print("rank:", A.rank())
        print("is_zero:", A.is_zero())
        #if not A.is_zero():
        #    print(A)

        print()



# ------------------------------------------------------------

# ------------------------ testing      ----------------------

def test():

    p = argv.get("p", 3)
    ring = element.FiniteField(p)

    space = Space(ring)
    f = space.parse("11. .11")
    assert -f == -1*f

    assert f.coequalizer(f).weak_eq(f.tgt.identity())

    U, V, W = Space(ring, 2), Space(ring, 3), Space(ring, 5)
    F = Space(ring, 2, 1) # fermionic space

    assert hash(U) is not None
    assert U+V == U+V
    assert U+V is U+V
    assert U+V+W is U+V+W
    assert len((U+V+W).items) == 3

    sUV = (U+V).get_swap()
    sVU = (V+U).get_swap()
    assert sVU*sUV == (U+V).identity()

    a = (U+U+U).get_swap([1,0,2])
    b = (U+U+U).get_swap([0,2,1])
    assert a*b*a == b*a*b

    # tensor --------------------

    sUV = (U@V).get_swap()
    sVU = (V@U).get_swap()
    assert sVU*sUV == (U@V).identity()

    s1 = (U@V@W).get_swap([1,0,2]) # V U W
    s2 = (V@U@W).get_swap([1,0,2]) # V U W
    assert s2*s1 == (U@V@W).identity()

    a = (U@U@U).get_swap([1,0,2])
    b = (U@U@U).get_swap([0,2,1])
    assert a*b*a == b*a*b

    s = (U@U@U).get_swap([1,2,0])
    assert s*s != s.src.identity()
    assert s*s*s == s.src.identity()

    s1 = (U@V@W).get_swap([1,2,0])
    s2 = s1.tgt.get_swap([1,2,0])
    s3 = s2.tgt.get_swap([1,2,0])
    assert s3*s2*s1 == s1.src.identity()

    sUU = (U@U).get_swap()
    assert sUU*sUU == (U@U).identity()

    sFF = (F@F).get_swap()
    assert sFF*sFF == (F@F).identity()

    A = Lin.rand(U, V, 1, 1)
    B = Lin.rand(U, V, 1, 1)
    assert A@B != B@A

    sUU = (U@U).get_swap()
    sVV = (V@V).get_swap()

    assert sUU*(A@B)*sVV == B@A

    # symmetric square
    for m in range(1, 5):
        U = Space(ring, m)
        UU = U@U
        I = UU.identity()
        s = UU.get_swap()
        f = I.coequalizer(s)
        assert f == f*s
        assert f.rank() == [1, 3, 6, 10][m-1]

    # alternating squares 
    for m in range(1, 5):
        #U = Space(ring, m, grade=1) # fermionic Space
        U = Space(ring, m)
        UU = U@U
        I = UU.identity()
        s = -UU.get_swap()
        f = I.coequalizer(s)
        assert f == f*s
        assert f.rank() == [0, 1, 3, 6][m-1]


def super_young(U, V, part):
    ring = U.ring

    n = sum(part)

    lookup = {}
    summands = []
    for idxs in cross([(0, 1)]*n):
        spaces = [[U, V][idx] for idx in idxs]
        space = reduce(matmul, spaces)
        lookup[idxs] = len(summands)
        summands.append(space)
        print(idxs, end=" ")
        print(space.n)

    src = reduce(add, summands)

    dsum = lambda U, V : U.direct_sum(V)

    G = Group.symmetric(n)
    hom = {}

    colim = src.identity()
    diagram = []
    for action in G:
        perm = tuple(action[i] for i in range(n))
        print(perm)
        sumswap = [None] * len(summands)
        fs = []
        for i, idxs in enumerate(cross([(0, 1)]*n)):
            jdxs = tuple(idxs[i] for i in perm)
            print(idxs, "-->", jdxs)
            sumswap[lookup[jdxs]] = i
            space = src.items[i]
            f = space.get_swap(perm)
            fs.append(f)
        print(sumswap)
        f = reduce(dsum, fs)
        print(f.src.name, "-->", f.tgt.name)
        g = f.tgt.get_swap(sumswap)
        print(g.src.name, "-->", g.tgt.name)
        assert f.tgt == g.src
        assert g.tgt == src
        #gf = g*f
        #diagram.append(gf)
        #colim = colim.coequalizer(gf)
        hom[action] = (g,f)
        print()

    #colim = src.identity().coequalizer(*diagram)
    young = Young(G, part)
    print(young)
#
#        horiz = None
#        for g in rowG:
#            P = action[g]
#            horiz = P if horiz is None else (horiz + P)
#
#        vert = None
#        for g in colG:
#            P = action[g]
#            s = g.sign()
#            P = ring.promote(s)*P
#            vert = P if vert is None else (vert + P)
#        A = horiz * vert


    rowG = young.get_rowperms()
    colG = young.get_colperms()
    print("rowG", len(rowG))
    print("colG", len(colG))

    horiz = None
    for action in rowG:
        (g,f) = hom[action]
        gf = g*f
        horiz = gf if horiz is None else (horiz + gf)
    assert horiz is not None

    vert = None
    for action in colG:
        sign = action.sign()*ring.one
        print(action, sign)
        g,f = hom[action]
        g = sign*g
        gf = g*f
        vert = gf if vert is None else (vert + gf)
    assert vert is not None
    P = horiz*vert

    print(P.rank())

#    print("diagram:", len(diagram))
#    colim = diagram[0].coequalizer(*diagram[1:])
#
#    assert colim.src == src
#    print(colim.shape)
#    print(colim)
#    print(colim.rank())
#
#    lhs = colim * diagram[0]
#    for f in diagram:
#        assert lhs == colim * f


def test_super():
    ring = element.Q

    U = Space(ring, 2, grade=0) # bosonic
    V = Space(ring, 2, grade=1) # fermionic

    lhs = (U@U@U).get_swap([0,2,1]).A
    rhs = (U@U@V).get_swap([0,2,1]).A
    assert eq(lhs, rhs)

    lhs = (U@V@V).get_swap([0,2,1]).A
    rhs = (U@U@U).get_swap([0,2,1]).A
    assert eq(lhs, -rhs)

    lhs = (V@V@V).get_swap([1,2,0]).A
    rhs = (U@U@U).get_swap([1,2,0]).A
    assert eq(lhs, rhs)

    dsum = lambda U, V : U.direct_sum(V)

    a = argv.get("a", 2)
    b = argv.get("b", 2)
    part = argv.get("part", (3,))

    U = Space(ring, a, grade=0, name="U") # bosonic 
    V = Space(ring, b, grade=1, name="V") # fermionic
    super_young(U, V, part)

    return


    for a in range(4):
      for b in range(4):
        U = Space(ring, a, grade=0) # bosonic 
        V = Space(ring, b, grade=1) # fermionic

        UU, UV, VU, VV = U@U, U@V, V@U, V@V

        src = UU + UV + VU + VV
        
        lhs = [
            UU.get_swap(), 
            (VU+UV).get_swap() * dsum(UV.get_swap(), VU.get_swap()),
            VV.get_swap()]
        lhs = reduce(dsum, lhs)

        rhs = src.identity()

        f = lhs.coequalizer(rhs)
        print("%3d"%f.rank(), end=" ")
      print()



def test_gf():
    ring = GF(4)
    x = ring.x
    space = Space(ring)

    m = argv.get("m", 3)
    n = argv.get("n", 4)

    one = ring.one
    A = elim.zeros(ring, m, n)
    for i in range(m):
      for j in range(n):
        A[i, j] = choice(ring.elements)

    H = Lin(Space(ring, m, 0), Space(ring, n, 1), A)
    #print(H)

    c = Chain([H])
    cc = c@c
    #for f in cc.lins:
    #    print(f)


def test_chain():
    p = argv.get("p", 3)
    ring = element.FiniteField(p)

    space = Space(ring)

    m = argv.get("m", 3)
    n = argv.get("n", 4)

    one = ring.one
    if argv.toric:
        V = Space(ring, n, 1)
        U = Space(ring, m, 0)

        A = Lin(U, V)
        for i in range(m):
            A[i, i] = one
            A[i, (i+1)%m] = -one
    elif argv.surface:
        V = Space(ring, m, 1)
        U = Space(ring, m-1, 0)

        A = Lin(U, V)
        for i in range(m-1):
            A[i, i] = one
            A[i, (i+1)%m] = -one
    else:
        V = Space(ring, n, 1)
        U = Space(ring, m, 0)

        A = Lin.rand(U, V)

    c = Chain([A])
    print(c)
    cc = c @ c
    print(cc)
    ccc = c@cc
    print(ccc)



def test_all():
    test()
    test_chain()



if __name__ == "__main__":


    fn = argv.next() or "test_all"
    fn = eval(fn)
    fn()

    print("OK")


