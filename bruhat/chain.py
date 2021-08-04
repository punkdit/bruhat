#!/usr/bin/env python3

"""
Graded vector spaces and linear maps over a ring, schur functors, ...
Use numpy object array's .
The space objects enforce type checking among the linear maps.

previous version: schur.py
"""

from random import shuffle, seed, choice
from functools import reduce
from operator import mul, add, matmul

import numpy

from bruhat.argv import argv
if argv.fast:
    from bruhat import _element as element
else:
    from bruhat import element

from bruhat import elim
from bruhat.elim import eq
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
    """
    vector Space over a field (aka a ring).
    These have a dimenion 'n', a 'grade', and a 'name' (for debugging).
    """
    def __init__(self, ring, n=0, grade=0, name="?"):
        assert isinstance(ring, element.Ring), ring.__class__
        assert type(n) is int
        assert grade is None or type(grade) is int, repr(grade)
        assert type(name) is str
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
        #if len(_items)==1:
        #    return _items[0]
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
        #print("send_into", self)
        for idx, lin in enumerate(lins):
            #print(lin)
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
        name = "@".join(item.name for item in items)
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
    """
    A _linear map between Space's,
    implemented using a 2d numpy array with object entries.
    """
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
        assert self.hom == other.hom, "%s != %s" % (self.hom, other.hom)
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
        assert other.tgt == self.src, "%s != %s" % (other.tgt, self.src)
        A = numpy.dot(self.A, other.A)
        return Lin(self.tgt, other.src, A)

    def __rshift__(self, other):
        return other * self

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

    def row_reduce(self):
        A = elim.row_reduce(self.ring, self.A, True)
        tgt = Space(self.ring, len(A))
        return Lin(tgt, self.src, A)

    def direct_sum(self, other):
        tgt = self.tgt + other.tgt
        src = self.src + other.src
        f = Lin.zeros(tgt, src)
        m, n = self.A.shape
        f.A[:m, :n] = self.A
        f.A[m:, n:] = other.A
        return f

    def _coequalizer(self, other, ref=None):
        assert self.src == other.src
        assert self.tgt == other.tgt
        ring = self.ring
        if ref is not None:
            A, B = elim.coequalizer(ring, self.A, other.A, ref.A)
            src = self.tgt
            tgt = Space(ring, len(A))
            a = Lin(tgt, src, A)
            b = Lin(ref.tgt, tgt, B)
            return a, b
        else:
            A = elim.coequalizer(ring, self.A, other.A)
            src = self.tgt
            tgt = Space(ring, len(A))
            return Lin(tgt, src, A)

    def coequalizer(self, *others, ref=None):
        if len(others)==1:
            return self._coequalizer(*others, ref=ref)
        assert ref is None, "not implemented..."
        colim = self.tgt.identity()
        prev = colim
        for other in others:
            assert prev.src == self.tgt
            assert self.hom == other.hom
            other = prev * other
            colim = (prev*self)._coequalizer(other)
            prev = colim * prev
        colim = prev
        assert colim.src == self.tgt
        return colim

    def pushout(self, other):
        assert self.src == other.src
        J, K = elim.pushout(self.ring, self.A, other.A)
        assert len(J) == len(K)
        tgt = Space(self.ring, len(J))
        J = Lin(tgt, self.tgt, J)
        K = Lin(tgt, other.tgt, K)
        return J, K

    def pullback(self, other):
        assert self.tgt == other.tgt
        J, K = elim.pullback(self.ring, self.A, other.A)
        n = J.shape[1]
        assert K.shape[1] == n
        src = Space(self.ring, n)
        J = Lin(self.src, src, J)
        K = Lin(other.src, src, K)
        return J, K


# ------------------------------------------------------------

# ------------------------              ----------------------

# 2-rig of matrices of Lin's ?

class Matrix(object):
    def __init__(self, tgts, srcs):
        self.tgts = list(tgts)
        self.srcs = list(srcs)


# ------------------------------------------------------------

# ------------------------ graded stuff ----------------------


class Seq(object):
    def __init__(self, lins):
        assert lins
        self.ring = lins[0].ring
        self.lins = list(lins)

    def __len__(self):
        return len(self.lins)

    def __getitem__(self, idx):
        return self.lins[idx]



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

    def get(self, idx):
        lins = self.lins
        n = len(lins)
        if idx < n:
            return lins[idx].src
        elif idx == n:
            return lins[n-1].tgt
        raise IndexError

    def get_grades(self):
        lins = self.lins
        if not lins:
            return [] # ?
        return [lin.src for lin in lins] + [lins[-1].tgt]

    def identity(chain):
        n = len(chain)
        lins = [chain.get(i).identity() for i in range(n+1)]
        return ChainMap(chain, chain, lins)

    def __str__(self):
        spaces = [lin.src for lin in self.lins] + [self.lins[-1].tgt]
        return "%s(%s)"%(self.__class__.__name__,
            "-->".join(str(space) for space in spaces))

    def __matmul__(self, other):
        assert isinstance(other, Chain)
        return MulChain(self, other)


class ChainMap(object):
    def __init__(self, tgt, src, lins=None, check=True):
        assert isinstance(tgt, Chain)
        assert isinstance(src, Chain)
        n = len(src)
        assert len(tgt) == n
        if lins is None:
            # zero map
            lins = [Lin(tgt.get(i), src.get(i)) for i in range(n+1)]
        assert len(lins) == n+1
        for i in range(n):
            assert src[i].src == lins[i].src
            assert tgt[i].src == lins[i].tgt
        if n:
            assert src[n-1].tgt == lins[n].src
            assert tgt[n-1].tgt == lins[n].tgt
        self.ring = tgt.ring
        self.tgt = tgt
        self.src = src
        self.hom = (tgt, src) # yes it's backwards, just like shape is.
        self.lins = list(lins)
        if check:
            self._check()

    def _check(self):
        tgt, src = self.hom
        lins = self.lins
        n = len(src)
        for i in range(n):
            assert tgt[i]*lins[i] == lins[i+1]*src[i], "not a chain map"

    def __str__(self):
        return "ChainMap(%s<---%s)"%(self.tgt, self.src)
    __repr__ = __str__

    def __len__(self):
        return len(self.lins)

    def __getitem__(self, idx):
        return self.lins[idx]

    def __eq__(self, other):
        assert self.hom == other.hom
        n = len(self)
        for i in range(n):
            if self[i] != other[i]:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        assert self.hom == other.hom
        n = len(self)
        lins = [self[i]+other[i] for i in range(n)]
        chain = ChainMap(self.tgt, self.src, lins)
        return chain

    def __sub__(self, other):
        assert self.hom == other.hom
        n = len(self)
        lins = [self[i]-other[i] for i in range(n)]
        chain = ChainMap(self.tgt, self.src, lins)
        return chain

    def __rmul__(self, r):
        r = self.ring.promote(r)
        lins = [r*lin for lin in self.lins]
        return ChainMap(self.tgt, self.src, lins)

    def coequalizer(self, other):
        assert self.hom == other.hom
        n = len(self)
        tgt = [] # construct new chain
        lins = [] # construct chain map
        idx = n-1 # working backwards down to idx=0
        while idx >= 0:
            if not lins:
                lin = self[idx].coequalizer(other[idx])
                lins.insert(0, lin)
            else:
                ref = lins[0] * self.tgt[idx]
                lin, f = self[idx].coequalizer(other[idx], ref=ref)
                lins.insert(0, lin)
                tgt.insert(0, f)
            idx -= 1
        tgt = Chain(tgt)
        cmap = ChainMap(tgt, self.tgt, lins)
        return cmap

    def cokernel(self):
        zero = ChainMap(self.tgt, self.src)
        return self.coequalizer(zero)



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

# ------------------------ testing      ----------------------


def test_young():
    # code ripped from qu.py

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



def test_young0():
    ring = element.Q

    n = argv.get("n", 3)
    part = argv.get("part", (1,1,1))
    assert sum(part) == n

    V = Space(ring, n)

    # build action of symmetric group on the space V
    items = list(range(n))
    gen1 = []
    gen2 = []
    for i in range(n-1):
        perm = dict((item, item) for item in items)
        perm[items[i]] = items[i+1]
        perm[items[i+1]] = items[i]
        g = Perm(perm, items)
        gen1.append(g)
        A = elim.zeros(ring, n, n)
        for k,v in perm.items():
            A[v,k] = ring.one
        lin = Lin(V, V, A)
        gen2.append(lin)

    perms = mulclose(gen1)
    G = Group(perms, items)

    #print(G)

    action = mulclose_hom(gen1, gen2)
    for g in G:
      for h in G:
        assert action[g*h] == action[g]*action[h] # check it's a group hom

    young = Young(G, part)
    rowG = young.get_rowperms()
    print("rowG", len(rowG))
    colG = young.get_colperms()
    print("colG", len(colG))

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

    print("part:", part)
    print(young)
    print("rank:", A.rank())
    print("is_zero:", A.is_zero())
    print(A)
    #if not A.is_zero():
    #    print(A)

    print()



def test(ring=element.Q):

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


dsum = lambda U, V : U.direct_sum(V)

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
        #print(idxs, end=" ")
        #print(space.n)

    src = reduce(add, summands)

    G = Group.symmetric(n)
    hom = {}

    colim = src.identity()
    for action in G:
        perm = tuple(action[i] for i in range(n))
        #print(perm)
        sumswap = [None] * len(summands)
        fs = []
        for i, idxs in enumerate(cross([(0, 1)]*n)):
            jdxs = tuple(idxs[i] for i in perm)
            #print(idxs, "-->", jdxs)
            sumswap[lookup[jdxs]] = i
            space = src.items[i]
            f = space.get_swap(perm)
            fs.append(f)
        #print(sumswap)
        f = reduce(dsum, fs)
        #print(f.src.name, "-->", f.tgt.name)
        g = f.tgt.get_swap(sumswap)
        #print(g.src.name, "-->", g.tgt.name)
        assert f.tgt == g.src
        assert g.tgt == src
        hom[action] = (g,f)
        #print()

    young = Young(G, part)
    #print("Young:")
    #print(young)

    rowperms = young.get_rowperms()
    colperms = young.get_colperms()
    #print("rowperms", len(rowperms))
    #print("colperms", len(colperms))

    horiz = None
    for action in rowperms:
        (g,f) = hom[action]
        gf = g*f
        horiz = gf if horiz is None else (horiz + gf)
    assert horiz is not None

    vert = None
    for action in colperms:
        sign = action.sign()*ring.one
        g,f = hom[action]
        g = sign*g
        gf = g*f
        vert = gf if vert is None else (vert + gf)
    assert vert is not None
    P = horiz*vert

    P = P.row_reduce()
    #print(P.rank())
    #print(P)

    #for a in src.items:
    #    print(a)

    even = src.identity()
    odd = src.identity()
    i = 0
    j = 0
    for space in src.items:
        j += space.n
        while i < j:
            if space.grade%2==0:
                odd[i, i] = ring.zero
            else:
                even[i, i] = ring.zero
            i += 1

    return (P*even).rank(), (P*odd).rank()




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

    n = argv.get("n", 4)
    part = argv.get("part", (2,))

    for a in range(n):
      for b in range(n):
        U = Space(ring, a, grade=0, name="U") # bosonic 
        V = Space(ring, b, grade=1, name="V") # fermionic
        evn, odd = super_young(U, V, part)
        print((evn, odd), end=" ", flush=True)
      print()


def test_symmetric_square():

    ring = element.Q

    result = [
        [0,  0,  1,  3 ],
        [1,  2,  4,  7 ],
        [3,  5,  8, 12 ],
        [6,  9, 13, 18 ]]

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
        assert f.rank() == result[a][b]
        #print("%3d"%f.rank(), end=" ")
      #print()


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
    #p = argv.get("p", 3)
    #ring = element.FiniteField(p)
    ring = element.Q

    space = Space(ring)

    m = argv.get("m", 3)
    n = argv.get("n", 4)

    one = ring.one
    if argv.toric:
        V = Space(ring, n, 1, "V")
        U = Space(ring, m, 0, "U")

        A = Lin(U, V)
        for i in range(m):
            A[i, i] = one
            A[i, (i+1)%m] = -one
    elif argv.surface:
        V = Space(ring, m, 1, "V")
        U = Space(ring, m-1, 0, "U")

        A = Lin(U, V)
        for i in range(m-1):
            A[i, i] = one
            A[i, (i+1)%m] = -one
    else:
        V = Space(ring, n, 1, "V")
        U = Space(ring, m, 0, "U")
        A = Lin.rand(U, V)

    c = Chain([A])
    print(c)

    cc = c @ c
    print(cc)
    #ccc = c@cc
    #print(ccc)


def test_chainmap():

    p = argv.get("p", 2)
    ring = element.FiniteField(p)
    #ring = element.Q

    one = ring.one

    space = Space(ring)

    def mk_ising(m, n):
        U = Space(ring, m, 0, "U")
        V = Space(ring, n, 1, "V")
    
        A = Lin(U, V) # U <--- V
        for i in range(m):
            A[i, i] = one
            if i+1<n:
                A[i, i+1] = -one
        return A

    A = mk_ising(3, 4)
    c = Chain([A])

    f = c.identity()
    zero = ChainMap(c, c)
    assert f != zero
    assert f == f
    assert f != 2*f
    assert f+f == 2*f
    assert f-f == zero

    B = mk_ising(2, 2)
    d = Chain([B])

    fn = Lin(A.src, B.src)
    for i in range(len(B.src)):
        fn[i, i] = one

    fm = Lin(A.tgt, B.tgt)
    for i in range(len(B.tgt)):
        fm[i, i] = one

    f = ChainMap(c, d, [fn, fm])

    # -----------------------------------

    m, n = 8, 10
    V1 = Space(ring, n, 1, "V1")
    V0 = Space(ring, m, 0, "V0")
    A = Lin.rand(V0, V1) # V0 <--- V1

    c = Chain([A])

    U0 = Space(ring, m, 0, "U0")
    f0 = Lin.rand(V0, U0)

    f1, B = A.pullback(f0)
    d = Chain([B])
    f = ChainMap(c, d, [f1, f0])
    
    # -----------------------------------
    # construct a chain map (and a chain) from a 
    # pullback of a grade zero map.

    m, n, p = 5, 6, 1
    V1 = Space(ring, n, 1, "V1")
    V0 = Space(ring, m, 0, "V0")
    #A = Lin.rand(V0, V1) # V0 <--- V1
    A = Lin(V0, V1)
    for i in range(m):
        A[i, i] = one
        A[i, (i+1)] = one

    a = Chain([A])
    #print("A:")
    #print(A)

    U0 = Space(ring, p, 0, "U0")
    f0 = Lin(V0, U0)
    for i in range(p):
        f0[i,i] = one

    f1, B = A.pullback(f0)
    b = Chain([B])
    #print("B:")
    #print(B)
    f = ChainMap(a, b, [f1, f0])
    #print(f0)
    #print(f1)

    g = ChainMap(a, b)
    h = f.coequalizer(g)
    #print(h)
    #print(h[0])
    #print(h[1])
    c = h.tgt
    C = c[0]
    #print(C.shape, C.rank())
    #print(C)
    
    # -----------------------------------
    # construct a 'puncture' of a repitition code 
    # as a cokernel 

    m, n, p = 5, 6, 1
    V1 = Space(ring, n, 1, "V1")
    V0 = Space(ring, m, 0, "V0")
    #A = Lin.rand(V0, V1) # V0 <--- V1
    A = Lin(V0, V1)
    for i in range(m):
        A[i, i] = one
        A[i, (i+1)] = one

    a = Chain([A])

    U1 = Space(ring, 1, 1, "U1")
    U0 = Space(ring, 2, 0, "U0")
    B = Lin(U0, U1)
    B[0, 0] = one
    B[1, 0] = one
    b = Chain([B])

    offset = 2 # where to puncture

    f0 = Lin(V0, U0)
    f0[0+offset,0] = one
    f0[1+offset,1] = one

    f1 = Lin(V1, U1)
    f1[offset+1,0] = one

    f = ChainMap(a, b, [f1, f0])
    h = f.cokernel()
    c = h.tgt

    # -----------------------------------
    # Here we puncture a bit from a parity check matrix
    # by using a cokernel. This deletes that column,
    # as well as all the rows with support therein.

    m, n, p = 5, 6, 1
    V1 = Space(ring, n, 1, "V1")
    V0 = Space(ring, m, 0, "V0")
    A = Lin.rand(V0, V1) # V0 <--- V1

    a = Chain([A])
    #print(A)
    col = 0
    rows = []
    for i in range(V0.n):
        if A[i, col]:
            rows.append(i)

    U1 = Space(ring, 1, 1, "U1") 
    U0 = Space(ring, len(rows), 0, "U0")
    B = Lin(U0, U1)
    for i in range(U0.n):
        B[i, 0] = one
    b = Chain([B])
    #print(B)

    f0 = Lin(V0, U0)
    for i, row in enumerate(rows):
        f0[row, i] = one
    #print(f0)

    f1 = Lin(V1, U1)
    f1[col,0] = one

    f = ChainMap(a, b, [f1, f0])
    h = f.cokernel()
    c = h.tgt

    #print(c[0])





def test_all():
    test(element.Q)

    test_super() # XXX broken
    test_symmetric_square()

    if not argv.fast:
        p = argv.get("p", 3)
        ring = element.FiniteField(p)
        test(ring)
        test_gf()
        test_chain()



if __name__ == "__main__":


    fn = argv.next() or "test_all"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK")


