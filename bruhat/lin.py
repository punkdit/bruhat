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
from bruhat import element

from bruhat import elim
from bruhat.elim import eq
from bruhat import solve

from bruhat.action import Perm, Group, mulclose, mulclose_hom
from bruhat.util import partitions, cross, allperms
from bruhat.smap import SMap


def shortstr(A):
    s = elim.shortstr(A)
    s = str(s)
    s = s.replace(" 0 ", " . ")
    return s


def _dot(ring, A, B):
    #C = numpy.dot(self.A, other.A) # argh...
    #print("dot", A.shape, B.shape)
    assert len(A.shape)==len(B.shape)==2
    m, n = A.shape
    assert B.shape[0] == n
    l = B.shape[1]
    C = elim.zeros(ring, m, l)
    for i in range(m):
      for j in range(l):
        for k in range(n):
          C[i,j] += A[i, k] * B[k, j]
    return C

dot = elim.dot # it's a tiny bit faster than _dot


def none_uniq(grades):
    if len(set(grades))==1:
        return grades[0]


none_add = lambda g0, g1 : g0+g1 if g0 is not None and g1 is not None else None
none_sub = lambda g0, g1 : g0-g1 if g0 is not None and g1 is not None else None
    
NO_NAME = "?"

class Space(object):
    """
    vector Space over a field (aka a ring).
    These have a dimenion 'n', a 'grade', and a 'name' (for debugging).
    """
    def __init__(self, ring, n=0, grade=0, name=NO_NAME, _dual=None):
        assert isinstance(ring, element.Type), ring.__class__
        assert type(n) is int
        assert grade is None or type(grade) is int, repr(grade)
        assert type(name) is str
        assert 0<=n
        self.ring = ring
        self.n = n
        self.grade = grade
        self.name = name
        self._dual = _dual
        #if self.n <= 1:
        #    self._dual = self # reasonable ?

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

    # XXX this is not consistent with Lin.__add__
    def __add__(self, other):
        #assert self.grade == other.grade
        return AddSpace(self.ring, self, other)

    def __matmul__(self, other):
        return MulSpace(self.ring, self, other)

    @property
    def dual(self):
        if self._dual is None:
            self._dual = DualSpace(self)
        return self._dual

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

    def asgrade(self, grade, name="?"):
        return Space(self.ring, self.n, grade, name)

    def get_normal(self, N=None, K=None, inverse=False, force=False):
        return self.identity()

    def get_add_slice(self, idx):
        if idx == 0:
            return slice(0, self.n)
        assert 0, "space %s not found in %s"%(idx, self)

    def get_slice(self, space):
        if space is self:
            return slice(0, self.n)
        assert 0, "space %s not found in %s"%(space, self)

    def get_add_items(self):
        return [self]

    def get_add_swap(self, perm=(0,)):
        assert tuple(perm) == (0,)
        return self.identity()

    def get_mul_swap(self, perm=(0,)):
        assert tuple(perm) == (0,)
        return self.identity()



class DualSpace(Space):
    def __init__(self, dual):
        assert type(dual) == Space, type(dual)
        Space.__init__(self, dual.ring, dual.n, -dual.grade, "~"+dual.name, dual)


# We record the bimonoidal structure (+,*) into AddSpace and MulSpace below.
# This is strictly _associative, ie., these are multi-arity operations.
# However, we don't make these two operations strictly unital,
# which may be a mistake as the unit's can be seen as nullary operations.

class AddSpace(Space):
    "direct sum of vector spaces"

    cache = {} # XXX use https://docs.python.org/3/library/weakref.html
    def __new__(cls, ring, *_items, N=None):
        assert isinstance(ring, element.Type), ring.__class__
        assert _items or N is not None
        if not _items:
            return N
        if len(_items)==1:
            return _items[0] # breaks MulChain ...
        items = []
        for item in _items:
            assert item.ring == ring
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
        #ring = items[0].ring
        #grade = items[0].grade
        grade = none_uniq([item.grade for item in items])
        n = sum(item.n for item in items)
        name = "("+"+".join(item.name for item in items)+")"
        Space.__init__(space, ring, n, grade, name)
        space.items = items
        cls.cache[key] = space
        assert len(items)>1
        return space

    def __init__(self, ring, *args, **kw):
        pass

    @classmethod
    def promote(cls, space):
        assert isinstance(space, Space)
        if isinstance(space, AddSpace):
            return space
        return AddSpace(space.ring, space)

    @property
    def dual(self):
        if self._dual is None:
            # distribute dual's
            self._dual = AddSpace(self.ring, *[space.dual for space in self.items])
        return self._dual

    def get_add_items(self):
        return list(self.items)

    def nullitor(self, inverse=False):
        "remove all zero (null) summands"
        items = self.items
        src = self
        tgt = [item for item in self.items if item.n]
        if len(tgt)>1:
            tgt = AddSpace(self.ring, *tgt)
        elif len(tgt):
            tgt = tgt[0]
        else:
            assert 0, "wah?"
        A = elim.identity(self.ring, self.n)
        if inverse:
            src, tgt = tgt, src # the same A works
        return Lin(tgt, src, A)

    def get_swap(self, perm=(1, 0)):
        perm = tuple(perm)
        items = self.items
        assert len(perm) == len(items), ("len(%s)!=%d"%(perm, len(items)))
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

    get_add_swap = get_swap

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

    def get_add_slice(self, idx):
        i = 0
        for jdx, item in enumerate(self.items):
            if jdx == idx:
                return slice(i, i+item.n)
            i += item.n
        assert 0, "index %d out of bounds"%(idx,)

    def get_slice(self, space):
        #assert 0, "only works for list of unique spaces!"
        found = None
        i = 0
        for idx, item in enumerate(self.items):
            if item == space:
                assert found is None, "ambiguous slice!"
                found = slice(i, i+item.n)
            i += item.n
        assert found is not None, "space %s not found in %s"%(space, self)
        return found

    def get_normal(self, N=None, K=None, inverse=False, force=False):
        # remove null summands, and recurse
        assert N is None or N.n == 0
        assert K is None or K.n == 1
        ring = self.ring
        spaces = list(self.items)
        assert len(spaces)>1, str(self)
        # depth-first recurse
        fs = [space.get_normal(N, K, force=force) for space in spaces]
        f = reduce(Lin.direct_sum, fs)
        assert isinstance(f.tgt, AddSpace)
        spaces = list(f.tgt.items)
        idx = 0
        while idx < len(spaces):
            if force and spaces[idx].n == 0 or spaces[idx]==N:
                spaces.pop(idx)
                tgt, src = AddSpace(ring, *spaces, N=N), f.tgt
                g = Lin(tgt, src, elim.identity(ring, tgt.n))
                f = g*f
                done = False
            else:
                idx += 1
        if inverse:
            f = f.transpose() # permutation matrix
        return f


class MulSpace(Space):
    "tensor product of vector spaces"

    cache = {} # XXX use https://docs.python.org/3/library/weakref.html
    def __new__(cls, ring, *_items, K=None):
        assert isinstance(ring, element.Type), ring.__class__
        assert _items or K is not None
        if not _items:
            return K
        if len(_items)==1:
            return _items[0] 
        items = []
        for item in _items:
            assert item.ring == ring
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
        #ring = items[0].ring
        #grade = sum(item.grade for item in items)
        grade = reduce(none_add, [item.grade for item in items])
        n = reduce(mul, [item.n for item in items], 1)
        name = "@".join(item.name for item in items)
        Space.__init__(space, ring, n, grade, name=name)
        space.items = items
        cls.cache[key] = space
        return space

    def __init__(self, ring, *args, **kw):
        pass

    @classmethod
    def promote(cls, space):
        assert isinstance(space, Space)
        if isinstance(space, MulSpace):
            return space
        return MulSpace(space.ring, space)

    @property
    def dual(self):
        if self._dual is None:
            # distribute dual's
            self._dual = MulSpace(self.ring, *[space.dual for space in self.items])
        return self._dual

    def get_normal(self, N=None, K=None, inverse=False, force=False):
        #print("\nMulSpace.get_normal", self.name)
        assert N is None or N.n == 0
        assert K is None or K.n == 1
        ring = self.ring

        # depth-first recurse
        spaces = list(self.items)
        assert len(spaces)>1, str(self)
        fs = [space.get_normal(N, K, force=force) for space in spaces]
        f = reduce(matmul, fs)
        assert isinstance(f.tgt, MulSpace)
        assert f.src == self
        spaces = list(f.tgt.items)

        #print("spaces:")
        #print([space.name for space in spaces])

        # first deal with N, K factors (*)
        idx = 0
        while idx < len(spaces):
            space = spaces[idx]
            if force and space.n == 0 or space == N:
                # kills everything
                tgt = N
                g = Lin.zero(tgt, f.tgt)
                f = g*f
                if inverse:
                    f = f.transpose() # permutation matrix
                return f # <-------------- return
            elif force and space.n == 1 or space == K:
                spaces.pop(idx)
            else:
                idx += 1

        #print("spaces:", [space.name for space in spaces])
        #print("self.items:", [space.name for space in self.items])

        if len(spaces) < len(f.tgt.items):
            tgt, src = MulSpace(ring, *spaces, K=K), f.tgt
            g = Lin(tgt, src, elim.identity(ring, tgt.n))
            f = g*f
            #print(tgt)

        #print("f.tgt", f.tgt)
        if spaces:
            gs = [space.get_normal(N, K, force=force) for space in spaces] # <--------- recurse
            spaces = [g.tgt for g in gs]
            g = reduce(matmul, gs)
            assert f.tgt == g.src
            f = g*f
            # go back to (*) ?

        # now distribute over AddSpace's
        for idx, space in enumerate(spaces):
            if isinstance(space, AddSpace):
                break
        else:
            if inverse:
                f = f.transpose() # permutation matrix
            return f # <---------------- return

        if idx+1 < len(spaces):
            head = spaces[:idx]
            lhs = spaces[idx]
            rhs = MulSpace(ring, *spaces[idx+1:], K=K)
            g = Lin.right_distributor(lhs, rhs)
            if head:
                head = MulSpace(ring, *head)
                g = head.identity() @ g
            f = g*f
        elif len(spaces) > 1:
            lhs = MulSpace(ring, *spaces[:idx], K=K)
            rhs = spaces[idx]
            #print("f =", f.homstr())
            assert f.tgt == lhs @ rhs
            #print("left_distributor", lhs, rhs)
            g = Lin.left_distributor(lhs, rhs)
            f = g*f

        g = f.tgt.get_normal(N, K, force=force) # <-------------- recurse
        f = g*f
        if inverse:
            f = f.transpose() # permutation matrix
        return f

    def unitor(self, inverse=False):
        "remove all tensor units"
        items = self.items
        src = self
        tgt = [item for item in self.items if item.n!=1]
        if len(tgt)>1:
            tgt = MulSpace(self.ring, *tgt)
        elif len(tgt):
            tgt = tgt[0]
        else:
            assert 0, "wah?"
        A = elim.identity(self.ring, self.n)
        if inverse:
            src, tgt = tgt, src # the same A works
        return Lin(tgt, src, A)

    def annihilator(self, inverse=False):
        "zero factors _annihilate the tensor"
        items = self.items
        zero = None
        for space in items:
            if space.n == 0:
                zero = space
        if zero is None:
            return self.identity()
        src = self
        tgt = zero
        if inverse:
            src, tgt = tgt, src # the same A works
        return Lin(tgt, src)

    def right_distributor(self, inverse=False):
        # These turn out to be identity matrices (on the nose).
        assert len(self.items) == 2, "not implemented: %s"%(self,)
        return Lin.right_distributor(*self.items, inverse=inverse)
#        src = self
#        left, right = self.items
#        assert isinstance(left, AddSpace)
#        tgt = AddSpace(self.ring, *[l@right for l in left.items])
#        if inverse:
#            tgt, src = src, tgt
#        A = elim.identity(self.ring, self.n)
#        return Lin(tgt, src, A)

    def left_distributor(self, inverse=False):
        # here we cheat: use right_distributor and swap's
        assert len(self.items) == 2, "not implemented: %s"%(self,)
        return Lin.left_distributor(*self.items, inverse=inverse)
#        S = self.get_swap((1, 0))
#        tgt = S.tgt
#        L = tgt.right_distributor()
#        assert isinstance(L.tgt, AddSpace)
#        swaps = []
#        for i, sumand in enumerate(L.tgt.items):
#            assert isinstance(sumand, MulSpace)
#            n = len(sumand.items)
#            assert n>=2
#            perm = [n-1]+list(range(n-1))
#            swap = sumand.get_swap(perm)
#            swaps.append(swap)
#        S1 = reduce(Lin.direct_sum, swaps)
#        R = S1*L*S
#        if inverse:
#            R = R.transpose() # it's a permutation matrix
#        return R

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

    get_mul_swap = get_swap


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
        if type(A) is list:
            A = elim.array(A)
        #for row in A:
        #    assert None not in row
        assert A.shape == (tgt.n, src.n), "%s != %s" % ( A.shape , (tgt.n, src.n) )
        self.hom = (tgt, src) # yes it's backwards, just like shape is.
        self.shape = A.shape
        self.A = A.copy()

    @classmethod
    def zero(cls, tgt, src):
        assert isinstance(tgt, Space), tgt
        assert isinstance(src, Space), src
        assert tgt.ring is src.ring
        ring = tgt.ring
        A = elim.zeros(ring, tgt.n, src.n)
        return Lin(tgt, src, A)

    @classmethod
    def rand(cls, tgt, src, a=1, b=1):
        assert isinstance(tgt, Space)
        assert isinstance(src, Space)
        assert tgt.ring is src.ring
        ring = tgt.ring
        A = elim.rand(ring, tgt.n, src.n, a, b)
        return Lin(tgt, src, A)

    @classmethod
    def iso(cls, tgt, src):
        assert isinstance(tgt, Space)
        assert isinstance(src, Space)
        assert tgt.ring is src.ring
        assert tgt.n == src.n
        ring = tgt.ring
        A = elim.identity(ring, tgt.n)
        return Lin(tgt, src, A)

    def is_zero(self):
        B = Lin.zero(self.tgt, self.src)
        return eq(self.A, B.A)

    def is_identity(self):
        assert self.tgt == self.src, "wah?"
        B = self.tgt.identity()
        return eq(self.A, B.A)

    def __str__(self):
        return shortstr(self.A)

    def __repr__(self):
        #s = str(self.A).replace("\n", " ")
        s = shortstr(self.A)
        return "Lin(%s, %s, %s)"%(self.tgt, self.src, s)

    def homstr(self):
        return "%s<---%s"%(self.tgt.name, self.src.name)

    def weak_eq(self, other):
        assert self.A.shape == other.A.shape
        return eq(self.A, other.A)

    def __eq__(self, other):
        #assert self.hom == other.hom, "%s != %s" % (self.hom, other.hom) # too strict ?
        if self.hom != other.hom:
            return True
        return eq(self.A, other.A)

    def __ne__(self, other):
        #assert self.hom == other.hom # too strict ...
        if self.hom != other.hom:
            return True
        return not eq(self.A, other.A)

    def __getitem__(self, idx):
        return self.A[idx]

    def __setitem__(self, idx, val):
        self.A[idx] = val

    # XXX should be direct sum here, consistent with Space.__add__ ?
    def __add__(self, other):
        assert self.hom == other.hom
        A = self.A + other.A
        return Lin(*self.hom, A)
    add = __add__

    def __sub__(self, other):
        assert self.hom == other.hom
        A = self.A - other.A
        return Lin(*self.hom, A)

    def __mul__(self, other):
        assert other.tgt == self.src, "%s != %s" % (other.tgt.name, self.src.name)
        A = dot(self.ring, self.A, other.A)
        return Lin(self.tgt, other.src, A)

    def __rshift__(self, other):
        return other * self

    def __rmul__(self, r):
        r = self.ring.promote(r)
        A = r*self.A
        return Lin(self.tgt, self.src, A)

    def __pow__(self, n):
        assert self.tgt == self.src, "not an endomorphism"
        assert n>=0
        if n==0:
            return self.tgt.identity()
        if n==1:
            return self
        op = self
        while n>1:
            op = op*self
            n -= 1
        return op

    def __neg__(self):
        A = -self.A
        return Lin(self.tgt, self.src, A)

    def __matmul__(self, other):
        src = self.src @ other.src
        tgt = self.tgt @ other.tgt
        ring = self.ring
        A, B = self.A, other.A
        if 0 in A.shape or 0 in B.shape:
            C = elim.zeros(ring, A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])
        else:
            #print("kron", A.shape, B.shape)
            C = numpy.kron(A, B)
            #print("\t", C.shape)
        return Lin(tgt, src, C)

    def sum(self, dim):
        A = self.A.sum(dim)
        assert len(A.shape) == 1
        K = Space(self.ring, 1, self.hom[dim].grade, "K")
        hom = list(self.hom)
        hom[dim] = K
        shape = list(self.shape)
        shape[dim] = 1
        A.shape = tuple(shape)
        return Lin(*hom, A)

    def intersect(self, other):
        assert self.src == other.src
        A = elim.intersect(self.A, other.A)
        tgt = Space(self.ring, len(A), self.tgt.grade)
        return Lin(tgt, self.src, A)

    @staticmethod
    def get_mulswap(ring, spaces, perm):
        # This is where strict _associative tensor product gets annoying
        assert len(spaces) == len(perm)
        if not len(spaces):
            return MulSpace(ring).identity()
        perms = []
        idx = 0
        for space in spaces:
            if isinstance(space, MulSpace):
                n = len(space.items)
                prm = tuple(i+idx for i in range(n))
                idx += n
            else:
                prm = (idx,)
                idx += 1
            perms.append(prm)
        perms = [perms[i] for i in perm]
        perm = reduce(add, perms, ())
        src = MulSpace(ring, *spaces)
        assert len(perm) == len(src.items)
        tgt = MulSpace(ring, *[src.items[i] for i in perm])
        f = src.get_swap(perm)
        assert f.tgt == tgt
        return f

    @staticmethod
    def right_distributor(lhs, rhs, inverse=False):
        "A@C+B@C <---- (A+B)@C"
        # These turn out to be identity matrices (on the nose).
        if not isinstance(lhs, AddSpace):
            # There's nothing to distribute
            return (lhs@rhs).identity()
        src = lhs @ rhs
        tgt = AddSpace(lhs.ring, *[l@rhs for l in lhs.items])
        if inverse:
            tgt, src = src, tgt
        A = elim.identity(src.ring, src.n)
        return Lin(tgt, src, A)

    @staticmethod
    def left_distributor(lhs, rhs, inverse=False):
        "A@B+A@C <---- A@(B+C)"
        # here we cheat: use right_distributor and swap's
        if not isinstance(rhs, AddSpace):
            # There's nothing to distribute
            return (lhs@rhs).identity()
        src = lhs @ rhs  # A@(B+C)
        S = Lin.get_mulswap(lhs.ring, [lhs, rhs], (1, 0)) # (B+C)@A <--- A@(B+C)
        assert S.src == src
        L = Lin.right_distributor(rhs, lhs)
        assert isinstance(L.tgt, AddSpace)
        offset = len(lhs.items) if isinstance(lhs, MulSpace) else 1
        swaps = []
        for i, sumand in enumerate(L.tgt.items):
            assert isinstance(sumand, MulSpace)
            n = len(sumand.items)
            assert offset<n
            perm = list(range(n))
            perm = perm[n-offset:]+perm[:n-offset]
            swap = sumand.get_swap(perm)
            swaps.append(swap)
        S1 = reduce(Lin.direct_sum, swaps)
        R = S1*L*S
        if inverse:
            R = R.transpose() # it's a permutation matrix
        return R

    @staticmethod
    def unit(K, U): # method of K ?
        assert isinstance(K, Space)
        assert isinstance(U, Space)
        ring = K.ring
        assert K.n == 1, "not tensor identity"
        tgt = U.dual @ U
        src = K
        #A = elim.zeros(ring, tgt.n, src.n)
        A = elim.identity(ring, U.n)
        A.shape = (tgt.n, src.n)
        lin = Lin(tgt, src, A)
        return lin

    @staticmethod
    def counit(K, U): # method of K ?
        assert isinstance(K, Space)
        assert isinstance(U, Space)
        ring = K.ring
        assert K.n == 1, "not tensor identity"
        tgt = K
        src = U @ U.dual
        A = elim.identity(ring, U.n)
        A.shape = (tgt.n, src.n)
        lin = Lin(tgt, src, A)
        return lin

    @staticmethod
    def diagonal(src):
        assert isinstance(src, Space)
        n = src.n
        tgt = src + src
        ring = src.ring
        one = ring.one
        A = elim.zeros(ring, tgt.n, src.n)
        for i in range(n):
           A[i, i] = one
           A[i+n, i] = one
        lin = Lin(tgt, src, A)
        return lin

    @staticmethod
    def codiagonal(tgt):
        lin = Lin.diagonal(tgt)
        return lin.transpose()

    @staticmethod
    def copy(src):
        assert isinstance(src, Space)
        n = src.n
        tgt = src @ src
        ring = src.ring
        one = ring.one
        A = elim.zeros(ring, tgt.n, src.n)
        for i in range(n):
         for j in range(n):
          for k in range(n):
            if i==j==k:
                A[i + n*j, k] = 1
        lin = Lin(tgt, src, A)
        return lin

    def transpose(self):
        src, tgt = self.tgt, self.src
        A = self.A.transpose()
        return Lin(tgt, src, A)

    def rank(self):
        return elim.rank(self.ring, self.A)

    def row_reduce(self):
        A = elim.row_reduce(self.ring, self.A, True)
        tgt = Space(self.ring, len(A))
        return Lin(tgt, self.src, A)

    def direct_sum(self, other):
        tgt = self.tgt + other.tgt
        src = self.src + other.src
        f = Lin.zero(tgt, src)
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

    def kernel(self, other=None):
        assert other is None, "todo"
        K = elim.kernel(self.ring, self.A)
        src = Space(self.ring, K.shape[1])
        K = Lin(self.src, src, K)
        return K

    def cokernel(self, other=None):
        assert other is None, "todo"
        K = elim.cokernel(self.ring, self.A)
        tgt = Space(self.ring, K.shape[0])
        lin = Lin(tgt, self.tgt, K)
        return lin



# ------------------------------------------------------------

# ------------------------ testing      ----------------------


def test_structure():
    ring = element.Q
    N = Space(ring, 0, 0, 'N') # null (direct_sum unit)
    K = Space(ring, 1, 0, 'K') # field (tensor unit)
    U = Space(ring, 2, 0, 'U')
    V = Space(ring, 3, 0, 'V')
    W = Space(ring, 5, 0, 'W')
    X = Space(ring, 2, 0, 'X')
    Y = Space(ring, 1, 0, 'Y')

    f = (V+N).nullitor()
    fi = (V+N).nullitor(inverse=True)
    assert f*fi == V.identity()
    assert fi*f == (V+N).identity()

    f = (V@K).unitor()
    fi = (V@K).unitor(inverse=True)
    assert f*fi == V.identity()
    assert fi*f == (V@K).identity()

    f = (V@N).annihilator()
    fi = (V@N).annihilator(inverse=True)
    assert f*fi == N.identity()
    assert fi*f == (V@N).identity()

    # These turn out to be identity matrices (on the nose).
    tgt, src = V@U + W@U + K@U, (V+W+K)@U
    f = src.right_distributor()
    assert f.tgt == tgt
    assert f.src == src
    fi = src.right_distributor(inverse=True)
    assert fi.tgt == src
    assert fi.src == tgt
    assert f*fi == tgt.identity()
    assert fi*f == src.identity()

    tgt, src = U@V@X + U@W@Y, U@(V@X+W@Y)
    f = src.left_distributor()
    assert f.homstr() == "(U@V@X+U@W@Y)<---U@(V@X+W@Y)"
    assert f.tgt == tgt
    fi = src.left_distributor(inverse=True)
    assert f*fi == tgt.identity()
    assert fi*f == src.identity()

    assert U.dual == U.dual
    assert U.dual.dual == U
    assert (U@V).dual == U.dual @ V.dual
    assert (U+V).dual == U.dual + V.dual

    cup = Lin.unit(K, W)
    cap = Lin.counit(K, W)
    assert cup.homstr() == "~W@W<---K"
    assert cap.homstr() == "K<---W@~W"

    swap = cup.tgt.get_swap()
    bubble = cap * swap * cup
    assert bubble[0,0] == W.n

    i_W = W.identity()
    r_W = (W@K).unitor() # right unitor
    l_W = (K@W).unitor() # left unitor

    lhs = l_W * (cap @ i_W) * (i_W @ cup)
    rhs = r_W
    assert lhs == rhs

    i_dW = W.dual.identity()
    r_dW = (W.dual@K).unitor() # right unitor
    l_dW = (K@W.dual).unitor() # left unitor
    lhs = r_dW * (i_dW @ cap) * (cup @ i_dW)
    rhs = l_dW
    assert lhs == rhs

    assert AddSpace(ring, W) == W
    assert MulSpace(ring, W) == W
    assert AddSpace(ring, N=N) == N
    assert MulSpace(ring, K=K) == K

    for (tgt, src) in [
            (U, K@U@K),
            (N, K@N@K),
            (N, N+N+N),
            (K+K, N+N+N+K+K),
            (U+U, (K+K)@U),
            (U@U+U@U, U@(K+K)@U),
            (U+U, K@(U+U)@K),
            (U@U, (U@K@K+U@N@N+V@K@N+V@N@K)@(U@K@K+U@N@N+V@K@N+V@N@K)),
            (U+V+U, (K@U@K+V+U)@K),
            (U@U, (K@U@U@K+K@U@U@N+K@U@U@N+N@U@U@K+N@U@U@N+N@U@U@N)@K),
        ]:
        #print("src =", src.name)
        try:
            f = src.get_normal(N, K)
        except AssertionError:
            print("get_normal %s failed" % (src.name))
            raise
        assert f.src == src, "%s should be %s<---%s"%(f.homstr(), tgt.name, src.name)
        assert f.tgt == tgt, "%s should be %s<---%s"%(f.homstr(), tgt.name, src.name)

        fi = src.get_normal(N, K, inverse=True)
        assert fi.src == tgt
        assert fi.tgt == src
        assert f*fi == tgt.identity()
        assert fi*f == src.identity()


def test_young():
    # code ripped from qu.py
    from bruhat.rep import Young

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
    from bruhat.rep import Young
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
    from bruhat.rep import Young
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


def chain_product(ring, fs):
    N = len(fs)

    grades = [[] for i in range(N+1)]
    #slookup = {}
    #tlookup = {}
    for items in cross([f.hom for f in fs]):
        space = MulSpace(ring, *items)
        grades[space.grade].append(space)
    grades = [AddSpace(ring, *spaces) for spaces in grades]
    #print([space.name for space in grades])

    srcs = [[] for i in range(N+1)]
    tgts = [[] for i in range(N+1)]
    items = [(f.src.identity(), f.tgt.identity()) for f in fs]
    for idx in range(N):
        _items = list(items)
        _items[idx] = (fs[idx],)
        for lins in cross(_items):
            lin = reduce(matmul, lins)
            srcs[lin.src.grade].append(lin)
            tgts[lin.tgt.grade].append(lin)
            #print([lin.homstr() for lin in lins])
            #print(lin.homstr())
            #print()
    
    #print([len(grade) for grade in srcs])

    chain = []
    for i in range(N):
        tgt = grades[i]
        src = grades[i+1]
        #print([space.name for space in tgt.get_add_items()])
        #print([space.name for space in src.get_add_items()])
        A = elim.zeros(ring, tgt.n, src.n)
        for lin in tgts[i]:
            #print("\t", lin.homstr())
            cols = src.get_slice(lin.src)
            rows = tgt.get_slice(lin.tgt)
            assert numpy.alltrue(A[rows, cols] == ring.zero)
            A[rows, cols] = lin.A
        lin = Lin(tgt, src, A)
        #print(shortstr(A))
        #print()
        chain.append(lin)

    return chain


def chain_perm(ring, fs, perm):
    #print("chain_perm")
            
    N = len(fs)
    assert N == len(perm)
    gs = [fs[i] for i in perm]

    def product(fs):
        grades = [[] for i in range(N+1)]
        for items in cross([f.hom for f in fs]):
            space = MulSpace(ring, *items)
            grades[space.grade].append(space)
        grades = [AddSpace(ring, *spaces) for spaces in grades]
        return grades

    sgrades = product(fs)
    tgrades = product(gs)

    cmap = []
    for i in range(N+1):
        tgt, src = tgrades[i], sgrades[i]
        #print(tgt.name, "<---", src.name)
        tgts = tgt.get_add_items()
        srcs = src.get_add_items()
        lins = [src.get_mul_swap(perm) for src in srcs]

        f = reduce(Lin.direct_sum, lins)
        #print(f.homstr())
        assert f.src == src

        swap = [None]*len(lins)
        for i, lin in enumerate(lins):
            assert tgts.count(lin.tgt) == 1 # careful !
            idx = tgts.index(lin.tgt)
            swap[idx] = i
        #print("swap:", swap)
        swap = f.tgt.get_add_swap(swap)
        f = swap * f

        #print(f.homstr())

        assert f.hom == (tgt, src)
        cmap.append(f)

        #print()


    return cmap


def test_chain_map(tgt, src, cmap):
    n = len(tgt)
    assert n == len(src) == len(cmap)-1
    #for i in range(n):
    #    print(tgt[i].homstr())
    #    print('\t', cmap[i].homstr())
    
    for i in range(n):
        lhs = cmap[i] * src[i]
        rhs = tgt[i] * cmap[i+1]
        assert lhs == rhs


def test_chain_symmetry(ring, lins):
    N = len(lins)
    src = chain_product(ring, lins)
    #print([lin.homstr() for lin in src])
    for perm in allperms(list(range(N))):
        #print("perm:", perm)
        tgt = chain_product(ring, [lins[i] for i in perm])
        #print([lin.homstr() for lin in tgt])
        cmap = chain_perm(ring, lins, perm)
        test_chain_map(tgt, src, cmap)


def cmap_product(ring, cmaps):
    for cmap in cmaps:
        assert len(cmap) == 2
        f0, f1 = cmap
        #print(f0.homstr(), f1.homstr())
    N = len(cmaps)

#    srcs = [[] for i in range(N+1)]
#    for items in cross([(cmap[0].src, cmap[1].src) for cmap in cmaps]):
#        space = MulSpace(ring, *items)
#        srcs[space.grade].append(space)
#    srcs = [AddSpace(ring, *spaces) for spaces in srcs]
#    #print([space.name for space in srcs])
#
#    tgts = [[] for i in range(N+1)]
#    for items in cross([(cmap[0].tgt, cmap[1].tgt) for cmap in cmaps]):
#        space = MulSpace(ring, *items)
#        tgts[space.grade].append(space)
#    tgts = [AddSpace(ring, *spaces) for spaces in tgts]
#    #print([space.name for space in tgts])

    linss = [[] for i in range(N+1)]
    for items in cross(cmaps):
        #print([f.homstr() for f in items])
        f = reduce(matmul, items)
        #print(f.homstr())
        assert f.src.grade == f.tgt.grade
        linss[f.src.grade].append(f)
    chain = [reduce(Lin.direct_sum, lins) for lins in linss]
    #print([lin.homstr() for lin in chain])
    return chain


def test_chain_product():

    p = 2
    ring = element.FiniteField(p)

    C1 = Space(ring, 1, 1, "C_1")
    C0 = Space(ring, 1, 0, "C_0")
    f = Lin(C0, C1, [[1]])

    chain = chain_product(ring, [f, f])

    D1 = Space(ring, 1, 1, "D_1")
    D0 = Space(ring, 1, 0, "D_0")
    g = Lin(D0, D1, [[1]])

    E1 = Space(ring, 1, 1, "E_1")
    E0 = Space(ring, 1, 0, "E_0")
    h = Lin(E0, E1, [[1]])

    test_chain_symmetry(ring, [f, g])
    test_chain_symmetry(ring, [f, g, h])
    test_chain_symmetry(ring, [f, f])
    test_chain_symmetry(ring, [f, f, h])
    test_chain_symmetry(ring, [f, f, f])


    C1 = Space(ring, 3, 1, "C_1")
    C0 = Space(ring, 2, 0, "C_0")
    f = Lin(C0, C1, [[1,1,0],[0,1,1]])

    chain = chain_product(ring, [f, f, f])
    #test_chain_symmetry(ring, [f, f, f])

    C1 = Space(ring, 3, 1, "C_1")
    C0 = Space(ring, 3, 0, "C_0")
    f = Lin(C0, C1, [[1,1,0],[0,1,1],[1,0,1]])

    D1 = Space(ring, 3, 1, "D_1")
    D0 = Space(ring, 3, 0, "D_0")
    g = Lin(D0, D1, [[1,1,0],[0,1,1],[1,0,1]])

    cmap = [
        Lin(D0, C0, [[0,1,0],[0,0,1],[1,0,0]]),
        Lin(D1, C1, [[0,1,0],[0,0,1],[1,0,0]]),
    ]

    test_chain_map([g], [f], cmap)

    ident = [
        Lin(D0, C0, [[1,0,0],[0,1,0],[0,0,1]]),
        Lin(D1, C1, [[1,0,0],[0,1,0],[0,0,1]]),
    ]
    test_chain_map([g], [f], ident)

    src = chain_product(ring, [f, f, f])
    tgt = chain_product(ring, [g, g, g])

    c3 = cmap_product(ring, [cmap, cmap, cmap])
    test_chain_map(tgt, src, c3)

    c3 = cmap_product(ring, [cmap, cmap, ident])
    test_chain_map(tgt, src, c3)


def test_kagome():

    p = 2
    ring = element.FiniteField(p)

    C1 = Space(ring, 3, 1, "C_1")
    C0 = Space(ring, 3, 0, "C_0")
    f = Lin(C0, C1, [[1,1,0],[0,1,1],[1,0,1]])

    chain = chain_product(ring, [f]*3)

    for bdy in chain:
        print(bdy.homstr())


def test_all():

    test_young()
    test_structure()
    test_super()
    test_symmetric_square()
    test(element.Q)

    p = argv.get("p", 3)
    ring = element.FiniteField(p)
    test(ring)



if __name__ == "__main__":


    fn = argv.next() or "test_all"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK")


