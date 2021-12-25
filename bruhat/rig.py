#!/usr/bin/env python3

"""
Kapranov-Voevodsky 2-vector spaces

Fail: too complicated, and i think it's also wrong.

"""

import operator
from functools import reduce
from random import randint

import numpy # for the matrix indexing and ops

from bruhat import element
from bruhat.element import GenericElement, Type, Keyed
Gen = GenericElement
from bruhat import elim
from bruhat.chain import Space, Lin, shortstr

#
# unwrapped has ob_ prefix
#

class Rig(Keyed, Type):
    def __init__(self, ob_zero, ob_one, ob_add, ob_mul):
        key = (ob_zero, ob_one, ob_add, ob_mul)
        Keyed.__init__(self, key)
        Type.__init__(self)
        self.zero = Gen(ob_zero, self)
        self.one = Gen(ob_one, self)
        self.ob_add = ob_add
        self.ob_mul = ob_mul

    def promote(self, ob):
        if isinstance(ob, Gen):
            assert ob.tp == self
            return ob
        gen = Gen(ob, self)
        return gen

    def add(self, a, b):
        value = self.ob_add(a.value, b.value)
        c = Gen(value, self)
        return c

    def mul(self, a, b):
        value = self.ob_mul(a.value, b.value)
        c = Gen(value, self)
        return c


class Lin2(object):
    "A 2-morphism in 2Vect"
    def __init__(self, tgt, src, A=None):
        assert isinstance(tgt, Lin)
        assert isinstance(src, Lin)
        assert tgt.ring is src.ring
        self.ring = tgt.ring
        self.tgt = tgt
        self.src = src
        shape = src.shape
        assert shape == tgt.shape
        rows, cols = shape
        if A is None:
            A = [[Lin(tgt[i,j].value, src[i,j].value) # note unwrapping here 
                for j in range(cols)] for i in range(rows)]
        A = elim.array(A)
        for j in range(cols):
            for i in range(rows):
                #print(i, j)
                #print( A[i, j].tgt )
                #print(tgt[i, j].value )
                assert A[i, j].tgt == tgt[i, j].value, (A[i, j].tgt, tgt[i, j].value) # unwrapped
                assert A[i, j].src == src[i, j].value # unwrapped
        self.A = A
        self.shape = shape

    def __str__(self):
        return str(self.A)

    def __getitem__(self, idx):
        #i, j = idx
        #rows, cols = self.shape
        #assert 0<=i<rows
        #assert 0<=j<cols
        return self.A[idx]

    def __mul__(left, right):
        "(vertical) composition of 2-morphism's"
        assert right.tgt == left.src
        A = left.A * right.A
        return Lin2(left.tgt, right.src, A)

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
        return Lin2(tgt, src, lins)

    def __eq__(self, other):
        assert self.tgt == other.tgt
        assert self.src == other.src
        return numpy.alltrue(self.A == other.A)

    def weak_eq(self, other):
        assert self.tgt.shape == other.tgt.shape
        assert self.src.shape == other.src.shape
        shape = self.tgt.shape 
        for i in range(shape[0]):
          for j in range(shape[1]):
            if not self[i,j].weak_eq(other[i,j]):
                return False
        return True

    @classmethod
    def identity(cls, ring, rig, lin):
        r"""
            |
            |
         V  *   U
            |
            |
           lin
        """
        assert lin.ring == rig
        src = tgt = lin
        V, U = lin.hom
        lins = [[
            Lin(
                lin[i,j].value, 
                lin[i,j].value, 
                elim.identity(ring, lin[i,j].value.n)) # unwrapped
            for j in range(U.n)] for i in range(V.n)]
        i_2 = Lin2(tgt, src, lins)
        return i_2

    @classmethod
    def left_unitor(cls, ring, rig, lin):
        r"""
            |
            |
        V   *   U
          . |
         .  |
           lin
        """
        assert lin.ring == rig
        tgt = lin
        V, U = lin.hom
        one = V.identity()
        lins = [[
            Lin(
                lin[i,j].value, 
                lin[i,j].value, 
                elim.identity(ring, lin[i,j].value.n)) # unwrapped
            for j in range(U.n)] for i in range(V.n)]
        i_2 = Lin2(tgt, src, lins)
        return i_2


    @classmethod
    def unit(cls, ring, rig, V1, V):
        r"""
            |   |
            | V1|
             \  /
              \/
               *
             V .  V
               .
        """
        assert isinstance(V1, Space)
        assert V1.ring == rig
        assert V1.n == 1
        assert isinstance(V, Space)
        assert V.ring == rig
        n = V.n
        src = V.identity()
        add = Lin(V1, V, [[rig.one for i in range(n)]])
        copy = Lin(V, V1, [[rig.one] for i in range(n)])
        tgt = copy * add
        
        arrays = [elim.zeros(ring, 1, 0), elim.identity(ring, 1)]
        lins = [[
            Lin(tgt[i,j].value, src[i,j].value, arrays[int(i==j)]) # unwrapped
            for j in range(n)] for i in range(n)]
        unit = Lin2(tgt, src, lins)
        return unit

    @classmethod
    def counit(cls, ring, rig, V1, V):
        r"""
               .
             V .  V
               *
              /\
             /  \
            | V1|
            |   |
        """
        assert isinstance(V1, Space)
        assert V1.ring == rig
        assert V1.n == 1
        assert isinstance(V, Space)
        assert V.ring == rig
        n = V.n
        add = Lin(V1, V, [[rig.one for i in range(n)]])
        copy = Lin(V, V1, [[rig.one] for i in range(n)])
        src = copy * add
        tgt = V.identity()
        
        arrays = [elim.zeros(ring, 0, 1), elim.identity(ring, 1)]
        lins = [[
            Lin(tgt[i,j].value, src[i,j].value, arrays[int(i==j)]) # unwrapped
            for j in range(n)] for i in range(n)]
        counit = Lin2(tgt, src, lins)
        return counit


def main():

    ring = element.Q
    N = Space(ring, 0, 0, 'N') # null (direct_sum unit)
    K = Space(ring, 1, 0, 'K') # field (tensor unit)

    # N, K are the zero and one of a rig.
    # We "wrap" these Space objects into a GenericElement.
    rig = Rig(N, K, Space.__add__, Space.__matmul__)

    # An Object is a Space over a rig
    V0 = Space(rig, 0, name='V0')
    V1 = Space(rig, 1, name='V1')
    #V2 = Space(rig, 2, name='V2')
    V2 = V1 + V1
    V3 = V1 + V1 + V1

    # 1-morphism's are Lin's over this rig
    Vi = V1.identity()
    #print(V2.identity())
    copy = Lin(V2, V1, [[rig.one],[rig.one]])
    delete = Lin(V0, V1)
    add = Lin(V1, V2, copy.A.transpose())
    zero = Lin(V1, V0)

    assert copy * zero == Lin(V2, V0)

    # 2-morphism's are Lin2's.
    # These are 2d arrays of Lin's.
    zero_2 = Lin2(Vi, Vi, [[Lin(K, K)]])

    # The left unitorator
    tgt_1, src_1 = add * zero.direct_sum(Vi), Vi
    tgt = tgt_1[0,0].value  # unwrapped
    src = src_1[0,0].value  # unwrapped
    assert tgt == K@N + K@K
    assert src == K

    lin = Lin(tgt, src, elim.identity(ring, 1))
    m_2 = Lin2(tgt_1, src_1, [[lin]])

    # The "bialgebrator" should be a 2-morphism from src to tgt
    a = copy.direct_sum(copy)
    b = a.tgt.get_swap((0, 2, 1, 3))
    c = add.direct_sum(add)
    tgt = (c*b*a)
    src = (copy * add)
    arrays = [elim.zeros(ring, 1, 1), elim.identity(ring, 1)]
    lins = [[
        Lin(tgt[i,j].value, src[i,j].value, arrays[int(i==j)]) # unwrapped
        for j in range(2)] for i in range(2)]
    bialgebrator = Lin2(tgt, src, lins)

    # ----------------------------

    def rand_1cell(ring, rig, m, n, maxdim=5, name="A"):
        A = elim.zeros(rig, m, n)
        for i in range(m):
          for j in range(n):
            dim = randint(0, 5)
            item = Space(ring, dim, 0, "%s_%d%d"%(name,i,j))
            A[i, j] = rig.promote(item)
        #print(reduce, add, [rig.one]*m, rig.zero)
        tgt = reduce(operator.add, [rig.one]*m)
        src = reduce(operator.add, [rig.one]*n)
        return Lin(tgt.value, src.value, A) # unwrapped

    V = V3
    Vi = V.identity()
    m_1 = rand_1cell(ring, rig, 3, 2)
    #print(type(Vi[0,0]))
    #print(type(m_1[0,0]))
    #print(type(Vi.tgt), type(m_1.tgt))
    print(m_1)

    def rand_2cell(ring, rig, m, n, maxim=5):
        tgt_1 = rand_1cell(ring, rig, m, n, name="tgt")
        src_1 = rand_1cell(ring, rig, m, n, name="src")
        lins = [[
            Lin(
                tgt_1[i,j].value, 
                src_1[i,j].value, 
                elim.rand(ring, tgt_1[i,j].value.n, 
                src_1[i,j].value.n))
            for j in range(n)]
            for i in range(m)]
        m_2 = Lin2(tgt, src, lins)
        return m_2

    m_2 = rand_2cell(ring, rig, 3, 2)
    print(m_2)

    # ----------------------------

    V = V1+V1+V1
    unit = Lin2.unit(ring, rig, V1, V)
    assert unit == Lin2.unit(ring, rig, V1, V)
    #print(unit)
    counit = Lin2.counit(ring, rig, V1, V)
    #print(counit)

    #print(counit * unit)
    assert unit == unit

    mid = unit * counit

    Vi = V.identity()
    Vii = Lin2.identity(ring, rig, Vi)

    lhs = (unit << Vii) * (Vii << counit)
    assert lhs == (unit*Vii) << (Vii*counit) # on the nose
    assert lhs == unit << counit # on the nose
    assert lhs.weak_eq(mid)

    rhs = (Vii << unit) * (counit<<Vii) 
    assert rhs == counit << unit # on the nose

    rhs = (Vii*counit) << (unit*Vii)
    assert rhs.weak_eq(unit << counit)
    assert rhs.weak_eq(mid)


if __name__ == "__main__":
    main()


