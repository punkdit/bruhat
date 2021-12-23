#!/usr/bin/env python3

"""
Kapranov-Voevodsky 2-vector spaces
"""


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
        assert src.shape == tgt.shape
        rows, cols = src.shape
        if A is None:
            A = [[Lin(tgt[i,j].value, src[i,j].value) # note unwrapping here 
                for j in range(cols)] for i in range(rows)]
        A = elim.array(A)
        for j in range(cols):
            for i in range(rows):
                #print(i, j)
                #print( A[i, j].tgt )
                #print(tgt[i, j].value )
                assert A[i, j].tgt == tgt[i, j].value # unwrapped
                assert A[i, j].src == src[i, j].value # unwrapped
        self.A = A

    def __str__(self):
        return str(self.A)



def test_rig(rig):
    one = rig.one
    zero = rig.zero
    assert one + zero == one


def test_structure():
    ring = element.Q
    N = Space(ring, 0, 0, 'N') # null (direct_sum unit)
    K = Space(ring, 1, 0, 'K') # field (tensor unit)
    U = Space(ring, 2, 0, 'U')
    V = Space(ring, 3, 0, 'V')
    W = Space(ring, 5, 0, 'W')

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

    tgt, src = U@V + U@W, U@(V+W)
    f = src.left_distributor()
    assert f.tgt == tgt
    fi = src.left_distributor(inverse=True)
    assert f*fi == tgt.identity()
    assert fi*f == src.identity()


def main():

    ring = element.Q
    N = Space(ring, 0, 0, 'N') # null (direct_sum unit)
    K = Space(ring, 1, 0, 'K') # field (tensor unit)

    # N, K are the zero and one of a rig.
    # We "wrap" these Space objects into a GenericElement.
    rig = Rig(N, K, Space.__add__, Space.__matmul__)

    # Objects are Space's over a rig
    V0 = Space(rig, 0, name='V0')
    V1 = Space(rig, 1, name='V1')
    #V2 = Space(rig, 2, name='V2')
    V2 = V1 + V1

    # 1-morphism's are Lin's over this rig
    i1 = V1.identity()
    copy = Lin(V2, V1, [[rig.one],[rig.one]])
    delete = Lin(V0, V1)
    add = Lin(V1, V2, copy.A.transpose())
    zero = Lin(V1, V0)

    assert copy * zero == Lin(V2, V0)

    # 2-morphism's are Lin2's.
    # These are 2d arrays of Lin's.
    zero_2 = Lin2(i1, i1, [[Lin(K, K)]])

    # The left unitorator
    tgt_1, src_1 = add * zero.direct_sum(i1), i1
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


if __name__ == "__main__":
    test_structure()
    main()


