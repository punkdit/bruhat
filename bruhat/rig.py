#!/usr/bin/env python3

"""
Kapranov-Voevodsky 2-vector spaces
"""


from bruhat import element
from bruhat.element import GenericElement, Type, Keyed
Gen = GenericElement
from bruhat import elim
from bruhat.chain import Space, Lin

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


def main():

    ring = element.Q
    N = Space(ring, 0, 0, 'N') # null
    K = Space(ring, 1, 0, 'K') # field
    rig = Rig(N, K, Space.__add__, Space.__matmul__)

    # Objects are Space's over a rig
    V0 = Space(rig, 0, name='V0')
    V1 = Space(rig, 1, name='V1')
    #V2 = Space(rig, 2, name='V2')
    V2 = V1 + V1

    i1 = V1.identity()

    #print(i1[0,0].value)
    #A = [[Lin(i1[0,0].value, i1[0,0].value )]]
    A = [[Lin(K, K)]]
    m = Lin2(i1, i1, A)

    copy = Lin(V2, V1, [[rig.one],[rig.one]])
    delete = Lin(V0, V1)
    add = Lin(V1, V2, copy.A.transpose())
    zero = Lin(V1, V0)

    #print(add * copy)
    #print(copy * add)
    assert copy * zero == Lin(V2, V0)
    
    tgt, src = add * zero.direct_sum(i1), i1

    print(tgt)
    print(src)

    m = Lin2(tgt, src)
    print(m)


if __name__ == "__main__":
    main()


