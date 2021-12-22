#!/usr/bin/env python3

"""

"""


from bruhat import element
from bruhat.element import GenericElement, Type, Keyed
Gen = GenericElement
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


def test_rig(rig):
    one = rig.one
    zero = rig.zero
    assert one + zero == one

def main():

    ring = element.Q
    rig = Rig(
        Space(ring, 0, 0, '0'),
        Space(ring, 1, 0, 'I'),
        Space.__add__,
        Space.__matmul__)

    V = Space(rig, 2)
    lin = Lin(V, V)
    I = V.identity()
    print(lin)
    print(I)
    print(I+I)
    print(I*I)
    print(I@I)
    


if __name__ == "__main__":
    main()


