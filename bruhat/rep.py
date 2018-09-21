#!/usr/bin/env python3

"""
category of reps of a finite group...

"""


from bruhat.action import Group, Perm
from bruhat.element import Type, Keyed, Element
from bruhat import element
from bruhat.vec import Space, Hom, Map

from argv import argv


class Cat(Keyed, Type):
    "tensor category of G-reps over a ring"

    def __init__(self, G, ring):
        Type.__init__(self) 
        Keyed.__init__(self, (G, ring))
        self.G = G
        self.ring = ring


class Rep(Element): # Object of the category
    """
        A G-rep over a ring.
        For each perm in the source Group G we map to a Hom(A, A)
    """

    def __init__(self, send_perms, space, cat):
        Element.__init__(self, cat)
        assert isinstance(space, Space)
        assert isinstance(cat, Cat)
        self.send_perms = dict(send_perms)
        self.space = space
        self.G = cat.G
        self.ring = cat.ring

    def check(self):
        G, space, send_perms = self.G, self.space, self.send_perms

        assert len(send_perms)==len(G.perms)
        for perm in G.perms:
            assert perm in send_perms
            f = send_perms[perm]
            hom = f.hom
            assert hom.src == space
            assert hom.tgt == space

        # Here we check that we have a homomorphism of groups.
        for g1 in G.perms:
            h1 = send_perms[g1]
            for g2 in G.perms:
                h2 = send_perms[g2]
                assert send_perms[g1*g2] == h1@h2

    #def __add__(self, other):

    def tensor(self, other):
        "tensor product"
        assert self.tp == other.tp
        space = self.space * other.space
        G = self.G
        send_perms = {}
        for g in G:
            a = self.send_perms[g]
            b = other.send_perms[g]
            rg = a*b
            send_perms[g] = rg
        return Rep(send_perms, space, self.tp)

    __mul__ = tensor

    def dump(self):
        for g in self.G:
            rg = self.send_perms[g]
            print(g)
            print(rg)
            print()


    @classmethod
    def perm_rep(cls, G, cat):
        
        ring = cat.ring
        one = ring.one
        basis = G.items
        V = Space(basis, ring)
    
        send_perms = {}
        hom = Hom(V, V)
        for g in G:
            rg = Map([((i, g[i]), one) for i in basis], hom)
            send_perms[g] = rg
        rep = cls(send_perms, V, cat)
        return rep


def test():

    ring = element.Z

    n = argv.get("n", 3)
    G = Group.symmetric(n)

    cat = Cat(G, ring)

    rep = Rep.perm_rep(G, cat)
    rep.check()

    r2 = rep * rep
    r2.check()

    r2.dump()


if __name__ == "__main__":

    test()






