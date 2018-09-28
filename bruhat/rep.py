#!/usr/bin/env python3

"""
category of reps of a finite group...

"""

import string

from bruhat.action import Group, Perm, conjugacy_subgroups
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

    def is_hom(self, rep1, rep2, f):
        G = self.G
        for g in G:
            g1 = rep1[g]
            g2 = rep2[g]
            if f*g1 != g2*f:
                return False
        return True


class Rep(Element): # Object of the category
    """
        A G-rep over a ring.
        For each perm in the source Group G we send_perms to a Map.
    """

    def __init__(self, send_perms, space, cat):
        Element.__init__(self, cat)
        assert isinstance(space, Space)
        assert isinstance(cat, Cat)
        self.send_perms = dict(send_perms)
        self.space = space
        self.G = cat.G
        self.ring = cat.ring

    def __eq__(self, other):
        assert isinstance(other, Rep)
        if self.tp != other.tp:
            return False
        if self.G != other.G:
            return False
        if self.space != other.space:
            return False
        return self.send_perms == other.send_perms

    def __ne__(self, other):
        return not self==other

    def __hash__(self, other):
        assert 0, "TODO"

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
                lhs = send_perms[g1*g2]
                rhs = h1*h2
                if lhs != rhs:
                    print("lhs =")
                    print(lhs.items)
                    print("rhs =")
                    print(rhs.items)
                    assert 0

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

    __matmul__ = tensor

    def dump(self):
        print("========= Rep ===========")
        print(self.space)
        for g in self.G:
            rg = self.send_perms[g]
            print(g)
            print(rg)
            print()

    def __getitem__(self, g):
        return self.send_perms[g]

    @classmethod
    def perm_rep(cls, tp):
        G = tp.G
        ring = tp.ring
        one = ring.one
        basis = G.items
        V = Space(basis, ring)
    
        send_perms = {}
        hom = Hom(V, V)
        for g in G:
            rg = Map([((g[i], i), one) for i in basis], hom)
            send_perms[g] = rg
        rep = cls(send_perms, V, tp)
        return rep

    @classmethod
    def mk_rep(cls, act, tp):
        G = tp.G
        ring = tp.ring
        one = ring.one
        gen = act.items
        V = Space(gen, ring)
    
        send_perms = {}
        hom = Hom(V, V)
        for g in G:
            ag = act[g]
            rg = Map([((ag[i], i), one) for i in gen], hom)
            #print(g, "--->")
            #print(rg)
            send_perms[g] = rg
        rep = cls(send_perms, V, tp)
        rep.check()
        return rep
    
    

def burnside(tp): # make it a method

    ring = tp.ring
    G = tp.G
    Hs = conjugacy_subgroups(G)
    names = {}

    letters = list(string.ascii_uppercase + string.ascii_lowercase)
    letters.remove("O")
    letters.remove("o")
    letters = letters + [l+"'" for l in letters] + [l+"''" for l in letters]
    assert len(letters) >= len(Hs)
    letters = letters[:len(Hs)]

    for i, H in enumerate(Hs):
        names[H] = "%s_0"%letters[i] # identity coset

    acts = {}
    for i, H in enumerate(Hs):
        cosets = G.left_cosets(H)
        assert len(G) == len(cosets)*len(H)

        items = [names[H]]
        letter = letters[i]
        assert H in cosets
        idx = 1
        for gH in cosets:
            if gH == H:
                continue
            name = "%s_%d"%(letter, idx)
            names[gH] = name
            items.append(name)
            idx += 1

        act = G.left_action(cosets)
        assert act.src is G
        act = act.rename(names, items)
        act.name = letter
        H.name = act.name
        acts[letter] = act
        assert len(act.components())==1 # transitive
        #print act.tgt.perms
        print("%s subgroup order = %d, number of cosets = %d, conjugates = %d" %(
            act.name, len(H), len(cosets), len(H.conjugates)))

    #print(list(names.values()))

    reps = {}
    for act in acts.values():
        rep = Rep.mk_rep(act, tp)
        reps[act.name] = rep

    arg = argv.next()
    assert "*" in arg
    left, right = arg.split("*")

    act0 = acts[left]
    act1 = acts[right]

    rep0 = reps[left]
    rep1 = reps[right]

    act2 = act0.pushout(act1)
    rep = Rep.mk_rep(act2, tp)

    U0 = Space(act0.items, ring)
    U1 = Space(act1.items, ring)
    UU = U0@U1
    print(UU)

    one = ring.one
    hom = Hom(U0, U1)
    for act in act2.components():

        items = [((y, x), one) for (x, y) in act.items]
        f = Map(items, hom)
        print(f)

        assert tp.is_hom(rep0, rep1, f)

        print("kernel:")
        g = f.kernel()
        print(g)
        print()

        print("cokernel:")
        g = f.cokernel()
        print(g)
        print()



def test():

    ring = element.Q

    n = argv.get("n", 3)
    G = Group.symmetric(n)

    tp = Cat(G, ring)
    burnside(tp)

    return


    rep = Rep.perm_rep(G, cat)
    rep.check()

    r2 = rep @ rep
    r2.check()

    r2.dump()


if __name__ == "__main__":

    test()






