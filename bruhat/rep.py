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
                rhs = h1@h2
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

    __mul__ = tensor

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

    acts = []
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
        acts.append(act)
        assert len(act.components())==1 # transitive
        #print act.tgt.perms
        print("%s subgroup order = %d, number of cosets = %d, conjugates = %d" %(
            act.name, len(H), len(cosets), len(H.conjugates)))

    print(list(names.values()))

    A, B, C, D = acts

    reps = []
    for act in acts:
        rep = Rep.mk_rep(act, tp)
        reps.append(rep)

    kA, kB, kC, kD = reps

    kCC = kC*kC
    #kCC.dump()

    CC = C.pushout(C)
    rep = Rep.mk_rep(CC, tp)
    #rep.dump()
    assert rep==kCC

    U = Space(C.items, ring)
    UU = U*U
    #print(UU)

    cup = U.cup
    ident = U.ident

    for act in CC.components():
        #rep = Rep.mk_rep(act, tp)
        #print(rep.space)
        space = Space(act.items, ring)
        print(space)
        f = space.inject_to(UU)
        #for g in G:
        #    fgf = f@rep[g]@f.transpose()
        ff = f@f.transpose()
        print(ff.hom)
        print(ff)
        r = (cup * ident) @ (ident * ff)
        #rep.dump()

    return

    for i in range(len(acts)):
      for j in range(i, len(acts)):
        A = acts[i]
        B = acts[j]
        C = A.pushout(B)
        #print(C.send_perms)

        for act in C.components():
        #    print(act.send_perms)
            send_perms = act.send_perms
            for g in G:
                h = send_perms[g]
                print(h.perm)



def test():

    ring = element.Q

    n = argv.get("n", 3)
    G = Group.symmetric(n)

    tp = Cat(G, ring)
    burnside(tp)

    return


    rep = Rep.perm_rep(G, cat)
    rep.check()

    r2 = rep * rep
    r2.check()

    r2.dump()


if __name__ == "__main__":

    test()






