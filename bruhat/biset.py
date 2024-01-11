#!/usr/bin/env python

"""
Given groups G, H,
A (G,H)-biset X has a left G action and a right H action
that commute w each other.

see also: adjoint.py for previous attempts
"""

from functools import lru_cache
cache = lru_cache(maxsize=None)

from bruhat.equ import Equ

from bruhat.action import Perm, Group, Action
from bruhat.argv import argv


def _trivial(items):
    I = Perm.identity(items)
    G = Group([I], items)
    send = {g:I for g in G}
    return G, send


class Biset(object):
    def __init__(self, lG, rG, l_send, r_send, items, check=True):
        assert isinstance(lG, Group)
        assert isinstance(rG, Group)
        assert isinstance(l_send, dict)
        assert isinstance(r_send, dict)
        self.lG = lG
        self.rG = rG
        self.l_send = dict(l_send)
        self.r_send = dict(r_send)
        self.items = list(items)
        if check:
            self.check()

    @property
    def data(self):
        return (self.lG, self.rG, self.l_send, self.r_send)

    def __str__(self):
        return "Biset(%s, %s, %s)"%(self.lG, self.rG, len(self.items))

    # Equality on-the-nose:
    def __eq__(self, other):
        assert isinstance(other, Biset)
        return self.data == other.data

    def __ne__(self, other):
        assert isinstance(other, Biset)
        return self.data != other.data

    # i am a set (list) of items
    def __getitem__(self, idx):
        return self.items[idx]

    # with a len'gth
    def __len__(self):
        return len(self.items)

    def __contains__(self, x):
        return x in self.items # ouch it's a list

# TODO
#    def __hash__(self):
#        send_perms = self.send_perms
#        send_perms = tuple((perm, send_perms[perm]) for perm in self.G)
#        return hash((self.G, send_perms))

    def check(self):
        lG, rG, l_send, r_send = self.data
        items = self.items

        for (G, send_perms, op) in [(lG, l_send, False), (rG, r_send, True)]:
            assert len(send_perms)==len(G.perms)
            for perm in G.perms:
                assert perm in send_perms
                perm = send_perms[perm]
                assert perm.items == items
    
            # Here we check that we have a homomorphism of groups.
            for g1 in G.perms:
              h1 = send_perms[g1]
              for g2 in G.perms:
                h2 = send_perms[g2]
                rhs = h2*h1 if op else h1*h2
                assert send_perms[g1*g2] == rhs

        for l in lG:
          for r in rG:
            lr = l_send[l] * r_send[r]
            rl = r_send[r] * l_send[l]
            assert lr == rl

    @classmethod
    def from_action(cls, action):
        "build Biset with a left action, and trivial right action"
        assert isinstance(action, Action)
        lG = action.G
        l_send = action.send_perms
        items = action.items
        #rG = Group.trivial()
        #I = Perm.identity(items)
        #r_send = {g:I for g in rG}
        rG, r_send = _trivial(items)
        return Biset(lG, rG, l_send, r_send, items)

    @classmethod
    def cayley(cls, G, lG, rG):
        "the left,right Cayley action of subgroups lG,rG on G"
        assert isinstance(G, Group)
        items = G.perms
        if lG is None:
            lG = G.trivial_subgroup()
        if rG is None:
            rG = G.trivial_subgroup()
        assert isinstance(lG, Group), type(lG)
        assert isinstance(rG, Group), type(rG)
        assert lG.items == rG.items == G.items
        l_send = {l:Perm({g : l*g for g in items}, items) for l in lG}
        r_send = {r:Perm({g : g*(~r) for g in items}, items) for r in rG}
        biset = Biset(lG, rG, l_send, r_send, items)
        return biset

    @classmethod
    def left_cayley(cls, G, H=None):
        "the left Cayley action of a subgroup on self"
        return cls.cayley(G, H, None)

    @classmethod
    def right_cayley(cls, G, H=None):
        "the right Cayley action of a subgroup on self"
        return cls.cayley(G, None, H)

    @classmethod
    def left_tautological(cls, G):
        assert isinstance(G, Group)
        items = G.items
        rG, r_send = _trivial(items)
        l_send = {g:g for g in G}
        return Biset(G, rG, l_send, r_send, items)

    #@cache
    def op(self):
        lG, rG, l_send, r_send = self.data

        lop_send = {~g:l_send[g] for g in lG}
        rop_send = {~g:r_send[g] for g in rG}
        return Biset(rG, lG, rop_send, lop_send, self.items)

    def __matmul__(X, Y):
        "horizontal composition"
        assert isinstance(Y, Biset)
        assert X.rG is Y.lG
        #print("__matmul__")
        G = X.rG
        XY = [(x, y) for x in X.items for y in Y.items]
        #print("X*Y =")
        #for xy in XY:
        #    print("\t", xy)
        lookup = {xy:Equ(xy) for xy in XY}
        for (x,y) in XY:
          for g in G:
            lhs = lookup[X.r_send[g][x], y]
            rhs = lookup[x, Y.l_send[g][y]]
            lhs.merge(rhs)
    
        equs = set(equ.top for equ in lookup.values())
        equs = list(equs)
        items = [tuple(equ.items) for equ in equs]
        #print("items =")
        #for item in items:
        #    print('\t', item)

        l_send = {}
        for l in X.lG:
            perm = {}
            for src in items:
                x, y = src[0]
                lx = l*x
                for tgt in items:
                    if (lx, y) in tgt:
                        break
                else:
                    assert 0
                perm[src] = tgt
            perm = Perm(perm, items)
            l_send[l] = perm

        #print("Y.items", Y.items)
        #print("Y.rG", Y.rG)

        r_send = {}
        for r in Y.rG:
            #print("\tr =", r)
            perm = {}
            for src in items:
                x, y = src[0]
                ry = r*y
                for tgt in items:
                    if (x, ry) in tgt:
                        break
                else:
                    assert 0
                perm[src] = tgt
            perm = Perm(perm, items)
            r_send[r] = perm

        biset = Biset(X.lG, Y.rG, l_send, r_send, items)
        return biset
        


def test():
    G = Group.symmetric(3)
    Hs = G.subgroups()

    for H in Hs:
        action = G.action_subgroup(H)
        x = Biset.from_action(action)
        xop = x.op()

    for H in Hs:
        if len(H) == 3:
            break

    print("H =", H)
    for g in H:
        print('\t', g)

    X = Biset.cayley(G, G, H)
    Y = Biset.left_tautological(H)
    XY = X@Y
    assert len(XY) == 6


if __name__ == "__main__":
    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "test"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds"%(time() - start_time))
    print()








