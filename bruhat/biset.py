#!/usr/bin/env python

"""
Given groups G, H,
A (G,H)-biset X has a left G action and a right H action
that commute w each other.

see also: adjoint.py, combinatorial.py for previous attempts
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
        self.l_send = dict(l_send) # send lg in lG --> Perm of items
        self.r_send = dict(r_send) # send rg in rG --> Perm of items
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

    def __lshift__(X, Y):
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

#    def __mul__(self, other):
#        assert self.lG == other.lG
#        assert self.rG == other.rG
#        items = []
#        for a1 in self.items:
#          for a2 in other.items:
#            items.append((a1, a2))
#        l_send = {}
#        r_send = {}
#        for (send, G) in [(l_send, self.lG), (r_send, self.rG)]:
#            for g in G:
#                perm = {}
#                h1 = self.l_send[g]
#                h2 = other.l_send[g]
#                for a1, a2 in items:
#                    perm[(a1, a2)] = h1(a1), h2(a2)
#                perm = Perm(perm, items)
#                send[g] = perm
#        return Action(self.G, send_perms, items)

    def to_action(self, check=True):
        lG, rG, l_send, r_send = self.data
        items = self.items
        G, l_proj, r_proj = lG.universal_product(rG)
        send_perms = {}
        for g in G:
            lg = l_proj[g]
            rg = r_proj[g]
            perm = l_send[lg] * ~r_send[rg]
            assert perm == ~r_send[rg] * l_send[lg]
            send_perms[g] = perm
        return Action(G, send_perms, self.items, check=check)

    def get_homs(self, other):
        assert self.lG is other.lG # strict !
        assert self.rG is other.rG # strict !
        A = self.to_action()
        B = other.to_action()
        B.G = A.G # HACK THIS
        for hom in A.get_homs(B):
            hom = Hom(self, other, hom.send_items)
            yield hom


class Hom(object):
    "Hom of Biset's"
    def __init__(self, src, tgt, send_items, check=True):
        assert isinstance(src, Biset)
        assert isinstance(tgt, Biset)
        assert src.lG == tgt.lG
        assert src.rG == tgt.rG
        #lG, rG = src.lG, src.rG
        self.src = src
        self.tgt = tgt
        self.lG = src.lG
        self.rG = src.rG
        self.send_items = dict(send_items)
        if check:
            self.check()

    def check(self):
        src = self.src
        tgt = self.tgt
        send_items = self.send_items
        for x,y in send_items.items():
            assert x in src
            assert y in tgt
        for x in src:
            assert x in send_items
        for g in self.lG:
            for item in src.items:
                left = send_items[src.l_send[g][item]]
                right = tgt.l_send[g][send_items[item]]
                if left != right:
                    print("item =", item)
                    print("g =", g, "src(g) =", src.l_send[g], "tgt(g) =", tgt.l_send[g])
                    print("send_items =", send_items)
                    print("%s != %s"%(left, right))
                    assert 0, "not a Hom of Action's"
        for g in self.rG:
            for item in src.items:
                left = send_items[src.r_send[g][item]]
                right = tgt.r_send[g][send_items[item]]
                assert left == right

    @classmethod
    def identity(cls, X):
        assert isinstance(X, Biset)
        send_items = {x:x for x in X}
        return Hom(X, X, send_items)

    def __str__(self):
        return "Hom(%s, %s, %s)"%(self.src, self.tgt, self.send_items)
    __repr__ = __str__

    def __eq__(self, other):
        assert self.lG is other.lG # strict 
        assert self.rG is other.rG # strict 
        assert self.src == other.src
        assert self.tgt == other.tgt
        return self.send_items == other.send_items

    def __ne__(self, other):
        assert self.lG is other.lG # strict 
        assert self.rG is other.rG # strict
        assert self.src == other.src
        assert self.tgt == other.tgt
        return self.send_items != other.send_items

    _hash = None
    def __hash__(self):
        if self._hash is not None:
            return self._hash
        pairs = list(self.send_items.items())
        pairs.sort(key = str) # canonical form
        pairs = tuple(pairs)
        self._hash = hash(pairs)
        return self._hash

    def compose(self, other):
        # other o self
        assert isinstance(other, Hom)
        assert self.tgt == other.src
        a = self.send_items
        b = other.send_items
        #send_items = [b[i] for i in a] # wut
        send_items = {i:b[a[i]] for i in a}
        return Hom(self.src, other.tgt, send_items)

    def __mul__(self, other):
        assert isinstance(other, Hom)
        return other.compose(self)





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
    XY = X<<Y
    assert len(XY) == 6

    f = Hom.identity(X)
    assert f*f == f

    A = XY.to_action()

    I = Hom.identity(XY)
    gen = []
    for g in XY.get_homs(XY):
        gen.append(g)
    assert len(gen) == 6
    assert gen.count(I) == 1
    for g in gen:
     for h in gen:
        assert g*h in gen


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








