#!/usr/bin/env python

"""
Groupoid's represented by Function's on hom sets.
"""

from time import time
start_time = time()

from bruhat.action import Group, Perm, mulclose_hom
from bruhat.argv import argv


class Function(object):
    def __init__(self, tgt, src, send={}):
        self.tgt = tgt
        self.src = src
        self.hom = (tgt, src)
        self.send = dict(send)
        self.inv = None

    # __eq__   == is
    # __hash__ == id(.)

    def __mul__(self, other):
        assert self.src is other.tgt # is or == ??
        function = self.send[other]
        return function

    def __setitem__(self, key, value):
        assert self.send.get(key, value) is value
        assert key.tgt is self.src # is or == ??
        assert value.tgt is self.tgt
        assert value.src is key.src
        self.send[key] = value

    def __len__(self):
        return len(self.send)

    def __invert__(self):
        return self.inv



class Groupoid(object):
    def __init__(self, functions, rank=None):
        obs = set()
        homs = {}
        for function in functions:
            if rank is None:
                rank = len(function)
            assert len(function) == rank
            homs.setdefault(function.hom, []).append(function)
            obs.add(function.src)
            obs.add(function.tgt)

        self.functions = list(functions)
        self.rank = rank
        self.gen = None
        self.homs = homs
        self.obs = obs
        self.check()

    def check(self):
        fs = self.functions
        for f in fs:
            src = (~f)*f
            tgt = f*~f
            assert ~(~f) == f
            assert f*src == f
            assert tgt*f == f
            for g in fs:
              if f.tgt != g.src:
                  continue
              assert g*tgt == g
              gf = g*f
              for h in fs:
                if h.src != g.tgt:
                  continue
                hgf = h*gf
                hg = h*g
                assert hgf == hg*f

    def __getitem__(self, idx):
        return self.functions[idx]

    def __len__(self):
        return len(self.functions)

    the_star = "*"
    @classmethod
    def from_group(cls, G):
        the_star = Groupoid.the_star
        lookup = {g : Function(the_star, the_star) for g in G}
        for g in G:
            a = lookup[g]
            a.inv = lookup[~g]
            for h in G:
                b = lookup[h]
                ab = lookup[g*h]
                a[b] = ab
        fs = [lookup[g] for g in G]
        S = Groupoid(fs, len(G))
        S.gen = [lookup[g] for g in G.gen]
        return S

    @classmethod
    def left_cayley(cls, G):
        fs = {}
        for src in G:
          for tgt in G:
            f = Function(tgt, src)
            fs[tgt, src] = f
        for a in G:
         for b in G:
          f = fs[b, a]
          f.inv = fs[a, b]
          for c in G:
            g = fs[c, b]
            g[f] = fs[c, a]
        
        S = Groupoid(list(fs.values()), len(G))
        return S
            
    @classmethod
    def left_action(cls, action):
        G = action.src
        items = action.items
        send_perms = action.send_perms
        T = Groupoid.from_group(G)
        hom = mulclose_hom(G.gen, T.gen)
        send_obs = {item : Groupoid.the_star for item in items}
        send_fns = {}
        homs = {}
        fs = []
        for src in items:
          for g in G:
            tgt = send_perms[g][src]
            function = Function(tgt, src)
            homs[tgt, src, g] = function
            send_fns[function] = hom[g]
            fs.append(function)
        for x in items:
          for g in G:
            y = send_perms[g][x]
            a = homs[y, x, g]
            a.inv = homs[x, y, ~g]
            for h in G:
                z = send_perms[h][y]
                b = homs[z, y, h]
                ba = homs[z, x, h*g]
                b[a] = ba
        S = Groupoid(fs)
        hom = Hom(T, S, send_obs, send_fns)
        return hom


class Hom(object):
    "A homomorphism of Groupoid's"
    def __init__(self, tgt, src, send_obs, send_fns):
        assert isinstance(tgt, Groupoid)
        assert isinstance(src, Groupoid)
        self.tgt = tgt
        self.src = src
        self.send_obs = dict(send_obs)
        self.send_fns = dict(send_fns)

        for (f,g) in send_fns.items():
            assert isinstance(f, Function)
            assert isinstance(g, Function)

        for f in src.functions:
          for g in src.functions:
            if g.tgt != f.src:
                continue
            fg = f*g
            assert send_fns[fg] == send_fns[f]*send_fns[g]

    def __getitem__(self, f):
        return self.send_fns[f]




def test():
    a, b = [0, 1]
    ia = Function(a, a)
    ib = Function(b, b)
    f = Function(b, a)
    g = Function(a, b)
    ia[g] = g
    ia[ia] = ia
    f[ia] = f
    f[g] = ib
    g[ib] = g
    g[f] = ia
    ib[f] = f
    ib[ib] = ib
    ia.inv = ia
    ib.inv = ib
    f.inv = g
    g.inv = f
    S = Groupoid([ia, ib, f, g], 2)
    assert ~f == g
    assert f*g == ib

    G = Group.symmetric([0, 1, 2])
    assert len(G) == 6

    S = Groupoid.from_group(G)

    S = Groupoid.left_cayley(G)

    Hs = list(G.subgroups())
    Hs.sort(key = len)
    for H in Hs:
        print(len(H), end='\n')
        action = G.action_subgroup(H)
        hom = Groupoid.left_action(action)
        T, S = hom.tgt, hom.src
        for g in T.gen:
            for x in S.obs:
                gen = [h for h in S.homs[x,x] if hom[h]==g]
                print(len(gen), end=" ")
            print()
            


if __name__ == "__main__":
    fn = argv.next() or "test"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))



