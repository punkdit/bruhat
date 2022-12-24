#!/usr/bin/env python

"""
Groupoid's represented by Function's on hom sets.
"""

from time import time
start_time = time()

from bruhat.action import Group, Perm
from bruhat.argv import argv


class Function(object):
    def __init__(self, tgt, src, send={}):
        self.tgt = tgt
        self.src = src
        self.hom = (tgt, src)
        self.send = dict(send)
        self.inv = None

    #def __eq__(self, other):
    #    return self.hom==other.hom and self.send==other.send

    def __mul__(self, other):
        assert self.src is other.tgt # is or == ??
        #function = {}
        #function = Function(self.tgt, other.src, function)
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
        #src, tgt = self.src, self.tgt
        #send = dict((v, k) for (k,v) in self.send.items())
        #return Function(tgt, src, send)



class Groupoid(object):
    def __init__(self, functions, rank=None):
        for function in functions:
            if rank is None:
                rank = len(function)
            assert len(function) == rank
        self.functions = list(functions)
        self.rank = rank
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
        homs = {}
        fs = []
        for src in items:
          for g in G:
            tgt = send_perms[g][src]
            function = Function(tgt, src)
            homs[tgt, src, g] = function
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
        return Groupoid(fs)


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
    S = Groupoid.left_cayley(G)

    for H in G.subgroups():
        print(len(H), end=' ')
        action = G.action_subgroup(H)
        S = Groupoid.left_action(action)
        print(len(S.functions))


if __name__ == "__main__":
    fn = argv.next() or "test"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))



