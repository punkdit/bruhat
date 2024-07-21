#!/usr/bin/env python
"""

Thinking of ordinal versus cardinal numbers as combinatorial species,
are they isomorphic or what?
"""

from bruhat.argv import argv
from bruhat.presh import Set, Function
from bruhat.util import all_perms, all_bijections
from bruhat.action import mulclose

def all_endo(X):
    n = len(X)
    Y = []
    for perm in all_perms(X):
        action = dict((X[i], perm[i]) for i in range(n))
        f = Function(X, X, action)
        Y.append(f)
    return Y


class Species(object):
    def get_ob(self, X):
        pass
    def get_func(self, f):
        pass
    def __call__(self, item):
        if isinstance(item, Set):
            return self.get_ob(item)
        if isinstance(item, Function):
            return self.get_func(item)
        assert 0

    def is_covar(F, X):
        I = Function.ident(X)
        FI = Function.ident(F(X))
        if F(I) != FI:
            return False
        for f in all_endo(X):
            Ff = F(f)
            for g in all_endo(X):
                Fg = F(g)
                Fgf = F(g*f)
                if Fg*Ff != Fgf:
                    return False
        return True

    def is_contravar(F, X):
        I = Function.ident(X)
        FI = Function.ident(F(X))
        if F(I) != FI:
            return False
        for f in all_endo(X):
            Ff = F(f)
            for g in all_endo(X):
                Fg = F(g)
                Fgf = F(g*f)
                if Ff*Fg != Fgf:
                    return False
        return True



class Cardinal(Species):
    def get_ob(self, X):
        Y = all_endo(X) # functions
        Y = Set(Y)
        return Y

    def get_func(self, f):
        X = f.src
        Y = f.tgt
        FX = self(X) # Function's
        FY = self(Y) # Function's
        #print(FY)
        action = {}
        fi = ~f
        for x in FX:
            fx = f*x*fi
            #print(fx)
            assert fx.src == fi.src
            assert fx in FY, "%s not in %s"%(fx, FY)
            action[x] = fx
        return Function(FX, FY, action)


class OpCardinal(Cardinal):
    def get_func(self, f):
        X = f.src
        Y = f.tgt
        FX = self(X) # Function's
        FY = self(Y) # Function's
        #print(FY)
        action = {}
        fi = ~f
        for y in FY:
            fy = fi*y*f
            assert fy in FY
            action[y] = fy
        return Function(FY, FX, action)




class Ordinal(Species):
    def get_ob(self, X):
        Y = list(all_perms(X)) # tuples
        Y = Set(Y)
        return Y

    def get_func(self, f):
        X = f.src
        Y = f.tgt
        GX = self(X) # Set of tuples
        GY = self(Y) # Set of tuples
        action = {}
        for x in GX:
            y = tuple(f(i) for i in x)
            assert y in GY
            action[x] = y
        return Function(GX, GY, action)


class OpOrdinal(Ordinal):
    def get_func(self, f):
        X = f.src
        GX = self(X)
        action = {}
        for x in GX:
            y = tuple(f(i) for i in x)
            assert y in GX
            action[y] = x
        return Function(GX, GX, action)


def find_nats(F, G, X):
    # search for a _natural transform F--->G, at X
    FX = F(X)
    GX = G(X)
    found = []
    for nat in all_bijections(FX, GX):
        nat = Function(FX, GX, nat)
        for f in all_endo(X):
            Ff = F(f)
            Gf = G(f)
            if Gf*nat != nat*Ff:
                break
        else:
            found.append(nat)
    return found



    

def main():

    F = Cardinal()
    G = Ordinal()

    X = Set("abc")
    Y = Set([1,2,3])

    FX = F(X)
    assert len(FX) == 6

    f = Function(X, Y, {"a":1, "b":2, "c":3})
    assert f*~f == Function.ident(Y)
    assert ~f*f == Function.ident(X)

    Ff = F(f)
    assert Ff.src == F(f.src)
    assert Ff.tgt == F(f.tgt)

    assert F.is_covar(X)

    GX = G(X)
    assert len(GX) == 6

    assert G.is_covar(X)

    assert not F.is_contravar(X)
    assert not G.is_contravar(X)

    assert not OpCardinal().is_covar(X)
    assert OpCardinal().is_contravar(X)

    assert not OpOrdinal().is_covar(X)
    assert OpOrdinal().is_contravar(X)

    natXs = find_nats(F, G, X)
    assert not natXs

    nats = find_nats(F, F, X)
    assert len(nats) == 2, len(nats)

    nats = find_nats(G, G, X)
    assert len(nats) == 6

    #X = Set("ab")
    #nats = find_nats(F, G, X)
    #assert len(nats) == 2


if __name__ == "__main__":
    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()



