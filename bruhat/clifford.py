#!/usr/bin/env python


"""
Implement Bruhat decomposition of the Clifford group
as described here:

Hadamard-free circuits expose the structure of the Clifford group
Sergey Bravyi, Dmitri Maslov
https://arxiv.org/abs/2003.09412

maybe ?

"""


from operator import mul, matmul
from functools import reduce, cache


from bruhat.action import mulclose
from bruhat.clifford_sage import Clifford
from bruhat.util import cross
from bruhat.argv import argv



class FreeGroup(object):
    def mul(self, g, h):
        word = Word(g.gens + h.gens, self)
        return word

free_group = FreeGroup()


class AbsGen(object):
    def __init__(self, name, parent=free_group):
        self.name = name
        self.parent = parent

    def __eq__(self, other):
        return self.name == other.name

    def __str__(self):
        return self.name
    __repr__ = __str__

    def __hash__(self):
        return hash(self.name)

    def __mul__(self, other):
        #return Word([self, other])
        #return Word(self.gens + other.gens)
        return self.parent.mul(self, other)

    def __pow__(self, n):
        assert n>=0
        if n==0:
            return Word([], self.parent)
        return reduce(self.parent.mul, [self]*n)


class Word(AbsGen):
    def __init__(self, gens=[], parent=free_group):
        self.gens = list(gens)
        name = "*".join(str(gen) for gen in gens) or "I"
        AbsGen.__init__(self, name, parent)

    def __len__(self):
        return len(self.gens)

    def promote(self, parent):
        return Word(self.gens, parent)

    def subs(self, src, tgt):
        n0, n1 = len(self), len(src)
        gens = self.gens
        src = src.gens
        for i in range(n0-n1+1):
            if gens[i:i+n1] == src:
                break
        else:
            return None
        gens = gens[:i] + tgt.gens + gens[i+n1:]
        return Word(gens, self.parent)


class FPGroup(FreeGroup):
    def __init__(self, gens, rels):
        self.rels = []
        I = Word([], self)
        for rel in rels:
            if isinstance(rel, Word):
                rel = rel.promote(self)
                self.rels.append((rel, I))
            else:
                a, b = rel
                a = a.promote(self)
                b = b.promote(self)
                self.rels.append((a, b))
        self.gens = [self.reduce(word.promote(self)) for word in gens]

    def reduce(self, g):
        rels = self.rels
        done = False
        while not done:
            done = True
            for rel in rels:
                k = g.subs(*rel)
                if k is not None:
                    g = k
                    done = False
        return g

    def mul(self, g, h):
        gh = FreeGroup.mul(self, g, h)
        gh = self.reduce(gh)
        return gh

    @property
    @cache
    def elements(self):
        return list(mulclose(self.gens))

    def __len__(self):
        return len(self.elements)

    def __getitem__(self, i):
        return self.elements[i]


def test_word():

    I = Word()
    U = Word("U")
    V = Word("V")

    assert U==U
    assert U*V != V*U
    assert I*U == U
    assert str(U*V) == "U*V"

    UU = U*U
    UUU = UU*U
    assert UUU == U*UU

    assert UUU.subs(UU, V) == V*U
    assert UU.subs(UU, I) == I

    G = FPGroup([U], [(U**7, I)])
    U = G.gens[0]

    assert len(mulclose([U])) == 7

    a = Word("a")
    b = Word("b")
    #G = FPGroup([a, b], [a**2, b**2, (a*b)**3]) # fail..
    G = FPGroup([a, b], [a**2, b**2, (a*b*a, b*a*b)])
    #for g in G:
    #    print(g)
    assert len(G) == 6




def test():

    test_word()
    test_normal_form()


def build_normal_form(lookup, factors):

    opss = [ mulclose([lookup[g] for g in factor]) for factor in factors ]
    print([len(f) for f in opss])

    found = set()
    for ops in cross(opss):
        g = reduce(mul, ops)
        found.add(g)
    print(len(found))




def test_normal_form():

    w = Word("w")
    X = Word("X")
    Z = Word("Z")

    cliff = Clifford(1)

    lookup = {
        w : cliff.w()**2,
        X : cliff.X(),
        Z : cliff.Z()
    }
    factors = [[w], [X], [Z]]
    build_normal_form(lookup, factors)

    cliff = Clifford(3)
    A = cliff.CX(0, 1)
    B = cliff.CX(1, 2)
    assert A*B != B*A

    A = cliff.CZ(0, 1)
    B = cliff.CZ(1, 2)
    assert A*B == B*A



    
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

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()





