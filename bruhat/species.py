#!/usr/bin/env python3

"""
Combinatorial species.

See chapter 4 here:
    COMBINATORIAL SPECIES AND LABELLED STRUCTURES 
    Brent Abraham Yorgey 
    https://www.cis.upenn.edu/~sweirich/papers/yorgey-thesis.pdf
"""

import string
letters = list(string.ascii_lowercase)

from bruhat.util import factorial, allperms, cross
from bruhat import series
from bruhat.series import Series, ring



class Set(object): # copied from bruhat.rel
    def __init__(self, items=[]):
        if type(items) is int:
            assert items <= len(letters)
            items = letters[:items]
        #items = [(x if type(x) is tuple else tuple(x)) for x in items]
        items = list(items)
        items.sort() # eeeeck: watch this !
        self.items = tuple(items)
        self.set_items = set(items)
        assert len(self.set_items) == len(self.items)
        #pairs = [(x, x) for x in items]
        #self.ident = Rel(pairs, self, self)

    def __str__(self):
        return "Set(%s)"%(str(list(self.items)))
    __repr__ = __str__

    @classmethod
    def promote(cls, items):
        if isinstance(items, Set):
            return items
        return Set(items)

    def __iter__(self):
        return iter(self.items)

    def __len__(self):
        return len(self.items)

    def __contains__(self, item):
        return item in self.set_items

    def __eq__(a, b):
        return a.items == b.items

    def __ne__(a, b):
        return a.items != b.items

    def __lt__(a, b):
        return a.items < b.items

    def __le__(a, b):
        return a.items < b.items

    def __hash__(self):
        return hash(self.items)

    def __add__(a, b):
        # categorical coproduct
        items = [(0, ai) for ai in a] + [(1, bi) for bi in b]
        return Set(items)

    def union(a, b):
        items = a.set_items.union(b.set_items)
        return Set(items)

    def __mul__(a, b):
        # categorical product
        items = [(ai, bi) for ai in a for bi in b]
        return Set(items)
        
    def all_partitions(self):
        # A partition is a Set of subsets of self
        # such that: disjoint, non-empty (?), covering.
        items = self.items
        if not items:
            return Set(empty) # ?
        #assert items
        if len(items) == 1:
            p = Set([self])
            return Set([p]) # only one partition here
        if len(items) == 2:
            p = Set([self])
            a, b = self
            q = Set([
        item = Set([items[0]])
        rest = Set(items[1:])
        for partition in rest.all_partitions():
            

    def all_parts2(self):
        "all the ways to break self into two subsets"
        items = self.items
        n = len(items)
    #    if n==0:
    #        yield items, items
    #        return
    #
    #    if n==1:
    #        yield [], items
    #        yield items, []
    #        return
    #
    #    if n==2:
    #        yield [], items
    #        yield [items[0]], [items[1]]
    #        yield [items[1]], [items[0]]
    #        yield items, []
    #        return
    
        bits = [(0, 1)]*n
        for idxs in cross(bits):
            left = Set([items[i] for i in range(n) if idxs[i]==0])
            right = Set([items[i] for i in range(n) if idxs[i]==1])
            yield left, right

assert len(list(Set(3).all_parts2())) == 8



empty = Set()

class Species(object):
    def __init__(self, f):
        self.apply = f

    def __getitem__(self, U):
        return self.apply(U)

    def __add__(F, G):
        return AddSpecies(F, G)

    def __mul__(F, G):
        "Cartesian/Hadamard product"
        return MulSpecies(F, G)

    def dot(F, G):
        "Partitional/Cauchy product"
        return DotSpecies(F, G)

    def __call__(F, G):
        return ComposeSpecies(F, G)

    def dirichlet(F, G):
        return DirichletSpecies(F, G)



class BinaryOpSpecies(Species):
    def __init__(self, F, G):
        self.F = F
        self.G = F


class AddSpecies(BinaryOpSpecies):
    def __getitem__(self, U):
        return self.F[U] + self.G[U]
    

class MulSpecies(BinaryOpSpecies):
    def __getitem__(self, U):
        return self.F[U] * self.G[U]


class DotSpecies(BinaryOpSpecies):
    def __getitem__(self, U):
        F = self.F
        G = self.G
        items = []
        for idx, (left, right) in enumerate(U.all_parts2()):
            for item in F[left] * G[right]:
                items.append((idx, item))
        return Set(items)


class ComposeSpecies(BinaryOpSpecies):
    def __getitem__(self, U):
        F = self.F
        G = self.G
        items = []
        # ...
        return Set(items)

    
class DirichletSpecies(BinaryOpSpecies):
    def __getitem__(self, U):
        F = self.F
        G = self.G
        items = []
        # ...
        return Set(items)


# ----------------------------------------------------

Zero = Species(lambda items : Set([]))
One = Species(lambda items : (Set([items]) if len(items)==0 else empty))
X = Species(lambda items : (Set([items]) if len(items)==1 else empty))
E = Species(lambda items : Set([items]))
Point = Species(lambda items : Set.promote(items))
List = Species(lambda items : Set(list(allperms(items))))


# ----------------------------------------------------
    
class GeneratingFunction(Series):
    def __init__(self, species):
        Series.__init__(self)
        self.species = species

class EGF(GeneratingFunction):
    def __getitem__(self, idx):
        U = Set(idx)
        FU = self.species[U]
        return ring.one*len(FU) / factorial(idx)

class OGF(GeneratingFunction):
    def __getitem__(self, idx):
        U = Set(idx)
        FU = self.species[U]
        return ring.one*len(FU)


# ----------------------------------------------------


def test():

    for n in range(0):
        U = Set(n)
        print("Zero:", Zero[U])
        print("X:", X[U])
        print("E:", E[U])
        print("Point:", Point[U])

    assert EGF(Zero).eq(series.zero)
    assert EGF(One).eq(series.one)
    assert EGF(X).eq(series.x)
    assert EGF(E).eq(series.exp)
    #print(EGF(E))

    f = EGF(List) # = 1/(1-x)
    for i in range(5):
        assert f[i] == 1

    assert EGF(X+X).eq(series.x + series.x)

    XX = X.dot(X)
    assert EGF(XX).eq(series.x * series.x)


if __name__ == "__main__":

    test()










