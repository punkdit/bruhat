#!/usr/bin/env python3

"""
Another attempt at action.py
Do everything with indexes.. might be more efficient than action.py
"""

import numpy

from bruhat.util import factorial
from bruhat.action import mulclose
from bruhat import algebraic
from bruhat.argv import argv



DEBUG = argv.get("debug", False)

class Perm(object):

    def __init__(self, perm):
        assert isinstance(perm, list) or isinstance(perm, numpy.ndarray)
        perm = numpy.array(perm)
        self.perm = perm.copy()
        self.rank = len(perm)
        if DEBUG:
            items = list(perm)
            items.sort()
            assert items == list(range(self.rank))

    @property
    def identity(self):
        return Perm(list(range(self.rank)))

    def __str__(self):
        return "Perm(%s)"%(self.perm,)
    __repr__ = __str__

    def __hash__(self):
        return hash(self.perm.tostring())

    def __eq__(self, other):
        return numpy.alltrue(self.perm == other.perm)

    def __ne__(self, other):
        return not numpy.alltrue(self.perm == other.perm)

    def __lt__(self, other):
        return self.perm.tostring() < other.perm.tostring()

    def __le__(self, other):
        return self.perm.tostring() <= other.perm.tostring()

    def __mul__(self, other):
        perm = self.perm[other.perm] # check this is the right order!
        return Perm(perm)

    def inverse(self):
        perm = self.perm
        q = [None]*self.rank
        for (src, tgt) in enumerate(perm):
            q[tgt] = src
        return Perm(q)

    def __getitem__(self, i):
        return self.perm[i]

    #def orbits(self):


class Group(object):
    def __init__(self, perms=None, gen=None):
        if perms is None:
            assert gen is not None
            perms = list(mulclose(gen))
        else:
            perms = list(perms)
        perms.sort()
        assert perms
        #gen = list(gen)
        #self.gen = gen
        self.perms = perms
        self.n = len(self.perms)
        self.rank = perms[0].rank
        self.lookup = dict((perm, idx) for (idx, perm) in enumerate(self.perms))
        if DEBUG:
            self.do_check()

    def do_check(self):
        lookup = self.lookup
        assert self.identity in lookup
        for g in self:
            for h in self:
                assert g*h in lookup
            assert g.inverse() in lookup
            assert g.inverse() * g == self.identity

    @property
    def identity(self):
        return Perm(list(range(self.rank)))

    def __str__(self):
        return "Group(order=%s)"%(self.n,)
    __repr__ = __str__

    def __len__(self):
        return self.n

    def __getitem__(self, idx):
        return self.perms[idx]

    def __eq__(self, other):
        return self.perms == other.perms

    @classmethod
    def trivial(cls, n):
        assert n>0
        perm = Perm(list(range(n)))
        return Group([perm])

    @classmethod
    def symmetric(cls, n):
        assert n>0
        gen = []
        items = list(range(n))
        for i in range(n-1):
            perm = list(items)
            perm[i] = i+1
            perm[i+1] = i
            gen.append(Perm(perm))
        G = Group(gen=gen)
        assert len(G) == factorial(n)
        return G

    @classmethod
    def alternating(cls, n):
        assert n>0
        gen = []
        items = list(range(n))
        for i in range(n-2):
            perm = list(items)
            perm[i] = i+1
            perm[i+1] = i+2
            perm[i+2] = i
            gen.append(Perm(perm))
        G = Group(gen=gen)
        assert len(G) == factorial(n)//2
        return G

    @classmethod
    def cyclic(cls, n):
        assert n>0
        perms = []
        for k in range(n):
            perm = [(i+k)%n for i in range(n)]
            perms.append(Perm(perm))
        G = Group(perms)
        assert len(G) == n
        return G

    @classmethod
    def from_action(cls, G, X):
        lookup = dict((v, idx) for (idx, v) in enumerate(X))
        perms = []
        for g in G:
            perm = [lookup[g*v] for v in X]
            perms.append(Perm(perm))
        G = Group(perms)
        return G

    def regular_action(self):
        "the left _regular gset: the cayley action"
        lookup = self.lookup
        perms = []
        for idx, p in enumerate(self):
            perm = []
            for jdx, q in enumerate(self):
                r = p*q # left action on self
                kdx = lookup[r]
                perm.append(kdx)
            perm = Perm(perm)
            perms.append(perm)
        tgt = Group(perms) # ARGH, shuffles the order of perms
        send_perms = [tgt.lookup[perm] for perm in perms]
        hom = Hom(self, tgt, send_perms)
        return hom

    @property
    def i(self):
        send_perms = list(range(self.n)) # ARGH should we use dict's for these?
        hom = Hom(self, self, send_perms)
        return hom
        
    def orbits(self):
        remain = set(range(self.rank))
        orbits = []
        while remain:
            i = iter(remain).__next__()
            remain.remove(i)
            orbit = set([i])
            for g in self:
                j = g[i]
                while j != i:
                    orbit.add(j)
                    j = g[j]
            remain.difference_update(orbit)
            orbit = list(orbit)
            orbit.sort()
            orbits.append(orbit)
        return orbits

    def get_sequence(G, n=5):
        yield len(G.orbits())
        i = j = G.i
        for count in range(n):
            j = i*j
            H = j.tgt
            yield len(H.orbits())


class Hom(object):
    def __init__(self, src, tgt, send_perms):
        assert isinstance(src, Group)
        assert isinstance(tgt, Group)
        self.src = src
        self.tgt = tgt
        self.send_perms = send_perms
        if DEBUG:
            self.do_check()

    def __eq__(self, other):
        assert self.src == other.src
        return self.tgt == other.tgt and self.send_perms == other.send_perms

    def __ne__(self, other):
        assert self.src == other.src
        return self.tgt != other.tgt or self.send_perms != other.send_perms

    def do_check(self):
        src = self.src
        tgt = self.tgt
        lookup = src.lookup
        send_perms = self.send_perms
        for idx, p in enumerate(src):
          for jdx, q in enumerate(src):
            kdx = lookup[p*q]
            lhs = tgt[send_perms[kdx]]
            rhs = tgt[send_perms[idx]] * tgt[send_perms[jdx]]
            assert lhs == rhs, (lhs, rhs)

    def __add__(left, right):
        assert left.src == right.src
        G = left.src
        send_left = left.send_perms
        send_right = right.send_perms
        perms = []
        lrank = left.tgt.rank
        rrank = right.tgt.rank
        for idx, p in enumerate(G):
            l = left.tgt[send_left[idx]]
            r = right.tgt[send_right[idx]]
            perm = numpy.array(range(lrank+rrank))
            perm[:lrank] = perm[l.perm]
            perm[lrank:] = perm[r.perm + lrank]
            #perm = perm.copy()
            perm = Perm(perm)
            perms.append(perm)
        tgt = Group(perms) # ARGH, shuffles the order of perms
        send_perms = [tgt.lookup[perm] for perm in perms]
        hom = Hom(G, tgt, send_perms)
        return hom

    def __mul__(left, right):
        assert left.src == right.src
        G = left.src
        send_left = left.send_perms
        send_right = right.send_perms
        perms = []
        lrank = left.tgt.rank
        rrank = right.tgt.rank
        for idx, p in enumerate(G):
            l = left.tgt[send_left[idx]]
            r = right.tgt[send_right[idx]]
            perm = numpy.array(range(lrank*rrank))
            perm.shape = (lrank, rrank)
            perm = perm[l.perm, :]
            perm = perm[:, r.perm]
            perm = perm.copy()
            perm.shape = lrank*rrank
            perm = Perm(perm)
            perms.append(perm)
        tgt = Group(perms) # ARGH, shuffles the order of perms
        send_perms = [tgt.lookup[perm] for perm in perms]
        hom = Hom(G, tgt, send_perms)
        return hom


def general_linear(n=3, p=2):
    G = algebraic.GL(n, p)
    v = numpy.array([0]*n)
    v[0] = 1
    v = algebraic.Op(v, p)
    X = set([])
    for g in G:
        X.add(g*v)
    X = list(X)
    G = Group.from_action(G, X)
    return G



def test():

    G = Group.trivial(1)
    H = Group.trivial(2)
    assert (G.i + G.i).tgt == H

    G = Group.symmetric(3)
    assert len(G) == 6
    assert G.orbits() == [[0, 1, 2]]
    
    hom = G.regular_action()
    R = hom.tgt
    assert R != G
    assert R == R.regular_action().tgt

    X = G.i
    XX = X*X
    assert XX.tgt.orbits() == [[0, 4, 8], [1, 2, 3, 5, 6, 7]]
    assert ((X+X).tgt.orbits()) == [[0, 1, 2], [3, 4, 5]]

    G = Group.cyclic(5)
    H = (G.i * G.i).tgt
    assert len(H.orbits()) == 5

    G = Group.alternating(4)
    H = (G.i * G.i).tgt
    #print(H.orbits())

    G = Group.alternating(5)
    #H = (G.i * G.i).tgt
    #assert len(H.orbits()) == 2

    #for G in [Group.symmetric(5), Group.alternating(5)]:

    G = Group.symmetric(3)
    assert list(G.get_sequence(5)) == [1, 2, 5, 14, 41, 122]

    G = general_linear(3, 2)
    assert len(G) == 168

    for size in G.get_sequence():
        print(size, end=" ", flush=True)
    print()



if __name__ == "__main__":

    test()


