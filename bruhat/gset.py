#!/usr/bin/env python3

"""
Another attempt at action.py
Here we do everything with indexes.. probably more efficient than action.py
"""

import numpy
scalar = numpy.int64

from bruhat.util import factorial
from bruhat.action import mulclose
from bruhat import algebraic
from bruhat.argv import argv



DEBUG = argv.get("debug", False)

class Perm(object):

    def __init__(self, perm):
        assert isinstance(perm, list) or isinstance(perm, numpy.ndarray)
        perm = numpy.array(perm, dtype=scalar)
        self.perm = perm.copy()
        self.rank = len(perm)
        if DEBUG:
            self.do_check()

    def do_check(self):
        items = list(self.perm)
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
        if isinstance(other, Perm):
            perm = self.perm[other.perm] # check this is the right order!
            result = Perm(perm)
        elif isinstance(other, (int, scalar)):
            result = self.perm[other]
        else:
            raise TypeError
        return result

    def inverse(self):
        perm = self.perm
        q = [None]*self.rank
        for (src, tgt) in enumerate(perm):
            q[tgt] = src
        return Perm(q)

    def __getitem__(self, i):
        return self.perm[i]

    #def orbits(self):

def compose(f, g):
    return [f[gi] for gi in g]


class Group(object):
    """
        A concrete group of permutations masquerading as an abstract group.
        This gets confusing once we start constructing GSet's below as
        actions of some abstract group on a concrete group. There is not much
        to distinguish these two kinds of group apart from your imagination.
    """
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
        assert len(self.lookup) == self.n
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
        return "Group(order=%s, rank=%s)"%(self.n, self.rank)
    __repr__ = __str__

    def longstr(self):
        return "Group(%s)"%(self.perms,)

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
    def dihedral(cls, n):
        assert n>0
        perms = []
        perm = [(i+1)%n for i in range(n)]
        perms.append(Perm(perm))
        perm = [(-i)%n for i in range(n)]
        perms.append(Perm(perm))
        G = Group(gen=perms)
        #assert len(G) == 2*n
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
        X = GSet(self, tgt, send_perms)
        return X

    @property
    def i(self):
        send_perms = list(range(self.n)) # ARGH should we use dict's for these?
        gset = GSet(self, self, send_perms)
        return gset
        
    def get_orbits(self):
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

    def get_components(self):
        orbits = self.get_orbits()
        Gs = []
        for orbit in orbits:
            G = Group.from_action(self, orbit)
            assert G.n <= self.n
            assert self.n % G.n == 0
            Gs.append(G)
        return Gs

    def get_sequence(G, n=5):
        yield len(G.get_orbits())
        i = j = G.i
        for count in range(n):
            j = GSet.mul(i, j).apex
            H = j.tgt
            yield len(H.get_orbits())


class GSet(object):
    """
        A morphism of concrete groups.
        This is (also) a G-set where the src is G.
    """
    def __init__(self, src, tgt, send_perms):
        assert isinstance(src, Group)
        assert isinstance(tgt, Group)
        self.src = src
        self.tgt = tgt
        self.send_perms = send_perms
        self.rank = tgt.rank
        if DEBUG:
            self.do_check()

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

    def __eq__(self, other):
        assert self.src == other.src
        return self.tgt == other.tgt and self.send_perms == other.send_perms

    def __ne__(self, other):
        assert self.src == other.src
        return self.tgt != other.tgt or self.send_perms != other.send_perms

    def __str__(self):
        #return "GSet(%s, %s, %s)"%(self.src.n, self.tgt.rank, self.send_perms)
        return "GSet(%s, %s, %s)"%(self.src, self.tgt, self.send_perms)
    __repr__ = __str__

    def compose(self, other):
        assert self.tgt == other.src
        a = self.send_perms
        b = other.send_perms
        send_perms = [b[i] for i in a]
        return GSet(self.src, other.tgt, send_perms)

    @classmethod
    def from_action(cls, G, X):
        lookup = dict((v, idx) for (idx, v) in enumerate(X))
        perms = []
        for g in G:
            perm = [lookup[g*v] for v in X]
            perms.append(Perm(perm))
        tgt = Group(perms)
        send_perms = [tgt.lookup[perm] for perm in perms]
        gset = cls(G, tgt, send_perms)
        return gset

    def get_orbits(self):
        src = self.src
        tgt = self.tgt
        send_perms = self.send_perms
        remain = set(range(tgt.rank))
        orbits = []
        while remain:
            i = iter(remain).__next__()
            remain.remove(i)
            orbit = set([i])
            for idx in send_perms:
                g = tgt[idx]
                j = g[i]
                while j != i:
                    orbit.add(j)
                    j = g[j]
            remain.difference_update(orbit)
            orbit = list(orbit)
            orbit.sort()
            orbits.append(orbit)
        return orbits

    def get_components(self):
        orbits = self.get_orbits()
        nats = []
        for orbit in orbits:
            gset = GSet.from_action(self.tgt, orbit)
            gset = self.compose(gset)
            send_items = list(orbit)
            nat = Nat(gset, self, send_items)
            nats.append(nat)
        return nats

    def add(left, right):
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
            perm = numpy.array(range(lrank+rrank), dtype=scalar)
            perm[:lrank] = perm[l.perm]
            perm[lrank:] = perm[r.perm + lrank]
            #perm = perm.copy()
            perm = Perm(perm)
            perms.append(perm)
        tgt = Group(perms) # ARGH, shuffles the order of perms
        send_perms = [tgt.lookup[perm] for perm in perms]
        gset = GSet(G, tgt, send_perms)
        send_items = [i for i in range(lrank)]
        left = Nat(left, gset, send_items)
        send_items = [i+lrank for i in range(rrank)]
        right = Nat(right, gset, send_items)
        return Cone(gset, [left, right], contra=True)

    def mul(left, right):
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
            perm = numpy.array(range(lrank*rrank), dtype=scalar)
            perm.shape = (lrank, rrank)
            perm = perm[l.perm, :]
            perm = perm[:, r.perm]
            perm = perm.copy()
            perm.shape = lrank*rrank
            perm = Perm(perm)
            perms.append(perm)
        tgt = Group(perms) # ARGH, shuffles the order of perms
        send_perms = [tgt.lookup[perm] for perm in perms]
        gset = GSet(G, tgt, send_perms)
        send_items = [i for i in range(lrank) for j in range(rrank)]
        left = Nat(gset, left, send_items)
        send_items = [j for i in range(lrank) for j in range(rrank)]
        right = Nat(gset, right, send_items)
        return Cone(gset, [left, right])


class Cone(object):
    def __init__(self, apex, legs, contra=False):
        self.apex = apex
        self.legs = list(legs)
        self.contra = contra

        for leg in legs:
            if contra:
                assert leg.tgt == apex
            else:
                assert leg.src == apex

    def __getitem__(self, idx):
        return self.legs[idx]



class Nat(object):
    """
        A morphism of G-sets. 
    """

    def __init__(self, src, tgt, send_items):
        assert isinstance(src, GSet)
        assert isinstance(tgt, GSet)
        G = src.src
        assert G == tgt.src
        self.src = src
        self.tgt = tgt
        self.G = G
        assert len(send_items) == src.tgt.rank
        self.send_items = list(send_items)
        if DEBUG:
            self.do_check()

    def do_check(self):
        src = self.src
        tgt = self.tgt
        send_items = self.send_items
        items = set(send_items)
        assert items.issubset(set(range(tgt.tgt.rank)))
        G = self.G
        for i in range(G.n):
            lg = src.tgt[src.send_perms[i]]
            left = compose(send_items, lg)
            rg = tgt.tgt[tgt.send_perms[i]]
            right = compose(rg, send_items)
            assert left == right

    def __str__(self):
        return "Nat(%s, %s, %s)"%(self.src, self.tgt, self.send_items)
    __repr__ = __str__


def general_linear(n=3, p=2):
    G = algebraic.GL(n, p)
    v = numpy.array([0]*n, dtype=scalar)
    v[0] = 1
    v = algebraic.Op(v, p)
    X = set([])
    for g in G:
        X.add(g*v)
    X = list(X)
    G = Group.from_action(G, X)
    return G



def test():

    add, mul = GSet.add, GSet.mul

    G = Group.trivial(1)
    H = Group.trivial(2)
    assert add(G.i, G.i).apex.tgt == H

    G = Group.symmetric(3)
    assert len(G) == 6
    assert G.get_orbits() == [[0, 1, 2]]
    
    gset = G.regular_action()
    R = gset.tgt
    assert R != G
    assert R == R.regular_action().tgt

    X = G.i
    XX = mul(X, X).apex
    X_X = add(X, X).apex
    assert XX.tgt.get_orbits() == [[0, 4, 8], [1, 2, 3, 5, 6, 7]]
    assert (X_X.tgt.get_orbits()) == [[0, 1, 2], [3, 4, 5]]

    G = Group.cyclic(3)
    H = add(G.i, G.i).apex
    for nat in H.get_components():
        pass

    G = Group.cyclic(5)
    H = mul(G.i, G.i).apex.tgt
    assert len(H.get_orbits()) == 5

    G = Group.alternating(5)
    H = mul(G.i, G.i).apex
    H = H.tgt
    Gs = H.get_components()
    assert len(Gs) == 2
    assert Gs[0].rank == 5
    assert Gs[1].rank == 20

    G = Group.alternating(5)
    H = mul(G.i, G.i).apex
    H = mul(G.i, H).apex
    items = [nat.src for nat in H.get_components()]
    A, B, C = items[0], items[1], items[4]
    assert A.rank == 5
    assert B.rank == 20
    assert C.rank == 60
    AB = mul(A, B).apex
    BC = mul(B, C).apex
    AB_C = mul(AB, C).apex
    A_BC = mul(A, BC).apex
    assert AB_C == A_BC # strict monoidal 

    G = Group.alternating(5)

    GG = mul(G.i, G.i).apex
    GG_G = mul(GG, G.i).apex
    G_GG = mul(G.i, GG).apex
    assert GG_G == G_GG # strict monoidal 

    G = Group.dihedral(5)

    G = Group.symmetric(3)
    assert list(G.get_sequence(5)) == [1, 2, 5, 14, 41, 122]

    G = general_linear(3, 2)
    assert len(G) == 168

    #for size in G.get_sequence():
    #    print(size, end=" ", flush=True)
    #print()


    print("OK")


if __name__ == "__main__":

    test()


