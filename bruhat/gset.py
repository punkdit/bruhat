#!/usr/bin/env python3

"""
Compute limits & colimits in the category of GSet's .

Based on bruhat/action.py, but 
here we do everything (set maps, permutations, etc.)
with indexes, which is probably more efficient than action.py.

It's also much more difficult using this than bruhat/action.py,
there's no typechecking on indexes. Arghhh.
We really should have a class for each set, and to typechecking
as in lin.py with the Space class etc...

see also: set.py ?

see also: relation.py
"""


from functools import reduce, cache
from operator import add, mul
from string import ascii_letters

import numpy
from numpy import all as alltrue
scalar = numpy.int64

from bruhat.algebraic import Algebraic, Matrix
from bruhat.util import factorial, cross, is_prime, choose
from bruhat.action import mulclose, mulclose_hom
from bruhat.equ import Equ
from bruhat.argv import argv


if argv.Z2:
    from bruhat import solve
    zeros = solve.zeros2
    array = solve.array2
    dot = solve.dot2
    rank = solve.rank
    nullity = solve.nullity
    shortstr = solve.shortstr

else:
    from bruhat import gelim
    zeros = gelim.zeros
    array = gelim.array
    dot = gelim.dot
    rank = gelim.rank
    nullity = gelim.nullity
    shortstr = gelim.shortstr


debug = argv.get("debug", False)


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



class Perm(object):

    def __init__(self, perm):
        assert isinstance(perm, list) or isinstance(perm, numpy.ndarray)
        perm = numpy.array(perm, dtype=scalar)
        self.perm = perm.copy()
        self.rank = len(perm)
        if debug:
            self.do_check()

    def do_check(self):
        items = list(self.perm)
        items.sort()
        assert items == list(range(self.rank))

    @property
    def identity(self):
        return Perm(list(range(self.rank)))

    def is_identity(self):
        for i, ii in enumerate(self.perm):
            if i != ii:
                return False
        return True

    def __str__(self):
        return "Perm(%s)"%([int(i) for i in self.perm])
    __repr__ = __str__

    def __hash__(self):
        return hash(self.perm.tobytes())

    def __eq__(self, other):
        if self.rank != other.rank:
            return False
        return self.perm.tobytes() == other.perm.tobytes()
        #return numpy.alltrue(self.perm == other.perm) # bottleneck for .subgroups
        #return numpy.allclose(self.perm, other.perm) # argh, much slower !

#    def __ne__(self, other):
#        if self.rank != other.rank:
#            return True
#        return not numpy.allclose(self.perm, other.perm)

    def __lt__(self, other):
        return self.perm.tobytes() < other.perm.tobytes()

    def __le__(self, other):
        return self.perm.tobytes() <= other.perm.tobytes()

    def __mul__(self, other):
        if isinstance(other, Perm):
            perm = self.perm[other.perm] # check this is the right order!
            result = Perm(perm)
        elif isinstance(other, (int, scalar)):
            result = self.perm[other]
        else:
            result = NotImplemented
        return result

    def __pow__(self, n):
        assert n>=0, "NotImplemented"
        if n==0:
            return self.identity
        g = self
        while n>1:
            g = self*g
            n -= 1
        return g

    def sign(self):
        perm = self.perm
        n = self.rank
        s = +1
        for i in range(n):
          for j in range(i+1,n):
            if perm[i]>perm[j]:
                s *= -1
        return s

    def cross(left, right):
        assert isinstance(right, Perm)
        m, n = left.rank, right.rank
        lperm, rperm = left.perm, right.perm
        perm = [l + r*m for r in rperm for l in lperm]
        return Perm(perm)

    def inverse(self):
        perm = self.perm
        q = [None]*self.rank
        for (src, tgt) in enumerate(perm):
            q[tgt] = src
        return Perm(q)
    __invert__ = inverse

    def __getitem__(self, i):
        return int(self.perm[i])

    def __len__(self):
        return len(self.perm)

    def fixed(self):
        perm = self.perm
        return [i for i in range(len(perm)) if perm[i]==i]

    def order(self):
        i = 1
        g = self
        while not g.is_identity():
            g = self*g
            i += 1
            #assert i <= len(self.items)+1 # how big can this get ??
        return i

    def get_orbits(self):
        perm = self.perm
        remain = set(int(i) for i in perm)
        orbits = []
        while remain:
            i = iter(remain).__next__()
            remain.remove(i)
            orbit = set([i])
            j = self[i]
            while j != i:
                orbit.add(j)
                j = self[j]
            remain.difference_update(orbit)
            orbit = list(orbit)
            orbit.sort()
            orbits.append(orbit)
        return orbits

    def cycles(self):
        remain = set(self.perm)
        cycles = []
        while remain:
            item = iter(remain).__next__()
            orbit = [item]
            #print(type(self), type(item), "__mul__")
            item1 = self[item]
            while item1 != item:
                orbit.append(item1)
                item1 = self[item1]
                assert len(orbit) <= self.rank
            assert orbit
            cycles.append(orbit)
            n = len(remain)
            for item in orbit:
                remain.remove(item)
            assert len(remain) < n
        return cycles

    def gap_fmt(perm):
        cs = perm.cycles()
        ss = []
        for c in cs:
            if len(c)==1:
                continue
            c = [int(i+1) for i in c]
            ss.append(tuple(c))
        return ''.join(str(s) for s in ss)


def gap_code(perms):
    s = [perm.gap_fmt() for perm in perms]
    s = "Group(%s);"%(', '.join(s))
    return s





def compose(f, g):
    return [f[gi] for gi in g]


class Coset(object):
    def __init__(self, perms):
        perms = list(perms)
        perms.sort() # yes !
        self.perms = perms
        self.data = numpy.array([perm.perm for perm in perms], dtype=scalar)
        self.n = len(self.perms)
        self.rank = perms[0].rank if perms else 0
        #self._hash = hash(tuple(self.perms))
        self._hash = hash(self.data.tobytes())

    def __str__(self):
        return "Coset(%d)"%(self.n,)
    __repr__ = __str__

    def __len__(self):
        return self.n

    def __getitem__(self, idx):
        return self.perms[idx]

    #def __contains__(self, g):
    #    return g in self.lookup

    def __eq__(self, other):
        assert isinstance(other, Coset)
        return alltrue(self.data == other.data)

    def __hash__(self):
        return self._hash

    def __lt__(self, other):
        return self.data.tobytes() < other.data.tobytes()

    def left_mul(self, g):
        assert g.rank == self.rank
        perms = [g*h for h in self]
        return Coset(perms)

    def left_mul_perm(self, g):
        assert g.rank == self.rank
        perms = [g*h for h in self]
        idxs = Perm([self.perms.index(perm) for perm in perms])
        return idxs

    def right_mul(self, g):
        assert g.rank == self.rank
        perms = [h*g for h in self]
        return Coset(perms)

    def intersect(self, other):
        perms = set(self.perms).intersection(set(other.perms))
        return Coset(perms)


class Group(object):
    """
        A concrete group of permutations masquerading as an abstract group.
        This gets confusing once we start constructing GSet's below as
        actions of some abstract group on a concrete group. There is not much
        to distinguish these two kinds of group apart from your imagination.
    """
    def __init__(self, perms=None, gens=None, items=None, verbose=False, build=True):
        if perms is None:
            assert gens
            assert isinstance(gens[0], Perm)
            self.gens = list(gens)
            if build:
                perms = list(mulclose(gens, verbose=verbose))
        else:
            perms = list(perms)
            self.gens = perms
        self.rank = int(self.gens[0].rank)
        self.identity = Perm(list(range(self.rank)))
        if not perms:
            return # <---------------- keep only gens <------ 
        assert perms
        canonical = list(perms)
        canonical.sort()
        self.canonical = canonical
        self.perms = perms
        self.n = len(self.perms)
        self.lookup = dict((perm, idx) for (idx, perm) in enumerate(self.perms))
        if items is None:
            items = list(range(self.rank))
        self.items = items
        self._str = str(self.canonical)
        self._hash = hash(self._str)
        self._subgroups = None # cache
        assert len(self.lookup) == self.n, (self.n, len(self.lookup))
        if debug:
            self.do_check()

    def do_check(self):
        lookup = self.lookup
        assert self.identity in lookup
        for g in self:
            for h in self:
                assert g*h in lookup
            assert g.inverse() in lookup
            assert g.inverse() * g == self.identity

    @classmethod
    def generate(cls, gens, verbose=False):
        return cls(None, gens)

    @cache
    def get_gens(self):
        perms = self.perms
        n = len(perms)
        for count in range(1,n):
          for gens in choose(perms, count):
            G = mulclose(gens)
            if len(G) == n:
                self.gens = list(gens)
                return self.gens

#    @property
#    def identity(self):
#        return Perm(list(range(self.rank)))

    def rankstr(self):
        return "Group(order=%s, rank=%s)"%(self.n, self.rank)

    def __str__(self):
        #return "Group(order=%s, rank=%s)"%(self.n, self.rank)
        return "Group(order=%s)"%(self.n,)
    __repr__ = __str__

    def longstr(self):
        return "Group(%s)"%(self.perms,)

    def __len__(self):
        return self.n

    def __getitem__(self, idx):
        return self.perms[idx]
        #return self.canonical[idx] # breaks test...

    def __contains__(self, g):
        return g in self.lookup

    def __eq__(self, other):
        return self.canonical == other.canonical

    def __ne__(self, other):
        return self.canonical != other.canonical

    def __lt__(self, other):
        return self.canonical < other.canonical

    def __hash__(self):
        return self._hash

    def __mul__(self, other):
        assert isinstance(other, Group)
        perms = [g.cross(h) for g in self for h in other]
        return Group(perms)

    @classmethod
    def trivial(cls, n):
        assert n>0
        perm = Perm(list(range(n)))
        return Group([perm])

    @classmethod
    def symmetric(cls, n):
        assert n>0
        gens = []
        items = list(range(n))
        for i in range(n-1):
            perm = list(items)
            perm[i] = i+1
            perm[i+1] = i
            gens.append(Perm(perm))
        G = Group(gens=gens)
        assert len(G) == factorial(n)
        return G

    @classmethod
    def alternating(cls, n):
        assert n>0
        gens = []
        items = list(range(n))
        for i in range(n-2):
            perm = list(items)
            perm[i] = i+1
            perm[i+1] = i+2
            perm[i+2] = i
            gens.append(Perm(perm))
        G = Group(gens=gens)
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
        G = Group(gens=perms)
        #assert len(G) == 2*n
        return G

    @classmethod
    def coxeter_bc(cls, n, check=False):
        Lookup = lambda perm, items : Perm([items.index(perm[item]) for item in items])
        items = [(i,0) for i in range(n)]
        items += [(i,1) for i in range(n)]
        perms = []
        for i in range(n-1):
            jtems = list(items)
            jtems[i:i+2] = jtems[i+1], jtems[i]
            jtems[n+i:n+i+2] = jtems[n+i+1], jtems[n+i]
            perms.append(Lookup(dict(zip(items, jtems)), items))
        jtems = list(items)
        jtems[n-1], jtems[2*n-1] = jtems[2*n-1], jtems[n-1]
        perms.append(Lookup(dict(zip(items, jtems)), items))
        G = Group.generate(perms)
        return G

    @classmethod
    def coxeter_d(cls, n, check=False):
        Lookup = lambda perm, items : Perm([items.index(perm[item]) for item in items])
        items = [(i,0) for i in range(n)]
        items += [(i,1) for i in range(n)]
        perms = []
        for i in range(n-2):
            jtems = list(items)
            jtems[i:i+3] = jtems[i+1], jtems[i+2], jtems[i]
            jtems[n+i:n+i+3] = jtems[n+i+1], jtems[n+i+2], jtems[n+i]
            perms.append(Lookup(dict(zip(items, jtems)), items))
        jtems = list(items)
        jtems[n-1], jtems[2*n-1] = jtems[2*n-1], jtems[n-1]
        perms.append(Lookup(dict(zip(items, jtems)), items))
        G = Group.generate(perms)
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

    @classmethod
    def from_table(cls, table):
        n = len(table)
        assert table.shape == (n,n)
        perms = [Perm(row) for row in table]
        G = Group(perms)
        return G

    @classmethod
    def from_concrete(cls, items):
        items = list(items)
        lookup = {item:i for (i,item) in enumerate(items)}
        G = Group([Perm([lookup[item*jtem] for jtem in items]) for item in items])
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
        tgt = Group(perms)
        #send_perms = list(range(self.n))
        #assert send_perms == [tgt.lookup[perm] for perm in perms]
        X = GSet(self, tgt)
        return X

    @property
    def i(self):
        #send_perms = list(range(self.n)) # ARGH should we use dict's for these?
        gset = GSet(self, self)
        return gset

    def trivial_gset(self, n=1):
        I = Perm(list(range(n)))
        tgt = Group([I])
        send_perms = [0]*len(self)
        X = GSet(self, tgt, send_perms)
        return X
        
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

    def get_atoms(self):
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

    def cyclic_subgroups(self, verbose=False):
        # find all cyclic subgroups
        I = self.identity
        trivial = Group([I])
        cyclic = set([trivial])
        for g0 in self:
            if g0==I:
                continue
            perms = [g0]
            g1 = g0
            while 1:
                g1 = g0*g1
                if g1==g0:
                    break
                perms.append(g1)
            group = Group(perms)
            assert len(group)>1
            cyclic.add(group)
        return cyclic

    def subgroups(self, verbose=False):
        if self._subgroups is not None:
            #return list(self._subgroups)
            for H in self._subgroups:
                yield H
            return
        I = self.identity
        trivial = Group([I])
        cyclic = self.cyclic_subgroups()
        #if verbose:
        #    print "Group.subgroups: cyclic:", len(cyclic)
        n = len(self) # order
        items = set(cyclic)
        items.add(trivial)
        items.add(self)
        for H in items:
            yield H
        if verbose:
            print("/", flush=True, end="")
        cyclic = list(cyclic)
        cyclic.sort(key = lambda G:(len(G), G._str)) # make this deterministic, (and fast..)
        bdy = list(G for G in cyclic if len(G)<n)
        while bdy:
            _bdy = list()
            # consider each group in bdy
            for G in bdy:
                if is_prime(n//len(G)):
                    # prime index
                    continue
                # enlarge by a cyclic subgroup
                for H in cyclic:
                    perms = set(G.perms+H.perms)
                    k = len(perms)
                    if k==n or k==len(G) or k==len(H):
                        continue
                    perms = mulclose(perms)
                    K = Group(perms)
                    if K not in items:
                        yield K
                        i = n//len(K)
                        _bdy.append(K)
                        items.add(K)
                        if verbose:
                            print('(%s)'%i, flush=True, end="")
                    #else:
                    #    print('/')
            bdy = _bdy
            #if verbose:
            #    print "items:", len(items)
            #    print "bdy:", len(bdy)
        if verbose:
            print()
        items = list(items)
        #items.sort(key = len)
        self._subgroups = items
        #return list(items)

    @cache
    def conjugacy_classes(G):
        remain = set(G)
        cgys = []
        while remain:
            g = iter(remain).__next__()
            remain.remove(g)
            cgy = [g]
            for h in G:
                k = h*g*(~h)
                if k in remain:
                    remain.remove(k)
                    cgy.append(k)
            cgy.sort()
            cgys.append(cgy)
        cgys.sort(key = lambda cgy:(len(cgy), cgy[0]))
        cgy_lookup = {}
        for (i,cgy) in enumerate(cgys):
            for g in cgy:
                cgy_lookup[g] = i
        G.cgy_lookup = cgy_lookup
        return cgys

    @cache
    def conjugacy_subgroups(G, Hs=None, sort=False, verbose=False):

        # Find all conjugacy classes of subgroups

        if Hs is None:
            Hs = list(G.subgroups(verbose=verbose))

        #print("conjugacy_subgroups: subgroups", len(Hs))

        equs = dict((H1, Equ(H1)) for H1 in Hs)
        for H1 in Hs:
            for g in G:
                if g in H1:
                    continue
                H2 = H1.conjugate(g)
                #print("%s == %s ? %s" % (H2, H1, H2==H1))
                if H2 == H1:
                    continue
                else:
                    #print len(H1), "~", len(H2)
                    if H2 not in equs:
                        equs[H2] = Equ(H2)
                    equs[H1].merge(equs[H2])

        # get equivalance classes
        equs = list(set(equ.top for equ in equs.values()))
        equs.sort(key = lambda equ : (-len(equ.items[0]), equ.items[0]._str))
        for equ in equs:
            #print "equ:", [len(H) for H in equ.items]
            for H in equ.items:
                H.conjugates = list(equ.items)
        #print "total:", len(equs)
        Hs = [equ.items[0] for equ in equs] # pick unique (up to conjugation)
        if sort:
            Hs.sort(key = lambda H : (-len(H), H._str))
        return Hs

    def conjugate(G, g):
        gi = ~g
        perms = [g*h*gi for h in G]
        return Group(perms)

    def is_conjugate_subgroup(G, H, J):
        if len(H) != len(J):
            return
        for g in G:
            gi = ~g
            perms = [g*h*gi for h in H]
            K = Group(perms)
            if K==J:
                return True
        return False

    def is_normal(G, H):
        for g in G:
            if H.conjugate(g) != H:
                return False
        return True

    def is_simple(G):
        for H in G.subgroups():
            if len(H)==1 or len(H)==len(G):
                continue
            if G.is_normal(H):
                return False
        return True

    def is_abelian(G):
        for g in G:
          for h in G:
            if g*h != h*g:
                return False
        return True

    def is_cyclic(G):
        for g in G:
            if G == Group.generate([g]):
                return True
        return False

    def cyclic_generator(G):
        for g in G:
            if G == Group.generate([g]):
                return g
        assert 0, "not cyclic"

    def intersect(G, H):
        G = set(G.lookup.keys())
        H = set(H.lookup.keys())
        GH = G.intersection(H)
        return Group(GH)

    def is_subgroup(G, H):
        assert isinstance(H, Group)
        for h in H:
            if h not in G:
                return False
        return True

    def action_subgroup(G, H):
        assert isinstance(H, Group)
        assert H.rank == G.rank
        assert G.is_subgroup(H)
        H0 = H
        H = Coset(H)
        cosets = set([H])
        remain = set(G)
        remain.difference_update(H)
        while remain:
            g = iter(remain).__next__()
            remain.remove(g)
            gH = H.left_mul(g)
            cosets.add(gH)
            remain.difference_update(gH)
        assert len(cosets) * len(H) == len(G)
        cosets = list(cosets)
        cosets.sort()
        assert cosets == G.left_cosets(H0)
        lookup = dict((gH, idx) for (idx, gH) in enumerate(cosets))
        perms = set()
        all_perms = []
        for h in G:
            perm = Perm([lookup[gH.left_mul(h)] for gH in cosets])
            perms.add(perm)
            all_perms.append(perm)
        H = Group(perms)
        send_perms = [H.lookup[perm] for perm in all_perms]
        return GSet(G, H, send_perms)

    def act_subgroup(G, H):
        assert isinstance(H, Group)
        assert H.rank == G.rank
        assert G.is_subgroup(H)
        H0 = H
        H = Coset(H)
        cosets = set([H])
        remain = set(G)
        remain.difference_update(H)
        while remain:
            g = iter(remain).__next__()
            remain.remove(g)
            gH = H.left_mul(g)
            cosets.add(gH)
            remain.difference_update(gH)
        assert len(cosets) * len(H) == len(G)
        cosets = list(cosets)
        cosets.sort()
        assert cosets == G.left_cosets(H0)
        lookup = dict((gH, idx) for (idx, gH) in enumerate(cosets))
        gens = []
        for h in G.gens:
            perm = Perm([lookup[gH.left_mul(h)] for gH in cosets])
            gens.append(perm)
        #H = Group(perms)
        #send_perms = [H.lookup[perm] for perm in all_perms]
        #return GSet(G, H, send_perms)
        return Group.generate(gens)

    def left_cosets(G, H):
        assert isinstance(H, Group)
        assert G.is_subgroup(H)
        H = Coset(H)
        cosets = set([H])
        for g in G:
            gH = H.left_mul(g)
            cosets.add(gH)
        cosets = list(cosets)
        cosets.sort()
        return cosets

    def distribution(G):
        dist = [[] for _ in range(len(G)+1)]
        for g in G:
            dist[g.order()].append(g)
        return dist

    def is_isomorphic(G, H):
        hom = G.isomorphism(H)
        return hom is not None
    isomorphic = is_isomorphic

    def isomorphism(G, H):
        "a not very clever search for a group isomorphism"
        assert isinstance(H, Group)
        if len(G) != len(H):
            return

        if G is H:
            return {g:g for g in G}

        if len(G) > 64:
            assert 0, "too big i think ... ?!?"

        left = G.distribution()
        right = H.distribution()
        l = [len(items) for items in left]
        r = [len(items) for items in right]
        if l != r:
            return
        #print(l)
        #print(r)

        pairs = [(g,h) for g in G for h in G]
        pair = None
        for (g,h) in pairs:
            K = mulclose([g,h])
            if len(K) < len(G):
                continue
            pair = (g,h)

        if pair is None:
            assert 0, "todo.. look for generating triples, etc."

        lg, lh = pair
        i,j = lg.order(), lh.order()
        for rg in right[i]:
          for rh in right[j]:
            hom = mulclose_hom( [lg, lh], [rg, rh] )
            #print(hom is not None, end=" ")
            if hom is None:
                continue
            assert len(hom) == len(G)
            if len(set(hom.values())) < len(G):
                return None
            return hom

    def count_hecke(self, other):
        assert isinstance(other, Group)
        lgen = self.gens
        rgen = other.gens
        assert len(lgen) == len(rgen)
        assert isinstance(lgen, list)
        assert isinstance(rgen, list)
        M = lgen[0].rank # rows
        N = rgen[0].rank # cols
    
        pairs = list(zip(lgen, rgen))
    
        remain = set(range(M)) # rows
        counts = []
        while remain:
            row = remain.pop()
    
            # now find the orbit of (row,0)
            bdy = [(row,0)]
            found = set(bdy)
            count = 1
            while bdy:
                if argv.verbose:
                    print("%d:%d"%(len(found),len(bdy)), end=" ", flush=True)
                _bdy = []
                while bdy:
                    i,j = bdy.pop()
                    for l,r in pairs:
                        tgt = l[i], r[j]
                        if tgt in found:
                            continue
                        found.add(tgt)
                        _bdy.append(tgt)
                        if tgt[1] == 0:
                            remain.remove(tgt[0])
                            count += 1
                bdy = _bdy
            if argv.verbose:
                print()
            #print("[%s]" % count, end="", flush=True)
            counts.append(count)
        #print()
    
        counts.sort()
        #print("orbits:", counts, len(counts))
        assert sum(counts) == M
    
        return len(counts)
        
    def count_hecke_injections(self, other, max_count=None):
        "this counts (i believe/hope) a single entry in the table of marks"
        assert isinstance(other, Group)
        lgen = self.gens
        rgen = other.gens
        assert len(lgen) == len(rgen) # require same parent ?
        assert isinstance(lgen, list)
        assert isinstance(rgen, list)
        M = lgen[0].rank # rows
        N = rgen[0].rank # cols
    
        pairs = list(zip(lgen, rgen))
    
        remain = set(range(M)) # rows
        counts = []
        while remain:
            row = remain.pop()
    
            # now find the orbit of (row,0)
            bdy = [(row,0)]
            found = set(bdy)
            count = 1
            while bdy and count == 1:
                if argv.verbose:
                    print("%d:%d"%(len(found),len(bdy)), end=" ", flush=True)
                _bdy = []
                while bdy and count == 1:
                    i,j = bdy.pop()
                    for l,r in pairs:
                        tgt = l[i], r[j]
                        if tgt in found:
                            continue
                        found.add(tgt)
                        _bdy.append(tgt)
                        if tgt[1] == 0:
                            if tgt[0] in remain:
                                remain.remove(tgt[0])
                            count += 1
                            break
                bdy = _bdy
            if argv.verbose:
                print()
            #print("[%s]" % count, end="", flush=True)
            if count == 1:
                counts.append(count)
                if max_count:
                    return 1
        return len(counts)




class GSet(object):
    """
        A morphism of concrete groups.
        This is (also) a GSet, where the src is G.
    """
    def __init__(self, src, tgt, send_perms=None):
        assert isinstance(src, Group)
        assert isinstance(tgt, Group)
        self.src = src
        self.G = src
        self.tgt = tgt
        if send_perms is None:
            assert len(src) == len(tgt)
            send_perms = list(range(len(src)))
        self.send_perms = send_perms
        self.rank = tgt.rank # the size of the set we act on
        self._sig = None
        if debug:
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

    def __eq__(self, other): # equality on the nose
        assert self.src == other.src
        return self.tgt == other.tgt and self.send_perms == other.send_perms

    def __ne__(self, other):
        assert self.src == other.src
        return self.tgt != other.tgt or self.send_perms != other.send_perms

    def __str__(self):
        #return "GSet(%s, %s, %s)"%(self.src.n, self.tgt.rank, self.send_perms)
        #return "GSet(%s, %s, %s)"%(self.src, self.tgt, self.send_perms)
        return "GSet(%s, rank=%s)"%(self.src, self.rank)
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

    def get_stabilizer(self, fix=0):
        src = self.src
        tgt = self.tgt
        send_perms = self.send_perms
        stab = []
        for idx, i in enumerate(send_perms):
            if tgt[i][fix] == fix:
                stab.append(src[idx])
        return Group(stab)

    def get_orbits(self):
        #print("get_orbits", self)
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

    def get_atoms(self):
        orbits = self.get_orbits()
        atoms = []
        for orbit in orbits:
            gset = GSet.from_action(self.tgt, orbit)
            gset = self.compose(gset)
            send_items = list(orbit)
            nat = Hom(gset, self, send_items)
            atoms.append(nat)
        return atoms

    def get_retract_linear(self): # random name ...
        #atoms = self.get_atoms()
        orbits = self.get_orbits()
        n = len(orbits)
        m = self.rank
        #print(orbits)
        #points = self.src.trivial_gset(n)
        proj = zeros(n, m)
        sect = zeros(m, n)
        for i, orbit in enumerate(orbits):
            for j in orbit:
                proj[i, j] = 1
            sect[j, i] = 1
        return proj, sect

    def get_retract(self): # random name ...
        #atoms = self.get_atoms()
        orbits = self.get_orbits()
        proj = [None]*self.rank
        sect = [None]*len(orbits)
        for i, orbit in enumerate(orbits):
            for j in orbit:
                proj[j] = i
            sect[i] = j
        return proj, sect

    def get_sequence(self, n=5):
        yield len(self.get_orbits())
        i = j = self
        for count in range(n):
            j = GSet.mul(i, j).apex
            yield len(j.get_orbits())

    def add(left, right, cone=None):
        assert left.src == right.src
        G = left.src
        send_left = left.send_perms
        send_right = right.send_perms
        perms = []
        lrank = left.rank
        rrank = right.rank
        for idx, p in enumerate(G):
            l = left.tgt[send_left[idx]]
            r = right.tgt[send_right[idx]]
            perm = numpy.array(range(lrank+rrank), dtype=scalar)
            perm[:lrank] = perm[l.perm]
            perm[lrank:] = perm[r.perm + lrank]
            #perm = perm.copy()
            perm = Perm(perm)
            perms.append(perm)
        tgt = Group(perms)
        gset = GSet(G, tgt)
        send_items = [i for i in range(lrank)]
        i_left = Hom(left, gset, send_items)
        send_items = [i+lrank for i in range(rrank)]
        i_right = Hom(right, gset, send_items)
        limit = Cone(gset, [i_left, i_right], contra=True) # _cocone
        if cone is None:
            return limit

        assert cone.contra
        assert cone[0].src is left
        assert cone[1].src is right

        l, r = cone
        apex = cone.apex
        send_items = [l.send_items[i] for i in range(l.rank)] + [r.send_items[i] for i in range(r.rank)]
        univ = Hom(gset, apex, send_items)
        return limit, univ

    def mul(left, right, cone=None):
        assert left.src == right.src
        G = left.src
        send_left = left.send_perms
        send_right = right.send_perms
        perms = []
        lookup = {}
        lrank = left.rank
        rrank = right.rank
        send_perms = []
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
            if perm not in lookup:
                lookup[perm] = len(perms)
                perms.append(perm)
            send_perms.append(lookup[perm])
        tgt = Group(perms)
        gset = GSet(G, tgt, send_perms)
        send_items = [i for i in range(lrank) for j in range(rrank)]
        p_left = Hom(gset, left, send_items)
        send_items = [j for i in range(lrank) for j in range(rrank)]
        p_right = Hom(gset, right, send_items)
        limit = Cone(gset, [p_left, p_right])
        if cone is None:
            return limit

        assert not cone.contra
        assert cone[0].tgt is left
        assert cone[1].tgt is right

        l, r = cone
        apex = cone.apex
        send_items = [l.send_items[i]*rrank + r.send_items[i] for i in range(apex.rank)]
        univ = Hom(apex, gset, send_items)
        return limit, univ

    def __mul__(left, right):
        cone = GSet.mul(left, right)
        return cone.apex

    def __add__(left, right):
        cone = GSet.add(left, right)
        return cone.apex

    def fixed_points(self, H):
        send_perms = self.send_perms
        fixed = set(range(self.rank))
        G = self.src
        tgt = self.tgt
        for g in H:
            idx = G.lookup[g]
            g = tgt[send_perms[idx]]
            fixed = fixed.intersection(g.fixed())
            if not fixed:
                break
        return fixed

    def signature(self, Hs=None):
        if Hs is None:
            if self._sig is not None:
                return self._sig
            Hs = self.src.subgroups()
            self._sig = tuple(len(self.fixed_points(H)) for H in Hs)
            return self._sig
        else:
            return tuple(len(self.fixed_points(H)) for H in Hs)

    def isomorphic(self, other):
        return self.signature() == other.signature()

    def get_identity(self):
        return Hom(self, self)



class Hom(object):
    """
        A morphism of GSet's.
    """

    def __init__(self, src, tgt, send_items=None):
        assert isinstance(src, GSet)
        assert isinstance(tgt, GSet)
        G = src.src
        assert G == tgt.src
        self.src = src
        self.tgt = tgt
        self.G = G
        if send_items is None:
            assert src is tgt
            send_items = list(range(src.tgt.rank))
        assert len(send_items) == src.tgt.rank
        self.rank = src.tgt.rank
        self.send_items = list(send_items)
        if debug:
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
        return "Hom(%s, %s, %s)"%(self.src, self.tgt, self.send_items)
    __repr__ = __str__

    def __eq__(self, other):
        assert self.G is other.G # too strict ?
        assert self.src == other.src
        assert self.tgt == other.tgt
        return self.send_items == other.send_items

    def __ne__(self, other):
        assert self.G is other.G # too strict ?
        assert self.src == other.src
        assert self.tgt == other.tgt
        return self.send_items != other.send_items

    def compose(self, other):
        # other o self
        assert self.tgt == other.src
        a = self.send_items
        b = other.send_items
        send_items = [b[i] for i in a]
        return Hom(self.src, other.tgt, send_items)

    def mul(f, g):
        assert isinstance(g, Hom)
        cone = GSet.mul(f.src, g.src)
        src = cone.apex
        cone = Cone(src, [cone[0].compose(f), cone[1].compose(g)])
        cone, univ = GSet.mul(f.tgt, g.tgt, cone)
        return univ
    __mul__ = mul



class Simplicial(object):
    """
        A simplicial object in some category
    """

    def __init__(self, X):
        I = X.get_identity()
        cone = Cone(X, [I, I])
        cone, diag = X.mul(X, cone)
        XX = cone.apex
        self.items = [X, XX]
        self.degenmaps = [[diag]] # diagonals
        self.facemaps = [[cone[1], cone[0]]] # projections
        if debug:
            self.check()

    def construct(self):
        items = self.items
        X = items[0]
        I = X.get_identity()
        X1 = items[-1]
        XX = X1*X
        facemaps = self.facemaps[-1]
        degenmaps = self.degenmaps[-1]
        dmaps = [f*I for f in degenmaps] + [I*degenmaps[-1]]
        fmaps = [f*I for f in facemaps] + [I*facemaps[-1]]
        for f in dmaps:
            assert f.src == X1
            assert f.tgt == XX
        for f in fmaps:
            assert f.src == XX
            assert f.tgt == X1
        self.items.append(XX)
        self.facemaps.append(fmaps)
        self.degenmaps.append(dmaps)
        if debug:
            self.check()

    def dump(self):
        items = self.items
        facemaps = self.facemaps
        degenmaps = self.degenmaps
        n = len(facemaps)
        for idx in range(n):
            print("face  [%d] --> [%d]" % (idx+1, idx))
            for d in facemaps[idx]:
                print('\t', d)
            print("degen [%d] --> [%d]" % (idx, idx+1))
            for s in degenmaps[idx]:
                print('\t', s)

    def check(self):
        "check that we satisfy the defining relations of a Simplicial object."
        items = self.items
        facemaps = self.facemaps
        degenmaps = self.degenmaps
        n = len(facemaps)
        assert n == len(degenmaps)
        assert n+1 == len(items)
        for idx in range(n):
            for d in facemaps[idx]:
                assert d.src == items[idx+1]
                assert d.tgt == items[idx]
            for s in degenmaps[idx]:
                assert s.src == items[idx]
                assert s.tgt == items[idx+1]
        idx = 1
        while idx < n:
            m = len(facemaps[idx])
            for i in range(m):
              for j in range(i+1, m):
                lhs = facemaps[idx][j].compose(facemaps[idx-1][i])
                rhs = facemaps[idx][i].compose(facemaps[idx-1][j-1])
                assert lhs==rhs
            m = len(degenmaps[idx-1])
            for i in range(m):
              for j in range(i, m):
                lhs = degenmaps[idx-1][j].compose(degenmaps[idx][i])
                rhs = degenmaps[idx-1][i].compose(degenmaps[idx][j+1])
                assert lhs==rhs
            for i, di in enumerate(facemaps[idx]):
              for j, sj in enumerate(degenmaps[idx]):
                lhs = sj.compose(di)
                if i<j:
                    rhs = facemaps[idx-1][i].compose(degenmaps[idx-1][j-1])
                elif i==j or i==j+1:
                    rhs = items[idx].get_identity()
                else: # i>j+1
                    rhs = facemaps[idx-1][i-1].compose(degenmaps[idx-1][j])
                assert lhs==rhs
            idx += 1

    def __getitem__(self, idx):
        items = self.items
        while idx >= len(items):
            self.construct()
        return items[idx]

    def get_bdy(self, idx):
        # facemaps[idx] go from idx+1 --> idx
        src = self[idx+1]
        tgt = self[idx]
        ds = self.facemaps[idx]
        A = zeros(tgt.rank, src.rank)
        sign = 1
        for d in ds:
            assert d.src == src
            assert d.tgt == tgt
            for i, j in enumerate(d.send_items):
                A[j, i] += sign
            sign = -sign
        A = array(A)
        return A

    def get_cobdy(self, idx): # ummm ???
        # degenmaps[idx] go from idx --> idx+1
        src = self[idx]
        tgt = self[idx+1]
        fs = self.degenmaps[idx]
        A = zeros(tgt.rank, src.rank)
        sign = 1
        for f in fs:
            assert f.src == src
            assert f.tgt == tgt
            for i, j in enumerate(f.send_items):
                A[j, i] += sign
            sign = -sign
        A = array(A)
        return A


class Orbiplex(Simplicial):
    "equivariant Simplicial object"
    def __init__(self, X):
        assert isinstance(X, GSet)
        self.s = Simplicial(X) 
        n = len(X.get_orbits())
        Y = X.G.trivial_gset(n)
        self.G = X.G
        self.items = [Y]
        self.degenmaps = [] # diagonals
        self.facemaps = [] # projections

    def construct(self):
        s = self.s
        items = self.items

        idx = len(items)
        Y0 = items[idx-1]

        X0 = s[idx-1]
        X1 = s[idx]

        n = len(X1.get_orbits())
        Y1 = self.G.trivial_gset(n)
        items.append(Y1)
        
        proj0, sect0 = X0.get_retract()
        proj1, sect1 = X1.get_retract()
        dmaps = []
        fmaps = []
        for f in s.facemaps[idx-1]: # idx --> idx-1
            send_items = [proj0[f.send_items[sect1[i]]] for i in range(Y1.rank)]
            f = Hom(Y1, Y0, send_items)
            fmaps.append(f)
        for d in s.degenmaps[idx-1]: # idx-1 --> idx
            send_items = [proj1[d.send_items[sect0[i]]] for i in range(Y0.rank)]
            d = Hom(Y0, Y1, send_items)
            dmaps.append(d)

        self.facemaps.append(fmaps)
        self.degenmaps.append(dmaps)
        if debug:
            self.check()


def GL(n=3, p=2):
    G = Algebraic.GL(n, p)
    v = numpy.array([0]*n, dtype=scalar)
    v[0] = 1
    v = Matrix(v, p)
    X = set([])
    for g in G:
        X.add(g*v)
    X = list(X)
    #for x in X:
    #    print(x)
    G = Group.from_action(G, X)
    return G


def Sp(nn, p=2):
    assert nn%2 == 0
    n = nn//2
    G = Algebraic.Sp(nn, p)
    v = [0]*nn
    v[0] = 1
    v = Matrix(v, p)
    items = set([])
    for g in G:
        gv = g*v
        assert type(gv) is Matrix
        if gv not in items:
            items.add(gv)
    items = list(items)
    items.sort(key = str)
    lookup = {}
    basis = []
    for i, v in enumerate(items):
        lookup[v] = i
        if str(v).count('1') == 1:
            basis.append(v)
    assert len(basis) == nn
    B = G.get_borel()
    G = Group.from_action(G, items)
    B = Group.from_action(B, items)
    print("G.rank =", G.rank)
    print("|G| =", len(G))
    print("|B| =", len(B))
    X = G.action_subgroup(B)
    print("|G/B| =", X.rank)
    assert len(X.get_orbits()) == 1
    for m in range(1, n+1):
        print()
        print("(%d %d)"%(n, m))
        gen = basis[:m]
        space = []
        for idxs in cross([(0,1)]*m):
            v = reduce(add, [idx*v for (idx,v) in zip(idxs,basis)])
            space.append(v)
        assert len(space) == 2**m
        P = []
        for g in G:
            found = []
            for v in gen:
                gv = items[g[lookup[v]]] # argh, this sucks
                #print(v, "-->", gv, gv in space)
                if gv not in space:
                    break
                #found.append("%s-->%s"%(v, gv))
            else:
                #print(found)
                P.append(g)
        print("|P| =", len(P))
        continue
        P = Group(P)
        Y = G.action_subgroup(P)
        print("|G/P| =", Y.rank)
        assert len(Y.get_orbits()) == 1
        XY = X*Y
        orbits = XY.get_orbits()
        print("B-P double cosets:", len(orbits))
    return G


def test_bruhat():
    print("test_bruhat")
    n = argv.get("n", 2)
    nn = 2*n
    G = Sp(nn)


def test_gset():

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

    s = Simplicial(X)
    s.construct()
    s.construct()
    s.dump()

    G = Group.cyclic(3)
    X = G.i
    I = X.get_identity()
    cone = Cone(X, [I, I], contra=True)
    cone, univ = add(X, X, cone)
    X_X = cone.apex
    for hom in X_X.get_atoms():
        assert hom.src.isomorphic(X)
        assert hom.compose(univ) == I

    G = Group.cyclic(5)
    H = mul(G.i, G.i).apex.tgt
    assert len(H.get_orbits()) == 5

    G = Group.alternating(5)
    H = mul(G.i, G.i).apex
    H = H.tgt
    Gs = H.get_atoms()
    assert len(Gs) == 2
    assert Gs[0].rank == 5
    assert Gs[1].rank == 20

    G = Group.alternating(5)
    H = mul(G.i, G.i).apex
    H = mul(G.i, H).apex
    items = [nat.src for nat in H.get_atoms()]
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

    G = Group.symmetric(3)
    assert list(G.get_sequence(5)) == [1, 2, 5, 14, 41, 122]

    G = GL(3, 2)
    assert len(G) == 168

    #for size in G.get_sequence():
    #    print(size, end=" ", flush=True)
    #print()

    print("OK")


def test_fano():
    G = GL(3, 2)
    assert len(G) == 168

    #G = Group.symmetric(3)

    Hs = list(G.subgroups(verbose=True))
    assert len(Hs) == 179

    Hs = G.conjugacy_subgroups(verbose=True)
    assert len(Hs) == 15

    print("Hs:", len(Hs))
    for H in Hs:
        print(H)

    return

    for g in G:
        if g.order() == 7:
            break
    print(g)

    H = Group.generate([g])
    assert len(H) == 7

    X = G.action_subgroup(H)
    print(X)

    XX = X*X
    print(XX)

    orbits = XX.get_orbits()
    print([len(o) for o in orbits])


def test_subgroups_only():
    G = Group.symmetric(5)
    print("|G| =", len(G))
    Hs = list(G.subgroups(verbose=True))
    print(len(Hs))


def test_subgroups():

    G = Group.cyclic(12)
    #print(G.cyclic_subgroups())
    assert len(G.cyclic_subgroups()) == 6 # 1, 2, 3, 4, 6, 12

    G = Group.symmetric(4)
    assert not G.is_simple()

    orders = []
    for H in G.subgroups():
        orders.append(len(H))
    assert len(orders) == 1+1+3+4+3+1+3+4+6+3+1
    assert orders.count(24) == 1
    assert orders.count(12) == 1
    assert orders.count(8) == 3
    assert orders.count(6) == 4
    assert orders.count(4) == 3+1+3
    assert orders.count(3) == 4
    assert orders.count(2) == 6+3
    assert orders.count(1) == 1

    G = Group.cyclic(6)
    for H in G.subgroups():
        assert G.is_normal(H)
    assert not G.is_simple()

    print("OK")



def test_left_right():

    n = argv.get("n", 3)

    G = Group.symmetric(n)

    print(G)

    act = G.regular_action()
    print(act)

    act.do_check()


    lookup = G.lookup
    perms = []
    for idx, p in enumerate(G):
        perm = []
        for jdx, q in enumerate(G):
            #r = p*q # left action on G
            r = q*p.inverse()
            kdx = lookup[r]
            perm.append(kdx)
        perm = Perm(perm)
        perms.append(perm)
    tgt = Group(perms)
    #send_perms = list(range(G.n))
    #assert send_perms == [tgt.lookup[perm] for perm in perms]
    X = GSet(G, tgt)
    X.do_check()


def test_homology():
    #X = Group.trivial(2).i
    #G = Group.dihedral(4)
    n = argv.get("n", 2)
    if argv.trivial:
        G = Group.trivial(n)
    elif argv.symmetric:
        G = Group.symmetric(n)
    elif argv.alternating:
        G = Group.alternating(n)
    elif argv.dihedral:
        G = Group.dihedral(n)
    elif argv.cyclic:
        G = Group.cyclic(n)
    else:
        return

    if argv.regular:
        X = G.regular_action()
    elif argv.double:
        X = G.i
        X = X + X
    else:
        X = G.i

    e = Orbiplex(X)
    e.construct()

#    for i in range(3):
#        proj, sect = e[i].get_retract_linear()
#        print(proj)
#        print(sect)
#        print(shortstr(dot(proj, sect))) # identity matrix
#        print()
#
#    return

#    for f in e.facemaps[0]:
#        print("proj:", f)
#    for f in e.degenmaps[0]:
#        print("diag:", f)

    degree = argv.get("degree", 4)
    bdys = [e.get_bdy(i) for i in range(degree)]
    #for i in range(degree):
    #    print(bdys[i])

    for i in range(degree-1):
        d0 = bdys[i] # 1 --> 0
        d1 = bdys[i+1] # 2 --> 1
        assert numpy.abs(dot(d0, d1)).sum() == 0
        kern = nullity(d0)
        im = rank(d1)
        print("betti %d: %d-%d = %d"%(i, kern, im, kern-im))

#
#    bdys = [e.get_cobdy(i) for i in range(5)]
#    for i in range(3):
#        d0 = bdys[i] # 1 <-- 0
#        #print(d0)
#        d1 = bdys[i+1] # 2 <-- 1
#        assert numpy.abs(dot(d1, d0)).sum() == 0
#        print("cobetti %d ="%i, rank(d0) - nullity(d1)) # == 0
    

def test_isomorphic():

    Gs = [
        Group.dihedral(6),
        Group.alternating(4),
        Group.symmetric(3) * Group.cyclic(2),
        Group.cyclic(2) * Group.symmetric(3),
        Group.cyclic(2) * Group.cyclic(3) * Group.cyclic(2),
        Group.cyclic(3) * Group.cyclic(2) * Group.cyclic(2),
    ]

    N = len(Gs)
    for i in range(N):
      for j in range(i+1, N):
        G,H = Gs[i], Gs[j]
        iso = G.isomorphic(H)
        #assert (iso is not None) == ((i,j) in [(0,2), (0,3), (2,3), (4,5)]) # XXX BROKEN
        assert (iso is not None) == (H.isomorphic(G) is not None)



def main():
    G = Group.dihedral(4)
    G = Group.symmetric(4)
    #G = Group.alternating(5)

    sigs = []
    for H in G.subgroups():
        X = G.action_subgroup(H)
        assert len(X.get_orbits()) == 1
        sigs.append(X.signature())
    sigs = set(sigs)
    sigs = list(sigs)
    sigs.sort()
    names = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    X = G.i
    #print(list(X.get_sequence()))

    mul = lambda a,b: GSet.mul(a, b).apex
    XX = mul(X, X)
    XXX = mul(XX, X)
    XXXX = mul(XXX, X)

    for A in [X, XX, XXX, XXXX]:
        Bs = [hom.src for hom in A.get_atoms()]
        items = [names[sigs.index(B.signature())] for B in Bs]
        print(len(items), "+".join(items), end=" --- ")
        for B in Bs:
            print(B.rank, end=" ")
        print()


def test_orbits():
    meth = argv.get("meth", "symmetric")
    fn = getattr(Group, meth)
    n1 = 7
    if meth != "symmetric":
        n1 = 6
    for n in range(2, n1):
        G = fn(n)
        X = G.i
        H = [g for g in G if g[0]==0]
        print("n=%s |G|=%d |H|=%d |X|=%d:"%(n, len(G), len(H), X.rank))
        print("\t\t\t\t", end="")
        Y = None
        for k in range(1, n+1):
            Y = X*Y if Y is not None else X
            print(len(Y.get_orbits()), end=" ", flush=True)
            if Y.rank > 40000:
                break
        print()
    print("done")


def test():
    test_isomorphic()
    test_gset()
    test_subgroups()
    test_homology()
    test_orbits()

    

if __name__ == "__main__":
    from time import time
    start_time = time()

    fn = argv.next() or "test"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds.\n"%(time() - start_time))


