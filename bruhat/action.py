#!/usr/bin/env python3

"""
Group actions.

see also: gset.py

"""

import sys
import string
from random import randint, shuffle
from functools import reduce, cache

from bruhat.util import factorial, all_subsets, write, uniqtuples, cross
from bruhat.equ import Equ, quotient_rep
from bruhat.argv import argv
from bruhat import isomorph
from bruhat.smap import SMap, tabulate

long = int


def mulclose_fast(gen, verbose=False, maxsize=None):
    els = set(gen)
    bdy = list(els)
    changed = True 
    while bdy:
        if verbose:
            print(len(els), end=" ", flush=True)
        _bdy = []
        for A in gen:
            for B in bdy:
                C = A*B  
                if C not in els: 
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        return els
        bdy = _bdy
    if verbose:
        print()
    return els 


mulclose = mulclose_fast

def mulclose_names(gen, names, verbose=False, maxsize=None):
    bdy = list(set(gen))
    assert len(names) == len(gen)
    names = dict((gen[i], (names[i],)) for i in range(len(gen)))
    changed = True 
    while bdy:
        _bdy = []
        for A in gen:
            for B in bdy:
                C = A*B  
                if C not in names: 
                    #els.add(C)
                    names[C] = names[A] + names[B]
                    _bdy.append(C)
                    if maxsize and len(names)>=maxsize:
                        return list(names)
        bdy = _bdy
    return names 



def mulclose_hom(gen1, gen2, verbose=False, maxsize=None):
    "build a group hom from generators: gen1 -> gen2"
    hom = {}
    assert len(gen1) == len(gen2)
    for i in range(len(gen1)):
        hom[gen1[i]] = gen2[i]
    bdy = list(gen1)
    changed = True 
    while bdy:
        #if verbose:
        #    print "mulclose:", len(hom)
        _bdy = []
        for A in gen1:
            for B in bdy:
                C1 = A*B  
                if C1 not in hom: 
                    hom[C1] = hom[A] * hom[B]
                    _bdy.append(C1)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
        bdy = _bdy
    return hom 


def identity(items):
    return dict((i, i) for i in items)



class Perm(object):

    """
    A permutation of a list of items.
    """
    def __init__(self, perm, items, word=''):
        #print "Perm.__init__", perm, items
        #if isinstance(perm, list):
        #    perm = tuple(perm)
        if perm and isinstance(perm, (list, tuple)) and isinstance(perm[0], (int, long)):
            perm = list(items[i] for i in perm)
        #print "\t", perm
        if not isinstance(perm, dict):
            perm = dict((perm[i], items[i]) for i in range(len(perm)))
        #print "\t", perm
        self.perm = perm # map item -> item
        #print "\t", perm
        set_items = set(items)
        self.set_items = set_items
        self.items = list(items)
        assert len(perm) == len(items), (perm, items)
        self.n = len(perm)
        self._str = None
        for key, value in perm.items():
            assert key in set_items, repr(key)
            assert value in set_items, repr(value)
        self.word = word
        self._str_cache = None
        self._hash_cache = None

    @classmethod
    def identity(cls, items, *args, **kw):
        n = len(items)
        perm = dict((item, item) for item in items)
        return Perm(perm, items, *args, **kw)

    def is_identity(self):
        for k, v in self.perm.items():
            if k != v:
                return False
        return True

    def rename(self, send_items, items):
        perm = {}
        assert len(items) == len(self.items)
        for item in self.items:
            jtem = send_items[item]
            assert jtem in items
            perm[jtem] = send_items[self.perm[item]]
        return Perm(perm, items)

    @classmethod
    def fromcycles(cls, cycles, items, *args, **kw):
        perm = {}
        for cycle in cycles:
            m = len(cycle)
            for i in range(m):
                perm[cycle[i]] = cycle[(i+1)%m]
        return Perm(perm, items, *args, **kw)

    def order(self):
        i = 1
        g = self
        while not g.is_identity():
            g = self*g
            i += 1
            #assert i <= len(self.items)+1 # how big can this get ??
        return i

    def restrict(self, items, *args, **kw):
        perm = dict((i, self.perm[i]) for i in items)
        return Perm(perm, items, *args, **kw)

    def fixes(self, items):
        items = set(items)
        for item in items:
            item = self(item)
            if item not in items:
                return False
        return True

    def intersection(self, other):
        assert self.items == other.items
        items = []
        for i in self.items:
            if self.perm[i] == other.perm[i]:
                items.append(i)
        return items

    def _X_str__(self):
        #return str(dict((i, self.perm[i]) for i in range(self.n)))
        #return str(dict((i, self.perm[i]) for i in range(self.n)))
        if self._str:
            return self._str
        perm = self.perm
        keys = perm.keys()
        keys.sort()
        items = ["%s:%s"%(key, perm[key]) for key in keys]
        s = "{%s}"%(', '.join(items))
        self._str = s
        return s

    def cycles(self):
        remain = set(self.set_items)
        cycles = []
        while remain:
            item = iter(remain).__next__()
            orbit = [item]
            #print(type(self), type(item), "__mul__")
            item1 = self(item)
            while item1 != item:
                orbit.append(item1)
                item1 = self(item1)
                assert len(orbit) <= len(self.items)
            assert orbit
            cycles.append(orbit)
            n = len(remain)
            for item in orbit:
                remain.remove(item)
            assert len(remain) < n
        return cycles

    def sign(self):
        s = 1
        for orbit in self.cycles():
            if len(orbit)%2 == 0:
                s *= -1
        return s

    def cycle_str(self):
        remain = set(self.set_items)
        s = []
#        print "__str__", self.perm, self.items
        while remain:
            item = iter(remain).__next__()
#            print "item:", item
            orbit = [item]
            item1 = self*item
#            print "item1:", item1
            while item1 != item:
                orbit.append(item1)
                item1 = self*item1
                assert len(orbit) <= len(self.items)
            s.append("(%s)"%(' '.join(str(item) for item in orbit)))
            assert orbit
            n = len(remain)
            for item in orbit:
                remain.remove(item)
            assert len(remain) < n
        #return "Perm(%s)"%(''.join(s))
#        print "__str__: end"
        return ''.join(s)
#    __repr__ = __str__

    def get_idxs(self):
        lookup = dict((v, k) for (k, v) in enumerate(self.items))
        idxs = [lookup[self.perm[i]] for i in self.items]
        return idxs

    def slowstr(self): # HOTSPOT
        if self._str_cache:
            return self._str_cache
        perm = self.perm
        items = self.items
        s = []
        for i, item in enumerate(items):
            j = items.index(perm[item]) # <--- oof
            s.append("%d:%d"%(i, j))
        s = "{%s}"%(', '.join(s))
        self._str_cache = s
        return s

    def str(self): # HOTSPOT
        if self._str_cache:
            return self._str_cache
        perm = self.perm
        items = self.items
        lookup = dict((v,k) for (k,v) in enumerate(items))
        s = []
        for i, item in enumerate(items):
            j = lookup[perm[item]]
            s.append("%d:%d"%(i, j))
        s = "{%s}"%(', '.join(s))
        self._str_cache = s
        return s

    def __str__(self):
        # this is the *index* action, not the real permutation action
        return "Perm(%s)"%self.str()
    __repr__ = __str__

    def __hash__(self):
        #return hash(str(self))
        if self._hash_cache is None:
            self._hash_cache = hash(str(self))
        return self._hash_cache

    """
    _hash = None
    def __hash__(self):
        if self._hash is None:
            self._hash = hash(str(self))
        return self._hash
    """

    def leftmul(self, action):
        perms = [self*perm for perm in action]
        action = Group(perms, action.items)
        return action

    def rightmul(self, action):
        perms = [perm*self for perm in action]
        action = Group(perms, action.items)
        return action

    def __mul__(self, other):
        # break this into three cases: 
        if isinstance(other, Group):
            perms = [self*perm for perm in other]
            action = Group(perms, self.items)
            return action # <----- return

        if not isinstance(other, Perm):
            item = self.perm[other]
            return item # <---- return

        assert self.items == other.items
        perm = {}
        #for item in self.items:
            #perm[item] = other.perm[self.perm[item]]
        for item in other.items:
            perm[item] = self.perm[other.perm[item]]
        return Perm(perm, self.items, self.word+other.word)

    # More strict than __mul__:
    def __call__(self, item):
        return self.perm[item]
    __getitem__ = __call__

    def __pow__(self, n):
        assert int(n)==n
        if n==0:
            return Perm.identity(self.items)
        if n<0:
            self = self.__invert__()
            n = -n
        g = self
        for i in range(n-1):
            g = self*g
        return g

    def __invert__(self):
        perm = {}
        for item in self.items:
            perm[self.perm[item]] = item
        return Perm(perm, self.items)
    inverse = __invert__

    def __eq__(self, other):
        assert self.items == other.items
        return self.perm == other.perm

    def __ne__(self, other):
        assert self.items == other.items
        return self.perm != other.perm

    def fixed(self):
        items = []
        for item in self.items:
            if self.perm[item] == item:
                items.append(item)
        return items

    def orbits(self):
        remain = set(self.items)
        orbits = []
        while remain:
            item = iter(remain).__next__()
            orbit = [item]
            item0 = item
            while 1:
                item = self.perm[item]
                if item == item0:
                    break
                orbit.append(item)
                assert len(orbit) <= len(remain)
            remain = remain.difference(orbit)
            orbits.append(orbit)
        return orbits

    def conjugacy_cls(self):
        "this is a partition of len(items)"
        orbits = self.orbits()
        sizes = [len(orbit) for orbit in orbits]
        sizes.sort() # uniq
        return tuple(sizes)

    def preserves_partition(self, part):
        for items in part:
            for item in items:
                jtem = self.perm[item]
                if jtem not in items:
                    return False
        return True


#class Species(object):
#    def __call__(self, group, items):
#        pass

#class PointedSpecies(Species):

class Item(object):
    def __init__(self, item, name=None):
        self.item = item
        if name is None:
            name = str(item)
        self.name = name # a canonical name

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return self.name
    def __repr__(self):
        return "I(%s)"%(self.name)
    def __eq__(self, other):
        return self.name == other.name
    def __ne__(self, other):
        return self.name != other.name


class TupleItem(Item):
    def __init__(self, items):
        Item.__init__(self, items)


class SetItem(Item):
    "a hashable unordered set"
    def __init__(self, items):
        items = list(items)
        #for i in range(len(items)):
        #    item = items[i]
        #    if isinstance(item, Item):
        #        pass
        #    elif isinstance(items[i], (int, long, str)):
        #        items[i] = Item(item)
        #    else:
        #        assert 0, repr(item)
        items.sort(key = lambda item : str(item))
        #items = tuple(items)
        Item.__init__(self, items)

    def __iter__(self):
        return iter(self.item)

    def __len__(self):
        return len(self.item)


def disjoint_union(items, _items):
    items = [(0, item) for item in items]
    _items = [(1, item) for item in _items]
    return items + _items
    

def all_functions(source, target):
    m, n = len(source), len(target)
    source = list(source)
    target = list(target)
    assert n**m < 1e8, "%d too big"%(n**m,)
    if m==0:
        yield {}
    elif n==0:
        return # no functions here
    elif m==1:
        for i in range(n):
            yield {source[0]: target[i]}
    else:
        for func in all_functions(source[1:], target):
            for i in range(n):
                _func = dict(func)
                _func[source[0]] = target[i]
                yield _func

assert len(list(all_functions('ab', 'abc'))) == 3**2
assert len(list(all_functions('abc', 'a'))) == 1
assert len(list(all_functions('a', 'abc'))) == 3


def __choose(items, k):
    "choose k elements"
    n = len(items)
    assert 0<=k<=n
    if k==0:
        yield [], items, [] # left, right, chosen
    elif k==n:
        yield items, [], [] # left, right, chosen
    elif k==1:
        for i in range(n):
            yield items[:i], items[i+1:], [items[i]]
    else:
        for left, right, chosen in __choose(items, k-1):
            n = len(right)
            for i in range(n):
                yield left + right[:i], right[i+1:], chosen+[right[i]]

def _choose(items, *ks):
    "yield a tuple "
    if len(ks)==0:
        yield items,
    elif len(ks)==1:
        k = ks[0]
        for left, right, chosen in __choose(items, k):
            yield items, chosen
    else:
        k = ks[0]
        for flag in _choose(items, *ks[:-1]):
            chosen = flag[-1]
            for chosen, _chosen in _choose(chosen, ks[-1]):
                yield flag + (_chosen,)


def choose(items, *ks):
    "choose k elements"
    items = list(items)
    _items = []
    #for left, right, chosen in _choose(items, k):
    #    _items.append((SetItem(left+right), SetItem(chosen)))
    for flag in _choose(items, *ks):
        flag = tuple(SetItem(item) for item in flag)
        _items.append(flag)
    return _items

items4 = list('abcd')
assert len(choose(items4, 0)) == 1
assert len(choose(items4, 1)) == 4
assert len(choose(items4, 2)) == 4*3//2
assert len(choose(items4, 3)) == 4
assert len(choose(items4, 4)) == 1

class Group(object):
    """
    A collection of Perm's.
    """

    def __init__(self, perms, items, check=False):
        perms = list(perms)
        self.perms = perms # ordered  ( see Group .str and .__hash__ )
        self.items = list(items)
        self.set_items = set(items) # unordered
        self.set_perms = set(perms) # unordered
        for perm in perms:
            assert isinstance(perm, Perm), type(perm)
            assert perm.items == self.items, (perm.items, items)
        self._str = None # cache
        self.conjugates = []

    def str(self):
        if not self._str:
            ss = [perm.str() for perm in self.perms]
            ss.sort() # <-- ordered 
            self._str = ''.join(ss)
        return self._str

    def __eq__(self, other): # HOTSPOT 
        if len(self.perms) != len(other.perms):
            return False
        return (self.set_items == other.set_items and self.set_perms == other.set_perms)

    def __ne__(self, other):
        if len(self.perms) != len(other.perms):
            return True
        return (self.set_items != other.set_items or self.set_perms != other.set_perms)

    def __hash__(self):
        return hash(self.str())

    def __contains__(self, g):
        assert g.items == self.items
        return g in self.set_perms

#    def __lt__(self, other):
#        return id(self) < id(other)

    @classmethod
    def generate(cls, gen, *args, **kw):
        items = kw.get("items")
        if items is not None:
            del kw["items"]
        elif gen:
            items = gen[0].items
        else:
            items = []
        if not gen:
            gen = [Perm.identity(items)]
        perms = list(mulclose(gen, *args))
        G = cls(perms, items, **kw)
        G.gen = gen
        G.gens = gen # deprecate this?
        return G

    def is_cyclic(self):
        n = len(self)
        for g in self:
            if g.order() == n:
                return True
        return False

    def cgy_cls(self):
        "conjugacy classes of elements of G"
        found = set() # map Perm to it's conjugacy class
        itemss = []
        for g in self:
            if g in found:
                continue
            items = set([g])
            itemss.append(items)
            found.add(g)
            for h in self:
                k = h*g*(~h)
                items.add(k)
                found.add(k)
        itemss.sort(key = lambda items : iter(items).__next__().order())
        return itemss

    @property
    def identity(self):
        p = Perm.identity(self.items)
        return p

    @classmethod
    def trivial(cls, items_or_n=1, check=False):
        "the trivial action on items fixes every item"
        if type(items_or_n) in (int, long):
            items = range(items_or_n)
        else:
            items = list(items_or_n)
        perm = Perm.identity(items)
        G = Group([perm], items, check=check)
        return G

    @classmethod
    def symmetric(cls, items_or_n, check=False):
        if type(items_or_n) in (int, long):
            items = range(items_or_n)
        else:
            items = list(items_or_n)
        perms = []
        n = len(items)
        for i in range(n-1):
            perm = dict((item, item) for item in items)
            perm[items[i]] = items[i+1]
            perm[items[i+1]] = items[i]
            perms.append(perm)
        if not perms:
            perms.append({items[0]:items[0]}) # trivial perm
        perms = [Perm(perm, items) for perm in perms]
        G = Group.generate(perms, check=check)
        assert len(G) == factorial(n)
        return G

    @classmethod
    def alternating(cls, items_or_n, check=False):
        if type(items_or_n) in (int, long):
            items = range(items_or_n)
        else:
            items = list(items_or_n)
        n = len(items)
        gen = []
        for i in range(n-2):
            perm = dict((j, j) for j in items)
            perm[items[i]] = items[i+1]
            perm[items[i+1]] = items[i+2]
            perm[items[i+2]] = items[i]
            perm = Perm(perm, items)
            gen.append(perm)
        perms = mulclose(gen)
        G = Group(perms, items)
        assert len(G)==factorial(n)//2
        return G

    @classmethod
    def cyclic(cls, items_or_n, check=False):
        if type(items_or_n) in (int, long):
            items = range(items_or_n)
        else:
            items = list(items_or_n)
        perms = []
        n = len(items)
        perms = [dict((items[i], items[(i+k)%n]) for i in range(n))
            for k in range(n)]
        assert len(perms) == n
        perms = [Perm(perm, items) for perm in perms]
        G = Group(perms, items, check=check)
        return G

    @classmethod
    def dihedral(cls, items_or_n, check=False):
        if type(items_or_n) in (int, long):
            items = range(items_or_n)
        else:
            items = list(items_or_n)
        perms = []
        n = len(items)
        perms = [
            dict((items[i], items[(i+1)%n]) for i in range(n)),
            dict((items[i], items[(-i)%n]) for i in range(n))]
        perms = [Perm(perm, items) for perm in perms]
        G = Group.generate(perms, check=check)
        assert len(G) == 2*n
        return G

    def __repr__(self):
        return "%s(%s, %s)"%(self.__class__.__name__, self.perms, self.items)

    def __str__(self):
        return "%s(%d, %d)"%(self.__class__.__name__, len(self.perms), len(self.items))

    def __len__(self):
        return len(self.perms)

    def __getitem__(self, idx):
        return self.perms[idx]
    
    def __mul__(self, other):
        if isinstance(other, Group):
            assert self.items == other.items
            perms = set(g*h for g in self for h in other)
            return Group(perms, self.items)
        elif isinstance(other, Perm):
            assert self.items == other.items
            perms = set(g*other for g in self)
            return Group(perms, self.items)
        raise TypeError

    def __add__(self, other):
        items = disjoint_union(self.items, other.items)
        perms = []
        for perm in self:
            _perm = {}
            for item in self.items:
                _perm[0, item] = 0, perm[item]
            for item in other.items:
                _perm[1, item] = 1, item # identity
            _perm = Perm(_perm, items)
            perms.append(_perm)
        for perm in other:
            _perm = {}
            for item in self.items:
                _perm[0, item] = 0, item # identity
            for item in other.items:
                _perm[1, item] = 1, perm[item]
            _perm = Perm(_perm, items)
            perms.append(_perm)
        perms = list(mulclose(perms))
        return Group(perms, items)

    def direct_product(self, other):
        items = [(i, j) for i in self.items for j in other.items]
        perms = []
        for g in self:
          for h in other:
            perm = {}
            for i, j in items:
                perm[i, j] = g[i], h[j]
            perm = Perm(perm, items)
            perms.append(perm)
        group = Group(perms, items)
        return group

    def fixed(self):
        fixed = {}
        for el in self.perms:
            fixed[el] = el.fixed()
        return fixed

    def stabilizer(self, *items):
        "subgroup that fixes every point in items"
        perms = []
        for g in self.perms:
            if len([item for item in items if g[item]==item])==len(items):
                perms.append(g)
        return Group(perms, self.items)

    def invariant(self, *items):
        "subgroup that maps set of items to itself"
        perms = []
        items = set(items)
        for g in self.perms:
            items1 = set(g[item] for item in items)
            if items1.issubset(items):
                perms.append(g)
        return Group(perms, self.items)

    def preserve_partition(self, part):
        H = [g for g in self if g.preserves_partition(part)]
        return Group(H, self.items)

    def orbit(self, item):
        items = set(g(item) for g in self.perms)
        return items

    def orbits(self):
        #print "orbits"
        #print self.perms
        #print self.items
        remain = set(self.set_items)
        orbits = []
        while remain:
            #print "remain:", remain
            item = iter(remain).__next__()
            orbit = set(g(item) for g in self.perms)
            #print "orbit:", orbit
            for item in orbit:
                remain.remove(item)
            orbits.append(orbit)
        return orbits

    def restrict(self, items):
        G = Group([perm.restrict(items) for perm in self.perms], items)
        return G

    def components(self):
        orbits = self.orbits()
        groups = [self.restrict(orbit) for orbit in orbits]
        return groups

    def shape_item(self): # HOTSPOT
        shape_item = []
        for item in self.items:
            shape = []
            for g in self:
                gitem = item
                n = 1
                while 1:
                    gitem = g(gitem)
                    if gitem == item:
                        break
                    n += 1
                shape.append(n)
            shape_item.append(tuple(shape))
        return shape_item

    def check(group):
        #print "check"
        orbits = group.orbits()
        assert sum(len(orbit) for orbit in orbits) == len(group.items)

        # orbit stabilizer theorem
        n = sum(len(group.stabilizer(item)) for item in group.items)
        assert n == len(group) * len(orbits)

        # Cauchy-Frobenius lemma
        assert n == sum(len(g.fixed()) for g in group)

    def subgroups_slow(self): # XXX SLOW
        "All subgroups, _acting on the same items."
        subs = set()
        I = self.identity
        items = self.items
        group = Group([I], self.items)
        subs = set([group, self])
        bdy = set()
        if len(self)>1:
            bdy.add(group)
        set_perms = self.set_perms
        while bdy:
            _bdy = set()
            for group in bdy:
                assert len(group)<len(self)
                # XXX do something clever with cosets...
                for perm in set_perms.difference(group.set_perms):
                    perms = mulclose(group.perms + [perm])
                    _group = Group(perms, items)
                    if _group not in subs:
                        _bdy.add(_group)
                        subs.add(_group)
            bdy = _bdy

        return subs

    def cyclic_subgroups(self, verbose=False):
        # find all cyclic subgroups
        I = self.identity
        trivial = Group([I], self.items)
        cyclic = set()
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
            group = Group(perms, self.items)
            assert len(group)>1
            cyclic.add(group)
        return cyclic

    @cache
    def subgroups(self, verbose=False):
        I = self.identity
        items = self.items
        trivial = Group([I], items)
        cyclic = self.cyclic_subgroups()
        #if verbose:
        #    print "Group.subgroups: cyclic:", len(cyclic)
        n = len(self) # order
        subs = set(cyclic)
        subs.add(trivial)
        subs.add(self)
        bdy = set(G for G in cyclic if len(G)<n)
        while bdy:
            _bdy = set()
            # consider each group in bdy
            for G in bdy:
                # enlarge by a cyclic subgroup
                for H in cyclic:
                    perms = set(G.perms+H.perms)
                    k = len(perms)
                    if k==n or k==len(G) or k==len(H):
                        continue
                    #perms = mulclose(perms)
                    perms = mulclose_fast(perms)
                    K = Group(perms, items)
                    if K not in subs:
                        _bdy.add(K)
                        subs.add(K)
                        if verbose:
                            write('.')
                    #else:
                    #    write('/')
            bdy = _bdy
            #if verbose:
            #    print "subs:", len(subs)
            #    print "bdy:", len(bdy)
        return subs

    def left_cosets(self, H=None):
        cosets = set()
        if H is not None:
            Hs = [H]
        else:
            Hs = self.subgroups()
        lookup = dict((g, g) for g in self) # remember canonical word 
        for action in Hs:
            for g in self:
                coset = Coset([lookup[g*h] for h in action], self.items)
                cosets.add(coset)
        return list(cosets)

    def right_cosets(self, H=None):
        cosets = set()
        if H is not None:
            Hs = [H]
        else:
            Hs = self.subgroups()
        lookup = dict((g, g) for g in self) # remember canonical word 
        for action in Hs:
            for g in self:
                coset = Coset([lookup[h*g] for h in action], self.items)
                cosets.add(coset)
        return list(cosets)

    def left_action(self, items, basepoint=None):
        send_perms = {}
        perms = []
        lookup = dict((item, item) for item in items)
        for g in self:
            perm = {}
            for item in items:
                perm[item] = lookup[g*item]
            perm = Perm(perm, items)
            perms.append(perm)
            send_perms[g] = perm
        f = quotient_rep(perms)
        send_perms = dict((k, f[v]) for (k, v) in send_perms.items())
        perms = list(f.values())
        hom = Action(self, send_perms, items, basepoint)
        return hom

    def action_subgroup(self, H):
        assert self.items == H.items
        assert self.is_subgroup(H)
        cosets = self.left_cosets(H)
        hom = self.left_action(cosets, H)
        return hom

    def tautological_action(self):
        send_perms = {g:g for g in self}
        action = Action(self, send_perms, self.items)
        return action

    def cayley_action(self, H=None):
        "the left Cayley action of a subgroup on self"
        if H is None:
            H = self
        items = self.perms
        send_perms = {}
        for g in H:
            perm = Perm({h : g*h for h in items}, items)
            send_perms[g] = perm
        action = Action(H, send_perms, items)
        return action

    def is_subgroup(self, H):
        assert H.items == self.items
        for g in H.perms:
            if not g in self.perms:
                return False
        return True

    def is_abelian(self):
        pairs = [(g, h) for g in self for h in self]
        shuffle(pairs)
        for g, h in pairs:
            if g*h != h*g:
                return False
        return True

    def regular_rep(self):
        items = range(len(self))
        lookup = dict((v,k) for (k,v) in enumerate(self.perms))
        perms = []
        for perm in self:
            _perm = {}
            for i in items:
                j = lookup[perm*self[i]]
                _perm[i] = j
            perms.append(Perm(_perm, items))
        return Group(perms, items)

    @classmethod
    def product(cls, H, J):
        perms = list(set(H.perms+J.perms))
        return cls.generate(perms)
    
    def conjugacy_subgroups(G, Hs=None):
    
        # Find all conjugacy classes of subgroups
    
        if Hs is None:
            Hs = G.subgroups()
        #print "subgroups:", len(Hs)
        #for H in Hs:
        #  for K in Hs:
        #    print (int(H.is_subgroup(K)) or "."),
        #  print
    
        equs = dict((H1, Equ(H1)) for H1 in Hs)
        for H1 in Hs:
            for g in G:
                if g in H1:
                    continue
                H2 = g * H1 * (~g) # conjugate
                if H2 == H1:
                    continue
                else:
                    #print len(H1), "~", len(H2)
                    if H2 not in equs:
                        equs[H2] = Equ(H2)
                    equs[H1].merge(equs[H2])
    
        # get equivalance classes
        equs = list(set(equ.top for equ in equs.values()))
        equs.sort(key = lambda equ : (-len(equ.items[0]), equ.items[0].str()))
        for equ in equs:
            #print "equ:", [len(H) for H in equ.items]
            for H in equ.items:
                H.conjugates = list(equ.items)
        #print "total:", len(equs)
        Hs = [equ.items[0] for equ in equs] # pick unique (up to conjugation)
        #Hs.sort(key = lambda H : (-len(H), H.str()))
        return Hs

    # | | | | | |                                       | | | | 
    # v v v v v v  probably should delete these methods v v v v 

    def choice(group, *ks):
        "choose k elements"
    
        items = group.items
        _items = choose(items, *ks)
        _group = []
        for g in group:
            perm = g.perm
            _perm = {}
            for flag in _items:
                _flag = tuple(SetItem(tuple(perm[item] for item in __items)) for __items in flag)
                _perm[flag] = _flag
            _g = Perm(_perm, _items)
            _group.append(_g)
        return Group(_group, _items)

    def is_hom(self, items, func):
        #print "is_hom", items, func
        # If we use a set here this could make a "subgroup" of self
        #group = set()
        group = []
        for g in self:
            perm = {} # try to send g to this guy
            for i in self.items:
                gi = g[i]
                f_i = func[i]
                f_gi = func[gi]
                j = perm.get(f_i)
                if j is None:
                    perm[f_i] = f_gi
                elif j != f_gi:
                    return # fail
            #group.add(Perm(perm, items))
            group.append(Perm(perm, items))
        group = Group(group, items)
        return group

    def all_homs(self, items, surjective=True):
        for func in all_functions(self.items, items):
            if surjective:
                if len(set(func.values())) < len(items):
                    continue
            group = self.is_hom(items, func)
            if group is not None:
                group.check()
                yield group, func

    def isomorphic(self, other):
        assert isinstance(other, Group)
        if len(self)!=len(other):
            return False
        n = len(self)
        for i in range(n):
          for j in range(n):
            g = self.perms[i] * self.perms[j]
            g = self.perms.index(g) # slow
            h = other.perms[i] * other.perms[j]
            h = other.perms.index(h) # slow
            if g!=h:
                return False
        return True

    def is_hom_iso(self, other, func):
        for i in range(len(self)):
            g = self[i]
            h = other[i]
            for i in self.items:
                gi = g[i]
                fi = func[i]
                fgi = func[gi]
                hfi = h(fi)
                if hfi != fgi:
                    return # fail
        return True

    def all_homs_iso(self, other, surjective=True):
        assert self.isomorphic(other)
        for func in all_functions(self.items, other.items):
            if surjective:
                if len(set(func.values())) < len(other.items):
                    continue
            if self.is_hom_iso(other, func):
                yield func


class Coset(Group):
    def intersect(self, other):
        assert self.items == other.items
        perms = self.set_perms.intersection(other.set_perms)
        return Coset(perms, self.items)
    intersection = intersect


def conjugacy_subgroups(G, Hs=None):
    Hs = G.conjugacy_subgroups(Hs)

    return Hs



class Action(object):
    """
        A Group _acting on a set, possibly with a basepoint.
        For each perm in the source Group G we map to a perm of items.
    """
    def __init__(self, G, send_perms, items, basepoint=None, check=False):
        assert isinstance(G, Group)
        self.G = G
        assert isinstance(send_perms, dict)
        self.send_perms = dict(send_perms)
        self.items = list(items)
        self.basepoint = basepoint
        self.repr = {}
        if basepoint is not None:
            assert basepoint in items
            for g in G:
                item = self.send_perms[g](basepoint)
                self.repr[item] = g
        if check:
            self.check()

    @property
    def src(self):
        assert 0, "use .G"
        return self.G

    # Equality on-the-nose:
    def __eq__(self, other):
        assert isinstance(other, Action)
        return (self.G==other.G and self.send_perms==other.send_perms)

    def __str__(self):
        return "Action(%s, %s)"%(self.G, len(self.items))

    def __ne__(self, other):
        assert isinstance(other, Action)
        return (self.G!=other.G or self.send_perms!=other.send_perms)

    def __hash__(self):
        send_perms = self.send_perms
        send_perms = tuple((perm, send_perms[perm]) for perm in self.G)
        return hash((self.G, send_perms))

    @classmethod
    def identity(cls, G, check=False):
        send_perms = dict((g, g) for g in G)
        return cls(G, send_perms, G.items, check=check)

    def check(self):
        G, items, send_perms = self.G, self.items, self.send_perms

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
            assert send_perms[g1*g2] == h1*h2

    def __call__(self, g):
        perm = self.send_perms[g]
        return perm

#    def __getitem__(self, g): # use __call__ ???
#        assert 0, "use __call__"
#        perm = self.send_perms[g]
#        return perm
#    # should __getitem__ index into .items ?? seems more useful/meaningful...

    # i am a set (list) of items
    def __getitem__(self, idx):
        return self.items[idx]

    # with a len'gth
    def __len__(self):
        return len(self.items)

    def __contains__(self, x):
        return x in self.items # ouch it's a list

    def get_repr(self, x):
        return self.repr[x]

    def rename(self, send_items, items):
        #G, items, send_perms = self.G, self.items, self.send_perms
        for item in self.items:
            assert item in send_items
        new_send = {}
        for g in self.G.perms:
            assert g in self.send_perms
            perm = self.send_perms[g]
            perm = perm.rename(send_items, items)
            new_send[g] = perm
        return Action(self.G, new_send, items)

    def coproduct(*Xs):
        self = Xs[0]
        for X in Xs:
            assert X.G is self.G
        items = [(i, x) for (i,X) in enumerate(Xs) for x in X]
        #print(items)
        send_perms = {}
        for g in self.G:
            perm = {}
            for (i,X) in enumerate(Xs):
                for x in X:
                    perm[i,x] = (i, X(g)[x])
            #print(perm)
            perm = Perm(perm, items)
            send_perms[g] = perm
        return Action(self.G, send_perms, items)
    __add__ = coproduct

    def product(self, other): # HOTSPOT
        assert self.G == other.G
        items = []
        for a1 in self.items:
          for a2 in other.items:
            items.append((a1, a2))
        send_perms = {}
        for g in self.G:
            perm = {}
            h1 = self.send_perms[g]
            h2 = other.send_perms[g]
            for a1, a2 in items:
                perm[(a1, a2)] = h1(a1), h2(a2)
            perm = Perm(perm, items)
            send_perms[g] = perm
        return Action(self.G, send_perms, items)
    __mul__ = product

    def hecke(self, other):
        import numpy
        assert self.G == other.G
        m, n = (len(self.items), len(other.items))
        marked = set((i, j) for i in range(m) for j in range(n))
        assert marked
        while len(marked):
            H = numpy.zeros((m, n), dtype=numpy.float64)
            i, j = iter(marked).__next__()
            marked.remove((i, j))
            H[i, j] = 1
            ai = self.items[i]
            aj = other.items[j]
            for g in self.G:
                bi = self.send_perms[g][ai]
                bj = other.send_perms[g][aj]
                ii = self.items.index(bi)
                jj = other.items.index(bj)
                if H[ii, jj] == 0:
                    H[ii, jj] = 1
                    marked.remove((ii, jj))
            yield H

    def orbits(self, G=None):
        #print "orbits"
        #print self.perms
        #print self.items
        if G is None:
            G = self.G
        send_perms = self.send_perms
        remain = set(self.items)
        orbits = []
        while remain:
            #print "remain:", remain
            item = iter(remain).__next__()
            orbit = set(send_perms[g](item) for g in G.perms)
            #print "orbit:", orbit
            for item in orbit:
                remain.remove(item)
            orbits.append(orbit)
        return orbits

    def get_components(self):
        G = self.G
        orbits = self.orbits()
        actions = []
        for orbit in orbits:
            send_perms = {}
            perms = []
            for perm in G:
                perm1 = self.send_perms[perm].restrict(orbit)
                send_perms[perm] = perm1
                perms.append(perm1)
            actions.append(Action(G, send_perms, orbit))
        return actions
    components = get_components # backward compat

    def get_atoms(self):
        G = self.G
        orbits = self.orbits()
        homs = []
        for orbit in orbits:
            send_perms = {}
            send_items = {item:item for item in orbit} # inclusion
            perms = []
            for perm in G:
                perm1 = self.send_perms[perm].restrict(orbit)
                send_perms[perm] = perm1
                perms.append(perm1)
            src = Action(G, send_perms, orbit)
            hom = Hom(src, self, send_items)
            homs.append(hom)
        return homs

    def _find_homs_atomic(X, Y):
        assert X.G is Y.G
        if not len(X.items):
            yield Hom(X, Y, {})
            return
        G = X.G
        #x = iter(X.items).__next__()
        x = X[0]
        for y in Y:
            send_items = {x:y}
            for g in G:
                gx = X(g)[x]
                gy = Y(g)[y]
                _gy = send_items.get(gx)
                if _gy is None:
                    send_items[gx] = gy
                elif _gy != gy:
                    break # not a Hom
            else:
                #print("Hom", send_items)
                yield Hom(X, Y, send_items)

    def find_homs(X, Y):
        assert X.G is Y.G
        Xs = X.get_components()
        Ys = Y.get_components()
        #print("find_homs")
        #print(Xs, Ys)
        for section in cross([Ys]*len(Xs)):
            homss = []
            for (Xi,Yi) in zip(Xs, section):
                homs = list(Xi._find_homs_atomic(Yi))
                homss.append(homs)
            for homs in cross(homss):
              send_items = {}
              for hom in homs:
                send_items.update(hom.send_items)
              #print(send_items)
              yield Hom(X, Y, send_items)

    def get_graph(self):
        graph = isomorph.Graph()
        G = self.G
        send_perms = self.send_perms
        items = self.items
        n = len(items)
        # fibres are all the same:
        fibres = [graph.add("fibre") for item in items]
        for i, perm in enumerate(G):
            lookup = dict((item, graph.add()) for item in items)
            for j in range(n):
                graph.join(fibres[j], lookup[items[j]])
            # each layer is unique and contains the "action of g"
            layer = graph.add("layer_%d"%i)
            for item in items:
                graph.join(lookup[item], layer)
                gitem = send_perms[perm](item)
                if gitem == item:
                    # don't need an edge here
                    continue
                graph.add_directed(lookup[item], lookup[gitem])
        return graph

    def get_shape(self):
        G = self.G
        send_perms = self.send_perms
        shape = []
        for perm in G:
            perm = send_perms[perm]
            shape.append(perm.conjugacy_cls())
        return shape

    def isomorphisms(self, other):
        """ yield all isomorphisms in the category of G-sets.
            Each isomorphism is a dict:item->item 

            self  maps G -> H1
            other maps G -> H2
            an isomorphism maps H1.items -> H2.items
        """
        assert isinstance(other, Action)
        assert self.G == other.G
        n = len(self.items)
        if n != len(other.items):
            return
        if self.get_shape() != other.get_shape():
            return
        graph0 = self.get_graph()
        graph1 = other.get_graph()
        for fn in isomorph.search(graph0, graph1):
        #for fn in isomorph.search_recursive(graph0, graph1):
            send_items = {}
            for i in range(n):
                send_items[self.items[i]] = other.items[fn[i]]
            yield send_items

    def check_isomorphism(self, other, send_items):
        assert isinstance(other, Action)
        for perm in self.G:
            perm1 = self.send_perms[perm]
            perm2 = other.send_perms[perm]
            for item1 in self.items:
                item2 = send_items[perm1[item1]]
                assert item2 == perm2[send_items[item1]]

    def slow_isomorphic(self, other, check=False):
        "is isomorphic in the category of G-sets"
        assert isinstance(other, Action)
        for send_items in self.isomorphisms(other):
            self.check_isomorphism(other, send_items)
            return True
        return False

    def fixed_points(self, H):
        send_perms = self.send_perms
        fixed = set(self.items)
        for g in H:
            g = send_perms[g]
            fixed = fixed.intersection(g.fixed())
            if not fixed:
                break
        return fixed

    def signature(self, Hs=None):
        if Hs is None:
            Hs = self.G.subgroups()
        return [len(self.fixed_points(H)) for H in Hs]

    def isomorphic(self, other):
        return self.signature() == other.signature()

    def _refute_isomorphism(self, other, Hs):
        "return True if it is impossible to find an isomorphism"
        assert isinstance(other, Action)
        assert self.G == other.G
        n = len(self.items)
        if n != len(other.items):
            return True

        n = len(self.G)
        for perm in self.G:
            perm1 = self.send_perms[perm]
            perm2 = other.send_perms[perm]
            if perm1.conjugacy_cls() != perm2.conjugacy_cls():
                return True

        # Not able to refute
        return False

    def refute_isomorphism(self, other, Hs):
        # see: http://math.stackexchange.com/a/1891096/360303
        ref = self._refute_isomorphism(other, Hs)
        if ref==True:
            assert self.signature(Hs) != other.signature(Hs)
        elif self.signature(Hs) != other.signature(Hs):
            ref = True
        return ref


class Hom(object):
    "Hom'omorphism of Action's"
    def __init__(self, src, tgt, send_items, check=True):
        assert isinstance(src, Action)
        assert isinstance(tgt, Action)
        G = src.G
        assert G == tgt.G
        self.src = src
        self.tgt = tgt
        self.G = G
        self.send_items = dict(send_items)
        if check:
            self.do_check()

    def do_check(self):
        src = self.src
        tgt = self.tgt
        send_items = self.send_items
        for x,y in send_items.items():
            assert x in src
            assert y in tgt
        for x in src:
            assert x in send_items
        for g in self.G:
            for item in src.items:
                left = send_items[src(g)[item]]
                right = tgt(g)[send_items[item]]
                if left != right:
                    print("item =", item)
                    print("g =", g, "src(g) =", src(g), "tgt(g) =", tgt(g))
                    print("send_items =", send_items)
                    print("%s != %s"%(left, right))
                    assert 0, "not a Hom of Action's"

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
        send_items = [b[i] for i in a]
        return Hom(self.src, other.tgt, send_items)

    def __mul__(self, other):
        assert isinstance(other, Hom)
        return other.compose(self)

#    def mul(f, g):
#        assert isinstance(g, Hom)
#        cone = Action.mul(f.src, g.src)
#        src = cone.apex
#        cone = Cone(src, [cone[0].compose(f), cone[1].compose(g)])
#        cone, univ = Action.mul(f.tgt, g.tgt, cone)
#        return univ
#    __mul__ = mul




r"""
Here we calculate the Burnside Ring corresponding to a particular G-set.

Finite kleinian geometry: 
(1) a set of "points" and a group G _acting on the set.
(2) (conjugacy classes of) subgroups of G are the "geometric properties" 

For example, a triangle.
This has symmetry group of order 6, and lattice
of subgroups:

      6
     / \
    /   \
    3    2
    \   /
     \ /
      1

We name the corresponding properties as:

      Frames
     / \
    /   \
Points  Orientations
    \   /
     \ /
      Nothing


More details:

We consider the group G = S_3 the permutation group on three things.
This is also the symmetry group of the triangle.

This group has 6 elements, and has 6 subgroups.
These subgroups have orders: 1, 2, 2, 2, 3 and 6.

Each subgroup preserves (stabilizes) something, which we
think of as a type of figure. 
The one element subgroup preserves any FRAME (which we call F)
Each 2 element subgroup preserves a POINT (P) 
The 3 element subgroup preserves any ORIENTATION (O)
The 6 element subgroup preserves any NOTHING (N)

The 2 element subgroups are all conjugates of
each other, so these three subgroups are all considered the "same" subgroup.
Conjugate subgroups correspond to
isomorphic objects in the category of G-sets.

Now, an actual thing of type F, P, O, or N corresponds to
a coset of the subgroup for that type.

Therefore, we have 6 things of type FRAME,
3 things of type POINT,
2 things of type ORIENTATION,
and 1 thing of type NOTHING.

"""


def test_action():

    for n in range(3, 6):
        items = range(n)

        G = Group.trivial(items, check=True)
        assert len(G.components()) == len(items)
    
        action = Action.identity(G)
        assert len(list(action.isomorphisms(action))) == factorial(n)
    
        G = Group.cyclic(items, check=True)
        assert len(G.components()) == 1
    
        action = Action.identity(G)
        assert len(list(action.isomorphisms(action))) == n
    
        H = Group.cyclic(list("abcdef"[:n]))
        send_perms = dict((G[i], H[i]) for i in range(len(G)))
        action1 = Action(G, send_perms, H.items)
        assert len(list(action.isomorphisms(action1))) == n

        G = Group.symmetric(items, check=True)
        assert len(G.components()) == 1


def get_P2():

    # Pauli group
    items = "+00 -00 +01 -01 +10 -10 +11 -11".split()
    #          0   1   2   3   4   5   6   7
    II = Perm((0, 1, 2, 3, 4, 5, 6, 7), items)
    XI = Perm((4, 5, 6, 7, 0, 1, 2, 3), items)
    IX = Perm((2, 3, 0, 1, 6, 7, 4, 5), items)
    ZI = Perm((0, 1, 2, 3, 5, 4, 7, 6), items)
    IZ = Perm((0, 1, 3, 2, 4, 5, 7, 6), items)
    w2 = Perm((1, 0, 3, 2, 5, 4, 7, 6), items)
    assert XI*XI==II
    assert ZI*ZI==II
    assert IX*IX==II
    assert IZ*IZ==II
    assert ZI*XI != XI*ZI
    assert ZI*XI == w2*XI*ZI
    assert ZI*XI*ZI*XI == w2
    assert IZ*IX != IX*IZ
    assert IZ*IX == w2*IX*IZ
    assert IZ*IX*IZ*IX == w2

    assert ZI*IX == IX*ZI

    P2 = Group.generate([XI, ZI, IX, IZ])
    assert len(P2)==32

    return P2


def main():

    n = argv.get("n", 3)
    items = range(n)


    if argv.dihedral:
        G = Group.dihedral(items, check=True)
        assert len(G.components()) == 1

    elif argv.cyclic:
        G = Group.cyclic(items, check=True)
        assert len(G.components()) == 1

    elif argv.symmetric or argv.sym:
        G = Group.symmetric(items, check=True)
        assert len(G.components()) == 1

    elif argv.alternating or argv.alt:
        G = Group.alternating(items, check=True)
        assert len(G.components()) == 1

    elif argv.A_2:
        G = Group.symmetric(range(3), check=True)

    elif argv.A_3:
        G = Group.symmetric(range(4), check=True)

    elif argv.A_4:
        G = Group.symmetric(range(5), check=True)

    elif argv.B_2:
        items = range(8)
        gen = [
            Perm({0:0, 1:2, 2:1, 3:3, 4:6, 5:7, 6:4, 7:5}, items), 
            Perm({0:1, 1:0, 2:3, 3:2, 4:4, 5:5, 6:7, 7:6}, items)]
        perms = mulclose(gen)
        G = Group(perms, items)
        assert len(G)==8

    elif argv.Q_8:
        items = range(8)
        gen = [
            Perm({0:1, 1:3, 3:6, 6:0, 2:5, 5:7, 7:4, 4:2}, items),
            Perm({0:2, 2:3, 3:7, 7:0, 1:4, 4:6, 6:5, 5:1}, items)]
        perms = mulclose(gen)
        G = Group(perms, items)
        assert len(G)==8

    elif argv.Q_8C_3:

        items = range(8)
        gen = [
            Perm({0:1, 1:3, 3:6, 6:0, 2:5, 5:7, 7:4, 4:2}, items),
            Perm({0:2, 2:3, 3:7, 7:0, 1:4, 4:6, 6:5, 5:1}, items)]
        perms = mulclose(gen)
        G = Group(perms, items)
        assert len(G)==8

        C3 = Group.cyclic(range(3))
        G = C3.direct_product(G)

    elif argv.P2:
        G = get_P2()

    elif argv.B_3:
        items = range(18)
        gen = [
            Perm({0:0, 1:2, 2:1, 3:3, 4:8, 5:9, 6:10, 7:11, 8:4, 9:5, 
                10:6, 11:7, 12:14, 13:15, 14:12, 15:13, 16:16, 17:17}, items),
            Perm({0:4, 1:5, 2:6, 3:7, 4:0, 5:1, 6:2, 7:3, 8:8, 9:10,
                10:9, 11:11, 12:12, 13:13, 14:16, 15:17, 16:14, 17:15}, items),
            Perm({0:0, 1:1, 2:2, 3:3, 4:5, 5:4, 6:7, 7:6, 8:9, 9:8,
                10:11, 11:10, 12:12, 13:13, 14:14, 15:15, 16:17, 17:16}, items)]
        perms = mulclose(gen)
        G = Group(perms, items)

    elif argv.D_4:
        items = range(24)
        gen = [
            Perm({0:0, 1:2, 2:1, 3:3, 4:12, 5:13, 6:14, 7:15, 8:16,
                9:17, 10:18, 11:19, 12:4, 13:5, 14:6, 15:7, 16:8, 17:9,
                18:10, 19:11, 20:20, 21:21, 22:22, 23:23}, items),
            Perm({0:4, 1:5, 2:6, 3:7, 4:0, 5:1, 6:2, 7:3, 8:8, 9:9,
                10:10, 11:11, 12:12, 13:14, 14:13, 15:15, 16:20, 17:21,
                18:22, 19:23, 20:16, 21:17, 22:18, 23:19}, items),
            Perm({0:0, 1:1, 2:2, 3:3, 4:8, 5:9, 6:10, 7:11, 8:4,
                9:5, 10:6, 11:7, 12:16, 13:17, 14:18, 15:19, 16:12, 17:13,
                18:14, 19:15, 20:20, 21:22, 22:21, 23:23}, items),
            Perm({0:0, 1:1, 2:2, 3:3, 4:9, 5:8, 6:11, 7:10, 8:5,
                9:4, 10:7, 11:6, 12:17, 13:16, 14:19, 15:18, 16:13, 17:12,
                18:15, 19:14, 20:23, 21:21, 22:22, 23:20}, items)
        ]
        perms = mulclose(gen)
        G = Group(perms, items)

    elif argv.projective:
        from bruhat import geometry
        g = geometry.projective(3)
        items = range(14)
        perms = [f for f in g.get_symmetry()]
        perms = [Perm(f, items) for f in perms]
        G = Group(perms, items)

    elif argv.fano:
        from bruhat import geometry
        g = geometry.fano()
        keys = list(g.items)
        keys.sort(key = str)
        #print keys
        lookup = dict((v, k) for (k, v) in enumerate(keys))
        points = keys[:7]
        lines = keys[7:]
        flags = []
        for l in lines:
            for p in l:
                flags.append((p, l))
        assert len(flags)==21
        perms = []
        for f in g.get_symmetry():
            perm = {}
            for i, flag in enumerate(flags):
                glag = keys[f[lookup[flag[0]]]], keys[f[lookup[flag[1]]]]
                perm[flag] = glag
            perm = Perm(perm, flags)
            perms.append(perm)
        G = Group(perms, flags)

    elif argv.schreier:
        schreier()
        return

    elif argv.mathieu:
        mathieu()
        return

    elif argv.desargues:
        desargues()
        return

    else:
        return

    check = argv.check

    if 0:
        shapes = list(G.shape_item())
        #print len(shapes), len(set(shapes))
        #for shape in shapes:
        #    print shape
        return

    if 0:
        action = Action.identity(G)
        #for iso in action.isomorphisms(action):
        #    print iso
    
        return

    #print "G:", len(G)

    if argv.test_projective:
        test_projective(G)

    if argv.conjugacy_subgroups:
        conjugacy_subgroups(G)
        return

    if argv.cgy_cls:
        itemss = G.cgy_cls()
        print("|G|=", len(G))
        for gs in itemss:
            gs = list(gs)
            print("Order:", gs[0].order(), "Size:", len(gs))

    if argv.burnside:
        burnside(G)
        return

    if argv.hecke:
        hecke(G)
        return

    if argv.orbiplex:

        #Gs = G.components()
        #G = Gs[0]
        #orbiplex(G)

        if argv.regular_rep:
            G = G.regular_rep()

        if argv.subgroups:
            Hs = conjugacy_subgroups(G)
        else:
            Hs = [G]

        for H in Hs:
            if len(H)==1:
                continue
            orbiplex(H)

    elif argv.subgroups or argv.cyclic_subgroups:
        if argv.cyclic_subgroups:
            Hs = G.cyclic_subgroups()
        else:
            Hs = G.subgroups()
        #print "subgroups:", len(Hs)
        Hs = list(Hs)
        Hs.sort(key = lambda H : len(H))
        for H in Hs:
            #if len(H)==len(G) or len(H)==1:
            #    continue
            #print "subgroup order=%d:"%len(H)
            perms = [perm.perm for perm in H]
            #print perms


def orbiplex(G, k=4):

    import numpy
    from gelim import zeros, dotx, rank, nullity

    #print "orbiplex: |G|=%d" % len(G)
    items = G.items

    nchains = {} # map tuple -> index
    dims = []
    bdys = []
    k = argv.get("k", k)
    for n in range(k):

        write("|C_%d| ="%n)
        #if n > len(items):
        #    break

        tpls = list(uniqtuples(items, n))
    
        perms = []
        for perm in G:
            _perm = {}
            for key in tpls:
                value = tuple(perm[i] for i in key)
                _perm[key] = value
            _perm = Perm(_perm, tpls)
            perms.append(_perm)
    
        G2 = Group(perms, tpls)

        orbits = list(G2.orbits())
        d = len(orbits)
        dims.append(d)
        write("%d,"%d)
        for idx, orbit in enumerate(orbits):
            for key in orbit:
                nchains[key] = idx

        if n==0:
            # skip first guy ?
            bdy = zeros(0, d)
            bdys.append(bdy)
            continue # <------- continue

        bdy = zeros(dims[-2], dims[-1])
        for idx, orbit in enumerate(orbits):
            key = iter(orbit).__next__()
            c = 1
            for m in range(n):
                # all the keys should be the same...
                dkey = key[:m] + key[m+1:]
                jdx = nchains[dkey]
                bdy[jdx, idx] += c
                c *= -1
        #print bdy
        bdys.append(bdy)

        if n>len(items):
            break

    print

    euler = 0
    c = -1
    for i, dim in enumerate(dims):
        euler += c*dim
        c *= -1
    #print "euler:", euler

    for i in range(len(bdys)-1):
        A = bdys[i]
        B = bdys[i+1]
        C = dotx(A, B)
        #print "C:"
        #print C
        assert numpy.abs(C).sum() == 0

        #   B      A
        #  ---> . ---> 
        im = rank(B)
        ker = nullity(A)
        hom = ker - im
        assert hom>=0
        #if hom > 0: # and B.shape[0]*B.shape[1]:
        #    print B.shape, A.shape
        #    print "im=%d, ker=%d" % (im, ker)
        #    print "H_%d = %d" % (i, hom)



def test_projective(G):

    Hs = G.cyclic_subgroups()

    #print [len(H) for H in Hs]
    H2s = [H for H in Hs if len(H)==2]
    H3s = [H for H in Hs if len(H)==3]

    pairs = list(zip(H2s, H3s))
    #print "pairs:", len(pairs)
    shuffle(pairs)
    for H2, H3 in pairs:
        H = Group.product(H2, H3)
        if len(H)!=6:
            write('.')
            continue
    else:
        assert 0

    #print
    #print H.is_abelian()
    #for perm in H:
    #    print perm.cycle_str()


def todot(graph, names={}, filename="graph.dot"):

    lookup = dict((node, i) for (i, node) in enumerate(graph.nodes()))
    f = open(filename, "w")
    if "Di" in graph.__class__.__name__:
        f.write("digraph the_graph\n")
        link = "->"
    else:
        f.write("graph the_graph\n")
        link = "--"
    f.write("{\n")
    for edge in graph.edges():
        i = lookup[edge[0]]
        j = lookup[edge[1]]
        i = names.get(edge[0], "n%s"%i)
        j = names.get(edge[1], "n%s"%j)
        f.write("    %s %s %s;\n" % (i, link, j))
    f.write("}\n")
    f.close()


def schreier():

    import networkx as nx

    # --------------- group ----------------
    n = argv.get("n", 5)
    check = argv.get("check", True)

    items = range(n)
    sigmas = []
    for i in range(n-1):
        perm = dict((item, item) for item in items)
        perm[items[i]] = items[i+1]
        perm[items[i+1]] = items[i]
        sigmas.append(perm)
    sigmas = [Perm(perm, items) for perm in sigmas]
    G = Group.generate(sigmas, check=check)

#    # Pauli group
#    items = range(4)
#    X = Perm((3, 1, 2, 0), items)
#    Z = Perm((2, 3, 0, 1), items)
#    I = Perm((0, 1, 2, 3), items)
#    w = Perm((3, 2, 1, 0), items)
#    assert X*X==I
#    assert Z*Z==I
#    assert Z*X != X*Z
#    assert Z*X*Z*X == w
#    G = Group.generate([X, Z])
#    assert len(G)==8

    # --------- subgroups ----------------
#    Hs = conjugacy_subgroups(G)

    s0 = sigmas[0]
    s1 = sigmas[1]
    s2 = sigmas[2]
    s3 = sigmas[3]
#    gen = [Perm(p, items) for p in 
#        [[1, 2, 0, 3, 4],
#         [1, 0, 2, 3, 4],
#         [0, 1, 2, 4, 3]]]
    gen = [s0, s1, s3]
    H = Group.generate(gen, check=check)
    assert len(H) == 12
    H = Group([g for g in H if g.sign()==1], items, check=check)
    assert len(H) == 6
    Hs = [H]

#    gen = [X*Z, w]
#    H = Group.generate(gen, check=check)
#    print "len(H):", len(H)
#    Hs = [H]

#    s0 = sigmas[0]
#    s1 = sigmas[1]
#    gen = [s0]
#    H = Group.generate(gen, check=check)
#    print "len(H):", len(H)
#    Hs = [H]

    # ---------- generators for schreier graph --------
#    print "\n"*5
#    gen = [
#        Perm([0,3,4,1,2], items), 
#        Perm([2,3,4,0,1], items), 
#        Perm([1,3,4,0,2], items)]
    gen = [ [3,0,4,1,2], [3,2,4,0,1], [1,3,4,0,2], ]
    gen = [ [1, 3, 4, 0, 2] ]
    gen = [ [3, 1, 4, 2, 0]]
    gen = [ [3, 4, 1, 0, 2]]
    gen = [ [4, 1, 3, 0, 2]]
    gen = [ [4, 3, 1, 2, 0]]
    gen = [ Perm(dict(enumerate(perm)), items) for perm in gen]
    gen = list(sigmas)
    #for g in gen:
    #    print g.get_idxs()
    #    assert g.sign()==-1 # Not
    #    print (g*g).get_idxs()
    #    print

#    _gen = []
#    for g in gen:
#        for _g in G:
#            if _g==g:
#                _gen.append(_g)
#                break
#    gen = _gen

#    #gen = [X, Z]
#    #gen = [s0, s1]
#
#    #return
#    gen = list(sigmas)

    cls = nx.MultiGraph
    for perm in gen:
        #print perm
        if not (perm*perm).is_identity():
            cls = nx.MultiDiGraph

    cosets = G.left_cosets(H)

    #print "H:", len(H),
    #print "cosets:", len(cosets)

    graph = cls()
    edges = set()
    hom = G.left_action(cosets)

#    for g in gen:
#      for h in gen:
#        if g is h:
#            continue
#        perm = hom.send_perms[g]
#        qerm = hom.send_perms[h]
#        print perm.intersection(qerm)
#        idxs = perm.get_idxs()
#        jdxs = qerm.get_idxs()
#        for i in range(len(idxs)):
#            if idxs[i]==jdxs[i]:
#                print i
#        #print len(idxs), len(set(idxs))
#    return

    names = {}
    for ii, i in enumerate(cosets):
        #print "coset:", i
        strs = []
        for g in i:
            idxs = g.get_idxs()
            strs.append(''.join(str(idx) for idx in idxs))
        strs.sort()
        #print "%2d"%ii, ' '.join(strs)
        #print "   ",
        graph.add_node(i, name=strs[0])
        names[i] = strs[0]
        for g in gen:
            j = hom.send_perms[g][i]
            #print "j:", j
            assert j in cosets
            if cls == nx.MultiGraph and (j, i) in edges:
                continue
            edges.add((i, j))
            #print "(%s -> %s)" % (cosets.index(i), cosets.index(j)),
            graph.add_edge(i, j)
        #print
    todot(graph, names)


def desargues():
    
    import networkx as nx

    # --------------- group ----------------
    n = 5
    check = argv.get("check", True)

    items = range(n)
    sigmas = []
    for i in range(n-1):
        perm = dict((item, item) for item in items)
        perm[items[i]] = items[i+1]
        perm[items[i+1]] = items[i]
        sigmas.append(perm)
    sigmas = [Perm(perm, items) for perm in sigmas]
    G = Group.generate(sigmas, check=check)

    configs = [(tuple(g.get_idxs()), parity) for g in G for parity in (+1, -1)]
    configs.sort()
    #print configs

    for g in G:
        #if g.get_idxs() != [3, 0, 4, 1, 2]:
        #if g.get_idxs() != [1, 3, 4, 0, 2]:
        if g.get_idxs() != [3, 2, 4, 0, 1]:
            continue
        #print g, g.sign()
        assert g.sign() == -1
        perm = {}
        for config in configs:
            _config = (tuple(g(i) for i in config[0]), g.sign()*config[1])
            #print '   ', config, "->", _config

    return

    # --------- subgroups ----------------
#    Hs = conjugacy_subgroups(G)

    s0 = sigmas[0]
    s1 = sigmas[1]
    s2 = sigmas[2]
    s3 = sigmas[3]
#    gen = [Perm(p, items) for p in 
#        [[1, 2, 0, 3, 4],
#         [1, 0, 2, 3, 4],
#         [0, 1, 2, 4, 3]]]
    gen = [s0, s1, s3]
    H = Group.generate(gen, check=check)
    assert len(H) == 12
    H = Group([g for g in H if g.sign()==1], items, check=check)
    assert len(H) == 6
    Hs = [H]

    # ---------- generators for schreier graph --------
#    print "\n"*5
#    gen = [
#        Perm([0,3,4,1,2], items), 
#        Perm([2,3,4,0,1], items), 
#        Perm([1,3,4,0,2], items)]
    gen = [ [0,3,4,1,2], [2,3,4,0,1], [1,3,4,0,2], ]
    gen = [ Perm(dict(enumerate(perm)), items) for perm in gen]
    #for g in gen:
        #print g.get_idxs()
        #assert g.sign()==-1 # Not
        #print (g*g).get_idxs()
        #print
    return

#    _gen = []
#    for g in gen:
#        for _g in G:
#            if _g==g:
#                _gen.append(_g)
#                break
#    gen = _gen

#    #gen = [X, Z]
#    #gen = [s0, s1]
#
#    #return
#    gen = list(sigmas)

    cls = nx.MultiGraph
    for perm in gen:
        #print perm
        if not (perm*perm).is_identity():
            cls = nx.MultiDiGraph

    cosets = G.left_cosets(H)

#    print "H:", len(H),
#    print "cosets:", len(cosets)
    graph = cls()
    edges = set()
    hom = G.left_action(cosets)
    for ii, i in enumerate(cosets):
        #print "coset:", i
        strs = []
        for g in i:
            idxs = g.get_idxs()
            strs.append(''.join(str(idx) for idx in idxs))
        strs.sort()
#        print "%2d"%ii, ' '.join(strs)
#        print "   ",
        graph.add_node(i, name=strs[0])
        for g in gen:
            j = hom.send_perms[g][i]
            #print "j:", j
            assert j in cosets
            if cls == nx.MultiGraph and (j, i) in edges:
                continue
            edges.add((i, j))
#            print "(%s -> %s)" % (cosets.index(i), cosets.index(j)),
            graph.add_edge(i, j)
#        print
    todot(graph)


def mathieu():

    # https://en.wikipedia.org/wiki/Mathieu_group
    # http://math.ucr.edu/home/baez/week234.html

    # Generate M_12
    # http://www.neverendingbooks.org/monsieur-mathieu
    items = range(12)
    "(1,2)(3,4)(5,8)(7,6)(9,12)(11,10)"
    "(0,1)(2,3)(4,7)(6,5)(8,11)(10,9)"
    a = Perm({0:1, 1:0, 2:3, 3:2, 4:7, 7:4, 6:5, 5:6, 8:11, 11:8, 10:9, 9:10}, items)
    b = Perm({0:1, 1:2, 2:0, 3:4, 4:5, 5:3, 6:6, 7:8, 8:9, 9:7, 10:10, 11:11}, items)
    perms = mulclose([a, b])
    M12 = Group(perms, items)
    # M12 sharply 5-transitive, simple sporadic
    assert len(M12)==95040

    # M11 sharply 4-transitive, simple sporadic
    perms = [a for a in perms if a(11)==11]
    M11 = Group(perms, items) # remove 11 ?
    assert len(M11)==7920

    # M10 sharply 3-transitive, = Alt(6) ?
    perms = [a for a in perms if a(10)==10]
    M10 = Group(perms, items) # remove 10 ?
    assert len(M10)==720

    # M9 sharply 2-transitive
    perms = [a for a in perms if a(9)==9]
    M9 = Group(perms, items) # remove 9 ?
    assert len(M9)==72

    # M8 sharply 1-transitive = Quaternion group
    perms = [a for a in perms if a(8)==8]
    M8 = Group(perms, items) # remove 8 ?
    assert len(M8)==8

    print("OK")


def hecke(G):

    import numpy
    from solve import shortstr

    Hs = conjugacy_subgroups(G)

    letters = list(string.ascii_uppercase + string.ascii_lowercase)
    letters = letters + [l+"'" for l in letters] + [l+"''" for l in letters]
    assert len(letters) >= len(Hs)
    letters = letters[:len(Hs)]

    homs = []
    for i, H in enumerate(Hs):
        cosets = G.left_cosets(H)
        assert len(G) == len(cosets)*len(H)

        hom = G.left_action(cosets)
        assert hom.src is G
        hom.name = letters[i]
        homs.append(hom)
        assert len(hom.components())==1 # transitive
        #print hom.tgt.perms
        #print "%s subgroup order = %d, number of cosets = %d" %(
        #    hom.name, len(H), len(cosets))

#    for i in range(len(homs)):
#      for j in range(i, len(homs)):
#        A = homs[i]
#        B = homs[j]
#        print "%s * %s" % (A.name, B.name)
#        for H in A.hecke(B):
#            print shortstr(H)
#            print "sum:", H.sum()
#            print

#        C = A.product(B)
#        assert C.src is G
#
#        #print C.items
#        print C.send_perms.keys()


r"""
\begin{array}{c|lcr}
n & \text{Left} & \text{Center} & \text{Right} \\
\hline
1 & 0.24 & 1 & 125 \\
2 & -1 & 189 & -8 \\
3 & -20 & 2000 & 1+10i
\end{array}
"""

class CFunc(object):
    """ class function on a group 
    """

    def __init__(self, G, func):
        self.G = G
        self.func = dict(func)

    def __str__(self):
        return "CFunc([%s])"%(', '.join(list(str(self.func[g]) for g in self.G)))
    __repr__ = __str__

    def dot(self, other, normalize=True):
        assert self.G is other.G
        x = 0
        for g in self.G:
            x += self[g] * other[~g]
        n = len(self.G)
        #print("x =", x)
        #assert x%n == 0, x
        if normalize:
            x = x/n
        if type(x) is float and x == int(x):
            x = int(x)
        return x

    def frobenius_schur(self):
        "The frobenius_schur indicator"
        # For a complex irreducible character,
        # 1: real, 0: complex, -1: quaternionic
        val = sum(self.func[g*g] for g in self.G)
        val = val/len(self.G)
        return val

    def __getitem__(self, g):
        return self.func[g]

#    def __len__(self):
#        return len(self.func)

    def __add__(self, other):
        assert self.G is other.G
        func = {}
        for g in self.G:
            func[g] = self[g] + other[g]
        return CFunc(self.G, func)

    def __sub__(self, other):
        assert self.G is other.G
        func = {}
        for g in self.G:
            func[g] = self[g] - other[g]
        return CFunc(self.G, func)

    def __mul__(self, other):
        assert self.G is other.G
        func = {}
        for g in self.G:
            func[g] = self[g] * other[g]
        return CFunc(self.G, func)

    def __rmul__(self, r):
        assert isinstance(r, int)
        func = {}
        for g in self.G:
            func[g] = r*self[g]
        return CFunc(self.G, func)

    @classmethod
    def from_action(cls, action):
        G = action.G
        func = {}
        for g in G:
            ga = action(g)
            func[g] = len(ga.fixed())
        return cls(G, func)

    @classmethod
    def show_table(cls, chis):
        G = chis[0].G
        itemss = G.cgy_cls()
        els = []
        for gs in itemss:
            gs = list(gs)
            g = gs[0]
            els.append(g)
            print("Order:", g.order(), "Size:", len(gs))
        for chi in chis:
            assert chi.G is G
            for g in els:
                print(str(chi[g]).rjust(3), end=" ")
            print()

    @classmethod
    def latex_table(cls, chis):
        G = chis[0].G
        itemss = G.cgy_cls()
        rows = [[r"\mathrm{class}"], [r"\mathrm{size}"], r"\hline"]
        els = []
        for gs in itemss:
            gs = list(gs)
            g = gs[0]
            els.append(g)
            rows[0].append(g.order())
            rows[1].append(len(gs))
        for idx, chi in enumerate(chis):
            assert chi.G is G
            side = "V_{%d}"%(idx+1)
            row = [side]+[chi[g] for g in els]
            rows.append(row)
        s = latex_nosep(rows, "r|" + "r"*len(itemss))
        print(s)


def latex_nosep(rows, desc):
    n = len(rows[0])
    lines = [("$$")]
    lines.append(r"\begin{array}{%s}"%(desc))
    for i, row in enumerate(rows):
        if type(row)==list:
            line = " & ".join(str(fld) for fld in row) + r" \\"
        else:
            line = str(row)
        lines.append(line)
    lines.append(r"\end{array}")
    lines.append("$$")
    s = "\n".join(lines)
    return s



def latex_dump(header, rows, sider=None, sep=True):
    header = list(header)
    n = len(header)
    lines = [("$$")]
    desc = "c|"*n if sep else "c"*n
    if sider:
        desc = "c|"+desc
        header.insert(0, "")
    lines.append(r"\begin{array}{|%s|}"%(desc))
    lines.append(r"\hline")
    lines.append(" & ".join(r"%s"%fld for fld in header) + r" \\")
    lines.append(r"\hline")
    for i, row in enumerate(rows):
        line = ''
        if sider:
            line = "%s & " % sider[i]
        line += " & ".join(str(fld) for fld in row) + r" \\"
        lines.append(line)
        if sep:
            lines.append(r"\hline")
    if not sep:
        lines.append(r"\hline")
    lines.append(r"\end{array}")
    lines.append("$$")
    s = "\n".join(lines)
    s = s.replace('||', '|')
    return s


def burnside(G, Hs=None):

    if argv.cyclic_subgroups:
        Hs = G.cyclic_subgroups()
        Hs = conjugacy_subgroups(G, Hs)

    if Hs is None:
        Hs = conjugacy_subgroups(G)

    letters = list(string.ascii_uppercase + string.ascii_lowercase)
    letters.remove("O")
    letters.remove("o")
    letters = letters + [l+"'" for l in letters] + [l+"''" for l in letters]
    assert len(letters) >= len(Hs)
    letters = letters[:len(Hs)]

    homs = []
    rows = []
    for i, H in enumerate(Hs):
        cosets = G.left_cosets(H)
        assert len(G) == len(cosets)*len(H)

        hom = G.left_action(cosets)
        assert hom.src is G
        hom.name = letters[i]
        H.name = hom.name
        homs.append(hom)
        assert len(hom.components())==1 # transitive
        #print hom.tgt.perms
        row = [hom.name, len(H), len(cosets), len(H.conjugates)]
        print("%s subgroup order = %d, number of cosets = %d, conjugates = %d" %
            tuple(row), end="")
        if H.is_cyclic():
            print(", cyclic")
            row.append(r"\checkmark")
        else:
            print()
            row.append(r"")
        rows.append(row)

    if argv.subgroups_only:
        return

    if argv.latex:
        print()
        print(r"\noindent Subgroups:")
        print()
        s = latex_dump(
            r'\text{subgroup} \text{order} \text{cosets} \text{conjugates} \text{cyclic}'.split(), rows)
        print(s)

    if argv.make_dot:
        arrows = []
        names = [H.name for H in Hs]
        parents = dict((name, []) for name in names)
        for H in Hs:
          for K in Hs:
            # Look for K a subgroup of H
            if len(K) >= len(H):
                continue
            for K1 in K.conjugates:
              if H.is_subgroup(K1):
                arrows.append((K.name, H.name))
                parents[K.name].append(H.name)
                break
        print("digraph")
        print("{")
        print("    rankdir = BT;")
        arrows = list(arrows)
        arrows.sort()
        for src, tgt in arrows:
            factor = False
            for p in parents[src]:
                if tgt in parents[p]:
                    break
            else:
                print("    %s -> %s;" % (src, tgt))
        print("}")
        #return

    if 0:
        # We don't need to do this again: isomorphic homs all
        # come from conjugate subgroups.
        f = quotient_rep(homs, Action.isomorphic)
        homs = list(set(f.values())) # uniq
        #print "homs:", len(homs)

    T = []
    for hom in homs:
        chi = CFunc.from_action(hom)
        print(hom.name, chi)
        T.append(chi)
    #CFunc.show_table(T)
    #return

    table = {}
    width = 0

    for i in range(len(homs)):
      for j in range(i, len(homs)):
        A = homs[i]
        B = homs[j]
        C = A.product(B)
        assert C.src is G
        write("%s*%s ="%(A.name, B.name))
        names = []
        for hom in C.components():
            assert hom.src is G
            name = '?'

            # We know it must be one of these possibilities:
            possible = [hom1 for hom1 in homs if not hom.refute_isomorphism(hom1, Hs)]
            assert possible
            if len(possible)==1:
                name = possible[0].name
                #write('.')

            else:
                assert 0
                return
                #write('/')
                for hom1 in possible:
                    if hom.isomorphic(hom1):
                        name = hom1.name
                        break

            names.append(name)
            if name == '?':
                assert 0

        #print len(C.components())
        uniq = set(names)
        counts = dict((name, 0) for name in uniq)
        for name in names:
            counts[name] += 1
        #print "+".join(names)
        names = list(uniq)
        names.sort()
        ss = []
        for name in names:
            c = counts[name]
            ss.append(name if c==1 else "%s*%s"%(c, name))
        s = '+'.join(ss)
        write(s+' ')
        width = max(width, len(s)+2)
        table[A.name, B.name] = s
        table[B.name, A.name] = s # commutative
      print

    space = 1
    table = dict((k, v.replace("*", "")) for (k, v) in table.items())

    rows = cols = [hom.name for hom in homs]

    if argv.latex:
        print()
        print("$$")
        print(latex_table(table, rows, cols, upper=argv.get("upper")))
        print("$$")

    print()
    s = str(tabulate(table, rows, cols, space))
    print(s)

    import zelim
    import numpy
    A = zelim.parse(s)
    print(repr(A))

    if argv.latex:
        print()
        print(r"\noindent Table of multiplicities:")
        print()
        s = latex_dump(cols, A, sider=cols, sep=False)
        print(s)

    L, B = zelim.zelim(A)
    if argv.latex:
        print()
        print(r"\noindent Upper triangular form:")
        print()
        s = latex_dump(cols, B, sep=False)
        s = s.replace("0 ", ". ")
        print(s)
    else:
        print("H~:")
        print(B)
#        print("UMU:")
#        UMU = numpy.dot(L, numpy.dot(A, L.transpose()))
#        print(UMU)

    print("rank:", len(B))

    # diagonal matrix
    LT = numpy.dot(L, T)
    for chi1 in LT:
      for chi2 in LT:
        print(chi1.dot(chi2), end=" ")
      print()

    print()
    print("Character table:")
    CFunc.show_table(LT)
    if argv.latex:
        print()
        print(r"\noindent Character table for image of $\beta$:")
        print()
        CFunc.latex_table(LT)

    m, n = B.shape
    items = []
    for row in B:
      for x in row:
        items.append(x)

    #if argv.latex:
    #    print("MatrixSpace(Z,%s,%s)(%s).hermite_form()" % (m, n, items))


def latex_table(table, rows, cols, upper=False):
    lines = []
    m, n = len(rows), len(cols)
    lines.append(r"\begin{array}{r|%s}"%('r'*n))
    lines.append(r"\times & %s \\" % (' & '.join(cols)))
    lines.append(r"\hline")
    for i in range(m):
        row = rows[i]
        if upper:
            line = [(table[row, col] if col>=row else " ") for col in cols]
        else:
            line = [table[row, col] for col in cols]
        line = r"%s & %s \\" % (row, ' & '.join(line))
        line = line.replace("*", '')
        lines.append(line)
    lines.append(r"\end{array}")
    s = '\n'.join(lines)
    s = s.replace("+", "\\thinplus ")
    return s


def is_hom(hom):
    for g in hom.keys():
      for h in hom.keys():
        k = g*h
        if k in hom:
            if hom[k] != hom[g]*hom[h]:
                return False
    return True


def close_hom(hom):
    hom = dict(hom)
    done = False
    gen = list(hom.keys())
    bdy = list(gen)
    while bdy:
        _bdy = []
        for g in gen:
          for h in bdy:
            k = g*h
            if k in hom:
                if hom[k] != hom[g]*hom[h]:
                    return None
            else:
                hom[k] = hom[g]*hom[h]
                _bdy.append(k)
        bdy = _bdy
    return hom


def find_all_homs(G, H, hom=None, remain=None):
    """
        Iterate through all homomorphisms from Group G to Group H.
        Yield each homomorphism as a dict : Perm -> Perm.
    """
    if hom is None:
        GI = G.identity
        hom = {GI : H.identity}
        #remain = [p for p in G if p != GI]

    for g in G:
        if g in hom:
            continue

        for h in H:
            hom[g] = h

            hom1 = close_hom(hom)
            if hom1 is None:
                pass

            elif len(hom1) == len(G):
                yield hom1

            else:
                for _hom in find_all_homs(G, H, hom1):
                    yield _hom



def test():

    items = list('abc')

    e = Perm.identity(items)

    g = Perm((1, 0, 2), items)
    assert g*g == e
    assert ~g == g

    h = Perm((1, 2, 0), items)
    assert h*h*h == e
    assert h*h != e
    assert h*h != g
    assert not (h*h == e)

    S3 = Group.generate([g, h])
    S3.check()
    assert len(S3) == 6

    assert len(list(S3.all_homs_iso(S3))) == 1

    #for g in S3:
    #    print g, g.fixed()

    items4 = list('abcd')
    g = Perm((1, 2, 3, 0), items4)
    h = Perm((1, 0, 2, 3), items4)
    S4 = Group.generate([g, h])
    assert len(S4)==24

    S4_22 = S4.choice(2)
    S4_22.check()

    # Pauli group
    X = Perm((3, 1, 2, 0), items4)
    Z = Perm((2, 3, 0, 1), items4)
    I = Perm((0, 1, 2, 3), items4)
    w = Perm((3, 2, 1, 0), items4)
    assert X*X==I
    assert Z*Z==I
    assert Z*X != X*Z
    assert Z*X*Z*X == w
    P1 = Group.generate([X, Z])
    assert len(P1)==8

    if 0:
        import numpy
        group = P1.square(P1)
        #print "orbits:", len(group.orbits())
        for orbit in group.orbits():
            #print orbit
            A = numpy.zeros((4, 4))
            for (i, j) in orbit:
                i, j = items4.index(i), items4.index(j)
                A[i, j] = 1
            #print A

    # Pauli group
    items = "+00 -00 +01 -01 +10 -10 +11 -11".split()
    #          0   1   2   3   4   5   6   7
    II = Perm((0, 1, 2, 3, 4, 5, 6, 7), items)
    XI = Perm((4, 5, 6, 7, 0, 1, 2, 3), items)
    IX = Perm((2, 3, 0, 1, 6, 7, 4, 5), items)
    ZI = Perm((0, 1, 2, 3, 5, 4, 7, 6), items)
    IZ = Perm((0, 1, 3, 2, 4, 5, 7, 6), items)
    w2 = Perm((1, 0, 3, 2, 5, 4, 7, 6), items)
    assert XI*XI==II
    assert ZI*ZI==II
    assert IX*IX==II
    assert IZ*IZ==II
    assert ZI*XI != XI*ZI
    assert ZI*XI == w2*XI*ZI
    assert ZI*XI*ZI*XI == w2
    assert IZ*IX != IX*IZ
    assert IZ*IX == w2*IX*IZ
    assert IZ*IX*IZ*IX == w2

    assert ZI*IX == IX*ZI

    P2 = Group.generate([XI, ZI, IX, IZ])
    assert len(P2)==32

    #homs = [f for f in find_all_homs(P1, P2, {I:II, w:w2})] # 152 solutions
    #assert len(homs)==152

    hom = {II:I} # 963 solutions
    hom = {II:I, w2:I} # 963 solutions
    hom = {II:I, w2:w, XI:X} # 0 solutions
    hom = {II:I, XI:X} # 64 solutions
    hom = {II:I, XI:X, IX:X} # 16 solutions
    hom = {II:I, w2:w} # 0 solutions
    count = 0
    for f in find_all_homs(P2, P1, hom):
        count += 1
    assert count==0, count

#    if 0:
#        #G1 = S4.choice(2)
#        G1 = S4.choice(3, 2, 1)
#        G2 = S4.choice(2)
#        group = G1.square(G2) # <----- FAIL
#        print "SxS:", len(group.items)
#        print "orbits:", len(group.orbits())
#        print
#
#        group = S4.choice(2).square()
#        print "##\n##"
#        print "SxS:", len(group.items)
#        print "orbits:", len(group.orbits())
#    
#        group = S4.choice(2, 1).square()
#        print "##\n#\n#"
#        print len(group)
#        print "SxS:", len(group.items)
#        print "orbits:", len(group.orbits())
#    
#        group = S4.choice(3, 2, 1).square()
#        print "#\n#\n#\n#"
#        print len(group)
#        print "SxS:", len(group.items)
#        print "orbits:", len(group.orbits())

    S4_211 = S4.choice(2, 1)
    assert len(S4_211.items)==12
    S4_211.check()

    assert S4.isomorphic(S4_22)

    assert len(list(S4_22.all_homs_iso(S4_22))) == 2
    assert len(list(S4.all_homs_iso(S4))) == 1

    #print len(list(S4.all_homs(list('ab'))))

    Z4 = Group.generate([Perm((1, 2, 3, 0), items4)])
    Z4.check()
    assert len(Z4)==4

    #print len(list(Z4.all_homs(list('abcd'))))

    assert len(list(Z4.all_homs(list('ab')))) == 2
    assert len(list(Z4.all_homs_iso(Z4))) == 4

    Z4_Z4 = Z4 + Z4
    assert len(Z4_Z4)==16
    Z4_Z4.check()
    assert len(Z4_Z4.orbits()) == 2

    group = Z4 * Z4
    group.check()
    #print len(group), len(group.orbits())

    #for g in Z4:
    #    print g, g.fixed()

    Z22 = Group.generate([Perm((1, 0, 3, 2), items4), Perm((2, 3, 0, 1), items4)])
    assert len(Z22)==4

    assert len(list(Z22.all_homs_iso(Z22))) == 4
    assert len(list(Z22.all_homs(list('ab')))) == 6
    #for group, func in Z22.all_homs(list('ab')):
    #    print func

    if 0:
        group = Z22.choice(2)
        #print "Z22.choice(2): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z22.choice(2, 1)
        #print "Z22.choice(2, 1): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z22.choice(3, 2, 1)
        #print "Z22.choice(3, 2, 1): orbits", [len(orbit) for orbit in group.orbits()]
    
        group = Z4.choice(2)
        #print "Z4.choice(2): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z4.choice(2, 1)
        #print "Z4.choice(2, 1): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z4.choice(3, 2, 1)
        #print "Z4.choice(3, 2, 1): orbits", [len(orbit) for orbit in group.orbits()]

#    print "fixed:",
#    for g in group:
#        #print g,
#        print len(g.fixed()),
#        #for item in g.fixed():
#        #    print "\t", item
#    print
#
#    print "orbits:"
##    for item in group.items:
##        print item, len(group.orbit(item)), ":", #repr(group.orbit(item)),
##        #for item in group.orbit(item):
##        #    print item,
##        print
#
#    print len(group.orbits())

    #print "OK"


def test_mul():

    n = argv.get("n", 5)
    items = range(n)

    for i in range(100):

        perms = []
        for count in range(3):
            perm = dict((i, i) for i in items)
            for _ in range(1):
                i = randint(0, n-1)
                j = randint(0, n-1)
                perm[i] = j
                perm[j] = i
            perms.append(Perm(perm, items))

        perms = mulclose_fast(perms)
        for g in perms:
          for h in perms:
            assert g*h in perms
        #write(".")

#test_mul()

def test_group():
    items = list('abcdefgh')
    mkperm = lambda seq : Perm([items.index(i) for i in seq], items)
    f = mkperm('ghabcdef')
    g = mkperm('dcbaefgh')

    G = Group.generate([f, g])
    #print(G)

def test_rowcol_perm():
    #rows, cols = 3, 7
    for rows in range(1, 30, 1):
      for cols in range(1, 30, 1):
        #items = [i + j*rows for i in range(rows) for j in range(cols)]
        items = [(i, j) for i in range(rows) for j in range(cols)]
        jtems = [(i, j) for j in range(cols) for i in range(rows)]
        perm = Perm(dict((items[idx], jtems[idx]) for idx in range(rows*cols)), items)
        #print(perm)
        r = (rows-1)//2 if rows%2 else rows//2
        c = (cols-1)//2 if cols%2 else cols//2
        lhs, rhs = perm.sign(), (-1)**(r*c)
        #print(lhs, rhs)
        assert lhs == rhs


if __name__ == "__main__":

    if argv.test:
        #test_action()
        #test()
        test_group()
        test_rowcol_perm()
        print("OK")

    elif argv.profile:
        import cProfile as profile
        profile.run("main()")
    else:
        name = argv.next() or "main"
        fn = eval(name)
        fn()

        print("OK\n")


