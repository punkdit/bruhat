#!/usr/bin/env python

import sys
import string
from random import randint, shuffle

from util import factorial, all_subsets, write
from argv import argv
import isomorph
from smap import SMap, tabulate


def mulclose(els, verbose=False, maxsize=None):
    els = set(els)
    changed = True
    while changed:
        if verbose:
            print "mulclose:", len(els)
        changed = False
        _els = list(els)
        for A in _els:
            for B in _els:
                C = A*B 
                if C not in els:
                    els.add(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
                    changed = True
    return els


def mulclose_fast(els, bdy=None):
    els = set(els)
    if bdy is None:
        bdy = [(i, j) for i in els for j in els]

    while bdy:
        _bdy = []
        for i, j in bdy:
            k = i*j
            if k not in els:
                _bdy.append((k, k))
                for kk in els:
                    _bdy.append((k, kk))
                    #_bdy.append((kk, k)) # i don't think we need this
                els.add(k)
        bdy = _bdy
    return els


def identity(items):
    return dict((i, i) for i in items)

class Set(object):
    def __init__(self, els):
        self.els = els

    def __str__(self):
        return str(self.els)


class Equ(object):
    "Equivalence class for an _equivalance relation"
    def __init__(self, item):
        self.items = [item]
        self.parent = None 

    def __str__(self):
        return "Equ(%s, top=%s)" % (
            self.items, self.top if self.parent else None)

    @property
    def top(self):
        top = self
        while top.parent:
            top = top.parent
        if top is not self:
            self.parent = top
        return top

#    def add(self, item):
#        assert item not in self.items
#        self.items.append(item)

    def merge(self, other):
        top = self
        while top.parent:
            top = top.parent
        while other.parent:
            other = other.parent
        if top is other:
            return
        other.parent = top
        top.items += other.items

    def eq(e1, e2):
        while e1.parent:
            e1 = e1.parent
        while e2.parent:
            e2 = e2.parent
        return e1 is e2


def quotient(items, relation=None):
    """ return a dict:item->items that sends each item 
        to the list of equivalant items 
        under the equivalance relation.
    """
    if relation is None:
        relation = lambda a,b : (a==b)
    equs = [Equ(item) for item in items]
    hom = {}
    n = len(items)
    for i in range(n):
        ei = equs[i]
        for j in range(i+1, n):
            ej = equs[j]
            if ei.eq(ej):
                continue
            if relation(items[i], items[j]):
                ei.merge(ej)
    for i in range(n):
        hom[items[i]] = list(equs[i].top.items)
    return hom

def quotient_rep(items, relation=None):
    """ return a dict:item->item that sends each item 
        to a representative equivalant item
        under the equivalance relation.
    """
    hom = quotient(items, relation)
    hom = dict((item, items[0]) for (item, items) in hom.items()) # confused ?
    return hom


assert quotient(range(7), (lambda i,j : (i-j)%3==0)) \
    == {0: [0, 3, 6], 1: [1, 4], 2: [2, 5], 3: [0, 3, 6], 4: [1, 4], 5: [2, 5], 6: [0, 3, 6]}

assert quotient_rep(range(7), (lambda i,j : (i-j)%3==0)) \
    == {0:0, 1:1, 2:2, 3:0, 4:1, 5:2, 6:0}


class Perm(object):

    """
    A permutation of a list of items.
    """
    def __init__(self, perm, items, word=''):
        #if isinstance(perm, list):
        #    perm = tuple(perm)
        if perm and isinstance(perm, (list, tuple)) and isinstance(perm[0], (int, long)):
            perm = list(items[i] for i in perm)
        if not isinstance(perm, dict):
            perm = dict((perm[i], items[i]) for i in range(len(perm)))
        self.perm = perm # map item -> item
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

    @classmethod
    def identity(cls, items, *args, **kw):
        n = len(items)
        perm = dict((item, item) for item in items)
        return Perm(perm, items, *args, **kw)

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

    def cycle_str(self):
        remain = set(self.set_items)
        s = []
#        print "__str__", self.perm, self.items
        while remain:
            item = iter(remain).next()
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

    def str(self): # HOTSPOT
        if self._str_cache:
            return self._str_cache
        perm = self.perm
        items = self.items
        s = []
        for i, item in enumerate(items):
            j = items.index(perm[item])
            s.append("%d:%d"%(i, j))
        s = "{%s}"%(', '.join(s))
        self._str_cache = s
        return s

    def __str__(self):
        return "Perm(%s)"%self.str()
    __repr__ = __str__

    def __hash__(self):
        return hash(str(self))

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
            item = iter(remain).next()
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
    assert n**m < 1e8, "too big"
    if m==0:
        yield {}
    elif n==0:
        return # no functions here
    elif m==1:
        for i in range(n):
            yield dict([(source[0], target[i])])
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
        self.perms = perms # ordered 
        self.items = list(items)
        self.set_items = set(items) # unordered
        self.set_perms = set(perms) # unordered
        for perm in perms:
            assert isinstance(perm, Perm), type(perm)
            assert perm.items == self.items, (perm.items, items)
        self._str = None # cache

    def str(self):
        if not self._str:
            ss = [perm.str() for perm in self.perms]
            ss.sort()
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

    @classmethod
    def generate(cls, perms, *args, **kw):
        items = perms[0].items
        perms = list(mulclose(perms, *args))
        return cls(perms, items, **kw)

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

    def __str__(self):
        return "Group(%s, %s)"%(self.perms, self.items)

    def __repr__(self):
        return "Group(%d, %d)"%(len(self.perms), len(self.items))

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
            item = iter(remain).next()
            orbit = set(g(item) for g in self.perms)
            #print "orbit:", orbit
            for item in orbit:
                remain.remove(item)
            orbits.append(orbit)
        return orbits

    def components(self):
        orbits = self.orbits()
        actions = [Group([perm.restrict(orbit) for perm in self.perms], orbit) for orbit in orbits]
        return actions

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
        "All subgroups, acting on the same items."
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

    def subgroups(self, verbose=False):
        I = self.identity
        items = self.items
        trivial = Group([I], items)
        cyclic = self.cyclic_subgroups()
        if verbose:
            print "Group.subgroups: cyclic:", len(cyclic)
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
            if verbose:
                print "subs:", len(subs)
                print "bdy:", len(bdy)
        return subs

    def left_cosets(self, H=None):
        cosets = set()
        if H is not None:
            Hs = [H]
        else:
            Hs = self.subgroups()
        for action in Hs:
            for g in self:
                coset = Group([g*h for h in action], self.items)
                cosets.add(coset)
        return list(cosets)

    def left_action(self, items):
        send_perms = {}
        perms = []
        for g in self:
            perm = {}
            for item in items:
                perm[item] = g*item
            perm = Perm(perm, items)
            perms.append(perm)
            send_perms[g] = perm
        f = quotient_rep(perms)
        send_perms = dict((k, f[v]) for (k, v) in send_perms.items())
        perms = list(f.values())
        hom = Action(self, send_perms, items)
        return hom

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

    @classmethod
    def product(cls, H, J):
        perms = list(set(H.perms+J.perms))
        return cls.generate(perms)

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


class Action(object):
    """
        A Group acting on a set.
        For each perm in the source Group G we map to a perm of items.
    """
    def __init__(self, G, send_perms, items, check=False):
        assert isinstance(G, Group)
        self.G = G
        assert isinstance(send_perms, dict)
        self.send_perms = dict(send_perms)
        self.items = list(items)
        if check:
            self.check()

    @property
    def src(self):
        return self.G

    # Equality on-the-nose:
    def __eq__(self, other):
        assert isinstance(other, Action)
        return (self.G==other.G and self.send_perms==other.send_perms)

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

    def pushout(self, other): # HOTSPOT
        assert self.src == other.src
        items = []
        perms = []
        for a1 in self.items:
          for a2 in other.items:
            items.append((a1, a2))
        send_perms = {}
        for g in self.src:
            perm = {}
            for a1, a2 in items:
                h1 = self.send_perms[g]
                h2 = other.send_perms[g]
                perm[(a1, a2)] = h1(a1), h2(a2)
            perm = Perm(perm, items)
            perms.append(perm)
            send_perms[g] = perm
        return Action(self.src, send_perms, items)

    def orbits(self):
        #print "orbits"
        #print self.perms
        #print self.items
        G = self.G
        send_perms = self.send_perms
        remain = set(self.items)
        orbits = []
        while remain:
            #print "remain:", remain
            item = iter(remain).next()
            orbit = set(send_perms[g](item) for g in G.perms)
            #print "orbit:", orbit
            for item in orbit:
                remain.remove(item)
            orbits.append(orbit)
        return orbits

    def components(self):
        src = self.src
        orbits = self.orbits()
        homs = []
        for orbit in orbits:
            send_perms = {}
            perms = []
            for perm in src:
                perm1 = self.send_perms[perm].restrict(orbit)
                send_perms[perm] = perm1
                perms.append(perm1)
            homs.append(Action(src, send_perms, orbit))
        return homs

    def get_graph(self):
        graph = isomorph.Graph()
        src = self.src
        send_perms = self.send_perms
        items = self.items
        n = len(items)
        # fibres are all the same:
        fibres = [graph.add("fibre") for item in items]
        for i, perm in enumerate(src):
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
        src = self.src
        send_perms = self.send_perms
        shape = []
        for perm in src:
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
        assert self.src == other.src
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
        for perm in self.src:
            perm1 = self.send_perms[perm]
            perm2 = other.send_perms[perm]
            for item1 in self.items:
                item2 = send_items[perm1[item1]]
                assert item2 == perm2[send_items[item1]]

    def isomorphic(self, other, check=False):
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

    def signature(self, Hs):
        return [len(self.fixed_points(H)) for H in Hs]

    def _refute_isomorphism(self, other, Hs):
        "return True if it is impossible to find an isomorphism"
        assert isinstance(other, Action)
        assert self.src == other.src
        n = len(self.items)
        if n != len(other.items):
            return True

        n = len(self.src)
        for perm in self.src:
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



r"""
Here we calculate the Burnside Ring corresponding to a particular G-set.

Finite kleinian geometry: 
(1) a set of "points" and a group G acting on the set.
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

    elif argv.B_2:
        items = range(8)
        gen = [
            Perm({0:0, 1:2, 2:1, 3:3, 4:6, 5:7, 6:4, 7:5}, items), 
            Perm({0:1, 1:0, 2:3, 3:2, 4:4, 5:5, 6:7, 7:6}, items)]
        perms = mulclose(gen)
        G = Group(perms, items)

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

    elif argv.projective:
        import geometry
        g = geometry.projective(3)
        items = range(14)
        perms = [f for f in g.get_symmetry()]
        perms = [Perm(f, items) for f in perms]
        G = Group(perms, items)

    else:
        return

    check = argv.check

    if 0:
        shapes = list(G.shape_item())
        print len(shapes), len(set(shapes))
        for shape in shapes:
            print shape
        return

    if 0:
        action = Action.identity(G)
        for iso in action.isomorphisms(action):
            print iso
    
        return

    print "G:", len(G)

    if argv.test_projective:
        test_projective(G)

    if argv.subgroups or argv.cyclic_subgroups:
        if argv.cyclic_subgroups:
            Hs = G.cyclic_subgroups()
        else:
            Hs = G.subgroups()
        print "subgroups:", len(Hs)
        Hs = list(Hs)
        Hs.sort(key = lambda H : len(H))
        for H in Hs:
            #if len(H)==len(G) or len(H)==1:
            #    continue
            print "subgroup order=%d:"%len(H)
            perms = [perm.perm for perm in H]
            print perms

    if argv.burnside:
        burnside(G)


def test_projective(G):

    Hs = G.cyclic_subgroups()

    print [len(H) for H in Hs]
    H2s = [H for H in Hs if len(H)==2]
    H3s = [H for H in Hs if len(H)==3]

    pairs = list(zip(H2s, H3s))
    print "pairs:", len(pairs)
    shuffle(pairs)
    for H2, H3 in pairs:
        H = Group.product(H2, H3)
        if len(H)!=6:
            write('.')
            continue
    else:
        assert 0

    print
    print H.is_abelian()
    for perm in H:
        print perm.cycle_str()


def burnside(G):

    # Find all conjugacy classes of subgroups
    Hs = G.subgroups()
    print "subgroups:", len(Hs)

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
                equs[H1].merge(equs[H2])

    # get equivalance classes
    equs = list(set(equ.top for equ in equs.values()))
    equs.sort(key = lambda equ : (-len(equ.items[0]), equ.items[0].str()))
    for equ in equs:
        print "equ:", [len(H) for H in equ.items]
    print "total:", len(equs)

    #return

    Hs = [equ.items[0] for equ in equs] # pick unique (up to conjugation)
    #Hs.sort(key = lambda H : (-len(H), H.str()))

    letters = list(string.uppercase + string.lowercase)
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
        print "%s subgroup size = %d, number of cosets = %d" %(hom.name, len(H), len(cosets))

    if 0:
        # We don't need to do this again: isomorphic homs all
        # come from conjugate subgroups.
        f = quotient_rep(homs, Action.isomorphic)
        homs = list(set(f.values())) # uniq
        print "homs:", len(homs)

    table = {}
    width = 0

    for i in range(len(homs)):
      for j in range(i, len(homs)):
        A = homs[i]
        B = homs[j]
        C = A.pushout(B)
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

    space = 2
    if argv.compact or 1:
        table = dict((k, v.replace("*", "")) for (k, v) in table.items())
        space = 1

    rows = cols = [hom.name for hom in homs]
    print
    print tabulate(table, rows, cols, space)
    print
    print latex_table(table, rows, cols)


r"""
\begin{array}{c|lcr}
n & \text{Left} & \text{Center} & \text{Right} \\
\hline
1 & 0.24 & 1 & 125 \\
2 & -1 & 189 & -8 \\
3 & -20 & 2000 & 1+10i
\end{array}
"""

def latex_table(table, rows, cols):
    lines = []
    m, n = len(rows), len(cols)
    lines.append(r"\begin{array}{r|%s}"%('r'*n))
    lines.append(r"* & %s \\" % (' & '.join(cols)))
    lines.append(r"\hline")
    for i in range(m):
        row = rows[i]
        line = r"%s & %s \\" % (row, ' & '.join(table[row, col] for col in cols))
        line = line.replace("*", '')
        lines.append(line)
    lines.append(r"\end{array}")
    s = '\n'.join(lines)
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
    for g in hom.keys():
      for h in hom.keys():
        k = g*h
        if k in hom:
            if hom[k] != hom[g]*hom[h]:
                return None
        else:
            hom[k] = hom[g]*hom[h]
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
        print "orbits:", len(group.orbits())
        for orbit in group.orbits():
            print orbit
            A = numpy.zeros((4, 4))
            for (i, j) in orbit:
                i, j = items4.index(i), items4.index(j)
                A[i, j] = 1
            print A

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

    if 0:
        #G1 = S4.choice(2)
        G1 = S4.choice(3, 2, 1)
        G2 = S4.choice(2)
        group = G1.square(G2) # <----- FAIL
        print "SxS:", len(group.items)
        print "orbits:", len(group.orbits())
        print

        group = S4.choice(2).square()
        print "##\n##"
        print "SxS:", len(group.items)
        print "orbits:", len(group.orbits())
    
        group = S4.choice(2, 1).square()
        print "##\n#\n#"
        print len(group)
        print "SxS:", len(group.items)
        print "orbits:", len(group.orbits())
    
        group = S4.choice(3, 2, 1).square()
        print "#\n#\n#\n#"
        print len(group)
        print "SxS:", len(group.items)
        print "orbits:", len(group.orbits())

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
        print "Z22.choice(2): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z22.choice(2, 1)
        print "Z22.choice(2, 1): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z22.choice(3, 2, 1)
        print "Z22.choice(3, 2, 1): orbits", [len(orbit) for orbit in group.orbits()]
    
        group = Z4.choice(2)
        print "Z4.choice(2): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z4.choice(2, 1)
        print "Z4.choice(2, 1): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z4.choice(3, 2, 1)
        print "Z4.choice(3, 2, 1): orbits", [len(orbit) for orbit in group.orbits()]

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

    print "OK"


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



if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")
    else:
        main()

    if argv.test:
        test_action()
        test()


