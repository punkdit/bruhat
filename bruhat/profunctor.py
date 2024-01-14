#!/usr/bin/env python3

"""

build very finite Category's and Functor's 
concretely, as how Group's are done concretely in action.py

previous work: biset.py, adjoint.py, etc. etc.

"""


from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)
from operator import add, mul
from string import ascii_letters

from bruhat import util
from bruhat.argv import argv

CHECK = argv.get("check", True)


def mulclose(gen, verbose=False, maxsize=None):
    els = set(gen)
    bdy = list(els)
    while bdy:
        if verbose:
            print(len(els), end=" ", flush=True)
        _bdy = []
        for A in gen:
            for B in bdy:
                if B.tgt != A.src:
                    continue
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




class Set(object):
    def __init__(self, items):
        items = list(items)
        #items.sort()
        self.items = items
        self.set_items = set(items)

    @classmethod
    def promote(cls, items):
        if isinstance(items, Set):
            return items
        return Set(items)

    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__, self.items,)
    __repr__ = __str__

    def __contains__(self, item):
        return item in self.set_items

    def __iter__(self):
        return iter(self.items)

    def __len__(self):
        return len(self.items)

    # __eq__ is '=='

    @property
    def i(self):
        send_items = {i:i for i in self.items}
        return Map(self, self, send_items)

    def include_into(self, other):
        tgt = other
        src = self
        send_items = {i:i for i in self.items}
        return Map(tgt, src, send_items)

    def all_maps(tgt, src):
        assert isinstance(src, Set)
        ms = [Map(tgt, src, m) for m in util.all_functions(src.items, tgt.items)]
        return Set(ms)

    def all_perms(self):
        ms = [Map(self, self, enumerate(m)) for m in util.allperms(self.items)]
        return Set(ms)

    def __add__(self, other):
        return AddSet(self, other)

    def __mul__(self, other):
        return MulSet(self, other)


# see: bruhat.lin for similar setup
class AddSet(Set):
    "disjoint union of Set's"

    cache = {} # XXX use https://docs.python.org/3/library/weakref.html
    def __new__(cls, *args, **kw):
        if args in cls.cache:
            return cls.cache[args]
        items = []
        for i,arg in enumerate(args):
            assert isinstance(arg, Set)
            items += [(x, i) for x in arg]
        self = object.__new__(cls)
        Set.__init__(self, items)
        cls.cache[args] = self
        return self

    def __init__(self, *args, **kw):
        pass


class MulSet(Set):
    "product of Set's"

    cache = {} # XXX use https://docs.python.org/3/library/weakref.html
    #@cache # fails
    def __new__(cls, *args, **kw):
        if args in cls.cache:
            return cls.cache[args]
        items = list(util.cross([arg.items for arg in args]))
        self = object.__new__(cls)
        Set.__init__(self, items)
        cls.cache[args] = self
        return self

    def __init__(self, *args, **kw):
        pass


# we use order (tgt, src) because __mul__ is right-to-left

class Map(object):
    "A function on sets tgt<--src"
    def __init__(self, tgt, src, send_items, check=CHECK):
        assert isinstance(tgt, Set)
        assert isinstance(src, Set)
        send_items = dict(send_items)
        self.key = tuple((i,send_items[i]) for i in src.items)
        self.data = src, tgt, send_items
        self.src = src
        self.tgt = tgt
        self.send_items = send_items
        if check:
            self.check()

    def __str__(self):
        return "Map(%s<--%s, %s)"%(self.tgt, self.src, self.key,)
    __repr__ = __str__

    def __eq__(self, other):
        assert self.src == other.src, "incomparable"
        assert self.tgt == other.tgt, "incomparable"
        return self.data == other.data

    def __hash__(self):
        return hash((self.src, self.tgt, self.key))

    def __getitem__(self, i):
        return self.send_items[i]

    def check(self):
        src, tgt, send_items = self.data
        for i in src:
            assert i in send_items
            assert send_items[i] in tgt

    def __mul__(self, other):
        assert self.src == other.tgt
        send_items = {i:self.send_items[other.send_items[i]] for i in other.src}
        return Map(self.tgt, other.src, send_items)

    def __invert__(self):
        src, tgt, send_items = self.data
        assert len(src) == len(tgt)
        send = {k:v for (v,k) in send_items.items()}
        assert len(send) == len(send_items)
        return Map(src, tgt, send)


class Category(object):
    """
    These are _concrete Category's: objects are (finite!) Set's
    and morphisms are Map's of Set's.
    """
    def __init__(self, obs, homs, op=None, check=CHECK):
        obs = Set.promote(obs)
        homs = dict(homs)
        self.obs = obs # Set({Set( ),...})
        self.homs = homs # {(tgt,src) : Set(of Map's)}
        self.data = obs, homs
        self._op = op
        if check:
            self.check()

    def check(self):
        obs, hom = self.data
        #print("Category.check")
        #print("\t", obs)
        for X in obs:
            assert isinstance(X, Set), "%s is not a Set"%(X,)
            assert (X,X) in hom
            assert X.i in hom[X,X]
            for Y in obs:
                assert (Y,X) in hom, list(hom.keys())
                Y_X = hom[Y,X]
                assert isinstance(Y_X, Set), type(Y_X)
                for f in Y_X:
                    assert isinstance(f, Map)
                    assert f.src == X
                    assert f.tgt == Y
                    
    def __str__(self):
        obs, homs = self.obs, self.homs
        names = {ob:ascii_letters[i] for (i,ob) in enumerate(obs)}
        homs = ', '.join(["%s<-%s-%s"%(names[t],len(ms),names[s])
            for (t,s),ms in homs.items()])
        obs = ', '.join(names[ob] for ob in obs)
        return "Category({%s}, {%s})"%(obs, homs)

    @cache
    def cayley(self):
        "the co-variant Cayley representation"
        map_obs = {x:[] for x in self.obs}
        for (t,s),hom in self.homs.items():
            map_obs[t] += list(hom)
        map_obs = {t:Set(map_obs[t]) for t in self.obs}
        obs = Set(map_obs.values())
        map_homs = {}
        homs = {}
        for (t,s),hom in self.homs.items():
            send = {f : Map(map_obs[t], map_obs[s], {g:f*g for g in map_obs[s]})
                for f in hom}
            tgt = Set(send.values())
            homs[map_obs[t],map_obs[s]] = tgt
            map_homs[t,s] = Map(tgt, hom, send)
        tgt = Category(obs, homs)
        map_obs = Map(obs, self.obs, map_obs)
        F = Functor(tgt, self, map_obs, map_homs)
        return F

    @cache
    def cayley_op(self):
        "the contra-variant Cayley representation"
        map_obs = {x:[] for x in self.obs}
        for (t,s),hom in self.homs.items():
            map_obs[s] += list(hom)
        map_obs = {s:Set(map_obs[s]) for s in self.obs}
        obs = Set(map_obs.values())
        map_homs = {}
        homs = {}
        for (t,s),hom in self.homs.items():
            send = {f : Map(map_obs[s], map_obs[t], {g:g*f for g in map_obs[t]})
                for f in hom}
            tgt = Set(send.values())
            homs[map_obs[s],map_obs[t]] = tgt
            map_homs[t,s] = Map(tgt, hom, send)
        tgt = Category(obs, homs, op=self)
        assert tgt._op is self
        map_obs = Map(obs, self.obs, map_obs)
        F = ContraFunctor(tgt, self, map_obs, map_homs)
        return F

    @property
    def op(self):
        if self._op is not None:
            return self._op
        cat = self.cayley_op().tgt
        self._op = cat
        return cat

    @classmethod
    def generate(cls, gen):
        fs = mulclose(gen)
        obs = set()
        for f in fs:
            obs.add(f.src)
            obs.add(f.tgt)
        # should we use a set() or a Set() ??
        homs = {(X,Y):set() for X in obs for Y in obs}
        for X in obs:
            homs[X, X].add(X.i) # identity
        for f in fs:
            key = f.tgt, f.src
            homs[key].add(f)
        homs = {key:Set(value) for (key,value) in homs.items()}
        return Category(obs, homs)

    @property
    def i(self):
        map_obs = self.obs.i
        homs = self.homs
        map_homs = {(t,s):homs[t,s].i for (t,s) in homs.keys()}
        F = Functor(self, self, map_obs, map_homs)
        return F

    def include_into(self, other):
        tgt = other
        src = self
        map_obs = self.obs.include_into(other.obs)
        map_homs = {(s,t):self.homs[s,t].include_into(other.homs[s,t]) 
            for (s,t) in self.homs.keys()}
        F = Functor(tgt, src, map_obs, map_homs)
        return F

    def __add__(self, other):
        return AddCategory(self, other)

    def __mul__(self, other):
        return MulCategory(self, other)



class AddCategory(Category):
    "coproduct of Category's"

    cache = {}
    def __new__(cls, *args, **kw):
        if args in cls.cache:
            return cls.cache[args]
        obs = AddSet(*[arg.obs for arg in args])
        print("AddCategory")
        print(obs)
        homs = {}
        for arg in args:
            print(arg)
            for (t,s),m in arg.homs.items():
                print(m)
        self = object.__new__(cls)
        Category.__init__(self, obs, homs)
        cls.cache[args] = self
        return self

    def __init__(self, *args, **kw):
        pass


class MulCategory(Category):
    "product of Category's"

    cache = {}
    def __new__(cls, *args, **kw):
        if args in cls.cache:
            return cls.cache[args]
        self = object.__new__(cls)
        Category.__init__(self, obs, homs)
        cls.cache[args] = self
        return self

    def __init__(self, *args, **kw):
        pass


class Functor(object):
    def __init__(self, tgt, src, map_obs, map_homs, op=None, check=CHECK):
        assert isinstance(tgt, Category)
        assert isinstance(src, Category)
        assert isinstance(map_obs, Map)
        assert map_obs.src == src.obs
        assert map_obs.tgt == tgt.obs
        self.tgt = tgt # Category
        self.src = src # Category
        self.map_obs = map_obs # Map
        self.map_homs = map_homs # {(t,s):Map}
        self._op = op
        self.data = (tgt, src, map_obs, map_homs,)
        if check:
            self.check()

    def check(self):
        (tgt, src, map_obs, map_homs, ) = self.data 
        assert map_obs.src == src.obs
        assert map_obs.tgt == tgt.obs
        for (t,s) in src.homs.keys():
            assert (t,s) in map_homs
            m = map_homs[t,s]
            assert isinstance(m, Map)
            assert m.src == src.homs[t,s]
            assert m.tgt == tgt.homs[map_obs[t], map_obs[s]]
        for x in src.obs:
          for y in src.obs:
            for z in src.obs:
              for f in src.homs[y,x]:
                for g in src.homs[z,y]:
                    assert map_homs[z,x][g*f] == map_homs[z,y][g]*map_homs[y,x][f]

    def __eq__(self, other):
        assert self.__class__ is other.__class__, "incomparable"
        assert self.tgt == other.tgt, "incomparable"
        assert self.src == other.src, "incomparable"
        return self.data == other.data

    def __hash__(self):
        return hash(self.data)

    def __str__(self):
        return "%s(\n\t%s\n\t<--\n\t%s)"%(self.__class__.__name__, self.tgt, self.src)

    def __mul__(self, other):
        assert self.__class__ is other.__class__
        assert self.src == other.tgt
        map_obs = Map(self.tgt.obs, other.src.obs,
            {x:self.map_obs[other.map_obs[x]] for x in other.src.obs})
        map_homs = {}
        for (t,s),m in other.map_homs.items():
            t1, s1 = other.map_obs[t], other.map_obs[s]
            m1 = self.map_homs[t1, s1]
            map_homs[t,s] = m1*m
        return self.__class__(self.tgt, other.src, map_obs, map_homs)

    def __invert__(self):
        (tgt, src, map_obs, map_homs, ) = self.data 
        map_homs = {(map_obs[t],map_obs[s]):~m for (t,s),m in map_homs.items()}
        map_obs = ~map_obs
        return self.__class__(src, tgt, map_obs, map_homs, )

    def __call__(self, item):
        return self.map_obs[item]

    @property
    def op(self):
        if self._op is not None:
            return self._op
        D, C = self.tgt, self.src
        Dop_D = D.cayley_op()
        Cop_C = C.cayley_op()
        #F = Dop_D * self * ~Cop_C # FAIL
        assert 0, "TODO"
        self._op = F
        return F


class ContraFunctor(Functor):
    def __invert__(self):
        (tgt, src, map_obs, map_homs, ) = self.data 
        map_homs = {(map_obs[t],map_obs[s]):~m for (s,t),m in map_homs.items()} # um.. seems to work
        map_obs = ~map_obs
        return self.__class__(src, tgt, map_obs, map_homs, )

    def check(self):
        (tgt, src, map_obs, map_homs, ) = self.data 
        assert map_obs.src == src.obs
        assert map_obs.tgt == tgt.obs
        for (t,s) in src.homs.keys():
            assert (t,s) in map_homs
            m = map_homs[t,s]
            assert isinstance(m, Map)
            assert m.src == src.homs[t,s]
            assert m.tgt == tgt.homs[map_obs[s], map_obs[t]] # twisted
        for x in src.obs:
          for y in src.obs:
            for z in src.obs:
              for f in src.homs[y,x]:
                for g in src.homs[z,y]:
                    # twisted
                    assert map_homs[z,x][g*f] == map_homs[y,x][f]*map_homs[z,y][g]


class Nat(object):
    "A natural transform of Functor's"
    def __init__(self, tgt, src, components, check=CHECK):
        assert isinstance(tgt, Functor)
        assert isinstance(src, Functor)
        assert (tgt.tgt,tgt.src) == (src.tgt,src.src)
        D, C = (tgt.tgt,tgt.src) 
        components = dict(components) # {X in src --> Map(tgt(X), src(X)) in D}
        self.tgt = tgt
        self.src = src
        self.components = components
        self.D = D
        self.C = C
        self.data = (tgt, src, components)
        if check:
            self.check()

    def check(self):
        (G, F, components) = self.data
        D, C = (F.tgt, F.src) 
        for x in C.obs:
            assert x in components
            m = components[x]
            assert m.tgt == G.map_obs[x]
            assert m.src == F.map_obs[x]
        for (y,x),f in C.items():
            Ff = F.map_homs[y,x][f]
            Gf = G.map_homs[y,x][f]
            mx = components[x]
            my = components[y]
            assert Gf*mx == my*Ff

    def __eq__(self, other):
        assert self.tgt == other.tgt, "incomparable"
        assert self.src == other.src, "incomparable"
        return self.components == other.components

    @cache
    def __hash__(self):
        components = self.components
        key = tuple((components[x],x) for x in self.src)
        return hash(key)

    def __mul__(self, other):
        "vertical _composition"
        assert isinstance(other, Nat)
        assert other.tgt == self.src
        D, C = self.D, self.C
        components = {x : self.components[x] * other.components[x] for x in C.obs}
        return Nat(self.tgt, other.src, components)

    def __lshift__(self, other):
        "horizontal _composition"
        assert isinstance(other, Nat)
        assert other.D == self.C
        # horizontal _composition is strictly associative
        TODO


def test():
    X = Set([0,1,2,3])
    Y = Set([0,1])
    assert Y != Set([0,1]) # once built, it is never the same

    assert len(X+Y) == 4+2
    assert len(X*Y) == 4*2
    assert len(Y*Y*Y) == 2**3
    assert Y*Y*Y is Y*Y*Y

    f = Map(Y, X, {0:0, 1:0, 2:1, 3:0})
    g = Map(Y, Y, {0:1, 1:0})
    assert g == Map(Y, Y, {0:1, 1:0})
    assert ~g == g
    h = Map(X, Y, {0:0, 1:2})
    assert X==X
    assert X!=Y
    assert f==f
    assert g != f*h
    assert f*X.i == f
    assert g*g == Y.i
    gf = g*f
    assert g*g*f == f

    assert hash(g*g) == hash(Y.i)

    fs = mulclose([f, g])
    assert set(fs) == set([f, g, g*f, g*g])

    C = Category.generate([f])
    assert len(C.obs) == 2
    D = Category.generate([f, g])
    DCi = C.include_into(D)
    E = Category.generate([f, g, h])
    EDi = D.include_into(E)
    ECi = C.include_into(E)

    Ci = C.i
    Di = D.i
    assert Di*DCi == DCi
    assert DCi*Ci == DCi
    assert ECi == EDi*DCi

    F = C.cayley()
    _C = F.tgt
    assert ~F*F == C.i
    assert F*~F == _C.i

    Fop = C.cayley_op()
    assert Fop is C.cayley_op()
    Cop = Fop.tgt
    assert Cop.op is C
    ~Fop

    Eop = E.op
    assert Eop.op is E

    ECi = C.include_into(E)


    A = Set(list(range(3)))
    B = Set('xy')
    M = Category([A], {(A,A):A.all_maps(A)})

    Mc = M.cayley()
    M.cayley_op()

    obs = [A, B]
    C = Category(obs, {(t,s):t.all_maps(s) for t in obs for s in obs})
    C.cayley()

    G = Category([A], {(A,A):A.all_perms()})
    print(G)
    print(G+G)
    print(G*G)



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

    print("OK: finished in %.3f seconds"%(time() - start_time))
    print()






