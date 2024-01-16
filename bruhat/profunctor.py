#!/usr/bin/env python3

"""
Here we build (very!) finite Category's and Functor's *concretely*,
as how Group/Action/Hom is done concretely in action.py

Previous work: biset.py, adjoint.py, combinatorial.py, action.py ...

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


invert_send = lambda send: dict((v,k) for (k,v) in send.items())


class Set(object):
    """
    These are the object's of our Category's and we don't use these for any
    other purpose (such as the Set of object's of a Category).
    """
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
        return "%s%s"%(self.__class__.__name__, self.set_items or "{}",)
    __repr__ = __str__

    def __contains__(self, item):
        return item in self.set_items

    def __iter__(self):
        return iter(self.items)

    def __len__(self):
        return len(self.items)

    # very important:
    # __eq__ is object identity !

    def __lt__(self, other):
        assert isinstance(other, Set)
        return id(self) < id(other)

    def issubset(self, other):
        assert isinstance(other, Set)
        return self.set_items.issubset(other.set_items)

    @property
    def i(self):
        send = {i:i for i in self.items}
        return Map(self, self, send)

    def include(tgt, src):
        send = {i:i for i in src.items}
        return Map(tgt, src, send)

    def all_maps(tgt, src):
        assert isinstance(src, Set)
        return {Map(tgt, src, m) for m in util.all_send(tgt.items, src.items)}

    def all_perms(self):
        items = self.items
        return {Map(self, self, zip(items, m)) for m in util.allperms(items)}

    def all_subsets(self):
        n = len(self.items)
        for idxs in util.all_subsets(n):
            items = [self.items[idx] for idx in idxs]
            yield Set(items)

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
    """
    A function on sets tgt<--src.
    These are the morphism's of our Category's and we don't use these for any
    higher purpose, to save on confusion and level-slips... we hope.
    """
    def __init__(self, tgt, src, send, check=CHECK):
        assert isinstance(tgt, Set), "Map tgt=%s not a Set"%(src,)
        assert isinstance(src, Set), "Map src=%s not a Set"%(src,)
        send = dict(send)
        self.key = tuple((i,send[i]) for i in src.items)
        self.data = tgt, src, send
        self.src = src
        self.tgt = tgt
        self.send = send
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

    def __contains__(self, i):
        return i in self.send

    def __getitem__(self, i):
        return self.send[i]

    def check(self):
        tgt, src, send = self.data
        for i in src:
            assert i in send
            assert send[i] in tgt, (tgt, src, send)

    def __mul__(self, other):
        assert self.src == other.tgt
        send = {i:self.send[other.send[i]] for i in other.src}
        return Map(self.tgt, other.src, send)

    def __invert__(self):
        tgt, src, send_items = self.data
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
        assert type(obs) in [list, set]
        obs = set(obs)
        names = {ob:ascii_letters[i] for (i,ob) in enumerate(obs)} # XXX more names...
        homs = dict(homs)
        self.obs = obs # [Set( ),...]
        #canonical_obs = list(obs)
        #canonical_obs.sort() # hmmm...
        #self.canonical_obs = canonical_obs
        self.homs = homs # {(tgt,src) : {Map,...,}}
        self.names = names
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
            assert type(hom[X,X]) is set, type(hom[X,X])
            assert X.i in hom[X,X]
            for Y in obs:
                assert (Y,X) in hom, list(hom.keys())
                Y_X = hom[Y,X]
                assert isinstance(Y_X, set), type(Y_X)
                for f in Y_X:
                    assert isinstance(f, Map)
                    assert f.src == X
                    assert f.tgt == Y
        for X in obs:
            for Y in obs:
                Y_X = hom[Y,X]
                for Z in obs:
                    Z_X = hom[Z,X]
                    Z_Y = hom[Z,Y]
                    for f in Y_X:
                      for g in Z_Y:
                        gf = g*f
                        assert gf in Z_X
                    
    def __str__(self):
        obs, homs, names = self.obs, self.homs, self.names
        homs = ', '.join([
            "%s<-%s-%s"%(names[t],len(ms),names[s])
            for (t,s),ms in homs.items() if len(ms)])
        obs = ', '.join(names[ob] for ob in obs)
        return "Category({%s}, {%s})"%(obs, homs)

    @cache
    def cayley(self):
        "the co-variant Cayley representation"
        send_obs = {x:[] for x in self.obs}
        for (t,s),hom in self.homs.items():
            send_obs[t] += list(hom)
        send_obs = {t:Set(send_obs[t]) for t in self.obs}
        obs = set(send_obs.values())
        send_homs = {}
        homs = {}
        for (t,s),hom in self.homs.items():
            send = {f : Map(send_obs[t], send_obs[s], {g:f*g for g in send_obs[s]})
                for f in hom}
            tgt = set(send.values())
            homs[send_obs[t],send_obs[s]] = tgt
            #send_homs[t,s] = Map(tgt, hom, send)
            send_homs[t,s] = send
        tgt = Category(obs, homs)
        #send_obs = Map(obs, self.obs, send_obs)
        F = Functor(tgt, self, send_obs, send_homs)
        return F

    @cache
    def cayley_op(self):
        "the contra-variant Cayley representation"
        send_obs = {x:[] for x in self.obs}
        for (t,s),hom in self.homs.items():
            send_obs[s] += list(hom)
        send_obs = {s:Set(send_obs[s]) for s in self.obs}
        obs = set(send_obs.values())
        send_homs = {}
        homs = {}
        for (t,s),hom in self.homs.items():
            send = {f : Map(send_obs[s], send_obs[t], {g:g*f for g in send_obs[t]})
                for f in hom}
            tgt = set(send.values())
            homs[send_obs[s],send_obs[t]] = tgt
            #send_homs[t,s] = Map(tgt, hom, send)
            send_homs[t,s] = send
        tgt = Category(obs, homs, op=self)
        assert tgt._op is self
        #send_obs = Map(obs, self.obs, send_obs)
        F = ContraFunctor(tgt, self, send_obs, send_homs)
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
        return Category(obs, homs)

    def include(tgt, src):
        send_obs = {x:x for x in src.obs}
        send_homs = {(s,t):{m:m for m in src.homs[s,t]} for (s,t) in src.homs.keys()}
        F = Functor(tgt, src, send_obs, send_homs)
        return F

    @property
    def i(self):
        return self.include(self)

    def __add__(self, other):
        return AddCategory(self, other)

    def __mul__(self, other): # the cartesian product
        return MulCategory(self, other)
    __matmul__ = __mul__ # the monoidal (tensor) product is the cartesian product

    @classmethod
    def full(cls, obs):
        return Category(obs, {(t,s):t.all_maps(s) for t in obs for s in obs})

    @classmethod
    def symmetric(cls, obs):
        "the symmetric groupoid"
        obs = [Set.promote(ob) for ob in obs]
        return Category(obs, {(t,s):(s.all_perms() if s==t else set())
            for s in obs for t in obs})

    @classmethod
    def poset(cls, ob):
        ob = Set.promote(ob)
        obs = list(ob.all_subsets())
        return Category(obs, {(t,s):({t.include(s)} if s.issubset(t) else set())
            for s in obs for t in obs})


class AddCategory(Category):
    "coproduct of Category's"

    cache = {}
    def __new__(cls, *args, **kw):
        if args in cls.cache:
            return cls.cache[args]
        lookup = {(i,ob):Set((x,i) for x in ob)
            for (i,arg) in enumerate(args) for ob in arg.obs}
        obs = set(lookup.values())
        #print("AddCategory")
        #print("obs =", obs)
        homs = {}
        for i,arg in enumerate(args):
            #print(arg, "i=", i)
            for (t,s),hom in arg.homs.items():
                #print(t, s)
                _hom = set()
                for m in hom:
                    #print("\t", m)
                    send = {(x,i):(m[x],i) for x in s}
                    m = Map(lookup[i,t], lookup[i,s], send)
                    #print("\t", m)
                    _hom.add(m)
                homs[lookup[i,t],lookup[i,s]] = _hom
        for s in obs:
          for t in obs:
            if (t,s) not in homs:
                homs[t,s] = set()
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
        #print("MulCategory")
        lookup = {obs:Set(util.cross([x.items for x in obs])) 
            for obs in util.cross([arg.obs for arg in args])}
        #print(len(lookup))
        obs = set(lookup.values())
        #print(obs)
        keys = list(lookup.keys())
        homs = {}
        for s in keys:
          for t in keys:
            #print(s,t)
            src = lookup[s]
            tgt = lookup[t]
            #print(src, tgt)
            hom = set()
            for ms in util.cross([arg.homs[t[i], s[i]] for i,arg in enumerate(args)]):
                #print('\t', [m.send for m in ms])
                send = {ks:tuple(m[k] for (m,k) in zip(ms,ks))
                    for ks in util.cross([m.send.keys() for m in ms])}
                #print('\t', send)
                hom.add(Map(tgt, src, send))
            homs[tgt, src] = hom
        #exit()
        self = object.__new__(cls)
        Category.__init__(self, obs, homs)
        cls.cache[args] = self
        return self

    def __init__(self, *args, **kw):
        pass


class Functor(object):
    def __init__(self, tgt, src, send_obs, send_homs, op=None, check=CHECK):
        assert isinstance(tgt, Category)
        assert isinstance(src, Category)
        assert isinstance(send_obs, dict)
        self.tgt = tgt # Category
        self.src = src # Category
        self.send_obs = send_obs # dict
        self.send_homs = send_homs # {(t,s):dict}
        self._op = op
        #obs = src.canonical_obs
        self.data = (tgt, src, send_obs, send_homs)
            #tuple((send_obs[ob],ob) for ob in obs), 
            #tuple((send_homs[s,t],s,t) for s in obs for t in obs))
        if check:
            self.check()

    def check(self):
        (tgt, src, send_obs, send_homs, ) = self.data 
        assert set(send_obs.keys()) == src.obs
        assert set(send_obs.values()).issubset(tgt.obs)
        for (t,s) in src.homs.keys():
            assert (t,s) in send_homs
            send = send_homs[t,s]
            assert isinstance(send, dict)
            assert set(send.keys()) == src.homs[t,s]
            assert set(send.values()).issubset(tgt.homs[send_obs[t],send_obs[s]])
        for x in src.obs:
          for y in src.obs:
            for z in src.obs:
              for f in src.homs[y,x]:
                for g in src.homs[z,y]:
                    #lhs,rhs = send_homs[z,x][g*f], send_homs[z,y][g]*send_homs[y,x][f]
                    #print(lhs, rhs)
                    assert send_homs[z,x][g*f] == send_homs[z,y][g]*send_homs[y,x][f]

    def dump(self):
        print(self.__class__.__name__)
        print(self.tgt)
        print(self.src)
        for ob in self.src.obs:
            print("\t%s<--%s"%(self(ob),ob)) 
        for t in self.src.obs:
          for s in self.src.obs:
            print('  ', t.items, "<---", s.items)
            for f in self.src.homs[t,s]:
                lhs, rhs = self(f).send, f.send
                print('    ', lhs, "<---", rhs)

    def __eq__(self, other):
        assert self.__class__ is other.__class__, "incomparable"
        assert self.tgt == other.tgt, "incomparable"
        assert self.src == other.src, "incomparable"
        return self.data == other.data

#    def __hash__(self):
#        return hash(self.data)

    def __str__(self):
        return "%s(\n\t%s\n\t<--\n\t%s)"%(self.__class__.__name__, self.tgt, self.src)

    def __mul__(self, other):
        "horizontal _composition of Functor's"
        assert self.__class__ is other.__class__
        assert self.src == other.tgt
        send_obs = {x:self.send_obs[other.send_obs[x]] for x in other.src.obs}
        send_homs = {}
        for (t,s),send in other.send_homs.items():
            t1, s1 = other.send_obs[t], other.send_obs[s]
            send1 = self.send_homs[t1, s1]
            #send_homs[t,s] = m1*m
            send_homs[t,s] = {m:send1[send[m]] for m in send}
        return self.__class__(self.tgt, other.src, send_obs, send_homs)
    __lshift__ = __mul__

    def __matmul__(self, other):
        "monoidal (cartesian) product of Functor's"
        if isinstance(other, Category):
            other = other.i
        assert isinstance(other, Functor)
        tgt = self.tgt @ other.tgt
        src = self.src @ other.src
        # TODO
        return Functor(tgt, src, send_obs, send_homs)

    def __invert__(self):
        (tgt, src, send_obs, send_homs, ) = self.data 
        send_homs = {(send_obs[t],send_obs[s]):invert_send(send) 
            for (t,s),send in send_homs.items()}
        send_obs = invert_send(send_obs)
        return self.__class__(src, tgt, send_obs, send_homs, )

    def __call__(self, item):
        if isinstance(item, Set):
            return self.send_obs[item]
        elif isinstance(item, Map):
            return self.send_homs[item.tgt, item.src][item]
        assert 0, item

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

    @property
    def i(self):
        components = {x:self(x).i for x in self.src.obs}
        return Nat(self, self, components)


class ContraFunctor(Functor):
    def __invert__(self):
        (tgt, src, send_obs, send_homs, ) = self.data 
        send_homs = {(send_obs[t],send_obs[s]):invert_send(send)
            for (s,t),send in send_homs.items()} # um.. seems to work
        send_obs = invert_send(send_obs)
        return self.__class__(src, tgt, send_obs, send_homs, )

    def check(self):
        (tgt, src, send_obs, send_homs, ) = self.data 
        #assert send_obs.src == src.obs
        #assert send_obs.tgt == tgt.obs
        assert set(send_obs.keys()) == src.obs
        assert set(send_obs.values()).issubset(tgt.obs)
        for (t,s) in src.homs.keys():
            assert (t,s) in send_homs
            #m = send_homs[t,s]
            #assert isinstance(m, Map)
            #assert m.src == src.homs[t,s]
            #assert m.tgt == tgt.homs[send_obs[s], send_obs[t]] # twisted
            send = send_homs[t,s]
            assert isinstance(send, dict)
            assert set(send.keys()) == src.homs[t,s]
            assert set(send.values()).issubset(tgt.homs[send_obs[s],send_obs[t]]) # twisted
        for x in src.obs:
          for y in src.obs:
            for z in src.obs:
              for f in src.homs[y,x]:
                for g in src.homs[z,y]:
                    # twisted
                    assert send_homs[z,x][g*f] == send_homs[y,x][f]*send_homs[z,y][g]


def _valid_functor_send(obs, send_obs, send_homs):
    for (t,s) in send_homs:
        send_ts = send_homs[t,s]
        for (t1,s1) in send_homs:
            if t!=s1:
                continue
            # t1 <--g-- s1==t <--f-- s
            send_t1s1 = send_homs[t1,s1]
            send_t1s = send_homs[t1,s]
            for f,f1 in send_ts.items():
              for g,g1 in send_t1s1.items():
                if send_t1s[g*f] != send_t1s1[g]*send_ts[f]:
                    return False
    return True


def all_functors(tgt, src): # warning: slow and stupid
    assert isinstance(tgt, Category)
    assert isinstance(src, Category)
    keys = list(src.homs.keys())
    for send_obs in util.all_send(tgt.obs, src.obs):
        #for send_homs in _all_send_homs(tgt, src, send_obs):
        #m_src = reduce(add, [list(hom) for hom in src.homs.values()])
        #for m in m_src:
            #print('\t', m)
        searchs = []
        for (t,s) in keys:
            _t, _s = send_obs[t], send_obs[s]
            search = util.all_send(tgt.homs[_t, _s], src.homs[t, s])
            searchs.append(list(search))
        if [] in searchs:
            continue
        #print(len(searchs))
        #for search in searchs:
            #print(search.__next__())
        for send_homs in util.cross(searchs):
            send_homs = {keys[i]:send for (i,send) in enumerate(send_homs)}
            #print('\t', send_homs)
            if _valid_functor_send(src.obs, send_obs, send_homs):
                yield Functor(tgt, src, send_obs, send_homs)


class Nat(object):
    """
    A natural transform of Functor's: tgt<===src is

    lhs <---- rhs
         tgt
          ^
          ||
          ||
          ||
         src
    lhs <---- rhs
    """
    def __init__(self, tgt, src, components, check=CHECK):
        assert isinstance(tgt, Functor)
        assert isinstance(src, Functor)
        assert (tgt.tgt,tgt.src) == (src.tgt,src.src)
        lhs, rhs = (tgt.tgt,tgt.src) 
        components = dict(components) # {X in src --> Map(tgt(X), src(X)) in lhs}
        self.tgt = tgt
        self.src = src
        self.components = components
        self.lhs = lhs
        self.rhs = rhs
        self.data = (tgt, src, components)
        if check:
            self.check()

    def check(self):
        (G, F, components) = self.data
        lhs, rhs = (F.tgt, F.src) 
        for x in rhs.obs:
            assert x in components
            m = components[x]
            assert m.tgt == G.send_obs[x]
            assert m.src == F.send_obs[x]
        for (y,x),hom in rhs.homs.items():
            for f in hom:
                Ff = F.send_homs[y,x][f]
                Gf = G.send_homs[y,x][f]
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

    def __str__(self):
        return "Nat(%s<--%s)"%(self.tgt, self.src)

    def __getitem__(self, x):
        return self.components[x]

    def __mul__(self, other):
        "vertical _composition"
        assert isinstance(other, Nat)
        assert other.tgt == self.src
        lhs, rhs = self.lhs, self.rhs
        components = {x : self.components[x] * other.components[x] for x in rhs.obs}
        return Nat(self.tgt, other.src, components)

    def __lshift__(self, other):
        "horizontal _composition"
        assert isinstance(other, Nat)
        assert other.lhs == self.rhs
        # horizontal _composition is strictly associative
        components = {}
        for x in other.src.src.obs:
            a = self[other.src(x)]
            b = self.tgt(other[x])
            m = b*a
            c = self.src(other[x])
            d = self[other.tgt(x)]
            assert m == d*c # check interchange law
            components[x] = m
        tgt = self.tgt << other.tgt
        src = self.src << other.src
        return Nat(tgt, src, components)


def _valid_nat(G, F, components):
    lhs, rhs = (F.tgt, F.src) 
    for (y,x),hom in rhs.homs.items():
        for f in hom:
            Ff = F.send_homs[y,x][f]
            Gf = G.send_homs[y,x][f]
            mx = components[x]
            my = components[y]
            if Gf*mx != my*Ff:
                return False
    return True

def all_nats(G, F): # warning: slow and stupid
    assert isinstance(G, Functor)
    assert isinstance(F, Functor)
    assert G.tgt == F.tgt
    assert G.src == F.src
    tgt, src = G.tgt, G.src
    obs = list(G.src.obs)
    n = len(obs)
    homs = []
    for x in obs:
        hom = tgt.homs[G(x), F(x)]
        if not hom:
            return
        homs.append(list(hom))
    for send in util.cross(homs):
        components = {obs[i]:send[i] for i in range(n)}
        if _valid_nat(G, F, components):
            yield Nat(G, F, components)


def test_nats(found):
    def search(u):
        for v in found:
            if v.tgt == u.tgt and v.src == u.src and u==v:
                return True
        return False
    for u in found:
      for v in found:
        if u.tgt is v.src:
            vu = v*u
            assert search(vu)
        vu = v<<u
        assert search(vu)


def test_category():
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
    DCi = D.include(C)
    E = Category.generate([f, g, h])
    EDi = E.include(D)
    ECi = E.include(C)

    #print("C =", C)
    #print("D =", D)
    assert len(list(all_functors(C,C))) == 3
    assert len(list(all_functors(C,D))) == 3
    assert len(list(all_functors(D,D))) == 9 # why 9?

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

    ECi = E.include(C)

    F = E.i
    print(F)
    Fi = F.i
    print(Fi)
    assert Fi*Fi == Fi
    assert Fi<<Fi == Fi

    #return

    A = Set(list(range(3)))
    B = Set('xy')
    M = Category([A], {(A,A):A.all_maps(A)})

    Mc = M.cayley()
    M.cayley_op()


    C = Category.full([A, B])
    print(C)
    C.cayley()
    print(C)

    G = Category([A], {(A,A):A.all_perms()})
    print(G)
    print(G+G)
    print(G+C)

    C = Category.full([Set('ab'), Set('xy')])
    D = Category.symmetric([Set('ab'), Set('xyz')])
    print("C =", C)
    print("D =", D)
    print("C*C =", C*C)
    print("C*D =", C*D)

    C = Category.poset('ab')
    assert len(C.obs) == 2**2
    sends = []
    M = []
    for f in all_functors(C,C):
        send = f.send_obs
        assert send not in sends
        sends.append(send)
        M.append(f)
    assert len(sends) == 36
    found = []
    for f in M:
      for g in M:
        for nat in all_nats(g, f):
            found.append(nat)
    assert len(found) == 400, len(found)
    #test_nats(found) # takes 40s !

    C = Category.poset('abc')
    assert len(C.obs) == 2**3

    # S3 has 6 automorphisms, 3 morphisms onto S3/A3 and 1 trivial == 10 total
    S3 = Category.symmetric(['abc'])
    M = list(all_functors(S3,S3))
    assert len(M) == 10 
    for f in M:
      for g in M:
        assert f*g in M
    found = []
    for f in M:
      for g in M:
        for nat in all_nats(g, f):
            found.append(nat)
    assert len(found) == 60
    test_nats(found)

    f = S3.i
    assert len(list(all_nats(f, f))) == 1
    f.dump()


def homotopy_pullback(G, H, K):
    from bruhat.action import Group, Action, Hom
    assert isinstance(G, Group)
    assert isinstance(H, Group)
    assert isinstance(K, Group)
    assert G.is_subgroup(H)
    assert G.is_subgroup(K)
    #lookup = {g:Set(G.items) for g in G}
    #obs = list(lookup.values())
    fans = {g:[] for g in G}
    _homs = {} # abstract homs
    for g1 in G:
      for g0 in G:
        # g0 ---> g1
        hom = []
        for h in H:
          for k in K:
            if k*g0 == g1*h:
                hom.append((h,k))
                fans[g1].append((h,k))
        _homs[g1,g0] = hom
        #print(len(hom), end=" ")
      #print()
    #print()
    # now we build the Cayley-Yoneda representation
    fans = {g:Set(fan) for (g,fan) in fans.items()}
    obs = list(fans.values())
    homs = {}
    for t in G:
      for s in G:
        tgt = fans[t]
        src = fans[s]
        hom = set()
        for (h,k) in _homs[t,s]:
            send = {(h0,k0):(h*h0,k*k0) for (h0,k0) in src}
            m = Map(tgt, src, send)
            hom.add(m)
        homs[fans[t],fans[s]] = hom
    cat = Category(obs, homs)
    return cat


def test_group():
    from bruhat.action import Group, Action, Hom
    n = 3
    G = Group.symmetric(n)
    X = Set(G.items)
    CG = Category([X], {(X,X):{Map(X,X,g.perm) for g in G}})
    acts = []
    Hs = list(G.subgroups())
    Hs.sort(key = len)
    print([len(H) for H in Hs])
    lookup = {}
    for H in Hs:
        act = G.action_subgroup(H)
        acts.append(act)
        CH = Category([X], {(X,X):{Map(X,X,g.perm) for g in H}})
        lookup[H] = CH

    pairs = {}
    for H in Hs:
      for K in Hs:
        if H.is_subgroup(K):
            i = lookup[H].include(lookup[K])
        #i = CG.include(CH)
        P = homotopy_pullback(G, H, K)

    for a in acts:
      for b in acts:
        homs = list(a.get_homs(b))
        print("%2s"%(len(homs) or '.'), end=' ')
      print()


def test_homotopy_pullback():
    from bruhat.action import Group, Action, Hom
    n = 4
    G = Group.symmetric(n)
    H = Group([g for g in G if g[0]==0], G.items) # S_3
    assert len(H) == 6
    assert G.is_subgroup(H)
    K = Group([g for g in G if g[0] in [0,1] and g[1] in [0,1]], G.items) # S_2xS_2
    assert len(K) == 4
    assert G.is_subgroup(K)

    C = homotopy_pullback(G, H, K)
    print(C)


def test():
    #test_category()
    #test_group()
    test_homotopy_pullback()


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






