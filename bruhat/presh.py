#!/usr/bin/env python3


"""
construct _presheaves on a finite site _category...
"""


class Cell(object):
    "abstract Cell in abstract category"
    def __init__(self, cat, desc=None, dim=0, src=None, tgt=None):
        assert isinstance(cat, ACat)
        try:
            hash(desc)
        except:
            raise TypeError("need hashable desc")
            
        self.cat = cat
        self.desc = desc or "%s(?)"%self.__class__.__name__
        self.dim = dim
        #self.src = src or self # ?
        #self.tgt = tgt or self # ?
        assert src is None or isinstance(src, Cell) and src.dim==dim-1
        assert tgt is None or isinstance(tgt, Cell) and tgt.dim==dim-1
        self.src = src
        self.tgt = tgt
        self.key = (self.__class__, desc, dim, src, tgt)

    #def __hash__(self):
    #    return hash((self.__class__, self.desc, self.dim, self.src, self.tgt))

    def __eq__(self, other):
        assert self.cat is other.cat
        assert self.dim == other.dim
        return self.key == other.key

    def __lt__(self, other):
        assert self.cat is other.cat
        assert self.dim == other.dim
        return self.key < other.key

    def __hash__(self):
        return hash(self.key) # ??

    def __str__(self):
        return self.desc

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__, self.desc)


class Object(Cell):
    "Object (0-cell) in abstract category"
    def __init__(self, cat, desc=None):
        Cell.__init__(self, cat, desc, 0)

    @classmethod
    def promote(cls, cat, desc):
        assert desc is not None
        if isinstance(desc, cls):
            return desc
        ob = cat.o_lookup.get(desc)
        if ob is None:
            ob = cls(cat, desc)
        return ob


class Morphism(Cell):
    "Morphism (1-cell) in abstract category"
    def __init__(self, cat, desc, src, tgt):
        src = Object.promote(cat, src)
        tgt = Object.promote(cat, tgt)
        Cell.__init__(self, cat, desc, 1, src, tgt)

    @classmethod
    def promote(cls, cat, f):
        assert f is not None
        if isinstance(f, cls):
            return f
        morph = cat.m_lookup.get(f)
        if morph is None:
            morph = cls(cat, *f)
        return morph

    def __mul__(g, f):
        "g.__mul__(f) is g*f : first do f, then f"
        if f.tgt != g.src:
            raise TypeError
        h = g.cat.comp[g, f]
        return h


#class Hom(object):
#    def __init__(self, cat, src, tgt, items=[]):
#        self.cat = cat
#        self.src = src
#        self.tgt = tgt
#        self.items = list(items)
#
#    def __str__(self):
#        return "Hom(%s)"%(self.items,)
#    __repr__ = __str__

class ACat(object):
    """
        An abstract category.
    """
    def __init__(self, obs=[], morphs=[], comp={}):
        self.o_lookup = {}
        self.m_lookup = {}

        assert len(set(obs)) == len(obs)
        assert len(set(morphs)) == len(morphs)

        o_lookup = dict((ob, Object.promote(self, ob)) for ob in obs)
        self.obs = set(o_lookup.values())
        self.o_lookup = o_lookup

        self.m_lookup = dict((f, Morphism.promote(self, f)) for f in morphs)
        self.morphs = set(self.m_lookup.values())

        self.ident = {}
        for ob in self.obs:
            desc = "1_"+ob.desc
            i = self.m_lookup.get(desc)
            if i is None:
                i = Morphism(self, desc, ob, ob)
            self.ident[ob] = i
            self.m_lookup[desc] = i
            self.morphs.add(i)

        self.comp = dict(((self.m_lookup[f], self.m_lookup[g]), self.m_lookup[h]) 
            for ((f, g), h) in comp.items())  # (f, g) -> h, ...

        for f in self.morphs:
            src = self.ident[f.src]
            tgt = self.ident[f.tgt]
            self.comp[f, src] = f
            self.comp[tgt, f] = f
        self.build()
        self.check()

    def __str__(self):
        return "ACat(%s, %s)"%(
            [ob.desc for ob in self.obs],
            [m.desc for m in self.morphs],
        )
    __repr__ = __str__

    def dump(self):
        for ob in self.obs:
            F = self.hom_into(ob)
            print("hom(_,%s):"%ob)
            print(F)
        print()

    def build(self):
        obs = self.obs
        morphs = self.morphs
        ident = self.ident
        comp = self.comp

    def check(self):
        obs = self.obs
        morphs = self.morphs
        ident = self.ident
        comp = self.comp

        assert len(set(obs)) == len(obs) # unique
        assert len(ident) == len(obs)

        for ob in obs:
            assert ob in ident
            i = ident[ob]
            f = ident[ob]
            assert f.src == ob
            assert f.tgt == ob

        for f in morphs:
            assert comp[f, ident[f.src]] == f
            assert comp[ident[f.tgt], f] == f

            for g in morphs:
                if f.tgt != g.src:
                    continue
                gf = comp[g, f] # algebraic order!
                assert gf.src == f.src
                assert gf.tgt == g.tgt
                for h in morphs:
                    if g.tgt != h.src:
                        continue
                    lhs = comp[h, gf]
                    rhs = comp[comp[h, g], f]
                    assert lhs == rhs # assoc

    def hom_into(self, tgt):
        obs = self.obs
        morphs = self.morphs
        ident = self.ident
        comp = self.comp

        
        # build the homs into tgt
        homs = dict((ob, []) for ob in obs)
        for f in morphs:
            if f.tgt == tgt:
                homs[f.src].append(f)
        homs = dict((ob, Set.promote(val)) for (ob, val) in homs.items())

        # now build the functions between these
        fs = {}
        for f in morphs:
            src = homs[f.tgt]
            tgt = homs[f.src]
            func = {}
            for x in src:
                func[x] = x*f
            func = Function(src, tgt, func)
            fs[f] = func

        tgt = CCat(homs.values(), fs.values())

        return Functor(self, tgt, homs, fs, covariant=False)


    def hom(self, src=None, tgt=None):
        assert src or tgt, "not implemented"

        if src is None:
            return self.hom_into(tgt)
        elif tgt is None:
            return self.hom_outof(src)
        else:
            assert 0, "not implemented"


class Functor(object):
    def __init__(self, src, tgt, action_obs, action_morphs, covariant=True):
        self.src = src
        self.tgt = tgt
        self.action_obs = dict(action_obs)
        self.action_morphs = dict(action_morphs)
        self.covariant = covariant
        self.check()

    def check(self):
        src = self.src
        tgt = self.tgt
        action_obs = self.action_obs
        action_morphs = self.action_morphs
        for ob in src.obs:
            assert action_obs[ob] in tgt.obs
        assert len(action_obs) == len(src.obs)
        for morph in src.morphs:
            assert action_morphs[morph] in tgt.morphs
        assert len(action_morphs) == len(src.morphs)

        if self.covariant:
            for f in src.morphs:
                assert action_obs[f.src] == action_morphs[f].src
                assert action_obs[f.tgt] == action_morphs[f].tgt

            for f in src.morphs:
              for g in src.morphs:
                if f.tgt != g.src:
                    continue
                h = g*f
                assert action_morphs[g]*action_morphs[f] == action_morphs[h]
        else:
            for f in src.morphs:
                assert action_obs[f.tgt] == action_morphs[f].src
                assert action_obs[f.src] == action_morphs[f].tgt

            for f in src.morphs:
              for g in src.morphs:
                if f.tgt != g.src:
                    continue
                h = g*f
                assert action_morphs[f]*action_morphs[g] == action_morphs[h]


    def __str__(self):
        return "%s(%s,\n    %s,\n    %s,\n    %s)"%(
            self.__class__.__name__, self.src, self.tgt, self.action_obs, self.action_morphs)
    __repr__ = __str__




class CCat(object):
    "A concrete category"

    def __init__(self, obs=[], morphs=[]):
        self.obs = [Set.promote(ob) for ob in obs]
        self.morphs = list(morphs)

    def __str__(self):
        return "CCat(%s)"%(self.obs,)
    __repr__ = __str__


class Set(object):
    "immutable, hashable set, needs comparable elements"
    def __init__(self, items):
        items = list(items)
        items.sort() # canonicalize
        self.set_items = set(items)
        self.items = tuple(items) # canonical form
        assert len(self.items) == len(self.set_items)
        #self._hash = hash(self.items)

    def __str__(self):
        return "Set%s"%(self.items,)
    __repr__ = __str__

    def __getitem__(self, idx):
        return self.items[idx]

    def __len__(self):
        return len(self.items)

    def __eq__(self, other):
        assert isinstance(other, Set)
        return self.items == other.items

    def __ne__(self, other):
        assert isinstance(other, Set)
        return self.items != other.items

    def __lt__(self, other):
        assert isinstance(other, Set)
        return self.items < other.items

    def __hash__(self):
        #return self._hash
        return hash(self.items)

    def __contains__(self, el):
        return el in self.set_items

    @classmethod
    def promote(cls, items):
        if isinstance(items, cls):
            return items
        return cls(items)


class Function(object):
    "a concrete morphism (in the category of sets)"
    def __init__(self, src, tgt, _action):
        self.src = Set.promote(src)
        self.tgt = Set.promote(tgt)
        action = {}
        for k, v in _action.items():
            assert k in self.src
            assert v in self.tgt
            action[k] = v
        assert len(_action) == len(src)
        self.action = action
        self.key = (src, tgt, action)

    def __str__(self):
        return "Function(%s)"%(self.action,)
    __repr__ = __str__

    def __mul__(g, f):
        assert isinstance(f, Function)
        if f.tgt != g.src:
            raise TypeError("not a composable pair")
        action = {}
        for k, v in f.action.items():
            v = g.action[v]
            action[k] = v
        return Function(f.src, g.tgt, action)

    @classmethod
    def ident(cls, src):
        action = dict((k, k) for k in src)
        return cls(src, src, action)

    def __ne__(self, other):
        return self.key != other.key

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        assert 0, "not implemented"


        



def main():

    # Example: Sets
    # one generic figure called P
    C = ACat(["P"])

    # Example: Bisets
    C = ACat(["X", "Y"])

    # Example: Bouqets
    # two generic figures: V for vertex, L for loop
    C = ACat(["V", "L"], [("v", "V", "L")])

    # Example: Graphs
    # vertex & arrow
    C = ACat(["V", "A"], [("s", "V", "A"), ("t", "V", "A")])

    # Example: Rgraphs
    # vertex & arrow
    s = ("s", "V", "A") # source
    t = ("t", "V", "A") # target
    l = ("l", "A", "V") # distinguished loop..
    sigma = ("sigma", "A", "A")
    tau = ("tau", "A", "A")
    C = ACat(["V", "A"], [s, t, l, sigma, tau], 
        {(l, s):"1_V", (l, t):"1_V", (s, l):sigma, (t, l):tau,
        (sigma, s):s, (sigma, t):s, (sigma, tau): sigma, (sigma, sigma):sigma,
        (tau, s):t, (tau, t):t, (tau, tau):tau, (tau, sigma):tau,
        (l, tau):l, (l, sigma):l,
    })
    #C.dump()

    # Example: mod-n dynamical system
    n = 4
    star = "*"
    morphs = [("1_*", star, star)] + [("f^%d"%i, star, star) for i in range(1, n)]
    assert len(morphs)==n
    comp = dict(((morphs[i], morphs[j]), morphs[(i+j)%n]) for i in range(n) for j in range(n))
    C = ACat([star], morphs, comp)
    #C.dump()

    # Example: height-n tree...

    print("OK")

if __name__ == "__main__":

    main()



