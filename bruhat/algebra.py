#!/usr/bin/env python3

"""
An algebra is a vector space with a bilinear multiplication.
"""


from bruhat.argv import argv
from bruhat.action import Group
from bruhat import element
from bruhat import elim
from bruhat.element import Keyed, Type, Element, GenericElement
from bruhat.vec import Space, Map, Hom
from bruhat.rep import Rep, Cat, tensor_rep, Young
from bruhat.util import partitions

# _based on bruhat.vec.Map:

class Algebra(Type):

    def __init__(self, gen, _struct, ring):
        """
            _struct : list of quads (a, b, c, x) with a, b, c in gen, x in ring.
        """
        self.ring = ring
        self.gen = tuple(gen)
        self.set_gen = set_gen = set(gen)
        Type.__init__(self)

        struct = {} # map (a, b) -> Vector
        uniq = set()
        for (a, b, c, val) in _struct:
            assert a in set_gen
            assert b in set_gen
            assert c in set_gen
            val = ring.promote(val)
            items = struct.setdefault((a, b), [])
            items.append((c, val))
            key = (a, b, c)
            assert key not in uniq, "duplicate struct entry: %s"%(key,)
            uniq.add(key) 
        for a in gen:
          for b in gen:
            k = (a, b)
            struct[k] = Vector(struct.get(k, []), self)
        self.struct = struct
        self.zero = Vector([], self)
        self.basis = [Vector([(g, ring.one)], self) for g in gen]

    # equality on the nose!
    def __eq__(self, other):
        return (self.ring == other.ring and 
            self.gen == other.gen and 
            self.struct == other.struct)

    def __ne__(self, other):
        return (self.ring != other.ring or 
            self.gen != other.gen or 
            self.struct != other.struct)

    def __getitem__(self, i):
        return self.gen[i]

    def __len__(self):
        return len(self.gen)

    def __contains__(self, x):
        return x in self.set_gen

    def __str__(self):
        return "%s(%s, %s)"%(
            self.__class__.__ne__, 
            str(list(self.gen)), self.ring)
    __repr__ = __str__

    def promote(self, g):
        if isinstance(g, Vector) and g.algebra is self:
            v = g
        elif g in self.set_gen:
            v = Vector([(g, self.ring.one)], self)
        else:
            v = None
        return v


class GModule(Type):
    def __init__(self, algebra, rep):
        self.algebra = algebra
        self.rep = rep
        self.hom = rep.hom
        Type.__init__(self)

    def action(self, g):
        "compute action of g on V"
        algebra = self.algebra
        rep = self.rep
        v = algebra.promote(g)
        space = rep.space
        hom = self.hom
        f = hom.zero_vector()
        for (g, val) in v.items:
            g = val * rep.action(g)
            f = f+g
        return f


class GroupAlgebra(Algebra):
    def __init__(self, G, ring):
        gen = list(G)
        #print(gen)
        struct = []
        for g in G:
          for h in G:
            struct.append((g, h, g*h, ring.one))
        Algebra.__init__(self, gen, struct, ring)
        self.G = G

    def extend(self, rep):
        assert rep.G is self.G
        return GModule(self, rep)


class Vector(Element):

    def __init__(self, _items, algebra):
        Element.__init__(self, algebra)
        ring = algebra.ring
        zero = ring.zero
        items = []
        keys = []
        for item in _items:
            i, _val = item
            val = ring.promote(_val)
            assert val is not None, (ring, _val)
            if val != zero: # make canonical. sparse
                items.append((i, val))
                keys.append(i)
            assert val.tp == ring
            assert i in algebra, "%s not found in %s" % (i, algebra)
        assert len(keys)==len(items), "duplicate key: %s" % (keys)
        #items.sort() # Make canonical. careful with this...
        self.items = tuple(items)
        self.map_items = dict(items)
        self.ring = ring
        self.algebra = algebra

    def __getitem__(self, k):
        assert k in self.algebra
        return self.map_items.get(k, self.ring.zero)

    def __str__(self):
        vals = [self[k] for k in self.algebra]
        return "[%s]" % (', '.join(str(val) for val in vals))
    __repr__ = __str__

    def to_array(self):
        vals = [self[k] for k in self.algebra]
        return vals

    def __eq__(a, b):
        return a.map_items == b.map_items

    def __ne__(a, b):
        return a.map_items != b.map_items

    def __hash__(self):
        assert 0, "needs canonical form"
        return hash(self.items)

    def __add__(a, b):
        algebra = a.algebra
        b = algebra.promote(b)
        assert algebra is b.algebra
        zero = a.ring.zero
        map_items = dict(a.items)
        for (i, u) in b.items:
            val = map_items.get(i, zero) + u
            map_items[i] = val
        items = list(map_items.items())
        return Vector(items, a.algebra)

    def __sub__(a, b):
        algebra = a.algebra
        b = algebra.promote(b)
        assert algebra is b.algebra
        zero = a.ring.zero
        map_items = dict(a.items)
        for (i, u) in b.items:
            val = map_items.get(i, zero) - u
            map_items[i] = val
        items = list(map_items.items())
        return Vector(items, a.algebra)

    def __neg__(a):
        items = [(i, -u) for (i, u) in a.items]
        return Vector(items, a.algebra)

    def __mul__(a, b):
        algebra = a.algebra
        b = algebra.promote(b)
        assert algebra is b.algebra
        c = algebra.zero
        for (x, ax) in a.items:
          for (y, by) in b.items:
            d = (ax*by) * algebra.struct[x, y]
            c = c + d
        return c

    def __rmul__(a, r):
        r = a.ring.promote(r)
        items = [(i, r*u) for (i, u) in a.items]
        return Vector(items, a.algebra)
        



def test():

    ring = element.Q

    n = argv.get("n", 3)
    G = Group.symmetric(n)

    algebra = GroupAlgebra(G, ring)
    v = algebra.zero
    for g in G:
        g = algebra.promote(g)
        v = v+g
    assert 2*v == v+v
    assert v-v == algebra.zero

    for g in algebra.basis:
        assert g*v == v # oh yeah !
        assert v*g == v # come get it!

    act = G.left_action(list(range(n)))
    #print(act)

    tp = Cat(G, ring)
    rep = Rep.mk_rep(act, tp)
    #rep.dump()

    module = algebra.extend(rep)
    #print(module.action(v))

    # ---------------------------------------

    space = Space(2, ring)
    rep = tensor_rep(G, space)
    #rep.check()

    module = algebra.extend(rep)
    #print(module.action(v))

    # ---------------------------------------

    d = len(space)

    # Build the young symmetrizers
    projs = []
    for part in partitions(n):
        t = Young(part)

        rowG = t.get_rowperms()
        colG = t.get_colperms()
        H = None
        for g in rowG:
            P = algebra.promote(g)
            H = P if H is None else (H + P)

        V = None
        for g in colG:
            P = algebra.promote(g)
            s = g.sign()
            P = s*P
            V = P if V is None else (V + P)
        A = V * H

        assert V*V == len(colG) * V
        assert H*H == len(rowG) * H
        print("part:", part)
        #print("H:")
        #print(module.action(H))
        #print("V:")
        #print(module.action(V))

        specht = []
        for v in algebra.basis:
            u = v*A
            specht.append(u.to_array())
        hom = Space(len(G), ring).endo_hom()
        specht = Map.from_array(specht, hom)
        specht = specht.transpose()
        #print(specht)
        specht = specht.image()
        dim = specht.shape[1]
        print("dim:", dim)

        if len(part) > d:
            continue

        P = module.action(A)
        P = P.transpose()
        print(P.rank())

#        break
#        print("H*V + V*H:")
#        print(H*V + V*H)
#        print(module.action(H*V + V*H))
#        projs.append(A)
#
#        break
#        #print(part)
#        #print(t)
#        #print(A)





if __name__ == "__main__":

    test()





