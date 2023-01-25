#!/usr/bin/env python3

"""
Implement the Cayley-Dickson _construction to
get complex numbers, quaternions, octonions, sedenions, etc.

copied from sedenions.py

Build the group G2(2) which is the automorphism group of
the integral octonions.
See book by Conway and Smith.

"""

import math, os
from random import choice
from functools import reduce
from operator import add, mul
from time import time
start_time = time()

from bruhat import element
from bruhat.util import choose, cross, all_perms
from bruhat.action import mulclose, mulclose_hom
from bruhat.gset import Perm, Group
from bruhat.argv import argv
from bruhat.isomorph import Point, Graph, search



class Number(object):
    def __init__(self, a):
        self.a = a
        self.shape = ()
        self.flat = (a,)

    def __str__(self):
        return str(self.a)

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__, self.a)

    def __hash__(self):
        assert 0
        print("__hash__", self.__class__.__name__)
        return hash(str(self))

    def promote(self, a):
        if not isinstance(a, Number):
            a = Number(a) # wrap it up
        return a

    def __add__(self, other):
        #assert 0
        assert self.__class__ is other.__class__
        return Number(self.a + other.a)

    def __sub__(self, other):
        #assert 0
        assert self.__class__ is other.__class__
        return Number(self.a - other.a)

    def __mul__(self, other):
        #assert 0
        assert self.__class__ is other.__class__
        return Number(self.a * other.a)

    def __rmul__(self, other):
        return Number(self.a * other)

    def __truediv__(self, other):
        return self.__class__(self.a / other)

    def __eq__(self, other):
        assert self.__class__ is other.__class__
        return self.a == other.a

    def __ne__(self, other):
        assert self.__class__ is other.__class__
        return self.a != other.a

    def __neg__(self):
        return Number(-self.a)

    def conj(self):
        return self

    def is_real(self):
        return True

    def is_zero(self):
        return self.a == 0

    def get_zero(self):
        return Number(0)


class Double(Number):
    def __init__(self, a, b):
        if not isinstance(a, Number):
            a = Number(a)
        if not isinstance(b, Number):
            b = Number(b)
        assert a.shape == b.shape
        self.shape = (a.shape, b.shape)
        assert isinstance(a, Number)
        self.a = a
        self.b = b
        self.flat = a.flat + b.flat

    def get_zero(self):
        return Double(self.a.get_zero(), self.b.get_zero())

    def promote(self, a):
        if not isinstance(a, Number):
            a = Number(a)
        if a.shape == self.shape:
            return a
        assert str(a.shape) in str(self.shape)
        a = self.a.promote(a)
        return Double(a, self.b.get_zero())

    def __repr__(self):
        #a, b = self.pair
        return "%s(%s, %s)"%(self.__class__.__name__, self.a, self.b)

    def __str__(self):
        return "(%s, %s)"%(self.a, self.b)

    def __lt__(self, other):
        return self.flat < other.flat

    def __add__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        a = self.a + other.a
        b = self.b + other.b
        return self.__class__(a, b)

    def __sub__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        a = self.a - other.a
        b = self.b - other.b
        return self.__class__(a, b)

    def __mul__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        a, b = self.a, self.b
        c, d = other.a, other.b
        x = self.__class__(a*c - d.conj()*b, d*a + b*c.conj())
        return x

    def __rmul__(self, other):
        return self.__class__(self.a * other, self.b * other)

    def __truediv__(self, other):
        return self.__class__(self.a / other, self.b / other)

    def __eq__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        return self.a == other.a and self.b == other.b

    def __ne__(self, other):
        other = self.promote(other)
        assert self.__class__ is other.__class__
        assert self.shape == other.shape
        return self.a != other.a or self.b != other.b

    def __hash__(self):
        return hash(str(self))

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.a, -self.b)

    def conj(self):
        return self.__class__(self.a.conj(), -self.b)

    def norm2(self):
        return self.conj() * self

    def is_real(self):
        return self.a.is_real() and self.b.is_zero()

    def is_zero(self):
        return self.a.is_zero() and self.b.is_zero()



def is_commutative(items):
    for a in items:
      for b in items:
        if a*b != b*a:
            return False
    return True


def is_anticommutative(items):
    for a in items:
      for b in items:
        if a!=b and a*b != -b*a:
            return False
    return True


def is_associative(items):
    for a in items:
      for b in items:
        for c in items:
            if a*(b*c) != (a*b)*c:
                return False
    return True


def is_alternative(items):
    for a in items:
      for b in items:
        if a*(b*b) != (a*b)*b:
            return False
    return True

def dot(a, b): 
    ab = reduce(add, [ai*bi for (ai,bi) in zip(a.flat,b.flat)])
    return ab



def get_geometry(imag):
    graph = Graph()
    N = len(imag)
    triples = set()
    cycles = []
    for idx in range(N):
        graph.add('p')
    for idx in range(N):
      for jdx in range(N):
        if idx==jdx:
            continue
        k = imag[idx]*imag[jdx]
        if k not in imag:
            continue
        kdx = imag.index(k)
        key = [idx, jdx, kdx]
        cycle = list(key)
        key.sort()
        key = tuple(key)
        if key in triples:
            continue
        triples.add(key)
        cycles.append(cycle)
        p = graph.add('l')
        graph.join(idx, p)
        graph.join(jdx, p)
        graph.join(kdx, p)

    return graph, cycles



def find_perms(imag):
    bag0, cycles = get_geometry(imag)
    bag1, cycles = get_geometry(imag)
    #print(cycles)

    N = 3
    struct = []
    for cycle in cycles:
        items = []
        for i in range(N):
            items.append(tuple(cycle[(i+j)%N] for j in range(N)))
        items.sort()
        #print items
        struct.append(items)
    struct.sort()

    for cycle in cycles:
        cycle = [bag0[i] for i in cycle]
        nbd = set(cycle[0].nbd)
        for point in cycle:
            nbd = nbd.intersection(point.nbd)
        assert len(nbd)==1

    #print(struct)

    perms = []
    count = 0
    total = 0
    for f in search(bag0, bag1):
        _struct = [[tuple(f[i] for i in cycle) for cycle in items] for items in struct ]
        for items in _struct:
            items.sort()
        _struct.sort()
        #print(_struct)
        if struct==_struct:
            count += 1
            #print("*")
        total += 1
        perms.append(f)

    return perms


class Frame(object):
    def __init__(self, els):
        els = list(els)
        assert len(els) == 7, len(els)
        els.sort()
        self.els = tuple(els)
        assert self.is_orthogonal()

    def __eq__(self, other):
        return self.els == other.els

    def __hash__(self):
        return hash(self.els)

    def __getitem__(self, idx):
        return self.els[idx]

    def __len__(self):
        return len(self.els)

    def is_orthogonal(self):
        for x in self:
          for y in self:
            if x is y:
                continue
            if dot(x, y) != 0:
                return False
        return True


def main():

    ring = element.Q
    one, zero = ring.one, ring.zero

    x = Number(2*one)
    y = Double(2*one, zero)
    assert x==y

    # ----------- double: complex --------------------

    one = Double(ring.one, ring.zero)
    i = Double(ring.zero, ring.one)
    zero = Double(ring.zero, ring.zero)
    assert i*i == -1

    cplex = [one, i]
    assert is_commutative(cplex)
    assert is_associative(cplex)

    # ----------- double: quaternions --------------------

    zero, one, i, j, k = [
        Double(zero, zero),
        Double(one, zero),
        Double(i, zero),
        Double(zero, one),
        Double(zero, i),
    ]

    for x in [i, j, k]:
        assert x*x == -1
        for y in [i, j, k]:
            if x==y:
                continue
            assert x*y == -y*x
    assert i*j == -j*i
    assert i*j*k == -1

    quaternions = [one, i, j, k]
    assert not is_commutative(quaternions)
    assert is_anticommutative(quaternions[1:])
    assert is_associative(quaternions)

    # ----------- double: octonions --------------------

    basis = [
        Double(one, zero), 
        Double(zero, one),
        Double(i, zero), 
        Double(j, zero), 
        Double(zero, i), 
        Double(k, zero),
        Double(zero, k),
        Double(zero, j), 
    ]
    zero = Double(zero, zero)

    obasis = [basis[i] for i in [0, 2, 3, 5, 1, 4, 7, 6]]

    fmt = "(" + 8*"%4s," + ")"
    #for u in obasis:
    #    print(fmt%u.flat)
    #return

    one = basis[0]
    iframe = basis[1:]
    assert one + one == 2*one
    for i in iframe:
        assert i*i == -one

    assert not is_commutative(basis)
    assert not is_associative(basis)
    assert is_anticommutative(basis[1:])
    assert is_alternative(basis)

    def get(desc, frame=iframe):
        x = zero
        sign = 1
        for a in desc:
            if a=='-':
                sign = -1
            elif a=='i':
                x = x+sign*one
                sign = +1
            else:
                a = int(a)
                x = x+sign*frame[a]
                sign = +1
        return x / 2
    i0, i1, i2, i3, i4, i5, i6 = iframe

#    found = 0
#    for iframe in all_perms(iframe):
#        i0, i1, i2, i3, i4, i5, i6 = iframe
#    
#        try:
#
#            #lhs = get("i026", iframe) * get("i013", iframe)
#            #rhs = get("0235", iframe)
#            #assert lhs == rhs
#
#            assert i1*i0 == i3
#            assert i0*i5 == i4
#            assert i2*i0 == i6
#            assert i1*i5 == i6
#
#            assert i1*i2 == i4
#            assert i2*i3 == i5
#            assert i3*i4 == i6
#            #assert i4*i5 == i0
#        
#            assert i2*i1 == -i4
#            assert i3*i2 == -i5
#            assert i4*i3 == -i6
#        
#            assert i2*i4 == i1
#            assert i3*i5 == i2
#            assert i4*i6 == i3
#            #assert i5*i0 == i4
#        
#            assert i4*i1 == i2
#            assert i5*i2 == i3
#            assert i6*i3 == i4
#            found += 1
#        except AssertionError:
#            pass
#        
#    assert found == 21, found

    assert i6 == iframe[6]

    # page 102: incorrect ?
    lhs = (one + i0 + i2 + i6) / 2
    assert lhs == get("i026")
    rhs = (one + i0 + i1 + i3) / 2
    assert rhs == get("i013")
    assert lhs*rhs == get("0156")

    assert get("0235") == (i0 + i2 + i3 + i5)/2
    assert get("02-35") == (i0 + i2 - i3 + i5)/2

    lhs = get('i356')
    rhs = get('-i356')
    assert lhs*lhs == rhs
    #assert get('i356') * get('0463') == get('-i461') # ???

    sgen = "0124 0235 0346 i450 0561 i602 i013 i356 146i 25i1 3612 4i23 5134 6245"
    sgen = sgen.split()
    gen = [get(a) for a in sgen]
    units = mulclose(gen) # just happens to work... but this is non-associative
    assert len(units) == 240
    units = set(units)
    assert one in units

    for i in basis:
        assert i in units
        assert -i in units

    # page 138
    sframe = [ -i0, +i1, +i2, -i3, +i4, -i5, -i6, ]
    hen = [get(a, sframe) for a in sgen]

    swap = mulclose_hom(gen, hen)
    assert len(swap) == len(units)

    #a = choice(list(units))
    #b = choice(list(units))
    #c = choice(list(units))
    #assert (a*b)*c == a*(b*c) # nope!

    # jframe, page 138
    j0 = get("-0124") 
    j1 = get("0-124") 
    j2 = get("01-24") 
    j3 = -i6
    j4 = get("012-4") 
    j5 = -i3
    j6 = -i5
    jframe = Frame([j0, j1, j2, j3, j4, j5, j6])

    assert j0 + j1 == i2 + i4
    assert j0 - j1 == i1 - i0
    
    #lhs, rhs = j1*j0 , (j2+j3-j4+j6)/2 # page 139
    lhs, rhs = j0*j1 , (j2-j3-j4-j6)/2
    assert lhs == rhs

    def find(lhs, frame=iframe):
        for jdxs in choose(list(range(7)), 4):
          js = [frame[jdx] for jdx in jdxs]
          for signs in cross([(-1, 1)]*4):
            x = reduce(add, [(sign*u)/2 for (sign, u) in zip(signs, js)])
            if lhs != x:
                continue
            #print("found:", jdxs, signs)
            decl = []
            for (jdx, sign) in zip(jdxs, signs):
                decl.append( str(jdx) if sign==1 else "-"+str(jdx) )
            decl = ''.join(decl)
            assert x == get(decl, frame)
            return decl

    for j0 in jframe:
        assert j0*j0 == -one
        for j1 in jframe:
            assert j0*j1 == (j1*j0).conj()
            #for j2 in jframe:
            #    assert (j0*j1)*j2 == j0*(j1*j2) # nope..

    # non-associative...
    #J = mulclose(jframe)
    #assert J.issubset(units)

    k0 = -i0
    k1 = get("-1263") 
    k2 = get("-2456")
    k3 = get("-1-2-63")
    k4 = get("-4135")
    k5 = get("-4-1-35")
    k6 = get("-2-4-56")
    kframe = Frame([k0, k1, k2, k3, k4, k5, k6])

    lhs, rhs = k0*k1, get("3-621")
    rhs = get("-1-2-36")
    assert lhs == rhs

    found = []
    for a in units:
        if dot(a, one) == 0:
            found.append(a)
    assert len(found) == 126

    found = []
    for a in units:
        if a==one:
            continue
        if a*a==one:
            continue
        if a*a*a==one:
            found.append(a)
    assert len(found) == 56

    lookup = {unit:idx for (idx, unit) in enumerate(units)}
    rlookup = {idx:unit for (unit, idx) in lookup.items()}

    n = len(units)
    swap = Perm([lookup[swap[rlookup[idx]]] for idx in range(n)])

    perms = []
    for x in found:
        xc = x.conj()
        perm = [None]*n
        for a in units:
            b = xc*a*x
            perm[lookup[a]] = lookup[b]
        assert None not in perm
        perm = Perm(perm)
        perms.append(perm)
    G = Group(None, perms)
    print("|G| =", len(G))

    perms = list(G) + [g*swap for g in G]
    G2 = Group(perms)
    print("|G2| =", len(G2))

    #for H in G.cyclic_subgroups():
    #for H in G.subgroups():
    #    print(len(H), end=" ", flush=True)
    #print()

    def preserves(g, frame):
        idxs = [lookup[u] for u in frame]
        for idx in idxs:
            if g[idx] not in idxs:
                return False
        return True

    def sends(g, iframe, jframe):
        idxs = [lookup[u] for u in iframe]
        jdxs = [lookup[u] for u in jframe]
        for idx in idxs:
            if g[idx] not in jdxs:
                return False
        return True

    if 0:
        # build kframe from an orbit
        k0 = -i0
        kdxs = set([lookup[k0]])
        for g in G:
            if not preserves(g, jframe):
                continue
            for kdx in list(kdxs):
                kdxs.add(g[kdx])
            if len(kdxs) == 7:
                break
    
        kframe = [rlookup[kdx] for kdx in kdxs]
        for k in kframe:
            print(find(k))

    kframe = "-13-4-5 1-236 2-456 -1345 -2-45-6 -1-2-36".split()
    kframe = Frame([-i0] + [get(desc) for desc in kframe])
    k0, k1, k2, k3, k4, k5, k6 = kframe

    stab = []
    for g in G:
        if sends(g, jframe, kframe):
            assert 0, "jframe -> kframe"

        lhs, rhs = preserves(g, jframe), preserves(g, kframe)
        assert lhs == rhs
        if lhs:
            stab.append(g)
    assert len(stab) == 21

    jframe = Frame(jframe)

    act = lambda g, u: rlookup[g[lookup[u]]]

    hom = {}
    for g in stab:
        lhs, rhs = act(g,j0), act(g,k0)
        assert lhs in jframe
        assert rhs in kframe
        #assert hom.get(lhs, rhs) == rhs
        hom[lhs] = rhs
    assert len(hom) == 7

    frames = set()
    for g in G:
        frame = Frame(act(g, k) for k in kframe)
        assert frame != jframe
        frames.add(frame)
    assert len(frames) == 6048 // 21 == 288

    # page 134
    stab = []
    for g in G:
        gi0 = act(g, i0)
        if gi0 == i0 or gi0 == -i0:
            stab.append(g)
    assert len(stab) == 96

    assert get("i450") in units

    # Build the Hurwitz units
    # page 56
    quat = [one, i4, i5, i0]
    found = set(quat+[-u for u in quat])
    for signs in cross([(-1, 1)]*4):
        x = reduce(add, [sign*u for (sign, u) in zip(signs, quat)])
        x = x/2
        assert x in units
        found.add(x)
    assert len(found) == 24

    # page 134
    stab = []
    for g in G:
        fix1 = set(act(g, i) for i in found)
        if fix1==found:
            stab.append(g)
    assert len(stab) == 96


if __name__ == "__main__":

    from argv import argv
    name = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%name)
    else:
        fn = eval(name)
        fn()

    print("OK: %.3f seconds\n"%(time() - start_time))






