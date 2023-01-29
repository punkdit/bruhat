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
from random import choice, random, randint
from functools import reduce
from operator import add, mul
from time import time
start_time = time()

from bruhat import element
from bruhat.util import choose, cross, all_perms
from bruhat.action import mulclose, mulclose_hom
from bruhat.gset import Perm, Group, Coset
from bruhat.argv import argv
from bruhat import isomorph # import Point, Graph, search


try:
    from huygens.namespace import *
    from huygens.pov import Mat
    from huygens import config
    #config(text="pdflatex", latex_header=r"""
    #\usepackage{amsmath}
    #\usepackage{amssymb}
    #""")
except ImportError:
    pass


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
    graph = isomorph.Graph()
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
    for f in isomorph.search(bag0, bag1):
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


def get_orbits(g, X):

    found = set()
    orbits = []
    for H in X:
        if H in found:
            continue
        orbit = [H]
        found.add(H)
        H1 = H
        while 1:
            H1 = H1.left_mul(g)
            if H1 == H:
                break
            assert H1 not in found
            orbit.append(H1)
            found.add(H1)
        orbits.append(orbit)
    return orbits
        

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

    if 0:
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

    # page 134: stabilize Gaussian integers
    H = []
    for g in G:
        gi0 = act(g, i0)
        if gi0 == i0 or gi0 == -i0:
            H.append(g)
    assert len(H) == 96

    # Build the Hurwitz units
    # page 56
    assert get("i450") in units
    quat = [one, i4, i5, i0]
    found = set(quat+[-u for u in quat])
    for signs in cross([(-1, 1)]*4):
        x = reduce(add, [sign*u for (sign, u) in zip(signs, quat)])
        x = x/2
        assert x in units
        found.add(x)
    assert len(found) == 24

    # page 134: stabilize Hurwitz integers
    K = []
    for g in G:
        fix1 = set(act(g, i) for i in found)
        if fix1==found:
            K.append(g)
    assert len(K) == 96

    HK = [g for g in H if g in K]
    print(len(HK))

    H = Group(H)
    K = Group(K)

    for g in G:
        if g.order() != 6:
            continue
        if g in H:
            continue
        if g in K:
            continue
        break

    #X = G.action_subgroup(H)
    #Y = G.action_subgroup(K)

    n = 63
    points = list(range(n))
    lines = [i+n for i in range(n)]

    X = G.left_cosets(H)
    lookup = {H:idx for (idx,H) in enumerate(X)}
    orbits = get_orbits(g, X)
    orbits.sort(key = len)
    print([len(orbit) for orbit in orbits])
    Xs = ([[lookup[H] for H in orbit] for orbit in orbits])

    Y = G.left_cosets(K)
    lookup = {H:idx for (idx,H) in enumerate(Y)}
    orbits = get_orbits(g, Y)
    orbits.sort(key = len)
    print([len(orbit) for orbit in orbits])
    Ys = ([[lookup[H]+n for H in orbit] for orbit in orbits])

    assert len(X) == len(Y) == n
    pairs = []
    for idx, x in enumerate(X):
      for jdx, y in enumerate(Y):
        xy = x.intersect(y)
        if len(xy):
            pairs.append((idx, jdx+n))
    print(pairs)

    geometry = Geometry(points, lines, pairs)
    #geometry.todot()

    for orbit in Xs:
        m = len(orbit)
        print('  ', end='')
        for i in range(m):
            print(geometry.dist[orbit[i], orbit[(i+1)%m]], end=' ')
        print()


class Geometry(object):
    def __init__(self, points, lines, pairs):
        nn = len(points)+len(lines)
        nbd = {i:[] for i in range(nn)}
        for (i,j) in pairs:
            nbd[i].append(j)
            nbd[j].append(i)
            assert (j,i) not in pairs
        nbd2 = {i:set() for i in range(nn)}
        for i in range(nn):
          for j in nbd[i]:
            nbd2[i].update(nbd[j])
          nbd2[i].remove(i)
        #print(nbd2)
    
        dist = {(i,i):0 for i in range(nn)}
        done = False
        while not done:
            done = True
            keys = list(dist.keys())
            for (i1, i2) in keys:
                d0 = dist[i1, i2]
                for i0 in nbd[i1]:
                    d = dist.get((i0, i2))
                    if d is None or d > d0+1:
                        dist[i0, i2] = d0+1
                        dist[i2, i0] = d0+1
                        done = False
        assert len(dist) == nn*nn
        assert max(dist.values()) == 6

        self.dist = dist
        self.nbd = nbd
        self.points = list(points)
        self.lines = list(lines)
        self.pairs = list(pairs)

    def get_geodesic(src, tgt):
        nbd = self.nbd
        tree = {src:src} # root
        bdy = [src]
        while tgt not in tree:
            assert bdy
            _bdy = []
            for idx in bdy:
                for jdx in nbd[idx]:
                    if jdx in tree:
                        continue
                    tree[jdx] = idx
                    _bdy.append(jdx)
            bdy = _bdy
        path = [tgt]
        while path[-1] != src:
            path.append(tree[path[-1]])
        return list(reversed(path))

    def render(self, layout):
        points = self.points
        lines = self.lines
        pairs = self.pairs
        cvs = Canvas()
    
        for (idx, jdx) in pairs:
            if idx in layout and jdx in layout:
                x0, y0 = layout[idx]
                x1, y1 = layout[jdx]
                cvs.stroke(path.line(x0, y0, x1, y1), [black.alpha(0.5)]+st_THICk)
    
        for idx in points:
            if idx in layout:
                x, y = layout[idx]
                cvs.stroke(path.circle(x, y, 0.2), [red.alpha(0.5)]+st_THICk)
        for idx in lines:
            if idx in layout:
                x, y = layout[idx]
                cvs.stroke(path.circle(x, y, 0.2), [blue.alpha(0.5)]+st_THICk)
    
        cvs.writePDFfile("G2_render.pdf")

    def todot(self, name="G2_graph.dot"):
        points = self.points
        lines = self.lines
        pairs = self.pairs
        f = open(name, 'w')
        out = lambda s : print(s, file=f)
        out("graph {")
        out("""
        node [
            shape = circle
            style = filled
            color = "#00000000"
            fillcolor = black
            width = 0.1
            height = 0.1
            label = ""
        ]
        edge [
            penwidth = 2.0
        ]
        """)
        for i in points:
            out("    %d [fillcolor=red];"%i)
        for i in lines:
            out("    %d [fillcolor=blue];"%(i))
        for (idx, jdx) in pairs:
            out("   %d -- %d;"%(idx, jdx))
        out("}")
        f.close()
    
        import os
        rval = os.system("neato G2_graph.dot -Tpdf > G2_graph.pdf") 
        assert rval == 0
        


def get_geometry():

    pairs = [(0, 63), (0, 87), (0, 119), (1, 64), (1, 97),
    (1, 111), (2, 65), (2, 90), (2, 106), (3, 66), (3, 76),
    (3, 96), (4, 67), (4, 113), (4, 117), (5, 68), (5, 104),
    (5, 118), (6, 69), (6, 103), (6, 121), (7, 70), (7, 88),
    (7, 105), (8, 71), (8, 109), (8, 110), (9, 72), (9, 82),
    (9, 125), (10, 63), (10, 75), (10, 100), (11, 72), (11,
    98), (11, 114), (12, 64), (12, 86), (12, 108), (13, 71),
    (13, 92), (13, 95), (14, 68), (14, 84), (14, 93), (15,
    70), (15, 83), (15, 120), (16, 67), (16, 78), (16, 99),
    (17, 66), (17, 102), (17, 123), (18, 73), (18, 101),
    (18, 124), (19, 74), (19, 94), (19, 122), (20, 69), (20,
    79), (20, 116), (21, 73), (21, 80), (21, 89), (22, 74),
    (22, 81), (22, 85), (23, 65), (23, 112), (23, 115), (24,
    75), (24, 103), (24, 108), (25, 64), (25, 71), (25, 77),
    (26, 81), (26, 105), (26, 114), (27, 82), (27, 104),
    (27, 109), (28, 80), (28, 110), (28, 117), (29, 81),
    (29, 102), (29, 118), (30, 78), (30, 92), (30, 124),
    (31, 85), (31, 87), (31, 101), (32, 84), (32, 95), (32,
    98), (33, 85), (33, 113), (33, 121), (34, 76), (34, 88),
    (34, 97), (35, 77), (35, 91), (35, 107), (36, 83), (36,
    86), (36, 123), (37, 79), (37, 111), (37, 119), (38,
    78), (38, 79), (38, 94), (39, 89), (39, 111), (39, 118),
    (40, 90), (40, 103), (40, 104), (41, 89), (41, 115),
    (41, 120), (42, 75), (42, 80), (42, 122), (43, 87), (43,
    96), (43, 109), (44, 66), (44, 69), (44, 91), (45, 93),
    (45, 112), (45, 116), (46, 92), (46, 100), (46, 102),
    (47, 65), (47, 74), (47, 77), (48, 86), (48, 93), (48,
    101), (49, 67), (49, 68), (49, 107), (50, 95), (50, 120),
    (50, 121), (51, 97), (51, 113), (51, 125), (52, 98),
    (52, 106), (52, 119), (53, 88), (53, 90), (53, 124),
    (54, 82), (54, 83), (54, 94), (55, 96), (55, 99), (55,
    115), (56, 99), (56, 108), (56, 114), (57, 100), (57,
    112), (57, 125), (58, 72), (58, 73), (58, 91), (59, 76),
    (59, 84), (59, 122), (60, 63), (60, 70), (60, 107), (61,
    105), (61, 110), (61, 116), (62, 106), (62, 117), (62, 123)]

    Xs = [
        [20, 61, 45], 
        [6, 26, 23, 44, 7, 57], 
        [8, 14, 37, 28, 48, 38], 
        [0, 4, 12, 19, 43, 49, 1, 42, 31, 16, 25, 59], 
        [2, 17, 15, 9, 40, 29, 41, 58, 53, 46, 50, 11], 
        [3, 60, 51, 24, 22, 55, 35, 34, 10, 33, 56, 47], 
        [5, 39, 21, 18, 30, 13, 32, 52, 62, 36, 54, 27],
    ]

    Ys = [
        [116], 
        [69, 105, 112], 
        [79, 110, 93], 
        [64, 122, 87, 67], 
        [72, 90, 102, 120], 
        [63, 113, 108, 74, 96, 107, 97, 75, 85, 99, 77, 76], 
        [65, 66, 70, 125, 103, 81, 115, 91, 88, 100, 121, 114], 
        [68, 111, 80, 101, 78, 71, 84, 119, 117, 86, 94, 109], 
        [73, 124, 92, 95, 98, 106, 123, 83, 82, 104, 118, 89],
    ]
    # interleaved Xs & Ys
    orbits = [
        [116], # blue
        [20, 61, 45], # red
        [69, 105, 112], # blue
        [6, 26, 23, 44, 7, 57], # red
        [79, 110, 93], # blue
        [8, 14, 37, 28, 48, 38], # red
        [64, 122, 87, 67], # blue
        [0, 4, 12, 19, 43, 49, 1, 42, 31, 16, 25, 59], # red
        [72, 90, 102, 120], # blue
        [2, 17, 15, 9, 40, 29, 41, 58, 53, 46, 50, 11], # red
        [63, 113, 108, 74, 96, 107, 97, 75, 85, 99, 77, 76], # blue
        [3, 60, 51, 24, 22, 55, 35, 34, 10, 33, 56, 47], # red
        [65, 66, 70, 125, 103, 81, 115, 91, 88, 100, 121, 114], # blue
        [68, 111, 80, 101, 78, 71, 84, 119, 117, 86, 94, 109], # blue
        [73, 124, 92, 95, 98, 106, 123, 83, 82, 104, 118, 89], # blue
        [5, 39, 21, 18, 30, 13, 32, 52, 62, 36, 54, 27], # red
    ] # cost: 1419

    # interleaved Xs & Ys
    orbits = [
        [116], # blue
        [20, 61, 45], # red
        [69, 105, 112], # blue
        [6, 26, 23, 44, 7, 57], # red
        [79, 110, 93], # blue
        [8, 14, 37, 28, 48, 38], # red
        [64, 122, 87, 67], # blue
        [0, 4, 12, 19, 43, 49, 1, 42, 31, 16, 25, 59], # red
        [72, 90, 102, 120], # blue
        [63, 113, 108, 74, 96, 107, 97, 75, 85, 99, 77, 76], # blue
        [3, 60, 51, 24, 22, 55, 35, 34, 10, 33, 56, 47], # red
        [68, 111, 80, 101, 78, 71, 84, 119, 117, 86, 94, 109], # blue
        [65, 66, 70, 125, 103, 81, 115, 91, 88, 100, 121, 114], # blue
        [2, 17, 15, 9, 40, 29, 41, 58, 53, 46, 50, 11], # red
        [73, 124, 92, 95, 98, 106, 123, 83, 82, 104, 118, 89], # blue
        [5, 39, 21, 18, 30, 13, 32, 52, 62, 36, 54, 27], # red
    ]  # cost: 1169

    return pairs, Xs, Ys, orbits


def make_autos():
    from pynauty import Graph, autgrp
    n = 63
    nn = 2*n
    pairs = get_geometry()[0]
    nbd = {i:[] for i in range(nn)}
    for (i,j) in pairs:
        nbd[i].append(j)
        nbd[j].append(i)
        assert (j,i) not in pairs

    graph = Graph(nn)
    for k,v in nbd.items():
        graph.connect_vertex(k, v)
    gens, sz1, sz2, orbits, norbits = autgrp(graph)
#    for g in gens:
#        print(g)
    print(len(gens))

    gens = [Perm(perm) for perm in gens]
    G = Group(None, gens)
    print(len(G))

    # element of G that has orbit sizes:
    # [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7] == [7]*18
    found = []
    for g in G:
        k = g.order()
        if k == 7:
            orbits = g.get_orbits()
            sizes = [len(o) for o in orbits]
            if sizes == [k]*len(sizes):
                found.append(g)
    print("found:", len(found))
    g = found[0]

    def make_graph(edges):
        nbd = {i:[] for i in range(len(orbits))}
        for (i,j) in edges:
            nbd[i].append(j)
        graph = isomorph.Graph()
        for i,o in enumerate(orbits):
            graph.add('p')
        for (i,j) in edges:
            graph.join(i, j)
        return graph

    graphs = []
    for g in found[:10]:
        orbits = g.get_orbits()
        lookup = {}
        for i, orbit in enumerate(orbits):
            for idx in orbit:
                lookup[idx] = i
            #print("%3d"%i, orbit)
    
        edges = set()
        for src, o in enumerate(orbits):
            for idx in o:
              for jdx in nbd[idx]:
                tgt = lookup[jdx]
                edges.add((src, tgt))
    
        graph0 = make_graph(edges)
        #graph1 = make_graph(edges)
        graphs.append(graph0)
    
    for i in range(10):
     for j in range(i+1, 10):
        autos = []
        for f in isomorph.search(graphs[i], graphs[j]):
            autos.append(f)
        #print(len(autos), end=' ', flush=True)
        assert len(autos)==6
    print()
    
    if 0:
        f = open("orbits.dot", 'w')
        output = lambda s : print(s, file=f)
        output("graph {")
        for (i,j) in edges:
            if i < j:
                output("  %d -- %d;" %(i, j))
        output("}")
        f.close()
    
        import os
        rval = os.system("circo orbits.dot -Tpdf > orbits.pdf") 
        assert rval == 0
    
        


def make_code():
    pairs = get_geometry()[0]

    from bruhat.solve import zeros2, shortstr, rank, row_reduce
    n = 63
    H = zeros2(n, n)
    for (i,j) in pairs:
        H[i,j-n] = 1
    H = 1-H
    print(shortstr(H))
    H1 = (row_reduce(H))
    print()
    print(shortstr(H1))
    print(H1.shape)



def optimize_graph_G2():
    pairs, Xs, Ys, orbits = get_geometry()

    # interleaved Xs & Ys
#    orbits = []
#    while Xs or Ys:
#        if Ys:
#            orbits.append(Ys.pop(0))
#        if Xs:
#            orbits.append(Xs.pop(0))
    #for orbit in orbits:
    #    print(orbit)

    def get_layout(thetas):
        layout = {}
        radius = 0.
        delta = 1.
        theta = 0.
        for theta, orbit in zip(thetas, orbits):
            for idx in orbit:
                x, y = radius * sin(theta), radius * cos(theta)
                layout[idx] = (x, y)
                theta += 2*pi / len(orbit)
            radius += delta
            #theta += pi / len(orbit) + pi / 48
        cost = 0.
        for (idx, jdx) in pairs:
            x0, y0 = layout[idx]
            x1, y1 = layout[jdx]
            cost += ((x1-x0)**2 + (y1-y0)**2)**0.5
        return cost, layout

    N = len(orbits)
    thetas = [random()*0*pi]*N
    best = None
    for trial in range(6):
      for i in range(2, N):
        for t in range(48):
            thetas[i] = t*2*pi/48
            cost, layout = get_layout(thetas)
            if best is None or best[0] > cost:
                best = (cost, layout, list(thetas))
                print("cost:", best[0])
        thetas = list(best[2])
    if 1:
        thetas = best[2]
        #thetas[4] += pi/12
        layout = get_layout(thetas)[1]
    elif 0:
        layout = best[1]
    else:
        thetas = best[2]
        for i in range(2, N):
            thetas[i] += i*pi/96
        layout = get_layout(thetas)[1]

    cvs = Canvas()

    n = 63

    for (idx, jdx) in pairs:
        x0, y0 = layout[idx]
        x1, y1 = layout[jdx]
        cvs.stroke(path.line(x0, y0, x1, y1))
    for idx in range(n):
        x, y = layout[idx]
        cvs.fill(path.circle(x, y, 0.1), [red])
        x, y = layout[n+idx]
        cvs.fill(path.circle(x, y, 0.1), [blue])

    cvs.writePDFfile("G2.pdf")
    


def graph_G2():
    pairs, Xs, Ys, orbits = get_geometry()

    n = 63
    nn = 2*63
    nbd = {i:[] for i in range(nn)}
    for (i,j) in pairs:
        nbd[i].append(j)
        nbd[j].append(i)
        assert (j,i) not in pairs
    nbd2 = {i:set() for i in range(nn)}
    for i in range(nn):
      for j in nbd[i]:
        nbd2[i].update(nbd[j])
      nbd2[i].remove(i)
    #print(nbd2)

    dist = {(i,i):0 for i in range(nn)}
    done = False
    while not done:
        done = True
        keys = list(dist.keys())
        for (i1, i2) in keys:
            d0 = dist[i1, i2]
            for i0 in nbd[i1]:
                d = dist.get((i0, i2))
                if d is None or d > d0+1:
                    dist[i0, i2] = d0+1
                    dist[i2, i0] = d0+1
                    done = False
    assert len(dist) == nn*nn
    assert max(dist.values()) == 6

    radius = 20.
    layout = {}

    bdy = [5, 27, 54, 36, 62, 52, 32, 13, 30, 18, 21, 39]
    thetas = [ (2*pi*i)/len(bdy) for i in range(len(bdy))]
    rank = {i : float(ii) for (ii,i) in enumerate(bdy)}
    
    for idx, theta in zip(bdy, thetas):
        layout[idx] = (radius*sin(theta), radius*cos(theta))

    parent = {i:[] for i in range(nn)}
    found = set(bdy)
    while bdy:
        _bdy = []
        for idx in bdy:
            for j in nbd[idx]:
                parent[j].append(idx)
                if j not in found:
                    found.add(j)
                    _bdy.append(j)
        bdy = _bdy
    
    #while len(layout) < nn:
    for trial in range(10):
        layer = []
        size = 1
        for idx in range(nn):
            if idx in layout:
                continue
            jdxs = [j for j in parent[idx] if layout.get(j) is not None]
            size = max(size, len(jdxs))

        for idx in range(nn):
            if idx in layout:
                continue
            jdxs = [j for j in parent[idx] if layout.get(j) is not None]
            assert len(jdxs) <= size
            if len(jdxs) == size:
                layer.append(idx)
        print(size, layer)
        for idx in layer:
            jdxs = [j for j in parent[idx] if layout.get(j) is not None]
            assert rank.get(idx) is None
            rank[idx] = sum(rank[j] for j in jdxs) / len(jdxs)
        #layer.sort(key = lambda idx : rank[idx])

        counter = 0
        for idx in layer:
            jdxs = parent[idx]
            jdxs = [jdx for jdx in jdxs if jdx in layout]
            jdxs.sort(key = lambda jdx : rank[jdx])
            pts = [layout[jdx] for jdx in jdxs]
            x = reduce(add, [p[0] for p in pts])/len(pts)
            y = reduce(add, [p[1] for p in pts])/len(pts)
            if abs(x) < 1e-4 and abs(y) < 1e-4:
                #jdxs.pop(randint(0, len(jdxs)-1))
                jdxs.pop(counter%len(jdxs))
                counter += 1
                pts = [layout[jdx] for jdx in jdxs]
                x = reduce(add, [p[0] for p in pts])/len(pts)
                y = reduce(add, [p[1] for p in pts])/len(pts)
            if len(jdxs) == 1:
                x *= 0.9
                y *= 0.9
            layout[idx] = x, y

    print(len(layout), "of", nn)
            
    for key in list(layout.keys()):
        x, y = layout[key]
        v = x + y*1.j
        u = 0.5 + 0.01*abs(v)*1.j
        v = u*v
        layout[key] = v.real, v.imag

    cvs = Canvas()

    for (idx, jdx) in pairs:
        if idx in layout and jdx in layout:
            x0, y0 = layout[idx]
            x1, y1 = layout[jdx]
            cvs.stroke(path.line(x0, y0, x1, y1), [black.alpha(0.5)]+st_THICk)

    for idx in range(n):
        if idx in layout:
            x, y = layout[idx]
            cvs.stroke(path.circle(x, y, 0.2), [red.alpha(0.5)]+st_THICk)
        if n+idx in layout:
            x, y = layout[n+idx]
            cvs.stroke(path.circle(x, y, 0.2), [blue.alpha(0.5)]+st_THICk)

    cvs.writePDFfile("G2_geometry.pdf")
    

def make_render():

    pairs, Xs, Ys, orbits = get_geometry()

    n = 63
    nn = 2*63
    nbd = {i:[] for i in range(nn)}
    for (i,j) in pairs:
        nbd[i].append(j)
        nbd[j].append(i)
        assert (j,i) not in pairs
    nbd2 = {i:set() for i in range(nn)}
    for i in range(nn):
      for j in nbd[i]:
        nbd2[i].update(nbd[j])
      nbd2[i].remove(i)
    #print(nbd2)

    dist = {(i,i):0 for i in range(nn)}
    done = False
    while not done:
        done = True
        keys = list(dist.keys())
        for (i1, i2) in keys:
            d0 = dist[i1, i2]
            for i0 in nbd[i1]:
                d = dist.get((i0, i2))
                if d is None or d > d0+1:
                    dist[i0, i2] = d0+1
                    dist[i2, i0] = d0+1
                    done = False
    assert len(dist) == nn*nn
    assert max(dist.values()) == 6

    def geodesic(src, tgt):
        tree = {src:src} # root
        bdy = [src]
        while tgt not in tree:
            assert bdy
            _bdy = []
            for idx in bdy:
                for jdx in nbd[idx]:
                    if jdx in tree:
                        continue
                    tree[jdx] = idx
                    _bdy.append(jdx)
            bdy = _bdy
        path = [tgt]
        while path[-1] != src:
            path.append(tree[path[-1]])
        return list(reversed(path))

    radius = 20
    layout = {}

    reds = list(range(n))
    blues = [i+n for i in range(n)]

    idx, jdx = 5, 36
    print(dist[idx, jdx])
    for kdx in nbd[idx]:
        print(dist[kdx, jdx])
    geo = geodesic(idx, jdx)
    for k in blues:
        if dist[k,jdx]==1 and k not in geo:
            break
    else:
        assert 0
    other = geodesic(k, idx)
    verts = geo + other[:-1]
    print(verts)

    assert len(verts) == 12

    for i, idx in enumerate(verts):
        theta = 2*i*pi/len(verts)
        x, y = radius*sin(theta), radius*cos(theta)
        layout[idx] = x, y

    for i in range(3):
        src, tgt = verts[2*i], verts[(2*i + 6)%12]
        for jdx in nbd[src]:
            if jdx not in layout:
                break
        else:
            assert 0
        geo = geodesic(jdx, tgt)[:-1]
        x0, y0 = layout[src]
        x1, y1 = layout[tgt]
        for i, idx in enumerate(geo):
            x = conv((i+1)/(len(geo)+1), x0, x1)
            y = conv((i+1)/(len(geo)+1), y0, y1)
            layout[idx] = x, y

    print("layout:", len(layout), "of", nn)
    cvs = Canvas()

    for (idx, jdx) in pairs:
        if idx in layout and jdx in layout:
            x0, y0 = layout[idx]
            x1, y1 = layout[jdx]
            cvs.stroke(path.line(x0, y0, x1, y1), [black.alpha(0.5)]+st_THICk)

    for idx in range(n):
        if idx in layout:
            x, y = layout[idx]
            cvs.fill(path.circle(x, y, 0.2), [red.alpha(1.0)]+st_THICk)
        if n+idx in layout:
            x, y = layout[n+idx]
            cvs.fill(path.circle(x, y, 0.2), [blue.alpha(1.0)]+st_THICk)

    cvs.writePDFfile("G2_render.pdf")


if __name__ == "__main__":

    from argv import argv
    name = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%name)
    else:
        print("%s()"%(name,))
        fn = eval(name)
        fn()

    print("OK: %.3f seconds\n"%(time() - start_time))






