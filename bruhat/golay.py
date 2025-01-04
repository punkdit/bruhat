#!/usr/bin/env python

"""
See book by Curtis (2024):
"The art of working with the Mathieu group M_24"

also:
"The Subgroups of M24, or How to Compute the Table of Marks of a Finite Group"
Gotz Pfeiffer

"""

from functools import reduce
from operator import add, mul

import numpy

import pynauty

from bruhat.argv import argv
from bruhat.solve import parse, dot2, enum2, array2, solve, zeros2, shortstr
from bruhat.gset import Perm, Group, Coset, mulclose, GL
from bruhat.smap import SMap

from huygens import config
config(text="pdflatex", latex_header=r"""
\usepackage{amsmath}
\usepackage{amssymb}
""")
from huygens.namespace import *

class Set: # copied from species.py
    def __init__(self, items=[]):
        if type(items) is int:
            assert items <= len(letters)
            items = letters[:items]
        #if items:
        #    tp = type(items[0])
        #    for item in items:
        #        assert type(item) is tp

        items = list(items)
        try:
            items.sort() 
        except TypeError:
            items.sort(key = str) 
        except AssertionError:
            items.sort(key = str) # wup

        self.items = tuple(items)
        self.set_items = set(items)
        assert len(self.set_items) == len(self.items), self.items

    def __str__(self):
        return "{%s}" % (', '.join(str(x) for x in self))

    def __repr__(self):
        return "Set(%s)"%(str(list(self.items)))

    @classmethod
    def promote(cls, items):
        if isinstance(items, Set):
            return items
        #if type(items) is GeneratorType:
        return Set(items)

    def __iter__(self):
        return iter(self.items)

    def __len__(self):
        return len(self.items)

    def __contains__(self, item):
        return item in self.set_items

    def __eq__(a, b):
        assert isinstance(b, Set), repr(b)
        return a.items == b.items

    def __ne__(a, b):
        assert isinstance(b, Set), repr(b)
        return a.items != b.items

    def __lt__(a, b):
        assert isinstance(b, Set), repr(b)
        return a.items < b.items

    def __le__(a, b):
        assert isinstance(b, Set), repr(b)
        return a.items < b.items

    def __hash__(self):
        return hash(self.items)

    def __add__(a, b):
        # categorical coproduct
        items = [(0, ai) for ai in a] + [(1, bi) for bi in b]
        return Set(items)

    def union(a, b):
        items = a.set_items.union(b.set_items)
        return Set(items)

    def __mul__(a, b):
        # categorical product
        items = [(ai, bi) for ai in a for bi in b]
        return Set(items)


class Struct:
    def __init__(self, n, partss):
        self.n = n
        self.partss
        

empty = Set()
star = Set(["*"])


class Vec:
    def __init__(self, *idxs):
        idxs = list(idxs)
        idxs.sort()
        v = zeros2(N)
        v[idxs] = 1
        self.v = v
        self.idxs = tuple(idxs)
        self.key = str(v)

    def __str__(self):
        return "Vec%s"%(self.idxs,)
    __repr__ = __str__

    def __eq__(self, other):
        return self.key == other.key

    def __hash__(self):
        return hash(self.key)

    def __add__(self, other):
        v = (self.v + other.v)%2
        return Vec(*[i for i in range(N) if v[i]])

    def __mul__(self, other):
        v = self.v * other.v
        return Vec(*[i for i in range(N) if v[i]])

    def __rmul__(self, g):
        idxs = [g[i] for i in self.idxs]
        return Vec(*idxs)

    def __lt__(self, other):
        return self.idxs < other.idxs

    #def __getitem__(self, i):
    #    return self.v[i]

    def __getitem__(self, i):
        return self.idxs[i]

    def __contains__(self, idx):
        return idx in self.idxs

    def sum(self):
        return self.v.sum()

    @classmethod
    def get_rect(self, idx, w, h, sep=0.1):
        dx = (w-2*sep)/6
        dy = h/4
        r, c = coords[idx]
        if c > 1:
            c += sep/dx
        if c > 4:
            c += sep/dx
        rect = (c*dx, -(r+1)*dy, dx, dy)
        return rect

    def render(self, w, h, sep=0.1, m=0.15, bg_fill=None):
        cvs = Canvas()
        p = path.rect(-m, -h-m, w + 2*m, h + 2*m)
        left = (self*bricks[0]).sum()
        if bg_fill is None:
            cl = {8:grey, 2:green.alpha(0.4), 4:blue.alpha(0.4), 0:red.alpha(0.4)}[left]
        else:
            cl = bg_fill
        cvs.fill(p, [cl])
        st = [black]
        for j, (r,c) in coords.items():
            p = path.rect(*self.get_rect(j, w, h, sep))
            if j in self:
                cvs.fill(p, st)
            else:
                cvs.fill(p, [white])
            cvs.stroke(p, [grey])
        return cvs



N = 24
infty = 23
table = [
    infty,  0, 11,  1, 22,  2,
        3, 19,  4, 20, 18, 10,
        6, 15, 16, 14,  8, 17,
        9,  5, 13, 21, 12,  7]
assert len(set(table)) == N

coords = {idx:(i//6, i%6) for (i,idx) in enumerate(table)}

terns = [
    [infty, 11, 22],
    [0, 20, 8],
    [15, 16, 18],
    [5, 4, 2],
    [9, 1, 17],
    [19, 14, 7],
    [6, 21, 12],
    [3, 13, 10]
]

items = reduce(add, terns)
assert len(set(items)) == 24


def render(vecs):
    print("render", len(vecs))

    # 759 == 3*11*23
    assert len(vecs) == 759

    cvs = Canvas()
    
    dx = dy = 0.5
    for i, v in enumerate(vecs):
        x = 0
        y = -i*2*dy
        for j, vj in enumerate(v):
            x = j*dx
            p = path.rect(x, y, dx, dy)
            if vj:
                cvs.fill(p)
            else:
                cvs.stroke(p)

    cvs.writePDFfile("golay.pdf")



bricks = [
    Vec(*(table[col + j + 6*row] for row in range(4) for j in [0,1]))
    for col in [0,2,4]
]


def mk_plot(octads):
    left = bricks[0]
    lookup = {o*left:[] for o in octads}
    for o in octads:
        lookup[o*left].append(o)
    print(len(lookup))
    keys = list(lookup.keys())
    keys.sort(key = lambda v:-v.sum())
    print(set(v.sum() for v in keys))

    for v in keys:
        lookup[v].sort()

    assert left in octads
    rows = [[left]] # weight 8

    items = [v for v in keys if v.sum() == 4] # weight 4
    idx = 0
    while idx < len(items):
        v = items[idx]
        w = left + v
        assert w in items
        rows.append(lookup[v] + lookup[w])
        items.remove(w)
        idx += 1

    items = [v for v in keys if v.sum() == 2] # weight 2
    for v in items:
        rows.append(lookup[v])

    items = [v for v in keys if v.sum() == 0] # weight 0
    for v in items:
        rows.append(lookup[v])

    w = 2.0
    h = 1.2
    cvs = Canvas()
    for i,row in enumerate(rows):
        y = -i * 1.7 * h
        for j,octad in enumerate(row):
            x = j*1.45*w
            cvs.insert(x, y, octad.render(w, h))

    cvs.writePDFfile("octads.pdf")

        
def mk_plot_1(octads):

    COLS = 16
    w = 2.0
    h = 1.2
    row = col = 0
    count = 0
    cvs = Canvas()
    for i, octad in enumerate(octads):
        x0 = col * 1.4 * w
        y0 = -row * 1.7 * h
        fg = octad.render(w, h)
        cvs.insert(x0, y0, fg)

        if count in [0, 280, 759-31]:
            row += 1
            col = 0
        else:
            col += 1
            if col == COLS:
                col = 0
                row += 1
        count += 1

    cvs.writePDFfile("octads.pdf")
    

# We use the same M24 that gap uses:
"""
gap> m24 := MathieuGroup(24);
Group([ (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23), (3,17,10,7,9)
  (4,13,14,19,5)(8,18,11,12,23)(15,20,22,21,16), (1,24)(2,23)(3,12)(4,16)(5,18)(6,10)(7,20)
  (8,14)(9,21)(11,17)(13,22)(15,19) ])
gap> (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23) in m24;
true
gap> (2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23) in m24;
false
gap> (24, 1)(12, 2)(4, 20)(5, 21)(7, 16)(17, 15)(10, 6)(14, 22) in m24;
true
gap> MaximalSubgroupClassReps(m24);
[ Group([ (1,13)(4,8)(6,17)(7,10)(12,19)(14,22)(15,16)(21,23), (1,24,22,10)(2,5,6,4)(3,14,9,12)(7,23,17,11)(8,15)(13,20) ]),
  Group([ (1,12)(2,11)(3,14)(4,16)(6,19)(8,24)(15,20)(18,23), (1,5,17,4)(2,7,6,22)(3,19,12,15)(8,9,21,14)(10,16)(11,20)(13,24)(18,23) ]),
  Group([ (1,9,10)(3,8,18)(4,21,6)(5,12,11)(7,19,22)(15,23,20), (1,20,2,10,3,9,16,14,19,22,23,7)(4,11,17,6,12,13)(5,21)(8,24,15,18) ]),
  Group([ (1,3,10,12,22,9,24,14)(2,17,4,23,6,7,5,11)(8,20,15,13)(19,21), (1,7)(2,8)(3,16)(4,22)(5,13)(6,10)(9,19)(11,15)(12,23)(14,17)(18,20)(21,24) ]),
  Group([ (1,4)(2,7)(3,8)(5,6)(9,13)(10,16)(11,19)(12,21)(14,20)(15,24)(17,22)(18,23), (1,2,22,8,10,23)(3,20,14,6,24,15,18,13,9,7,19,11)(4,12,17,16)(5,21) ]),
  Group([ (1,16,10,13,9,4,5,19,2,7,17,22,11,21)(3,23,24,15,14,20,18)(6,12), (1,16,13,15,4,17)(2,7,5,10,11,14)(3,22)(6,8)(9,18,24)(19,21,23) ]),
  Group([ (1,4)(2,7)(3,8)(5,6)(9,13)(10,16)(11,19)(12,21)(14,20)(15,24)(17,22)(18,23), (1,17,11,12)(2,10,5,9)(3,16,13,23)(4,7,8,20)(6,19)(14,22)(15,21)(18,24) ]),
  Group([ (1,3)(2,13)(4,7)(5,6)(8,14)(9,20)(10,17)(11,16)(12,21)(15,22)(18,19)(23,24), (1,13)(2,11)(3,7)(4,22)(5,21)(6,10)(8,12)(9,16)(14,20)(15,18)(17,23)(19,24), (1,11)(2,12)(3,8)(4,10)(5,19)(6,9)(7,21)(13,22)(14,23)(15,18)(16,17)(20,24) ]),
  Group([ (1,15)(2,19)(3,22)(4,12)(5,8)(6,14)(7,21)(9,17)(10,20)(11,23)(13,18)(16,24), (1,23,21,15)(2,8,13,20)(3,12,24,22)(4,6,16,19)(5,17,9,10)(7,14,11,18) ]) ]
"""

def m24_gen():
    g = Perm([(i+1)%infty for i in range(23)] + [infty])
    swaps = [
        (infty, 0), (11,1), (3, 19), (4, 20), 
        (6, 15), (16, 14), (9, 5), (13, 21)]

    #gapstr = (str([(i+1,j+1) for i,j in swaps]).replace(", (", "("))

    idxs = list(range(N))
    for (i,j) in swaps:
        idxs[i], idxs[j] = idxs[j], idxs[i]
    h = Perm(idxs)
    gen = [g,h]
    return gen


def main():
    gen = m24_gen()
    #print(gen)
    return

    #M24 = mulclose([g,h], verbose=True) # nah

    octad = Vec(infty,19,15,5,11,1,22,2)
    octads = {octad}
    bdy = list(octads)
    while bdy:
        _bdy = []
        for g in gen:
            for octad in bdy:
                o = g * octad
                if o not in octads:
                    octads.add(o)
                    _bdy.append(o)
        bdy = _bdy
    assert len(octads) == 759
    octads = list(octads)
    octads.sort()
    #for octad in octads[:10]:
        #print(octad)

    def find(*items):
        o = None
        for octad in octads:
            for i in items:
                if i not in octad:
                    break
            else:
                assert o is None
                o = octad
        return o

    the_octad = Vec(1, 2, 3, 4, 5, 8, 11, 13)
    assert find(1,2,3,4,5) == the_octad

    # idx:(row,col)
    coords = {idx:(i//6, i%6) for (i,idx) in enumerate(table)}
    def vstr(vec):
        smap = SMap()
        for i in range(N):
            smap[coords[i]] = '*' if i in vec.idxs else '.'
        return smap

    the_octad = bricks[0]

    def o_key(v):
        vs = [brick*v for brick in bricks]
        ss = [v.sum() for v in vs]
        return tuple(zip(ss, vs))

    octads.sort(key = o_key, reverse=True)

    #mk_plot(octads)

    # the miracle octad generator
    str_mog = """
    11...1.1 11..1.1. 1.1.1.1. 1.1..1.1 1..1..11 1..111.. 
             11.1.1.. 111.1... 1..11..1 1.11.1.. 1....111 
    1..1.1.1 111...1. 11.1...1 1..1.11. 1...11.1 1.11...1
    1..11.1. 1.1.11.. .1.111.. 1111.... 11...11. 11..1..1
    1.1..11. 1...1.11 1.111... 11..11.. 111....1 11.1..1.
    1.1.1..1 1.11..1. 1...111. 11....11 11.11... 111..1..
    """.strip().split()
    assert len(str_mog) == 35

    mog = [parse(v) for v in str_mog]
    items = set()
    for v in list(mog):
        assert v.shape == (1,8)
        u = 1 - v
        items.add(shortstr(v))
        items.add(shortstr(u))
    assert len(items) == 70
    items = list(items)
    items.sort()
    from bruhat.util import choose
    idxs = list(range(8))
    for sub in choose(idxs, 4):
        v = ['.']*8
        for i in sub:
            v[i] = '1'
        assert ''.join(v) in items

    idxs = [infty, 0, 3, 19, 6, 15, 9, 5]
    idx_mog = [Vec(*[idxs[i] for i in range(8) if v[i]=='1']) for v in str_mog]
    mog = {Vec(*[idxs[i] for i in range(8) if v[i]=='1']):[] for v in str_mog}

    #lookup = {u : [] for u in mog}
    for octad in octads:
        #if octad == the_octad:
        #    continue
        u = octad * the_octad
        if u.sum() != 4:
            continue
        v = the_octad + u
        if u in mog:
            mog[u].append(octad)
        elif v in mog:
            #mog[v].append(octad)
            pass # don't need these
        else:
            assert 0

    for items in mog.values():
        assert len(items) == 4

    w, h = 3.2, 1.8
    bg_fill = (0.5*blue + 0.6*white) + 0.1*green
    print(bg_fill)
    cvs = Canvas()
    #for idx,v in enumerate(idx_mog):
    pos = {}
    for count in range(36):
        row = count // 6
        col = count % 6
        x,y = 1.5*w*col, -1.7*h*row
        if count < 6:
            idx = count
        elif count == 6:
            the_pos = x, y
            continue # see below
        else:
            idx = count - 1
        pos[idx] = x, y
        v = idx_mog[idx]
        cvs.insert(x, y, v.render(w, h, bg_fill=bg_fill))

    dx = w/6
    dy = h/4

    bg_fill = 0.2*red + 0.8*white
    print(bg_fill)
    x, y = the_pos
    cvs.insert(x, y, Vec().render(w, h, bg_fill=bg_fill))
    #cvs.stroke(path.rect(x, y-h, w, h))
    #cvs.text(x, y-h, r"$\infty$")
    for idx,val in coords.items():
        r, c = val
        s = str(idx) if idx != 23 else r"$\infty$"
        x1, y1 = (x + c*dx + 0.08, y - dy*(r+1) + 0.1)
        if len(s) == 1:
            x1 += 0.15
        #rect = Vec.get_rect(idx, w, h)
        #x1, y1 = rect[:2]
        cvs.text(x1, y1, s)

    sep = 0.1
    dx = (w-2*sep)/6
    dy = h/4

    mask = bricks[1] + bricks[2]
    cls = [black, white, blue, green]
    for (u,items) in mog.items():
        idx = idx_mog.index(u)
        #print(idx, pos[idx])
        x, y = pos[idx]
        assert len(items) == 4
        items.sort(key = lambda v:11 not in v)
        for i,v in enumerate(items):
            v = v*mask
            #print("\t", v.idxs)
            for idx in v.idxs:
                r = Vec.get_rect(idx, w, h)
                x1, y1 = r[:2]
                #cvs.text(x+x1, y+y1, str(idx))
                mx, my = x+x1+0.5*dx, y+y1+0.5*dy
                if i==0:
                    cvs.fill(path.rect(x+x1, y+y1, 0.96*dx, 0.96*dy))
                elif i==2:
                    cvs.fill(path.circle(mx, my, 0.25*dx))
                elif i==3:
                    cvs.stroke(path.circle(mx, my, 0.25*dx))


    cvs.writePDFfile("mog.pdf")


def get_hecke():

    def pstr(pair):
        smap = SMap()
        smap[0,0] = vstr(pair[0])
        smap[0,8] = vstr(pair[1])
        return smap

    pairs = [(a,b) for a in octads for b in octads]
    print("pairs:", len(pairs))
    orbits = []
    remain = set(pairs)
    found = set()
    while remain:
        pair = iter(remain).__next__()
        remain.remove(pair)
        found.add(pair)
        #print(pstr(pair), "OK\n")
        a, b = pair
        print(a*b, '\n')
        orbit = bdy = [pair]
        while bdy:
            _bdy = []
            for pair in bdy:
              for g in gen:
                other = (g*pair[0], g*pair[1])
                if other not in found:
                    _bdy.append(other)
                    found.add(other)
                    remain.remove(other)
            orbit += _bdy
            bdy = _bdy
        orbits.append(orbit)
        if len(orbit)==1:
            pair = orbit[0]
            print(pstr(pair), '\n')
            for g in gen:
                other = (g*pair[0], g*pair[1])
                print(pstr(other), '\n')
            return
        print(len(orbit), len(orbit)//len(octads))
    print()

    print([len(orbit) for orbit in orbits])

    # 

def from_gapstr(s, n):
    items = list(range(n))
    flds = s.replace("(", " ").replace(")", "")
    flds = flds.split()
    flds = [[int(i)-1 for i in fld.split(',')] for fld in flds]
    #print("from_gapstr", flds)

    perm = {i:i for i in range(n)}
    for cycle in flds:
        m = len(cycle)
        for i in range(m):
            perm[cycle[i]] = cycle[(i+1)%m]
    #print(perm)
    assert len(perm) == n
    assert len(set(perm.values())) == n
    perm = [perm[i] for i in items]
    return Perm(perm)
        
        
def slow_find_orbit(gen, H):
    found = {H}
    bdy = list(found)
    while bdy:
        _bdy = []
        for g in gen:
            for H in bdy:
                K = H.left_mul(g)
                if K in found:
                    continue
                found.add(K)
                _bdy.append(K)
        bdy = _bdy
        print(len(found))
    return found


def small_index_find_orbit(H, gen):
    found = [H.identity]
    bdy = list(found)
    while bdy:
        _bdy = []
        for g in gen:
            for h in bdy:
                j = g*h
                jinv = ~j
                for f in found:
                    if jinv*f in H:
                        break
                else:
                    found.append(j)
                    _bdy.append(j)
        bdy = _bdy
        print(len(found))
    return found



def main_maximal():

    desc = """
    (1,13)(4,8)(6,17)(7,10)(12,19)(14,22)(15,16)(21,23), (1,24,22,10)(2,5,6,4)(3,14,9,12)(7,23,17,11)(8,15)(13,20)
    (1,12)(2,11)(3,14)(4,16)(6,19)(8,24)(15,20)(18,23), (1,5,17,4)(2,7,6,22)(3,19,12,15)(8,9,21,14)(10,16)(11,20)(13,24)(18,23)
    (1,9,10)(3,8,18)(4,21,6)(5,12,11)(7,19,22)(15,23,20), (1,20,2,10,3,9,16,14,19,22,23,7)(4,11,17,6,12,13)(5,21)(8,24,15,18)
    (1,3,10,12,22,9,24,14)(2,17,4,23,6,7,5,11)(8,20,15,13)(19,21), (1,7)(2,8)(3,16)(4,22)(5,13)(6,10)(9,19)(11,15)(12,23)(14,17)(18,20)(21,24)
    (1,4)(2,7)(3,8)(5,6)(9,13)(10,16)(11,19)(12,21)(14,20)(15,24)(17,22)(18,23), (1,2,22,8,10,23)(3,20,14,6,24,15,18,13,9,7,19,11)(4,12,17,16)(5,21)
    (1,16,10,13,9,4,5,19,2,7,17,22,11,21)(3,23,24,15,14,20,18)(6,12), (1,16,13,15,4,17)(2,7,5,10,11,14)(3,22)(6,8)(9,18,24)(19,21,23)
    (1,4)(2,7)(3,8)(5,6)(9,13)(10,16)(11,19)(12,21)(14,20)(15,24)(17,22)(18,23), (1,17,11,12)(2,10,5,9)(3,16,13,23)(4,7,8,20)(6,19)(14,22)(15,21)(18,24)
    (1,3)(2,13)(4,7)(5,6)(8,14)(9,20)(10,17)(11,16)(12,21)(15,22)(18,19)(23,24), (1,13)(2,11)(3,7)(4,22)(5,21)(6,10)(8,12)(9,16)(14,20)(15,18)(17,23)(19,24), (1,11)(2,12)(3,8)(4,10)(5,19)(6,9)(7,21)(13,22)(14,23)(15,18)(16,17)(20,24)
    (1,15)(2,19)(3,22)(4,12)(5,8)(6,14)(7,21)(9,17)(10,20)(11,23)(13,18)(16,24), (1,23,21,15)(2,8,13,20)(3,12,24,22)(4,6,16,19)(5,17,9,10)(7,14,11,18)
    """.strip().split("\n")

    names = "m23 n22 ea8 n12 e3s6 n21 nea l23 l7".split()
    gens = {}
    for idx,line in enumerate(desc):
        line = line.strip()
        gen = line.split(", ")
        gen = [from_gapstr(s, 24) for s in gen]
        #G = Group(gen)
        #G = mulclose(gen, verbose=True) 
        gens[names[idx]] = gen

    octern = Group.generate(gens["l7"])
    assert len(octern) == 168

    H = octern
    #H = Group.generate(gens["n21"])
    print(len(H))

    gen = m24_gen()

    H = Coset(H)
    found = {H[0]}
    bdy = list(found)
    while bdy:
        _bdy = []
        for g in gen:
            for h in bdy:
                gh = g*h
                ghH = H.left_mul(gh)
                j = ghH[0]
                if j not in found:
                    found.add(j)
                    _bdy.append(j)
        bdy = _bdy
        print(len(found))

    return



    

def test_fano():
    G = GL(3, 2)
    assert len(G) == 168

    #Hs = list(G.subgroups(verbose=True))
    #print(len(Hs))

    found = []
    for g in G:
        if g.order() == 7:
            found.append(g)
            print(g)
    print(len(found))
    return

    H = Group.generate([g])
    assert len(H) == 7

    X = G.action_subgroup(H)
    print(X)

    XX = X*X
    print(XX)

    orbits = XX.get_orbits()
    print([len(o) for o in orbits])




def main_golay():
    H = parse("""
    1...........11...111.1.1
    .1...........11...111.11
    ..1.........1111.11.1...
    ...1.........1111.11.1..
    ....1.........1111.11.1.
    .....1......11.11..11..1
    ......1......11.11..11.1
    .......1......11.11..111
    ........1...11.111...11.
    .........1..1.1.1..1.111
    ..........1.1..1..11111.
    ...........11...111.1.11
    """) # Golay code
    
    
    # RM [[16,6,4]]
    H1 = parse("""
    11111111........
    ....11111111....
    ........11111111
    11..11..11..11..
    .11..11..11..11.
    """)
    
    
    m, n = H.shape
    wenum = {i:[] for i in range(n+1)}
    span = []
    for u in enum2(m):
        v = dot2(u, H)
        d = v.sum()
        wenum[d].append(v)
        if d:
            span.append(v)
    
    print([len(wenum[i]) for i in range(n+1)])
    
    #for d in range(1, n+1):
    #    if wenum[d]:
    #        break
    #for v in wenum[d]:
    #    print(shortstr(v))

    #render(wenum[8])

    N = len(span)
    graph = pynauty.Graph(N+n)
    
    colours = {d:[] for d in range(n+1)}
    for idx, v in enumerate(span):
        #print(idx, v, v.sum())
        d = v.sum()
        colours[d].append(idx)
        for i in range(n):
            if v[i]:
                graph.connect_vertex(N+i, idx)
    
    labels = []
    for d in range(n+1):
        if colours[d]:
            labels.append(colours[d])
    labels.append(list(range(N, N+n)))
    print([len(lbl) for lbl in labels])
    
    labels = [set(l) for l in labels]
    
    #fix = labels[1].pop()
    #labels.insert(0, {fix})
    
    items = []
    for l in labels:
        items += list(l)
    items.sort()
    #print(items, N+n)
    assert items == list(range(N+n))
    print(N+n, "vertices")
    
    
    graph.set_vertex_coloring(labels)
    
    aut = pynauty.autgrp(graph)
    print(len(aut))
    gen = aut[0]
    order = int(aut[1])
    print("autos:", order)
    #for perm in gen:
    #    print(perm)
    
    for perm in gen:
        f = [perm[i] - N for i in range(N, N+n)]
        print(f)





if __name__ == "__main__":
    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds"%(time() - start_time))
    print()






