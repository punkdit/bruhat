#!/usr/bin/env python

"""
See book by Curtis (2024):
"The art of working with the mathieu group M_24"

"""

import numpy

import pynauty

from bruhat.argv import argv
from bruhat.solve import parse, dot2, enum2, array2, solve, zeros2
from bruhat.gset import Perm, Group, mulclose
from bruhat.smap import SMap

from huygens.namespace import *

N = 24

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


class Vec(object):
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
        return self.keys[i]
    def __contains__(self, idx):
        return idx in self.idxs


def main():

    inf = 23
    g = Perm([(i+1)%inf for i in range(23)] + [inf])
    swaps = [
        (inf, 0), (11,1), (3, 19), (4, 20), 
        (6, 15), (16, 14), (9, 5), (13, 21)]
    idxs = list(range(N))
    for (i,j) in swaps:
        idxs[i], idxs[j] = idxs[j], idxs[i]
    h = Perm(idxs)
    gen = [g,h]

    #M24 = mulclose([g,h], verbose=True) # nah

    octad = Vec(inf,19,15,5,11,1,22,2)
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

    table = [
        inf, 0, 11, 1, 22, 2,
        3, 19, 4, 20, 18, 10,
        6, 15, 16, 14, 8, 17,
        9, 5, 13, 21, 12, 7]
    assert len(set(table)) == N

    bricks = [
        Vec(*(table[col + j + 6*row] for row in range(4) for j in [0,1]))
        for col in [0,2,4]
    ]
    print(bricks)
        
    # idx:(row,col)
    coords = {idx:(i//6, i%6) for (i,idx) in enumerate(table)}
    def vstr(vec):
        smap = SMap()
        for i in range(N):
            smap[coords[i]] = '*' if i in vec.idxs else '.'
        return smap

    print(vstr(bricks[0]))

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

    return

    cvs = Canvas()
    dx = dy = 0.3
    for i,o in enumerate(octads):
        row = i//23
        col = i%23
        x0 = col * 8 * dx
        y0 = -row * 6 * dy
        st = [black]
        #if o == the_octad:
        #    print(o)
        #    st = [orange]
        for j, (r,c) in coords.items():
            p = path.rect(x0+c*dx, y0-r*dy, dx, dy)
            cvs.stroke(p)
            if j in o:
                cvs.fill(p, st)

    cvs.writePDFfile("octads.pdf")
    

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







