#!/usr/bin/env python
"""
Use group theory to build surface codes, color codes, etc.
"""

import string, os
from random import randint, choice, random
from time import sleep, time
from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)
from operator import matmul
from string import ascii_lowercase

import numpy
from numpy import alltrue, zeros, dot

from bruhat.util import cross
from bruhat import solve 
from bruhat.solve import (
    array2, zeros2, shortstr, dot2, solve2, linear_independent, row_reduce, find_kernel,
    span, intersect, rank, enum2, shortstrx, identity2, eq2, pseudo_inverse)
from bruhat.action import Perm, Group, Coset, mulclose, close_hom, is_hom
from bruhat.todd_coxeter import Schreier
from bruhat.argv import argv
from bruhat.smap import SMap
from bruhat import lins_db


class Geometry(object):
    "A geometry specified by a Coxeter reflection group"
    def __init__(self, orders, lins_idx=0, build=True):
        ngens = len(orders)+1
        self.gens = {c:(i,) for (i,c) in enumerate(ascii_lowercase[:ngens])}
        orders = tuple(orders)
        assert orders in lins_db.db, (list(lins_db.db.keys()))
        #print("oeqc.Geometry.__init__:", len(lins_db.db[orders]))
        rels = lins_db.db[orders][lins_idx]
        rels = lins_db.parse(rels, **self.gens)
        self.orders = orders
        self.ngens = ngens
        self.rels = rels
        self.dim = len(orders)
        if build:
            self.G = self.build_group()

    def build_graph(self, figure=None, hgens=None):
        gens = [(i,) for i in range(5)]
        ngens = self.ngens
        rels = [gens[i]*2 for i in range(ngens)]
        orders = self.orders
        for i in range(ngens-1):
            order = orders[i]
            if order is not None:
                rels.append( (gens[i]+gens[i+1])*order )
            for j in range(i+2, ngens):
                rels.append( (gens[i]+gens[j])*2 )
        rels = rels + self.rels
        #print(rels)
        graph = Schreier(ngens, rels)
        if figure is not None:
            assert len(figure) == len(gens)
            gens = [gens[i] for i, fig in enumerate(figure) if fig] or [G.identity]
            graph.build(gens)
        elif hgens is not None:
            graph.build(hgens)
        else:
            graph.build()
        return graph

    def build_group(self):
        graph = self.build_graph()
        G = graph.get_group()
        return G

    def get_cosets(self, figure):
        G = self.G
        gens = G.gens
        assert len(figure) == len(gens)
        gens = [gens[i] for i, fig in enumerate(figure) if fig] or [G.identity]
        #print("gens:", gens)
        H = Group.generate(gens)
        #pairs = G.left_cosets(H)
        cosets = G.left_cosets(H)
        return cosets


def test():
    start_idx = argv.get("start_idx", 1)
    stop_idx = argv.get("stop_idx", None)
    key = argv.get("key", (4,3,4))
    index = argv.get("index", 1000)
    if key not in lins_db.db:
        print("lins_db.build_db...", end='', flush=True)
        lins_db.build_db(key, index)
        print(" done")
    n = len(lins_db.db[key])
    idx = start_idx
    while idx < n and (stop_idx is None or idx < stop_idx):
        main(key, idx)
        idx += 1


def main(key, idx):
    print("idx =", idx)
    N = len(key)+1
    geometry = Geometry(key, idx, False)
    total = geometry.build_graph()
    #total = total.compress()
    words = total.get_words()
    n = len(total)
    #print("|G| =", n)

    a, b, c, d = [(i,) for i in range(N)]

    colours = "red blue darkgreen purple".split()

    figs = """
    0000 0001 0010 0011 
    0100 0101 0110 0111 
    1000 1001 1010 1011 
    1100 1101 1110 1111 
    """.strip().split()

    gens = [a,b,c,d]
    graphs = []
    for fig in figs:
        hgens = [gens[i] for i in range(N) if fig[i]=='1']
        graph = geometry.build_graph(hgens=hgens).compress()
        graph.fig = fig
        #if len(graph) > 1:
        graphs.append(graph)
    graphs.sort(key = len)
    k = len(graphs)
    names = ascii_lowercase[:k]
    for g,name in zip(graphs, names):
        g.name = name
        #print(name, len(g), g.fig)
    
    As = {}
    for i in range(k):
      g = graphs[i]
      for j in range(k):
        h = graphs[j]
        A = zeros2(len(g), len(h))
        for w in words:
            idx = g.follow_path(0,w)
            jdx = h.follow_path(0,w)
            A[idx,jdx] = 1
        key = (g.name, h.name)
        #print("%s.%s(%d,%d)"%(g.name, h.name, A.sum(0)[0], A.sum(1)[0]), end=" ")
        As[g.fig, h.fig] = A
      #print()

    if 0:
        edges.get_dot(colours, "edges_434")
        faces = geometry.build_graph(hgens = [a, b, d])
        faces.get_dot(colours, "faces_434")
        bodis = geometry.build_graph(hgens = [a, b, c])
        print(len(verts), len(edges), len(faces), len(bodis))

    #from bruhat.qcode import QCode
    from qumba.csscode import CSSCode, distance_z3_css

    Hx = numpy.concatenate((
        As["0111", "0000"],
        As["1011", "0000"],
        As["1101", "0000"],
        As["1110", "0000"],
    ))
    #print(shortstr(Hx))
    #print(Hx.shape)
    #print(Hx.sum(0), Hx.sum(1))

    mx, n = Hx.shape
    Hz = numpy.concatenate((
        As["0011", "0000"],
        As["0110", "0000"],
        As["0101", "0000"],
        As["1001", "0000"],
        As["1010", "0000"],
        As["1100", "0000"],
    ))
    #print("bodis:", Hx.shape)
    #print("faces:", Hz.shape)
    Hx = linear_independent(Hx)
    Hz = linear_independent(Hz)
    #print("bodis:", Hx.shape)
    #print("faces:", Hz.shape)
    #code = QCode.build_css(Hx, Hz, check=True)

    C = dot2(Hx, Hz.transpose())
    if C.sum():
        #print("non-commutative")
        return

    code = CSSCode(Hx=Hx, Hz=Hz, check=True)
    if code.k == 0:
        return

    if argv.distance:
        d_x, d_z = distance_z3_css(code, code.n>40)
        if d_x <= 2 or d_z <= 2:
            return
        print("distance:", d_x, d_z)
    print(code)
    if argv.show:
        print("Hx:")
        print(shortstr(Hx))
        print("Hz:")
        print(shortstr(Hz))

    ops = dump_transverse(code.Hx, code.Lx, 2)

    for op in ops:
        if '1' not in str(op):
            continue
        #print()
        #print(op)
        op = op%2
        #for (i,x) in enumerate(op):
        #    if x>1:
        #        op[i] = 1
        #print(op)
        for (key,A) in As.items():
            if key[1] != "0000" or key[0] == "0000":
                continue
            print(key[0], shortstr(dot2(A, op)))


def dump_transverse(Hx, Lx, t=3):
    import CSSLO
    SX,LX,SZ,LZ = CSSLO.CSSCode(Hx, Lx)
    N = 1<<t
    zList, qList, V, K_M = CSSLO.comm_method(SX, LX, SZ, t, compact=True, debug=False)
    for z,q in zip(zList,qList):
        #print(z, q)
        print(CSSLO.CP2Str(2*q,V,N),"=>",CSSLO.z2Str(z,N))
    print()
    return zList

    
    
if __name__ == "__main__":

    from time import time

    start_time = time()


    profile = argv.profile
    name = argv.next()
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        test()


    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")


