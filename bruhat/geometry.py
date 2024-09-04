#!/usr/bin/env python
"""
Use group theory to build surface codes, color codes, etc.

https://en.wikipedia.org/wiki/Bitruncated_cubic_honeycomb

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
    span, intersect, rank, enum2, shortstrx, identity2, eq2, pseudo_inverse, parse)
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
    if argv.idx:
        main(key, argv.idx)
        return
    print("key =", key)
    n = len(lins_db.db[key])
    idx = start_idx
    while idx < n and (stop_idx is None or idx < stop_idx):
        main(key, idx)
        idx += 1


def main(key, index):
    N = len(key)+1
    geometry = Geometry(key, index, False)
    total = geometry.build_graph()
    #total = total.compress()
    words = total.get_words()
    n = len(total)
    print("idx = %d, |G| = %d"%(index, n))

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
        #graph.dump()
        assert n % len(graph) == 0, len(graph)
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
        #print(g.fig, h.fig)
        #print(shortstr(A), A.shape)
        As[g.fig, h.fig] = A
      #print()

    if 0:
        edges.get_dot(colours, "edges_434")
        faces = geometry.build_graph(hgens = [a, b, d])
        faces.get_dot(colours, "faces_434")
        bodis = geometry.build_graph(hgens = [a, b, c])
        print(len(verts), len(edges), len(faces), len(bodis))

    keys = list(As.keys())
    keys.sort()


    if argv.show_incidence:
        for fig in figs:
          for gig in figs:
            B = dot2(As[fig, "0000"], As[gig, "0000"].transpose())
            print(fig, gig, end=" ")
            if B.sum() == 0: # very interesting !
                print(".", end="")
            else:
                print("X", end="")
            print()
        print()
        #return

    #from bruhat.qcode import QCode
    from qumba.csscode import CSSCode, distance_z3_css, distance_lower_bound_z3

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

    if 0:
        assert Hx[:, n-1].sum() == 0
        assert Hz[:, n-1].sum() == 0
    
        Hx = Hx[:, :n-1]
        Hz = Hz[:, :n-1]

    C = dot2(Hx, Hz.transpose())
    if C.sum():
        print("/", end="", flush=True)
        return

    code = CSSCode(Hx=Hx, Hz=Hz, check=True)
    #print(code.longstr())
    #print(shortstr(code.Lx))
    #print("-"*code.n)
    #print(shortstr(code.Lz))
    #print()
    if code.k == 0:
        return

    if distance_lower_bound_z3(code.Hz, code.Lz, 2) is not None:
        return
    if distance_lower_bound_z3(code.Hx, code.Lx, 2) is not None:
        return

    #print()
    #print("idx =", index)
    print("[[%d, %d, >2]]"%(code.n, code.k))
    if argv.distance:
        if argv.dual:
            d_z, d_x = distance_z3_css(code.get_dual(), True)
        else:
            d_x, d_z = distance_z3_css(code, True)
        print("[[%d, %d, (%d,%d)]]"%(code.n, code.k, d_x, d_z))

    if argv.show:
        print("Hx:")
        print(shortstr(Hx))
        print("Hz:")
        print(shortstr(Hz))

    #print("Hx", Hx.sum(0), Hx.sum(1))
    #print("Hz", Hz.sum(0), Hz.sum(1))

    dump_transverse(code.Hx, code.Lx, 3, argv.show_all)
    #print("dual:")
    #dump_transverse(code.Hz, code.Lz, 3)
    return

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


def dump_transverse(Hx, Lx, t=3, show_all=False):
    import CSSLO
    SX,LX,SZ,LZ = CSSLO.CSSCode(Hx, Lx)
    N = 1<<t
    zList, qList, V, K_M = CSSLO.comm_method(SX, LX, SZ, t, compact=True, debug=False)
    for z,q in zip(zList,qList):
        #print(z, q)
        lhs, rhs = CSSLO.CP2Str(2*q,V,N), CSSLO.z2Str(z,N)
        if "CCZ" in lhs or "T" in lhs or show_all:
            print(lhs, "=>", rhs)
    print()
    return zList


def test_colour():
    # https://arxiv.org/pdf/1304.3709

    s = """
    IIIIIIIXXXXXXXX
    IIIXXXXIIIIXXXX
    IXXIIXXIIXXIIXX
    XIXIXIXIXIXIXIX
    """
    Hx = parse(s.replace("I","0").replace("X","1"))
    print(Hx)

    s = """
    .......ZZZZZZZZ
    ...ZZZZ....ZZZZ
    .ZZ..ZZ..ZZ..ZZ
    Z.Z.Z.Z.Z.Z.Z.Z
    """
    Hz = parse(s.replace(".","0").replace("Z","1"))
    print(Hz)
    ops = list(Hz)
    for h in Hz:
      for j in Hz:
        if (h*j).sum() == 4:
            ops.append(h*j)
    ops = array2(ops)
    print(shortstr(ops), ops.shape)
    Hz = ops
    Hz = linear_independent(Hz)
    print("Hz:", len(Hz))
    
    from qumba.csscode import CSSCode, distance_z3_css
    code = CSSCode(Hx=Hx, Hz=Hz, check=True)
    print(code)
    #assert distance_z3_css(code) == (3,3)
    print( distance_z3_css(code) )

    dump_transverse(code.Hx, code.Lx, 3)
    ops = dump_transverse(code.Hx, code.Lx, 2)
    print(ops)

    op = ops[0]
    print(dot2(Hz, op))

    for i,x in enumerate(op):
        if x!=2:
            op[i] = 0
    print(op)
    print(dot2(Hx, op))


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


