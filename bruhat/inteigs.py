#!/usr/bin/env python2

"""
search for integer eigenvectors of graphs.
"""

from __future__ import print_function

import sys, os
from random import gauss, random, seed, shuffle, choice
from math import sin, cos, pi

#import numpy


import networkx as nx
from networkx.generators import small, classic


from equiv import get_autos
from visualize import pos_circ
from argv import argv


if argv.seed is not None:
    seed(argv.seed)


def search(nodes, nbd, eigval, vals=[1, -1, 2, -2, 3, -3], vec=None, verbose=False):
    assert nodes

    if vec is None:
        node = nodes[0]

        for val in vals:
            if val >= 0:
                vec = {node : val}
                for vec in search(nodes, nbd, eigval, vals, vec, verbose): # recurse
                    if sum(abs(x) for x in vec.values()) != 0:
                        yield vec
        return

    vec = dict(vec)

    indent = "  "*len(vec)
    if verbose:
        print(indent+"==============")
        print(indent+"vec:", vec)

    satisfied = set()
    bdy = set(vec.keys())
    for i in vec:
        for j in nbd[i]:
            if j not in vec:
                break
        else:
            x = sum(vec[j] for j in nbd[i])
            if x != vec[i]*eigval:
                if verbose:
                    print(indent+"Not an evec")
                return 
            satisfied.add(i)
            bdy.remove(i)

    if verbose:
        print(indent+"satisfied:", satisfied)
        print(indent+"bdy:", bdy)

    if not bdy:
        yield vec
        return

    i = list(bdy)[0] # choose a component to satisfy
    
    remain = [j for j in nbd[i] if j not in vec]
    if len(remain)==1:
        # there's only one value that might work
        x = eigval * vec[i]
        k = None
        for j in nbd[i]:
            if j in vec:
                x -= vec[j]
            else:
                assert k is None
                k = j
        if verbose:
            print(indent+"force: vec[%d] = %d"%(k, x))
        if x not in vals:
            return
        vec[k] = x
        for vec1 in search(nodes, nbd, eigval, vals, vec, verbose): # <--- recurse
            yield vec1

    else:
        for j in nbd[i]:
            if j not in vec:
                break
        else:
            assert 0

        for val in vals:
            vec[j] = val
            for vec1 in search(nodes, nbd, eigval, vals, vec, verbose): # <--- recurse
                yield vec1
            del vec[j]
    

"""
bull_graph chvatal_graph complete_graph cubical_graph
cycle_graph desargues_graph diamond_graph dodecahedral_graph
empty_graph frucht_graph heawood_graph house_graph house_x_graph
icosahedral_graph krackhardt_kite_graph make_small_graph
make_small_undirected_graph moebius_kantor_graph octahedral_graph
pappus_graph path_graph petersen_graph sedgewick_maze_graph
tetrahedral_graph truncated_cube_graph truncated_tetrahedron_graph
tutte_graph
"""

def draw_graph(graph, vec, pts=None, name="output"):
    import pyx 
    from pyx import path, deco, trafo, style, text, color, deformer
    from pyx.color import rgb, cmyk
    from pyx.color import rgbfromhexstring as rgbhex

    black = rgb(0., 0., 0.) 
    blue = rgb(0., 0., 0.8)
    lred = rgb(1., 0.4, 0.4)
    red = rgb(1., 0.0, 0.0)
    green = rgb(0., 1.0, 0.0)
    lgreen = rgb(0.4, 1.0, 0.4)
    dgreen = rgb(0.0, 0.4, 0.0)
    white = rgb(1., 1., 1.) 

    if pts is None:
        pts = pos_circ(graph)

    directed = isinstance(graph, (nx.DiGraph, nx.MultiDiGraph))

    #W = 10.
    #H = 10.

    R = 1.5 
    r = 0.2 

    c = pyx.canvas.canvas()

    for edge in graph.edges():
        src, tgt = edge

        x0, y0 = pts[src]
        x1, y1 = pts[tgt]

        color = black
        c.stroke(path.line(R*x0, R*y0, R*x1, R*y1), [color])

    for node in graph.nodes():
        x, y = pts[node]

        val = vec[node]
        assert int(val) == val

        #p = path.circle(R*x, R*y, 0.2)
        #c.fill(p, [white])

        p = path.circle(R*x, R*y, 0.10)
        c.fill(p, [black])

        if val > 0:
            color = dgreen
        elif val < 0:
            color = red
            val = -val
        else:
            continue

        for i in range(val):
            p = path.circle(R*x, R*y, (i+2)*0.15)
            c.stroke(p, [color, style.linewidth.THick])

    c.writePDFfile(name)


def strvec(vec):
    if vec is None:
        return "None"
    n = max(list(vec.keys()))
    items = [vec.get(i) for i in range(n+1)]
    return str(items)

def tutte_coxeter():
    edges = [
        (0, 1), (0, 17), (0, 29), (1, 2), (1, 22), 
        (2, 3), (2, 9), (3, 4), (3, 26), (4, 5), 
        (4, 13), (5, 6), (5, 18), (6, 7), (6, 23), 
        (7, 8), (7, 28), (8, 9), (8, 15), (9, 10), 
        (10, 11), (10, 19), (11, 12), (11, 24), (12, 13), 
        (12, 29), (13, 14), (14, 15), (14, 21), (15, 16), 
        (16, 17), (16, 25), (17, 18), (18, 19), (19, 20), 
        (20, 21), (20, 27), (21, 22), (22, 23), (23, 24), 
        (24, 25), (25, 26), (26, 27), (27, 28), (28, 29)]
    g = nx.Graph()
    for edge in edges:
        g.add_edge(*edge)
    return g


def cubic_twisted():
    edges = []
    for i in range(10):
        edges.append((i, (i+1)%10))
    edges.extend([(0, 5), (1, 3), (2, 6), (4, 8), (7, 9)])
    g = nx.Graph()
    for edge in edges:
        g.add_edge(*edge)
    return g


def g11_graph():
    edges = []
    for i in range(10):
        edges.append((i, (i+1)%10))
    edges.extend([(0, 5), (1, 8), (2, 9), (3, 6), (4, 7)])
    g = nx.Graph()
    for edge in edges:
        g.add_edge(*edge)
    return g


def g10_graph():
    tree = [
        (0, 1), (0, 2), (0, 3), 
        (1, 4), (1, 5), (2, 6), (2, 7), (3, 8), (3, 9)]
    edges = list(tree) + [(a+10, b+10) for (a, b) in tree]
    edges.extend([
        (4, 14), (4, 16), (5, 17), (5, 18), (6, 14), (6, 18),
        (7, 15), (7, 19), (8, 15), (8, 16), (9, 17), (9, 19)])
    g = nx.Graph()
    for edge in edges:
        g.add_edge(*edge)
    return g


def main():

    name = argv.next()

    builder = getattr(small, name, None)

    if builder is None:
        builder = getattr(classic, name, None)

    if builder is None:
        builder = globals().get(name)

    if builder is None:
        print(name, "not found")
        return

    n = argv.n
    if n is not None:
        graph = builder(n)
    else:
        graph = builder()

    print("|nodes|=%d, edges=|%d|"%(len(graph.nodes()), len(graph.edges())))

    eigval = argv.get("eigval", 0)

    nodes = list(graph.nodes())
    edges = list(graph.edges())
    nbd = dict((i, []) for i in nodes)
    for i, j in edges:
        nbd[i].append(j)
        nbd[j].append(i)
    print(edges)
    print(nbd)

    verbose = argv.verbose
    G = get_autos(graph)
    print("|G|", len(G))

    valrange = argv.get("valrange")
    if valrange:
        vals = range(-valrange, 0) + range(1, valrange+1)
    else:
        vals = argv.get("vals", [1,-1,2,-2,3,-3])

    for vec in search(nodes, nbd, eigval, vals, verbose=verbose):
        print("vec = ", strvec(vec))
        H = get_stab(G, vec)
        print("H:",len(H))
        J = get_stab(G, vec, -1)
        print("J:",len(J))
        if argv.draw:
            draw_graph(graph, vec)
            break
        if not argv.alleigs:
            break

def get_stab(G, vec, rep=1):
    H = []
    for g in G:
        u = {}
        for i,x in vec.items():
            u[g[i]] = rep*x
        if u==vec:
            H.append(g)
    return H



if __name__ == "__main__":

    main()



