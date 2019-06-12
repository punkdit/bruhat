#!/usr/bin/env python2

"""
generate pdf visualizations of small symmetric networkx graphs.
Use the automorphisms of the graph to layout nodes in
concentric circles.
"""

from __future__ import print_function

import sys, os
from random import gauss, random, seed, shuffle, choice
from math import sin, cos, pi

#import numpy


import networkx as nx
from networkx.generators import small, classic


from bruhat.equiv import get_autos
from bruhat.argv import argv


if argv.seed is not None:
    seed(argv.seed)


# whoop
node_colors = {}
edge_colors = {}

def draw_graph(graph, pts=None, name="output"):
    import pyx 
    from pyx import path, deco, trafo, style, text, color, deformer
    from pyx.color import rgb, cmyk
    from pyx.color import rgbfromhexstring as rgbhex

    black = rgb(0., 0., 0.) 
    blue = rgb(0., 0., 0.8)
    lred = rgb(1., 0.4, 0.4)
    red = rgb(1., 0.0, 0.0)
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
        if edge in edge_colors:
            color = eval(edge_colors[edge])
        c.stroke(path.line(R*x0, R*y0, R*x1, R*y1), [color])

    for node in graph.nodes():
        x, y = pts[node]
        p = path.circle(R*x, R*y, r)

        color = white
        if node in node_colors:
            color = eval(node_colors[node])
        c.fill(p, [color])
        c.stroke(p, [black])

    c.writePDFfile(name)
    

def pos_rand(graph):

    nodes = list(graph.nodes())

    # try random arrangements, save the best.
    best_score = 1e100
    trials = argv.get("trials", 1000)
    for trial in range(trials):

        pts = {}
        for key in nodes:
            pts[key] = 10*random(), 10*random()

        score = metric(graph, pts)
        #print("score:", score)
        if score < best_score:
            best_pts = pts
            best_score = score

    metric(graph, best_pts, verbose=True)
    print("best_score:", best_score)
    return best_pts


def metric(graph, pts, verbose=False):
    "evaluate graph layout: lower score is better"

    #global node_colors # debug hack
    #node_colors = {}

    score = 0.

    # start with edge lengths
    for a, b in graph.edges():
        x0, y0 = pts[a]
        x1, y1 = pts[b]
        score += ((x1-x0)**2 + (y1-y0)**2)**0.5

#    if not verbose:
#        return score

    EPSILON = 1e-3

    radius = 0.2 # radius of nodes
    for a, b in graph.edges():
        xa, ya = pts[a]
        for c in graph.nodes():
            if c==a or c==b:
                continue
            # does c come close to the edge a--b ?
            xb, yb = pts[b] # mutates!
            xc, yc = pts[c]

            # work relative to node a
            xc -= xa
            yc -= ya
            xb -= xa
            yb -= ya

            rc = xc**2 + yc**2 
            rb = xb**2 + yb**2

            if rc > rb - EPSILON:
                continue # too far away

            dot = xc*xb + yc*yb
            if dot < EPSILON:
                continue # wrong way

            dist2 = (rc - (dot**2/rb)) + EPSILON
            if dist2 < 0.:
                print(dist2)
                assert 0
            if dist2**0.5 < radius:
                #if verbose:
                #print("a:", pts[a])
                #print("b:", pts[b])
                #print("c:", pts[c])
                #print("c:", xc, yc)
                #print("b:", xb, yb)
                #print(rc, rb)
                #print(a, b, c, "dist", dist)
                score += 10
                #node_colors[c] = "red"
                #edge_colors[a, b] = "red"
                #return score

    return score


def pos_circ(graph):

    G = get_autos(graph)

    print("autos:", len(G))
    if len(G)==1:
        print("trivial autos...")
        return None

    if len(G.orbits()) == 1:
        print("transitive graph")

    orbit_size = argv.orbit_size # look for orbit of this size

    best_gs = None
    best_orbit = []
    descs = []
    for g in G:
        if g.is_identity():
            continue
        #print(g)
        orbits = g.orbits()
        sizes = [len(orbit) for orbit in orbits]
        sizes.sort()
        #print(sizes)
        descs.append(tuple(sizes))
        if orbit_size is not None:
            if orbit_size not in sizes:
                continue
        #if set(sizes) != set([2]):
        #    continue
        for orbit in orbits:
            if len(orbit) > len(best_orbit):
                best_orbit = orbit
                best_gs = [g]
            elif len(orbit) == len(best_orbit):
                best_gs.append(g)

    if argv.desc:
        descs.sort()
        for desc in descs:
            print(desc)

    if not best_gs:
        print("no orbits found")
        return

    # print some stuff...
    strs = set()
    for g in best_gs:
        orbits = g.orbits()
        sizes = [len(orbit) for orbit in orbits]
        sizes.sort()
        strs.add(str(sizes))
    for s in strs:
        print("orbit sizes:", s)
    print("best_gs:", len(best_gs))

    best_pts = None
    best_score = 999999.*len(graph.edges())

    # now just try many random arrangements, save the best.
    trials = argv.get("trials", 1000)
    for trial in range(trials):

        g = choice(best_gs)
        orbits = list(g.orbits())
        shuffle(orbits)
        #orbits.sort(key = lambda o:len(o))
        #thetas = [0.]*len(orbits)
        thetas = [2*pi*random() for _ in orbits]
        thetas[0] = 0.
        pts = pos_orbits(graph, orbits, thetas)
        score = metric(graph, pts)
        #print("score:", score)
        if score < best_score:
            best_pts = pts
            best_score = score

    metric(graph, best_pts, verbose=True)
    print("best_score:", best_score)
    return best_pts


def pos_orbits(graph, orbits, thetas):

    assert len(thetas) == len(orbits)
    dR = 1.5
    R = 0.
    if len(orbits[0])>1:
        R = dR

    pts = {} # here we go
    for idx, orbit in enumerate(orbits):
        #theta = 2*pi*random()
        theta = thetas[idx]
        dtheta = 2*pi/len(orbit)
        for node in orbit:
            pts[node] = R*sin(theta), R*cos(theta)
            theta += dtheta
        R += dR

    #pts = pos_rand(graph)
    return pts

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

def main():

    name = argv.next()

    builder = getattr(small, name, None)

    if builder is None:
        builder = getattr(classic, name, None)

    if builder is None:
        print(name, "not found")
        return

    n = argv.n
    if n is not None:
        graph = builder(n)
    else:
        graph = builder()

    print("|nodes|=%d, edges=|%d|"%(len(graph.nodes()), len(graph.edges())))

    print(list(graph.nodes()))
    print(list(graph.edges()))

    if argv.pos_rand:
        pts = pos_rand(graph)
    elif argv.pos_circ:
        pts = pos_circ(graph)
    else:
        pts = pos_circ(graph)

    output = argv.get("output", "output")
    if pts is not None:
        draw_graph(graph, pts, output)
    else:
        print("failed to draw")



if __name__ == "__main__":

    main()



