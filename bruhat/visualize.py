#!/usr/bin/env python2

"""
generate pdf visualizations of small symmetric networkx graphs.
Use the automorphisms of the graph to layout nodes in
concentric circles.
"""

from __future__ import print_function

import sys, os
from random import gauss, random, seed, shuffle
from math import sin, cos, pi

#import numpy


import networkx as nx
from networkx.generators import small, classic


from equiv import get_autos
from argv import argv


if argv.seed is not None:
    seed(argv.seed)



def draw_graph(graph, pts, name):
    import pyx 
    from pyx import path, deco, trafo, style, text, color, deformer
    from pyx.color import rgb, cmyk
    from pyx.color import rgbfromhexstring as rgbhex

    black = rgb(0., 0., 0.) 
    blue = rgb(0., 0., 0.8)
    lred = rgb(1., 0.4, 0.4)
    white = rgb(1., 1., 1.) 


    #W = 10.
    #H = 10.

    R = 1.5
    r = 0.2

    c = pyx.canvas.canvas()

    for edge in graph.edges():
        src, tgt = edge

        x0, y0 = pts[src]
        x1, y1 = pts[tgt]

        c.stroke(path.line(R*x0, R*y0, R*x1, R*y1))

    for node in graph.nodes():
        x, y = pts[node]
        p = path.circle(R*x, R*y, r)
        c.fill(p, [white])
        c.stroke(p, [black])

    c.writePDFfile(name)
    

def pos_rand(graph):

    nodes = list(graph.nodes())

    pts = {}
    for key in nodes:
        pts[key] = 10*random(), 10*random()

    return pts


def metric(graph, pts):
    "sum of edge lengths"
    score = 0.
    for a, b in graph.edges():
        x0, y0 = pts[a]
        x1, y1 = pts[b]
        score += ((x1-x0)**2 + (y1-y0)**2)**0.5
    return score


def pos_circ(graph):

    G = get_autos(graph)

    print("autos:", len(G))

    if len(G.orbits()) == 1:
        print("transitive graph")

    best_g = None
    best_orbit = []
    for g in G:
        if g.is_identity():
            continue
        #print(g)
        for orbit in g.orbits():
            if len(orbit) > len(best_orbit):
                best_orbit = orbit
                best_g = g
        #print(g.orbits())

    orbits = best_g.orbits()
    print(orbits)


    best_pts = None
    best_score = 999999.*len(graph.edges())

    for trial in range(1000):

        shuffle(orbits)
        orbits.sort(key = lambda o:len(o))
        #thetas = [0.]*len(orbits)
        thetas = [2*pi*random() for _ in orbits]
        pts = pos_orbits(graph, orbits, thetas)
        score = metric(graph, pts)
        #print("score:", score)
        if score < best_score:
            best_pts = pts
            best_score = score

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
make_small_undirected_graph moebius_kantor_graph nx octahedral_graph
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

    draw_graph(graph, pts, "output")



if __name__ == "__main__":

    main()



