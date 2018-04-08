#!/usr/bin/env python3

import sys, os
from random import gauss, random, seed

#import numpy


import networkx as nx
from networkx.generators import small, classic

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

    R = 10.
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


    nodes = list(graph.nodes())

    pts = {}
    for key in nodes:
        pts[key] = random(), random()

    draw_graph(graph, pts, "output")



if __name__ == "__main__":

    main()



