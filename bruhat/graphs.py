#!/usr/bin/env python3

import sys, os

from math import sin, cos, pi

import numpy

from action import Perm, Group
from argv import argv
import isomorph


def mkgraph(nodes, edges):
    idxs = list(nodes)

    points = [isomorph.Point("", ii) for ii, i in enumerate(idxs)]
    lookup = {}
    for ii, point in enumerate(points):
        lookup[point.idx] = point

    for (i, j) in edges:
        lookup[i].nbd.append(lookup[j])
        lookup[j].nbd.append(lookup[i])

    return isomorph.Graph(points)


def get_autos(nodes, edges):
    graph0 = mkgraph(nodes, edges)
    graph1 = mkgraph(nodes, edges)

    items = [point.idx for point in graph0]
    perms = []
    for f in isomorph.search(graph0, graph1):
        perm = Perm(f, items)
        perms.append(perm)
        #print perm
    G = Group(perms, items)
    return G




class Graph(object):

    def __init__(self, nodes, edges, layout, G=None):
        self.nodes = list(nodes)
        self.edges = list(edges)
        self.layout = dict(layout)
        self.G = G

    def get_adjacency(self):
        nodes = self.nodes
        edges = self.edges
        n = len(nodes)
        #A = numpy.zeros((n, n), dtype=numpy.int)
        A = numpy.zeros((n, n))
        for (i, j) in edges:
            A[i, j] = 1
            A[j, i] = 1
        return A

    def get_metric(self):
        nodes = self.nodes
        edges = self.edges
        layout = self.layout
        n = len(nodes)
        #A = numpy.zeros((n, n), dtype=numpy.int)
        A = numpy.zeros((n, n))
        for i in range(n):
            x0, y0 = layout[i]
            for j in range(i+1, n):
                x1, y1 = layout[j]
                r = (x0-x1)**2 + (y0-y1)**2
                A[i, j] = r
                A[j, i] = r
        return A

    def permute(self, g):
        perm = g.perm
        nodes = list(self.nodes) # same nodes
        #edges = [(perm[i], perm[j]) for (i, j) in self.edges]
        edges = list(self.edges) # same edges
        layout = dict((i, self.layout[perm[i]]) for i in nodes) # move the layout!
        return Graph(nodes, edges, layout)

    def get_autos(self):
        "find the graph autos that act by isometries (on the layout)"
        nodes = self.nodes
        edges = self.edges
        G = get_autos(nodes, edges)
        A0 = self.get_metric()
        perms = []
        for g in G:
            graph = self.permute(g)
            A1 = graph.get_metric()
            if numpy.allclose(A0, A1):
                perms.append(g)
        G = Group(perms, g.items)
        return G

    def draw(self, vec=None, name="output"):
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

        nodes = self.nodes
        edges = self.edges
        layout = self.layout
    
        #W = 10.
        #H = 10.
    
        R = 1.5 
        r = 0.2 
    
        c = pyx.canvas.canvas()
    
        for edge in edges:
            src, tgt = edge
    
            x0, y0 = layout[src]
            x1, y1 = layout[tgt]
    
            color = black
            c.stroke(path.line(R*x0, R*y0, R*x1, R*y1), [color])
    
        for node in nodes:
            x, y = layout[node]
    
            #p = path.circle(R*x, R*y, 0.2)
            #c.fill(p, [white])
    
            p = path.circle(R*x, R*y, 0.10)
            c.fill(p, [black])
    
            val = 0
            if vec is not None:
                val = vec[node]
                assert int(val) == val
    
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



def petersen_graph():
    R = 2.0
    r = 1.0
    nodes = range(10)
    layout = {}
    w = 2*pi/6
    for i in range(6):
        layout[i] = R*cos(w*i), R*sin(w*i)
    w = 2*pi/3
    for i in range(6, 9):
        layout[i] = r*sin(w*i), -r*cos(w*i)
    layout[9] = (0., 0.)
    edges = [(i, (i+1)%6) for i in range(6)]
    for i in range(3):
        edges.append((6+i, 9))
        edges.append((6+i, 2*i))
        edges.append((6+i, (2*i+3)%6))

    if 0:
        # TODO: should be able to construct G from all graph
        # automorphisms. Choose the ones that preserve the distance between nodes.
    
        perm = {}
        for i in range(6):
            perm[i] = (i+2)%6
        for i in range(3):
            perm[i+6] = (i+1)%3 + 6
        perm[9] = 9
        g = Perm(perm, nodes)
    
        perm = {0:3, 1:2, 2:1, 3:0, 4:5, 5:4, 6:6, 7:8, 8:7, 9:9}
        h = Perm(perm, nodes)
    
        G = Group.generate([g, h])
        assert len(G) == 6

    yield Graph(nodes, edges, layout)

    # ------------

    edges = [(i, (i+1)%5) for i in range(5)]
    for i in range(5):
        edges.append((i, i+5))
        edges.append((i+5, (i+2)%5+5))

    w = 2*pi / 5
    layout = {}
    for i in range(5):
        layout[i] = R*sin(w*i), R*cos(w*i)
        layout[i+5] = r*sin(w*i), r*cos(w*i)

    yield Graph(nodes, edges, layout)


def cubical_graph():
    R = 2.0
    r = 1.0

    nodes = range(8)

    edges = [(i, (i+1)%4) for i in range(4)]
    for i in range(4):
        edges.append((i, i+4))
        edges.append((i+4, (i+1)%4+4))

    w = 2*pi / 4
    layout = {}
    for i in range(4):
        theta = w*(i+0.5)
        layout[i] = R*sin(theta), R*cos(theta)
        layout[i+4] = r*sin(theta), r*cos(theta)

    yield Graph(nodes, edges, layout)


def cycle_graph(n=6):
    R = 2.0

    nodes = range(n)
    edges = []
    for i in range(n):
        edges.append((i, (i+1)%n))

    layout = {}
    w = 2 * pi / n
    for i in range(n):
        layout[i] = R*sin(w*i), R*cos(w*i)

    yield Graph(nodes, edges, layout)



def complete_graph(n=4):
    R = 2.0

    nodes = range(n)
    edges = []
    for i in range(n):
      for j in range(i+1, n):
        edges.append((i, j))

    layout = {}
    w = 2 * pi / n
    for i in range(n):
        layout[i] = R*sin(w*i), R*cos(w*i)

    yield Graph(nodes, edges, layout)

    if n%2==0:
        layout = {}
        layout[0] = 0., 0.
        w = 2 * pi / (n-1)
        for i in range(1, n):
            layout[i] = R*sin(w*i), R*cos(w*i)
    
        yield Graph(nodes, edges, layout)




def complete_bipartite_graph(n=3, m=None):
    R = 2.0
    if m is None:
        m = n

    nodes = range(n+m)
    edges = []
    for i in range(n):
      for j in range(m):
        edges.append((i, n+j))

    assert n==m
    layout = {}
    w = pi / n
    for i in range(n):
        layout[i] = R*sin(2*w*i), R*cos(2*w*i)
        layout[i+n] = R*sin(2*w*i+w), R*cos(2*w*i+w)

    yield Graph(nodes, edges, layout)


# TODO: find closed path that visits all nodes with autos that act by rotation.
# layout in a circle.

def main():
    name = argv.next()
    fn = eval(name)

    arg = argv.next()
    if arg is not None:
        items = fn(int(arg))
    else:
        items = fn()

    for layout in items:
    
        G = layout.get_autos()
        print("|G| =", len(G))
        for g in G:
            print(g)

    layout.draw()

    A = layout.get_adjacency()
    print(A)
    vals, vecs = numpy.linalg.eigh(A)
    print(vals)


if __name__ == "__main__":
    main()




