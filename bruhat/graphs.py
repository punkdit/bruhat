#!/usr/bin/env python3

import sys, os

from math import sin, cos, pi

import numpy

from bruhat.action import Perm, Group
from bruhat.argv import argv
from bruhat import isomorph, inteigs


if 0:
    import pyx 
    from pyx import path, deco, trafo, style, text, color, deformer
    from pyx.color import rgb, cmyk
    from pyx.color import rgbfromhexstring as rgbhex

if 1:
    from huygens import front as pyx
    from huygens.front import *
    rgb = color.rgb

    black = rgb(0., 0., 0.) 
    blue = rgb(0., 0., 0.8)
    lred = rgb(1., 0.4, 0.4)
    red = rgb(1., 0.0, 0.0)
    green = rgb(0., 1.0, 0.0)
    lgreen = rgb(0.4, 1.0, 0.4)
    dgreen = rgb(0.0, 0.4, 0.0)
    white = rgb(1., 1., 1.) 

    center = [text.halign.boxcenter, text.valign.middle]

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



def zelim(A, verbose=False):

    A = A.copy()
    n = len(A)
    assert A.shape == (n, n)

    row = 0

    while row+1 < n:

        #print A
        #print "row = %d, col = %d" % (row, col)
        #assert A[row, :col].sum() == 0
        while row<n and A[row].sum() == 0:
            row += 1
    
        if row==n:
            break

        col = 0 
        while A[row, col] == 0:
            col += 1

        val0 = A[row, col]
        assert val0
        for row1 in range(row + 1, n): 
            val1 = A[row1, col]
            if val1 == 0:
                continue
            assert val1 % val0 == 0
            r = val1 // val0
            A[row1] -= r*A[row]
            assert A.min() >= 0
        #print
        #print shortstr(A)
    
        row += 1
    
    A = numpy.array([row for row in A if row.sum()])
    m, n = A.shape
    assert n==len(header)

#    if verbose:
#        print header
#        print A
#        
#        print "$$"
#        print latex_table(A, ' '*m, header)
#        print "$$"

    assert A.min() >= 0

    return A


def rank(A):
    A = zelim(A)
    idx = 0
    while idx<len(A):
        if A[idx].sum() == 0:
            break
        idx += 1
    return idx




class Graph(object):

    def __init__(self, nodes, edges, layout, fudge={}):
        self.nodes = list(nodes)
        assert set(nodes) == set(range(len(nodes)))
        self.edges = list(edges)
        nbd = dict((i, []) for i in nodes)
        for i, j in edges:
            nbd[i].append(j)
            nbd[j].append(i)
        self.nbd = nbd
        self.layout = dict(layout)
        i = self.nodes[0]
        self.valence = len(nbd[i])

        # displace some nodes by a "small" amount
        fudge = dict(fudge)
        for i in nodes:
            if fudge.get(i) is None:
                fudge[i] = 0., 0.
        self.fudge = fudge

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

    _G = None
    def get_autos(self):
        "find the graph autos"
        if self._G is not None:
            return self._G
        nodes = self.nodes
        edges = self.edges
        G = get_autos(nodes, edges)
        self._G = G
        return G

    def get_isoms(self):
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

    def get_inteigs(self, eigval, nnz=None):
        valss = []
        for r in range(1, 4):
            if argv.nonzero:
                vals = list(range(-r, 0)) + list(range(1, r+1))
            else:
                vals = list(range(-r, r+1))
            valss.append(vals)
        for vals in valss:
            found = False
            for vs in inteigs.search(self.nodes, self.nbd, eigval, vals):
                vec = Vec(self, vs)
                if nnz is not None:
                    if vec.get_nnz()==nnz:
                        yield vec
                else:
                    yield vec
                found = True
            if found:
                break

    def get_dimension(self, vecs):
        nodes = self.nodes
        n = len(nodes)
        A = numpy.zeros((n, n), dtype=numpy.int)
        for i in range(n):
            for vec in vecs:
                A[i, vec[i]] = 1
        fixfix
        dim = rank(A)
        return dim

    def get_orbit(self, vec):
        orbit = [vec]
        for g in self.get_autos():
            u = vec.act(g)
            if u not in orbit:
                orbit.append(u)
        return orbit

    def draw(self, vec=None, name=None, c=None, trafo=[]):
        nodes = self.nodes
        edges = self.edges
    
        R = 1.5 
        r = 0.2 
    
        if c is None:
            c = pyx.canvas.canvas()

        # fudge the layout
        layout = dict(self.layout)
        for i in nodes:
            dx, dy = self.fudge[i]
            x, y = layout[i]
            layout[i] = x+dx, y+dy
    
        for edge in edges:
            src, tgt = edge
    
            x0, y0 = layout[src]
            x1, y1 = layout[tgt]
    
            color = black

            delta = self.fudge.get(edge)
            if delta is None:
                c.stroke(path.line(R*x0, R*y0, R*x1, R*y1), [color]+trafo)
            else:
                dx, dy = delta
                #x2, y2 = (x0+x1)/2. + dx, (y0+y1)/2. + dy
                x2, y2 = (2*x0+x1)/3. + dx, (2*y0+y1)/3. + dy
                x3, y3 = (x0+2*x1)/3. + dx, (y0+2*y1)/3. + dy
                ps = [
                    path.moveto(R*x0, R*y0), 
                    path.curveto(R*x2, R*y2, R*x3, R*y3, R*x1, R*y1)]
                c.stroke(path.path(*ps), [color]+trafo)
            # path.curveto(x1, y1, x2, y2, x3, y3)
            # path.rcurveto(dx1, dy1, dx2, dy2, dx3, dy3)
    
        for node in nodes:
            x, y = layout[node]
    
            p = path.circle(R*x, R*y, 0.2)
            c.fill(p, [white]+trafo)
    
            p = path.circle(R*x, R*y, 0.10)
            c.fill(p, [black]+trafo)

            if argv.show_nodes:
                c.text(R*x+0.1, R*y+0.1, node, trafo)
    
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
                p = path.circle(R*x, R*y, (i+2)*0.10)
                c.stroke(p, [color, style.linewidth.THick]+trafo)
    
        if name is not None:
            print("Graph.draw: writePDFfile %s" % name)
            c.writePDFfile(name)
            c.writeSVGfile(name)
        return c


class Vec(object):
    def __init__(self, graph, vs):
        self.graph = graph
        self.nodes = graph.nodes
        self.nodes.sort()
        self.vs = dict(vs)
        assert set(self.nodes) == set(vs.keys())
    
    def __rmul__(self, r):
        w = {}
        for idx in self.nodes:
            w[idx] = r*self[idx]
        return Vec(self.graph, w)
    
    def __mul__(self, v):
        w = {}
        for idx in self.nodes:
            w[idx] = self[idx] * v[idx]
        return Vec(self.graph, w)
    
    
    def __add__(self, v):
        w = {}
        for idx in self.nodes:
            w[idx] = self[idx] + v[idx]
        return Vec(self.graph, w)

    def __getitem__(self, key):
        return self.vs[key]

    def __iter__(self):
        return iter(self.vs)

    def __hash__(self):
        vs = self.vs
        key = tuple((i, vs[i]) for i in self.nodes)
        return hash(key)

    def __eq__(self, other):
        for idx in self.nodes:
            if self.vs[idx] != other.vs[idx]:
                return False
        return True
    
    def __ne__(self, other):
        for idx in self.nodes:
            if self.vs[idx] != other.vs[idx]:
                return True
        return False

    def get_nnz(self):
        "non-zero entries"
        i = 0
        for v in self.vs.values():
            if v != 0:
                i += 1
        return i
    
    def act(self, g):
        u = dict((i, self.vs[g[i]]) for i in self.nodes)
        return Vec(self.graph, u)



def petersen_graph(R=2.0, r=1.0):
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


def double_cycle_graph(n=6, R=2.0, r=1.0):
    nodes = range(2*n)

    edges = [(i, (i+1)%n) for i in range(n)]
    for i in range(n):
        edges.append((i, i+n))
        edges.append((i+n, (i+1)%n+n))

    w = 2*pi / n
    layout = {}
    for i in range(n):
        theta = w*(i+0.5) 
        layout[i] = R*sin(theta), R*cos(theta)
        layout[i+n] = r*sin(theta), r*cos(theta)

    yield Graph(nodes, edges, layout)


def cubical_graph(*args, **kw):
    return double_cycle_graph(4, *args, **kw)


def cycle_graph(n=6, R=2.0):
    nodes = range(n)
    edges = []
    for i in range(n):
        edges.append((i, (i+1)%n))

    layout = {}
    w = 2 * pi / n
    for i in range(n):
        layout[i] = R*sin(w*i), R*cos(w*i)

    yield Graph(nodes, edges, layout)

    layout = {}
    w = 2 * pi / n
    for i in range(n):
        layout[i] = R*sin(w*i+0.5*w), R*cos(w*i+0.5*w)

    yield Graph(nodes, edges, layout)



def complete_graph(n=4, R=1.5):
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




def complete_bipartite_graph(n=3, m=None, R=1.5):
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



def cubic_twisted_graph(R=2.0):
    nodes = range(10)
    edges = []
    for i in range(10):
        edges.append((i, (i+1)%10))
    edges.extend([(0, 5), (1, 3), (2, 6), (4, 8), (7, 9)])

    layout = {}
    w = 2 * pi / 10
    for i in range(10):
        layout[i] = R*sin(w*i), R*cos(w*i)

    yield Graph(nodes, edges, layout)


def cayley_A4_graph():
    nodes = range(12)
    edges = [(0, 1), (1, 2), (2, 0),
        (1, 4), (2, 5), (0, 3),
        (3, 6), (3, 9), (4, 7), (4, 10), (5, 8), (5, 11),
        (6, 9), (9, 7), (7, 10), (10, 8), (8, 11), (11, 6)]
    layout = {}
    
    r0 = 0.7
    r1 = 1.5
    r2 = 2.4
    w = 2 * pi / 3
    d = 0.13
    for i in range(3):
        layout[i] = r0*sin(w*i), r0*cos(w*i)
        layout[i+3] = r1*sin(w*i), r1*cos(w*i)
        layout[i+6] = r2*sin(w*(i-d)), r2*cos(w*(i-d))
        layout[i+9] = r2*sin(w*(i+d)), r2*cos(w*(i+d))

    yield Graph(nodes, edges, layout)


def g13_graph():
    nodes = range(24)

    edges = []
    for i in range(3):
        a, b = (2*i+1, (2*i+2)%6)
        edges.append((a, b))
        for j in range(1, 4):
            a, b = (2*i, (2*i+1)%6)
            edges.append((a + 6*j, b + 6*j))

    edges.extend([
        (0, 10), (1, 9), (2, 18), (3, 23), (4, 14), (5, 13)])
    edges.extend([
        (7, 20), (8, 19), (6, 17), (11, 12), (15, 22), (16, 21)])

    for i in range(3):
        for j in range(4):
            edges.append((i+6*j, i+3+6*j))

    assert len(set(edges)) == len(edges)
    for i, j in edges:
        assert (j, i) not in edges

    layout = {}
    w = 2*pi/6
    r0 = 1.0
    r1 = 3.0
    for i in range(6):
        theta = w*(i-0.5)
        layout[i] = r0*sin(theta), r0*cos(theta)

    w1 = 2*pi/3
    for j in range(3):
        x, y = r1*sin(-w1*j), r1*cos(-w1*j)
        for i in range(6):
            theta = w*(i-0.5)
            layout[i + 6*(j+1)] = x + r0*sin(theta), y + r0*cos(theta)

    fudge = {}

    fudge[(6, 17)] = -0.5, 0.
    fudge[(16, 21)] = 0., -0.5
    fudge[(7, 20)] = 0.5, 0.

    yield Graph(nodes, edges, layout, fudge)


def g11_graph():
    nodes = range(10)
    edges = []
    for i in range(10):
        edges.append((i, (i+1)%10))
    edges.extend([(0, 5), (1, 8), (2, 9), (3, 6), (4, 7)])
    R = 2.0
    layout = {}
    w = pi / 10
    for i in range(10):
        layout[i] = R*sin(2*w*i), R*cos(2*w*i)
    yield Graph(nodes, edges, layout)


def g10_graph():
    nodes = range(20)
    tree = [
        (0, 1), (0, 2), (0, 3), 
        (1, 4), (1, 5), (2, 6), (2, 7), (3, 8), (3, 9)]
    edges = list(tree) + [(a+10, b+10) for (a, b) in tree]
    edges.extend([
        (4, 14), (4, 16), (5, 17), (5, 18), (6, 14), (6, 18),
        (7, 15), (7, 19), (8, 15), (8, 16), (9, 17), (9, 19)])
    w = 1.
    h = 1.
    layout = {}
    layout[0] = (0., 2.5*h)
    idx = 1
    for i in range(3):
        layout[idx] = (w, (i+1.5)*h)
        idx += 1
    for i in range(6):
        layout[idx] = (2*w, i*h)
        idx += 1
    for i in range(10):
        x, y = layout[i]
        layout[i+10] = 6*w - x, y
    yield Graph(nodes, edges, layout)


def fano_graph():
    "incidence geometry of the fano plane"
    nodes = range(14)
    edges = []
    for i in range(14):
        edges.append((i, (i+1)%14))
    for i in range(7):
        edges.append((2*i, (2*i+9)%14))

    r0 = 1.
    w = 2*pi / 14
    layout = {}
    for i in range(14):
        theta = w*i
        layout[i] = r0*sin(theta), r0*cos(theta)

    yield Graph(nodes, edges, layout)
    

def tutte_coxeter_graph():
    nodes = range(30)
    edges = []
    for i in range(10):
        edges.append((i, (i+3)%10))

        edges.append((i, i+10))
        if i<5:
            edges.append((i+10, (i+5)%10+10))
        edges.append((i+10, i+20))
        edges.append((i+20, (i+1)%10+20))

    layout = {}
    r0 = 1.5
    r1 = r0 + 0.6
    r2 = r1 + 0.6
    w = 2*pi / 10
    for i in range(10):
        theta = w*i
        layout[i] = r0*sin(theta + 0.5*w), r0*cos(theta + 0.5*w)
        layout[i+10] = r1*sin(theta), r1*cos(theta)
        layout[i+20] = r2*sin(theta), r2*cos(theta)

    for i, j in edges:
        assert (j, i) not in edges
    yield Graph(nodes, edges, layout)


def desargues_graph():
    nodes = range(20)
    edges = []
    for i in range(6):
        edges.append((i, (i+1)%6))
        edges.append((i, i+6))
        edges.append((i+6, i+12))
        edges.append((i+12, (i+1)%6+12))
    for i in range(3):
        edges.append((2*i+6, 18))
        edges.append((2*i+7, 19))
    layout = {}

    r0 = 1.
    r1 = r0 + 0.6
    r2 = r1 + 0.6
    w = 2*pi / 6
    for i in range(6):
        theta = w*i
        layout[i] = r0*sin(theta), r0*cos(theta)
        layout[i+6] = r1*sin(theta), r1*cos(theta)
        layout[i+12] = r2*sin(theta), r2*cos(theta)

    layout[18] = 0., 0.
    layout[19] = 0., 0.

    fudge = {}
    fudge[18] = -0.3*r0, 0.
    fudge[19] = +0.3*r0, 0.
    rz = 0.6
    for i in range(3):
        theta = 2*(i+0.5)*pi/3
        edge = (2*i+6, 18)
        fudge[edge] = rz*sin(theta), rz*cos(theta)
        theta = 2*(i+0.0)*pi/3
        edge = (2*i+7, 19)
        fudge[edge] = -rz*sin(theta), -rz*cos(theta)

    yield Graph(nodes, edges, layout, fudge) # hexagonal layout

    r0 = 1.
    r1 = r0 + 1.
    r2 = r1 + 1.

    edges = []
    for i in range(10):
        edges.append((i, (i+3)%10))
        edges.append((i, i+10))
        edges.append((i+10, (i+1)%10+10))

    w = 2*pi / 10
    for i in range(10):
        theta = w*i
        layout[i] = r0*sin(theta), r0*cos(theta)
        layout[i+10] = r1*sin(theta), r1*cos(theta)

    yield Graph(nodes, edges, layout) # decagonal layout

    r0 = 3.

    edges = []
    for i in range(20):
        edges.append((i, (i+1)%20))
    for i in range(20):
        if i%4==0:
            edges.append((i, (i+9)%20))
        if i%4==1:
            edges.append((i, (i-9)%20))
        if i%4==2:
            edges.append((i, (i+5)%20))
        if i%4==3:
            edges.append((i, (i-5)%20))

    w = 2*pi / 20
    for i in range(20):
        theta = w*(i-0.5)
        layout[i] = r0*sin(theta), r0*cos(theta)

    yield Graph(nodes, edges, layout) # 20-sided polygonal layout


# TODO: find closed path that visits all nodes with autos that act by rotation.
# layout in a circle.

def main():
    graph_name = argv.next()
    fn = eval(graph_name)

    n = argv.get("n")
    if n is not None:
        items = fn(int(n))
        graph_name = graph_name + "_" + str(n)
    else:
        items = fn()

    graphs = list(items)

    layout = argv.get("layout", 0)

    graph = graphs[layout]
    
    if argv.autos:
        G = graph.get_autos()
        print("|autos| =", len(G))
        #for g in G:
        #    print(g)
        #break

    if argv.draw:
        graph.draw(name=graph_name)
        return

    autos = graph.get_autos()
    isoms = graph.get_isoms()
    print("|isoms| =", len(isoms))

    A = graph.get_adjacency()
    #print(A)
    vals, vecs = numpy.linalg.eigh(A)
    print("evals:", vals)

    eigval = argv.get("eigval", vals[0])

    name = argv.get("name", "output")
    graph.draw(name=name)

    for vec in graph.get_inteigs(eigval):
        #print(vec)
        count = 0
        orbit = [vec]
        for g in autos:
            u = vec.act(g)
            if u==vec:
                count += 1
            else:
                orbit.append(u)
        #dim = graph.get_dimension(orbit)
        #print("count=%d, dim=%d" % (count, dim))
        print(count, end=" ")

        if argv.draw:
            graph.draw(vec, name)
            break
        if count==36:
            graph.draw(vec, name)
            break

    print()



def render():
    graphs = list(desargues_graph())

    hexa = graphs[0]
    deca = graphs[1]

    graph = graphs[0]
    graphs[0].draw(None, "desargues_0")
    graphs[1].draw(None, "desargues_1")

    G = graph.get_autos()
    print("|G| =", len(G))

    isoms = graph.get_isoms()
    print("isoms =", len(isoms))

    if 0:
        top = graph.valence
        vecs = list(graph.get_inteigs(-top))
        assert len(vecs)==1
        nvec = vecs[0]

        for pvec in graph.get_inteigs(2):
            break
        vec = (pvec * nvec)
        graph.draw(vec, "desargues_1_n2")
    
        for pvec in graph.get_inteigs(1):
            break
        vec = (pvec * nvec)
        graph.draw(vec, "desargues_1_n1")
    
    if 0:
        for val in [3, 2, 1]:
            for pvec in graph.get_inteigs(val):
                break
            graph.draw(pvec, "desargues_0_%d"%val)
        
            nvec = iter(graph.get_inteigs(-val)).__next__()
            graph.draw(nvec, "desargues_0_n%d"%val)

    if 0:
      for val in [2, 1]:
        nnz = None
        if val==1:
            nnz = 8
        vec = iter(hexa.get_inteigs(val, nnz)).__next__()
        orbit = hexa.get_orbit(vec)
        print("orbit:", len(orbit))
        for i, pvec in enumerate(orbit):
            name = "desargues_hexa_v%d_i%d"%(val, i)
            print('%s.svg<img src="%s.svg"><br>' % (name, name))
            hexa.draw(pvec, name)
        
        vec = iter(deca.get_inteigs(val, nnz)).__next__()
        orbit = deca.get_orbit(vec)
        print("orbit:", len(orbit))
        for i, pvec in enumerate(orbit):
            name = "desargues_deca_v%d_i%d"%(val, i)
            print('%s.svg<img src="%s.svg"><br>' % (name, name))
            deca.draw(pvec, name)
    
        vec = iter(hexa.get_inteigs(-val, nnz)).__next__()
        orbit = hexa.get_orbit(vec)
        print("orbit:", len(orbit))
        for i, pvec in enumerate(orbit):
            name = "desargues_hexa_vn%d_i%d"%(val, i)
            print('%s.svg<img src="%s.svg"><br>' % (name, name))
            hexa.draw(pvec, name)
        
        vec = iter(deca.get_inteigs(-val, nnz)).__next__()
        orbit = deca.get_orbit(vec)
        print("orbit:", len(orbit))
        for i, pvec in enumerate(orbit):
            name = "desargues_deca_vn%d_i%d"%(val, i)
            print('%s.svg<img src="%s.svg"><br>' % (name, name))
            deca.draw(pvec, name)
    
    if 1:
        vec = iter(deca.get_inteigs(2)).__next__()
        orbit = deca.get_orbit(vec)
        print("orbit:", len(orbit))

        # show an equation here
        #svecs = [orbit[i] for i in [5, 7, 19, 15, 3]]
        #coefs = [1, -1, 1, -1, 1] # sum with these coefs to get zero
        svecs = [orbit[i] for i in [5, 19, 3, 7, 15]]
        symbols = "++=+ "
    
        c = pyx.canvas.canvas()
        x, y = 0., 0.
        r = 0.43
        dx = 8.0 * r
        #u = dict((i, 0) for i in deca.nodes)
        for idx in range(5):
            vec = svecs[idx]
            #vec = coefs[idx]*vec
            #u = u+vec
            deca.draw(vec, None, c, [trafo.scale(r, r), trafo.translate(x, y)])
            c.text(x+0.5*dx, y, symbols[idx], center)
            x += dx
    
        # should be zero
        #deca.draw(u, None, c, [trafo.scale(r, r), trafo.translate(x, y)])
    
        name = "output/desargues_deca_v2_4_1"
        c.writeSVGfile(name)
        c.writePDFfile(name)



if __name__ == "__main__":
    if argv.render:
        render()
    else:
        main()


