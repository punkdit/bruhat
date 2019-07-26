#!/usr/bin/env python3

"""
Here we demonstrate the galois connection between
spectra and symmetry groups of graphs.
"""

from __future__ import print_function

import sys, os

import os, sys

import numpy
import numpy.linalg
import networkx as nx
from networkx.generators import small, classic

from bruhat import isomorph
from bruhat.action import Perm, Group
#from solve import shortstr

from bruhat.argv import argv 


def mkgraph(nxgraph):
    idxs = list(nxgraph.nodes())
    #print idxs
    #print list(nxgraph.edges())

    points = [isomorph.Point("", ii) for ii, i in enumerate(idxs)]
    lookup = {}
    for ii, point in enumerate(points):
        lookup[point.idx] = point

    for (i, j) in nxgraph.edges():
        lookup[i].nbd.append(lookup[j])
        lookup[j].nbd.append(lookup[i])

    return isomorph.Graph(points)


def get_autos(nxgraph):
    graph0 = mkgraph(nxgraph)
    graph1 = mkgraph(nxgraph)

    items = [point.idx for point in graph0]
    perms = []
    for f in isomorph.search(graph0, graph1):
        perm = Perm(f, items)
        perms.append(perm)
        #print perm
    G = Group(perms, items)
    return G


def build_orbigraph(nxgraph, G):

    A = nx.adjacency_matrix(nxgraph)
    A = A.todense()
    #print A

    n = len(A)
    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    count = 0
    fs = set()
    for g in G:
        f = [None]*n
        for i, j in g.perm.items():
            assert i<n and j<n
            f[i] = j
            graph.add_edge(i, j)
        f = tuple(f)
        if f in fs:
            #write('/')
            continue # <---- continue
        fs.add(f)
        #write('.')
        count += 1
    #print
    #print "isomorphisms:", count

    equs = nx.connected_components(graph)
    equs = list(equs)
    m = len(equs)

    #print "components:", m

    P = numpy.zeros((n, m))
    Q = numpy.zeros((m, n))
    for i, equ in enumerate(equs):
        for j in equ:
            P[j, i] = 1
        Q[i, j] = 1

    #print shortstr(P)
    #print shortstr(Q)

    A = numpy.dot(Q, numpy.dot(A, P))

    return A


def divide(source, G):

    #print list(source.edges())
    for i, j in source.edges():
        assert i<=j

    nodes = nx.Graph()
    for i in source.nodes():
        nodes.add_node(i)

    edges = nx.Graph()
    for (i, j) in source.edges():
        edges.add_node((i, j))

    for g in G:
        for i in source.nodes():
            nodes.add_edge(i, g[i])
        for e in source.edges():
            i, j = g(e[0]), g(e[1])
            if i>j:
                i, j = j, i
            edges.add_edge(e, (i, j))

    nodes = nx.connected_components(nodes)
    nodes = list(nodes)

    print("nodes:", len(nodes), nodes)
    lookup = {}
    for idx, node in enumerate(nodes):
        for i in node:
            assert lookup.get(i) is None
            lookup[i] = idx

    edges = nx.connected_components(edges)
    edges = list(edges)

    print("edges:", len(edges), edges)

    n = len(nodes)
    A = numpy.zeros((n, n))

    for edge in edges:
        edge = list(edge)[0]
        i = lookup[edge[0]]
        j = lookup[edge[1]]
        A[i, j] += 1
        if i!=j:
            A[j, i] += 1

    #print A
    #print
    return A


def cayley(gen):

    G = Group.generate(gen)
    lookup = dict((g, i) for (i, g) in enumerate(G))
    graph = nx.Graph() # undirected graph, must have gen closed under inverse
    for i, g in enumerate(G):
        graph.add_node(i)

    #for g in gen:
    #    assert (g*g).is_identity()
    for g in gen:
        assert g.inverse() in gen # undirected graph!

    edges = set()
    for g in G:
        for h in gen:
            hg = h*g
            if (hg, g) in edges:
                continue
            edges.add((g, hg))
    for g, h in edges:
        graph.add_edge(lookup[g], lookup[h])

    return graph




def main():

#    for name in """desargues_graph petersen_graph 
#        dodecahedral_graph icosahedral_graph""".split():

    name = argv.next()
    if name is None:
        return

    if name == "cayley":

        n = argv.get("n", 3)
        items = range(n)
        gen = []
        for i in range(n-1):
            perm = dict((item, item) for item in items)
            perm[items[i]] = items[i+1]
            perm[items[i+1]] = items[i]
            gen.append(perm)
        gen = [Perm(perm, items) for perm in gen]
        for g in gen:
            print(g)

        graph = cayley(gen)

    elif name == "transpositions":

        n = argv.get("n", 3)
        items = range(n)
        gen = []
        for i in range(n-1):
            for j in range(i+1, n):
                perm = dict((item, item) for item in items)
                perm[items[i]] = items[j]
                perm[items[j]] = items[i]
                gen.append(perm)
        gen = [Perm(perm, items) for perm in gen]

        print("gen:", len(gen))
        graph = cayley(gen)

    elif name == "cycles":

        n = argv.get("n", 3)
        items = range(n)
        gen = []
        for i in range(n):
          for j in range(i+1, n):
            for k in range(j+1, n):
                perm = dict((item, item) for item in items)
                perm[items[i]] = items[j]
                perm[items[j]] = items[k]
                perm[items[k]] = items[i]
                gen.append(perm)
                perm = dict((item, item) for item in items)
                perm[items[i]] = items[k]
                perm[items[k]] = items[j]
                perm[items[j]] = items[i]
                gen.append(perm)
        gen = [Perm(perm, items) for perm in gen]

        print("gen:", len(gen))
        graph = cayley(gen)

    else:
    
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

    if argv.visualize:
        from visualize import draw_graph
        draw_graph(graph)
        return

    if argv.spec:
        print("spec:",)
        spec = nx.adjacency_spectrum(graph)
        spec = [x.real for x in spec]
        spec.sort()
        print(' '.join("%5f"%x for x in spec))
        return

    G = get_autos(graph)
    print("|G|=%d"%len(G))

    if argv.all_subgroups:
        Hs = list(G.subgroups())
        Hs.sort(key = lambda H : len(H))
    else:
        Hs = [Group.generate([G.identity])]

    for H in Hs:
        print("|H|=%d"%len(H))

        A = build_orbigraph(graph, H)
        print(A)

        spec = numpy.linalg.eigvals(A)

        spec = [x.real for x in spec]
        spec.sort(reverse=True)
        s = ' '.join(("%5f"%x).rjust(9) for x in spec)
        s = s.replace(".000000", "")
        print("spec:", s)
        print
    return




if __name__ == "__main__":

    main()


