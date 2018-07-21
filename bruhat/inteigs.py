#!/usr/bin/env python2

"""
Eigenvectors of graphs.
"""

from __future__ import print_function

import sys, os
from random import gauss, random, seed, shuffle, choice
from math import sin, cos, pi

#import numpy


import networkx as nx
from networkx.generators import small, classic


from equiv import get_autos
from argv import argv


if argv.seed is not None:
    seed(argv.seed)


def search(nodes, nbd, eigval, vals=[1, -1, 2, -2, 3, -3], vec=None, verbose=False):
    assert nodes

    if vec is None:
        node = nodes[0]

        found = False
        for val in vals:
            if val >= 0:
                vec = {node : val}
                for vec in search(nodes, nbd, eigval, vals, vec, verbose): # recurse
                    yield vec
                    found = True
            #if found:
            #    break
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


def strvec(vec):
    if vec is None:
        return "None"
    n = max(list(vec.keys()))
    items = [vec.get(i) for i in range(n)]
    return str(items)


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

    vals = argv.get("vals", [1,-1,2,-2,3,-3])

    for vec in search(nodes, nbd, eigval, vals, verbose=verbose):
        print("vec = ", strvec(vec))
        H = get_stab(G, vec)
        print("H:",len(H))
        J = get_stab(G, vec, -1)
        print("J:",len(J))
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



