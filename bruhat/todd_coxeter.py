#!/usr/bin/env python3

"""
Build groups from presentations,
and Schreier coset graphs using Todd-Coxeter algorithm.
See also bruhat.o_todd_coxeter for an object oriented version.
See also bruhat.hyperbolic


Implementation from: 
https://math.berkeley.edu/~kmill/notes/todd_coxeter.html

See also:
https://arxiv.org/pdf/2012.09271.pdf p10.

"""

import sys
from random import randint, seed, choice
from time import sleep
from functools import reduce
from operator import mul

import numpy
from numpy import concatenate

#from bruhat.gset import Perm, Group, Coset, mulclose # FAIL
from bruhat.action import Perm, Group, Coset, mulclose, close_hom
from bruhat.util import cross
from bruhat.smap import SMap
from bruhat.argv import argv

    
def cycle(perm):
    # The permutations can be written in cycle notation fairly easily.
    parts = []
    for i in range(len(perm)):
        part = [str(i+1)]
        k = perm[i]
        while k != i:
            if k < i: break
            part.append(str(k+1))
            k = perm[k]
        else:
            parts.append(" ".join(part))
    return "("+")(".join(parts)+")"
    

class Schreier(object):
    """ 
    A Schreier coset graph coming from a group acting on the cosets
    of a subgroup.
    """
    DEBUG = False

    def __init__(self, ngens, rels):
        self.labels = []
        self.neighbors = []
        self.ngens = ngens
        self.rels = list(rels)

        for rel in rels:
            for gen in rel:
                assert type(gen) is int, repr(gen)
                assert 0 <= gen < ngens, repr(rel)

        #The labels variable is a list of _numbers, with the property
        #that labels[i] <= i. This is a union-find data structure
        #for keeping track of the vertex quotients for the Schreier
        #graph so far. The find function _looks up the current-lowest
        #label for a particular labeled vertex.
        
        #For vertices which have not been removed, the neighbors
        #data structure contains a list of all the neighboring
        #vertices, with None for non-existent neighbors — non-existent
        #neighbors are presumed to be a portion of the free group’s
        #universal Schreier graph, namely a tree. Each entry is
        #another vertex label. Vertex labels may be out of date,
        #so find must be run on them to ensure validity.
        #
        #The to_visit variable keeps track of which vertices have
        #been marked. Any vertex whose label is less than to_visit
        #is considered marked.
        
    def get_label(self, c):
        labels = self.labels
        assert c is not None
        c2 = labels[c]
        if c == c2:
            return c
        else:
            c2 = self.get_label(c2) # recurse
            labels[c] = c2
            return c2

#    def get_label(self, c):
#        labels = self.labels
#        c1 = labels[c]
#        while 1:
#            c2 = labels[c1]
#            if c2 == c1:
#                break
#            c1 = c2
#        labels[c] = c1
#        return c1
    
    # The identification procedure is called unify. 

    """
    I think I see the bug -- in the unify routine, when it's copying over and
    unifying the neighbors list, it immediately calls unify(n1, n2).  This is
    not correct, since this might result in c1 or c2 no longer satisfying c1 ==
    find(c1) or c2 == find(c2), and as a result the neighbors will either get
    stored into an old neighbors list or read from an old neighbors list.  I
    think the fix is to create a to_unify = [] array right before the loop,
    replace unify(n1, n2) with to_unify.append((n1, n2)), and then after the
    loop have a second loop for n1, n2 in to_unify: unify(n1, n2).
    """
    
    def unify_recursive(self, c1, c2):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        c1 = self.get_label(c1)
        c2 = self.get_label(c2)
        if c1 == c2:
            return
        c1, c2 = min(c1, c2), max(c1, c2)
        if self.DEBUG:
            print("unify:", c2, "-->", c1)
        labels[c2] = c1
        to_unify = []
        for d in range(ngens):
            n1 = neighbors[c1][d]
            n2 = neighbors[c2][d]
            if n1 == None:
                neighbors[c1][d] = n2
            elif n2 != None:
                #self.unify(n1, n2) # recurse
                to_unify.append((n1, n2))
        for n1, n2 in to_unify:
            self.unify_recursive(n1, n2) # recurse
    
    def unify(self, c1, c2):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        to_unify = [(c1, c2)]
        while to_unify:
            c1, c2 = to_unify.pop()
            c1 = self.get_label(c1)
            c2 = self.get_label(c2)
            if c1 == c2:
                continue
            c1, c2 = min(c1, c2), max(c1, c2)
            if self.DEBUG:
                print("unify:", c2, "-->", c1)
            labels[c2] = c1
            for d in range(ngens):
                n1 = neighbors[c1][d]
                n2 = neighbors[c2][d]
                if n1 == None:
                    neighbors[c1][d] = n2
                elif n2 != None:
                    to_unify.append((n1, n2))

    #unify = unify_recursive
    
    #    It takes two vertices, and makes the one with the greater
    #    label refer to the lesser one with labels[c2] = c1. Then,
    #    it moves all of the neighbors of the deleted vertex by
    #    recursively calling unify. In case the neighbor being
    #    replaced is non-existent, which is the base case of the
    #    recursion, we just record it as a neighbor. There is
    #    no need to record this identification on the other end
    #    of the edge because the other end is referring to c2,
    #    which is now c1.
    
    def add_vertex(self):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        c = len(labels)
        labels.append(c)
        neighbors.append(ngens*[None])
        return c

    # For following paths, we have the following two functions: 
    
    def follow_step(self, c, d):
        labels = self.labels
        neighbors = self.neighbors
        c = self.get_label(c)
        ns = neighbors[c]
        if ns[d] == None:
            ns[d] = self.add_vertex()
        return self.get_label(ns[d])
    
    def follow_path(self, idx, word):
        idx = self.get_label(idx)
        for jdx in reversed(word):
            idx = self.follow_step(idx, jdx)
        return idx
    
    #    The first takes a vertex and finds the neighbor in the
    #    d direction, creating a new vertex in that direction
    #    if needed, and the second follows a list of directions
    #    to find the end of a path. The follow function creates
    #    new neighbors with the add_vertex() function.

    def dump(self):
        labels = self.labels
        neighbors = self.neighbors
        print(labels)
        for idx, row in enumerate(neighbors):
            jdx = self.get_label(idx)
            if idx == jdx:
                print('\t', idx, row)
        print("-"*70)

    def __len__(self):
        return len([k for k in enumerate(self.labels) if k[0]==k[1]])
    
    #    Finally, there is the core algorithm: 
    def build(self, hgens=[], maxsize=None):
        labels = self.labels
        neighbors = self.neighbors
        rels = self.rels

        for rel in hgens:
            for gen in rel:
                assert 0 <= gen < self.ngens, repr(rel)

        start = self.add_vertex()
    
        if self.DEBUG:
            self.dump()
        
        for hgen in hgens:
            self.unify(self.follow_path(start, hgen), start)
        
        if self.DEBUG:
            self.dump()

        to_visit = 0
        while to_visit < len(labels):
            if self.DEBUG:
                print("to_visit =", to_visit)
            c = self.get_label(to_visit)
            if c == to_visit:
                for rel in rels:
                    self.unify(self.follow_path(c, rel), c)
                if self.DEBUG:
                    self.dump()
            to_visit += 1
            if maxsize and len(neighbors) > maxsize:
                return None
        return self
        
    #    It creates the start vertex, adds all of the relations
    #    for H as relators at this basepoint, and then for each
    #    unmarked vertex, adds all of the relations from rels.
    #    Notice how unify is being used to turn paths into relators.
    
    def get_gens(self):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        cosets = [] # act on these
        for idx, row in enumerate(neighbors):
            jdx = self.get_label(idx)
            if idx == jdx:
                cosets.append(idx)
        n = len(cosets)
        #print("cosets:", n)
        lookup = dict((v,k) for (k,v) in enumerate(cosets))
        gens = []
        items = list(range(n))
        for i in range(ngens):
            perm = {}
            for idx in cosets:
                nbd = neighbors[idx]
                assert nbd[i] is not None, nbd
                src, tgt = lookup[idx], lookup[self.get_label(nbd[i])]
                assert perm.get(src) is None
                perm[src] = tgt
            perm = Perm(perm, items)
            gens.append(perm)
        return gens

    def get_group(self):
        gens = self.get_gens()
        G = Group.generate(gens)
        return G

    def get_words(self):
        cosets = set([c for i, c in enumerate(self.labels) if i == c])
        ngens = self.ngens
        words = [()]
        found = set([0])
        bdy = list(words)
        while bdy:
            _bdy = []
            for word in bdy:
              for gen in range(ngens):
                w = (gen,) + word # left-multiply 
                idx = self.follow_path(0, w)
                if idx not in found:
                    found.add(idx)
                    _bdy.append(w)
                    words.append(w)
            bdy = _bdy
        return words

    def get_poincare(self, ring, q):
        "the poincare polynomial"
        dists = [len(w) for w in self.get_words()]
        N = max(dists)
        p = ring.zero
        for i in range(N+1):
            p += dists.count(i) * q**i
        return p

    def find_words(self):
        bdy = [0]
        words = {0:()}
        while bdy:
            _bdy = []
            for idx in bdy:
                word = words[idx]
                nbd = self.neighbors[idx]
                for j, jdx in enumerate(nbd):
                    if jdx is None:
                        continue
                    jdx = self.labels[jdx]
                    if jdx not in words:
                        _bdy.append(jdx)
                        words[jdx] = (j,)+word
            bdy = _bdy
        return words
    
    def show(self):
        #    After this, the data structures contain the Schreier
        #    graph for G/H .  This can be interpreted as a permutation
        #    representation, for instance using
        
        cosets = [c for i, c in enumerate(self.labels) if i == c]
        
        perms = [[cosets.index(self.follow_step(c, d)) for i, c in enumerate(cosets)]
                 for d in range(self.ngens)]
        
        #    to enumerate the cosets (which are vertices which have
        #    the least label in their equivalency classes), and then
        #    taking edges of the Schreier graph for each generator
        #    to construct a permutation.
        #    
    
        for d in range(self.ngens):
            print("g%d ="%d, cycle(perms[d]))

    def get_dot(self, colours, name=None):
        cosets = [c for i, c in enumerate(self.labels) if i == c]
        perms = [[cosets.index(self.follow_step(c, d)) for i, c in enumerate(cosets)]
                 for d in range(self.ngens)]
        lines = []
        lines.append("graph the_graph")
        lines.append("{")
        lines.append("""
        node [
		shape = circle
		style = filled
		color = "#00000000"
                fillcolor = black
                width = 0.1
                height = 0.1
                label = ""
	]
        """)
        n = len(cosets)
        for cl, perm in zip(colours, perms):
            for i, j in enumerate(perm):
                if i > j:
                    continue
                if cl is not None:
                    lines.append('    %d -- %d [color="%s"];'%(i, j, cl))
        lines.append("}")
        s = '\n'.join(lines)
        if name is None:
            return s
        f = open(name+".dot", 'w')
        print(s, file=f)
        f.close()
        print("Schreier.get_dot: wrote", name)
        import os
        os.popen("dot -Tpdf %s.dot > %s.pdf" % (name,name))

    # -----------------------------------------------------------
    # constructors

    @staticmethod
    def make_reflection(N, links, build=True):
        #print("make_reflection", links)
        rels = []
        for i in range(N):
            rels.append((i, i))
            for j in range(i+1, N):
                m = links.get((i, j), 2)
                rels.append( (i, j) * m)
        #print("make_reflection", rels)
        graph = Schreier(N, rels)
        if build:
            graph.build()
        return graph

    @staticmethod
    def make_A(N):
        "A series coxeter group"
        links = {}
        for i in range(N-1):
            links[(i, i+1)] = 3
        return Schreier.make_reflection(N, links)

    @staticmethod
    def make_B(N):
        "B/C series coxeter group"
        links = {}
        for i in range(N-1):
            links[(i, i+1)] = 4 if i==N-2 else 3
        return Schreier.make_reflection(N, links)
    make_C = make_B

    @staticmethod
    def make_D(N):
        links = {}
        for i in range(N-2):
            links[(i, i+1)] = 3
        links[(N-3, N-1)] = 3
        return Schreier.make_reflection(N, links)

#    @staticmethod
#    def make_E(N):
#        "E series coxeter group"
#        links = {}
#        for i in range(N-2):
#            links[(i, i+1)] = 3
#        links[2, N-1] = 3
#        return Schreier.make_reflection(N, links)

    @staticmethod
    def make_E(N):
        "E series coxeter group"
        links = {}
        for i in range(1, N-1):
            links[(i, i+1)] = 3
        links[0, 3] = 3
        return Schreier.make_reflection(N, links)

    @staticmethod
    def make_F(N=4):
        "F series coxeter group"
        links = {}
        assert N==4
        for i in range(N):
            links[(i, i+1)] = 4 if i==1 else 3
        return Schreier.make_reflection(N, links)

    @staticmethod
    def make_G(N=2):
        "G series coxeter group"
        assert N==2
        links = {(0,1):6}
        return Schreier.make_reflection(N, links)

    
def test():
    ngens = 4 # include the inverses
    rels = [
        (1, 0), # a^-1a
        (3, 2), # b^-1b
        (0, 0, 0), #a^3
        (2, 2), # b^2
        (0, 2, 0, 2) # abab
    ]
    graph = Schreier(ngens, rels)
    graph.build()
    #print(len(graph))
    #graph.show()
    assert len(graph) == 6 # S_3

    ngens = 2
    a, ai, b, bi = range(2*ngens)
    rels = [ (ai, a), (bi, b), ]
    rels += [ (a,a), (b,b), (a,b)*3 ]
    graph = Schreier(2*ngens, rels)
    graph.build()
    assert len(graph) == 6 # S_3
    #print(len(graph.labels))
    
    ngens = 3
    a, ai, b, bi, c, ci = range(2*ngens)
    rels = [ (ai, a), (bi, b), (c,ci)]
    rels += [ (a,a), (b,b), (c,c), (a,c)*2, (a,b)*3, (b,c)*3 ]
    graph = Schreier(2*ngens, rels)
    graph.build()
    assert len(graph) == 24 # S_4
    #print(len(graph.labels))

    ngens = 4
    a, ai, b, bi, c, ci, d, di = range(2*ngens)
    rels = [ (ai, a), (bi, b), (c,ci), (d,di)]
    rels += [ (a,a), (b,b), (c,c), (d,d), (a,c)*2, (a,d)*2, (b,d)*2, (a,b)*3, (b,c)*4, (c,d)*3 ]
    graph = Schreier(2*ngens, rels)
    graph.build()
    assert len(graph) == 1152 # F_4
    #print(len(graph.labels))
    
    if argv.slow:
        ngens = 4
        a, ai, b, bi, c, ci, d, di = range(2*ngens)
        rels = [ (ai, a), (bi, b), (c,ci), (d,di)]
        rels += [ (a,a), (b,b), (c,c), (d,d), (a,c)*2, (a,d)*2, (b,d)*2, (a,b)*5, (b,c)*3, (c,d)*3 ]
        graph = Schreier(2*ngens, rels)
        graph.build()
        assert len(graph) == 14400 # H_4
        print(len(graph.labels))
    
    #return

    # Klein Quartic
    rels = [ (ai, a), (bi, b), (a,)*3, (b,)*7, (a,b)*2 ]
    rels += [ (a,bi,bi)*4 ]
    graph = Schreier(2*ngens, rels)
    graph.build()
    assert len(graph) == 168, len(graph)

    ngens = 2
    a, ai, b, bi = range(2*ngens)
    rels = [ (ai, a), (bi, b), (a,)*5, (b,b), (b,a)*5]
    #rels += [ (a,a,a,b)*4 ] # 360
    rels += [ (a,a,a,b,a,b,a)*2 ]
    graph = Schreier(2*ngens, rels)
    graph.build()
    #print(len(graph))
    assert len(graph) == 80

    # Bring's curve reflection group
    ngens = 3
    a, ai, b, bi, c, ci = range(2*ngens)
    rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*5, (b,c)*5, (a,c)*2,
    ]
    a1 = (b,a)
    b1 = (a,c)
    rels += [ (3*a1+b1)*3 ]
    graph = Schreier(2*ngens, rels)
    graph.build()
    assert len(graph) == 120 # == 12 * 10
    G = graph.get_group()
    assert len(G) == 120

    #graph = Schreier(6, [
    #    (ai, a), (bi, b), (ci, c), 
    #    (a,a), (b,)*5, (c,)*3, (a,b)*3, (a,c)*3, (b,b,c)*2])
    #graph.build()
    #print(len(graph))

    ngens = 8
    rels = [(1, 0), (3, 2), (5, 4), (0, 0), (2, 2), (4, 4),
        (6, 7), (7, 4, 0, 4), (4, 6, 4, 6, 4, 6, 4, 6), (0, 4,
        0, 4, 0, 4, 0, 4), (2, 4, 2, 4, 2, 4, 2, 4), (0, 2, 0,
        2, 0, 2, 0, 2, 0, 2), (0, 4, 0, 4, 0, 4, 0, 4), (2, 4,
        2, 4), (2, 0, 2, 0, 2, 0, 0, 6, 2, 0, 2, 0, 2, 0, 0,
        6, 2, 0, 2, 0, 2, 0, 0, 6)]
    graph = Schreier(ngens, rels)
    graph.build()

    assert len(graph) == 240


def test_coxeter():

    for make, N, size in [
        (Schreier.make_A, 4, 120),
        (Schreier.make_B, 3, 48),
        (Schreier.make_B, 4, 384),
        (Schreier.make_B, 5, 3840),
        (Schreier.make_D, 4, 192),
        (Schreier.make_D, 5, 1920),
        (Schreier.make_E, 6, 51840), # slow...
        #(Schreier.make_E, 7, 2903040), # too big
        #(Schreier.make_E, 8, 696729600), # way too big
        (Schreier.make_F, 4, 1152),
        (Schreier.make_G, 2, 12),
    ]:
        #print()
        print(make, N, size)
        graph = make(N)
        assert len(graph) == size, ("%s != %s"%(len(graph), size))

    if 0:
        print("go...")
        N = 6
        graph = Schreier.make_E(N)
        rels = list(graph.rels)
        a = N
        b = N+1
        rels.append((a, a))
        rels.append((b, b))
        for i in [0, 2, 3, 4]:
            rels.append((a, i)*2) # commutes
        for i in [0, 1, 3, 5]:
            rels.append((b, i)*2) # commutes
    
        rels.append((a, 1, a, 5))
        rels.append((b, 2, b, 4))
        rels.append((a, b)*2)
    
        graph = Schreier(N+2, rels)
        graph.build()
        print(len(graph))
    
        return

    # fold D_4 fishtail to get B_4
    N = 4
    graph = Schreier.make_D(N)
    rels = list(graph.rels)
    # add another reflection
    rels.append((N, N))
    for i in range(N-2):
        rels.append((i, N)*2)
    rels.append((N, N-2, N, N-1))
    graph = Schreier(N+1, rels)
    graph.build()
    assert len(graph) == 384, len(graph)

    # fold again 
    N = 5
    rels = list(graph.rels)
    # add another reflection
    rels.append((N, N))
    for j in [1, 3]:
        rels.append((i, N)*2)
    rels.append((N, 0, 2, N-1))
    graph = Schreier(N+1, rels)
    graph.build()
    assert len(graph) == 12, len(graph) # G_2 !!
    G = graph.get_group()
    gens = G.gens
    #for a in G:
    #  for b in G:
    #    print( (a*b).order(), end = " " )
    #  print()
    
    # fold D_5 fishtail to get B_5
    N = 5
    graph = Schreier.make_D(N)
    rels = list(graph.rels)
    # add another reflection
    rels.append((N, N))
    for i in range(N-2):
        rels.append((i, N)*2)
    rels.append((N, N-2, N, N-1))
    graph = Schreier(N+1, rels)
    graph.build()
    assert len(graph) == 3840


def test_F4():

    ngens = 4
    a, b, c, d = range(ngens)
    rels = [(a, a), (b, b), (c,c), (d,d)]
    rels += [(a,a), (b,b), (c,c), (d,d), (a,c)*2, (a,d)*2, (b,d)*2, (a,b)*3, (b,c)*4, (c,d)*3]
    graph = Schreier(ngens, rels)
    graph.build()
    assert len(graph) == 1152 # F_4
    
    gens = [(a,), (b,), (c,), (d,)]

    colours = "#FF0000 #3333B2 #009900 #CCCCCC".split()
    colours[-1] = None

    for i in range(ngens):
        graph = Schreier(ngens, rels)
        gs = list(gens)
        gs.pop(i)
        graph.build(gs)
        cs = list(colours)
        #cs[i] = None
        s = graph.get_dot(cs)
        name = "schreier_F4_%d"%(i,)

        f = open(name+".dot", 'w')
        print(s, file=f)
        f.close()
        print("wrote", name)
        import os
        os.popen("dot -Tpdf %s.dot > %s.pdf" % (name,name))
    

if __name__ == "__main__":

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
        test_coxeter()

    print("OK\n")

        
