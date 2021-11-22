#!/usr/bin/env python3

from bruhat.gset import Perm, Group
from bruhat.argv import argv

# Todd-Coxeter algorithm: compute finite group from generators and relations.
# Implementation from: 
# https://math.berkeley.edu/~kmill/notes/todd_coxeter.html

# See also:
# https://arxiv.org/pdf/2012.09271.pdf p10.

    
def cycle(perm):
    # The permutations can be written in cycle notation fairly
    # easily. One way is with
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
    of a normal subgroup.
    """
    DEBUG = False

    def __init__(self, ngens, rels, hgens=[]):
        self.labels = []
        self.neighbors = []
        self.ngens = ngens
        self.rels = list(rels)
        self.hgens = list(hgens)

        #The labels variable is a list of _numbers, with the property
        #that labels[i] <= i. This is a union-find data structure
        #for keeping track of the vertex quotients for the Schreier
        #graph so far. The find function looks up the current-lowest
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
        
        #We use 0 and 1 for a and a−1, and 2 and 3 for b and b−1.
        #For simplicity of the core algorithm, we introduce a−1a
        #and b−1b as explicit relations.

    def get_label(self, c):
        labels = self.labels
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
    
    def unify(self, c1, c2):
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
        for d in range(ngens):
            n1 = neighbors[c1][d]
            n2 = neighbors[c2][d]
            if n1 == None:
                neighbors[c1][d] = n2
            elif n2 != None:
                self.unify(n1, n2) # recurse
    
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
    
    def follow_path(self, c, ds):
        c = self.get_label(c)
        for d in reversed(ds):
            c = self.follow_step(c, d)
        return c
    
    #    The first takes a vertex and finds the neighbor in the
    #    d direction, creating a new vertex in that direction
    #    if needed, and the second follows a list of directions
    #    to find the end of a path. The follow function creates
    #    new neighbors with the add_vertex() function.

#    # Uses bruhat.action version of Group's and Perm's:
#    def get_group(self):
#        labels = self.labels
#        neighbors = self.neighbors
#        ngens = self.ngens
#        items = [] # act on these
#        for idx, row in enumerate(neighbors):
#            jdx = self.get_label(idx)
#            if idx == jdx:
#                items.append(idx)
#        gens = []
#        for i in range(ngens):
#            perm = {}
#            for idx in items:
#                nbd = neighbors[idx]
#                perm[idx] = self.get_label(nbd[i])
#            perm = Perm(perm, items)
#            gens.append(perm)
#        G = Group.generate(gens)
#        return G
    
    def get_group(self):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        items = [] # act on these
        for idx, row in enumerate(neighbors):
            jdx = self.get_label(idx)
            if idx == jdx:
                items.append(idx)
        n = len(items)
        lookup = dict((v,k) for (k,v) in enumerate(items))
        gens = []
        for i in range(ngens):
            perm = [None]*n
            for idx in items:
                nbd = neighbors[idx]
                src, tgt = lookup[idx], lookup[self.get_label(nbd[i])]
                assert perm[src] is None
                perm[src] = tgt
            assert None not in perm
            perm = Perm(perm)
            gens.append(perm)
        G = Group(None, gens)
        return G
    
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
    def build(self):
        labels = self.labels
        neighbors = self.neighbors
        rels = self.rels
        hgens = self.hgens

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
                for rel in rels + hgens:
                    self.unify(self.follow_path(c, rel), c)
                if self.DEBUG:
                    self.dump()
            to_visit += 1
            #assert len(neighbors)  < 60
        
    #    It creates the start vertex, adds all of the relations
    #    for H as relators at this basepoint, and then for each
    #    unmarked vertex, adds all of the relations from rels.
    #    Notice how unify is being used to turn paths into relators.
    #    
    
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

    
def test():
    ngens = 4 # include the inverses
    rels = [
        (1, 0), # a^-1a
        (3, 2), # b^-1b
        (0, 0, 0), #a^3
        (2, 2), # b^2
        (0, 2, 0, 2) # abab
    ]
    hgens = [
        #(2,), # b
    ]
    graph = Schreier(ngens, rels, hgens)
    graph.build()
    #print(len(graph))
    #graph.show()
    assert len(graph) == 6 # S_3

    ngens = 2
    a, ai, b, bi = range(2*ngens)
    rels = [ (ai, a), (bi, b), ]
    hgens = [ (a,a), (b,b), (a,b)*3 ]
    graph = Schreier(2*ngens, rels, hgens)
    graph.build()
    assert len(graph) == 6 # S_3
    
    ngens = 3
    a, ai, b, bi, c, ci = range(2*ngens)
    rels = [ (ai, a), (bi, b), (c,ci)]
    hgens = [ (a,a), (b,b), (c,c), (a,c)*2, (a,b)*3, (b,c)*3 ]
    graph = Schreier(2*ngens, rels, hgens)
    graph.build()
    assert len(graph) == 24 # S_4
    
    ngens = 4
    a, ai, b, bi, c, ci, d, di = range(2*ngens)
    rels = [ (ai, a), (bi, b), (c,ci), (d,di)]
    hgens = [ (a,a), (b,b), (c,c), (d,d), (a,c)*2, (a,d)*2, (b,d)*2, (a,b)*3, (b,c)*4, (c,d)*3 ]
    graph = Schreier(2*ngens, rels, hgens)
    graph.build()
    assert len(graph) == 1152 # F_4
    
    if argv.slow:
        ngens = 4
        a, ai, b, bi, c, ci, d, di = range(2*ngens)
        rels = [ (ai, a), (bi, b), (c,ci), (d,di)]
        hgens = [ (a,a), (b,b), (c,c), (d,d), (a,c)*2, (a,d)*2, (b,d)*2, (a,b)*5, (b,c)*3, (c,d)*3 ]
        graph = Schreier(2*ngens, rels, hgens)
        graph.build()
        assert len(graph) == 14400 # H_4
    

    # Klein Quartic
    rels = [ (ai, a), (bi, b), (a,)*3, (b,)*7, (a,b)*2 ]
    hgens = [ (a,bi,bi)*4 ]
    graph = Schreier(2*ngens, rels, hgens)
    graph.build()
    assert len(graph) == 168, len(graph)

    ngens = 2
    a, ai, b, bi = range(2*ngens)
    rels = [ (ai, a), (bi, b), (a,)*5, (b,b), (b,a)*5]
    #hgens = [ (a,a,a,b)*4 ] # 360
    hgens = [ (a,a,a,b,a,b,a)*2 ]
    graph = Schreier(2*ngens, rels, hgens)
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
    hgens = [ (3*a1+b1)*3 ]
    graph = Schreier(2*ngens, rels, hgens)
    graph.build()
    assert len(graph) == 120 # == 12 * 10
    G = graph.get_group()
    assert len(G) == 120


def make_bring():
    # Bring's curve rotation group
    ngens = 2
    a, ai, b, bi = range(2*ngens)
    rels = [ (ai, a), (bi, b), (a,)*5, (b,b), (b,a)*5]
    hgens = [ (a,a,a,b)*3 ]
    graph = Schreier(2*ngens, rels, hgens)
    #graph.DEBUG = True
    graph.build()
    assert len(graph) == 60 # == 12 * 5
    G = graph.get_group()
    #burnside(G)
    assert len(G.gens) == 4
    ai, a, bi, b = G.gens
    assert a==~ai
    assert b==~bi

    assert a.order() == 5
    assert b.order() == 2
    H = Group.generate([a])

    faces = G.left_cosets(H)
    assert len(faces) == 12
    #hom = G.left_action(faces)

    ab = a*b
    K = Group.generate([ab])
    vertices = G.left_cosets(K)
    assert len(vertices) == 12

    J = Group.generate([b])
    edges = G.left_cosets(J)
    assert len(edges) == 30

    from bruhat.solve import shortstr, zeros2, dot2
    from numpy import alltrue, zeros, dot

    def get_adj(left, right):
        A = zeros2((len(left), len(right)))
        for i, l in enumerate(left):
          for j, r in enumerate(right):
            lr = l.intersect(r)
            A[i, j] = len(lr)
        return A

    Hz = get_adj(faces, edges)
    Hxt = get_adj(edges, vertices)
    print(shortstr(Hz))
    print()
    print(shortstr(Hxt))
    print()

    assert alltrue(dot2(Hz, Hxt)==0)

    # ------------ _oriented version ---------------
    
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
    hgens = [ (3*a1+b1)*3 ]
    graph = Schreier(2*ngens, rels, hgens)
    graph.build()
    assert len(graph) == 120 # == 12 * 10

    # --------------- Group ------------------
    G = graph.get_group()
    assert len(G) == 120

    a, ai, b, bi, c, ci = G.gens
    a1 = b*a
    b1 = a*c

    assert a1.order() == 5
    assert b1.order() == 2

    bc = b*c
    ba = b*a
    assert bc.order() == 5
    assert ba.order() == 5
    L = Group.generate([a1, b1])
    assert len(L) == 60, len(L)
    orients = G.left_cosets(L)
    assert len(orients) == 2
    #L, gL = orients
    orients = list(orients)

    H = Group.generate([a, b])
    assert len([g for g in H if g in L]) == 5
    faces = G.left_cosets(H)
    assert len(faces) == 12 # unoriented faces
    #hom = G.left_action(faces)
    o_faces = [face.intersect(orients[0]) for face in faces] # oriented faces
    r_faces = [face.intersect(orients[1]) for face in faces] # reverse oriented faces

    K = Group.generate([b, c])
    vertices = G.left_cosets(K)
    assert len(vertices) == 12 # unoriented vertices
    o_vertices = [vertex.intersect(orients[0]) for vertex in vertices] # oriented vertices

    J = Group.generate([a, c])
    u_edges = G.left_cosets(J)
    assert len(u_edges) == 30 # unoriented edges

    J = Group.generate([c])
    edges = G.left_cosets(J)
    assert len(edges) == 60, len(edges) # oriented edges ?

    pairs = {}
    for e in u_edges:
        pairs[e] = []
        for e1 in edges:
            if e.intersect(e1):
                pairs[e].append(e1)
        assert len(pairs[e]) == 2

    def shortstr(A):
        s = str(A)
        s = s.replace("0", '.')
        return s

    Hz = zeros((len(o_faces), len(u_edges)), dtype=int)
    for i, l in enumerate(o_faces):
      for j, e in enumerate(u_edges):
        le = l.intersect(e)
        if not le:
            continue
        e1, e2 = pairs[e]
        if e1.intersect(le):
            Hz[i, j] = 1
        elif e2.intersect(le):
            Hz[i, j] = -1
        else:
            assert 0

    Hxt = zeros((len(u_edges), len(o_vertices)), dtype=int)
    for i, e in enumerate(u_edges):
      for j, v in enumerate(o_vertices):
        ev = e.intersect(v)
        if not ev:
            continue
        e1, e2 = pairs[e]
        if e1.intersect(ev):
            Hxt[i, j] = 1
        elif e2.intersect(ev):
            Hxt[i, j] = -1
        else:
            assert 0

    print(shortstr(Hz))
    print()
    print(shortstr(Hxt))
    print()

    assert alltrue(dot(Hz, Hxt)==0)



if __name__ == "__main__":


    if argv.profile:
        import cProfile as profile
        profile.run("make_bring()")

    else:
        test()
        make_bring()

    print("OK\n")

        
