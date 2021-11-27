#!/usr/bin/env python3

import sys
from random import randint, seed, choice


from bruhat.gset import Perm, Group, Coset, mulclose
from bruhat.argv import argv

#sys.setrecursionlimit(8000)

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
                assert 0 <= gen < ngens, repr(rel)

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
            #assert len(neighbors)  < 60
            if maxsize and len(neighbors) > maxsize:
                return False
        return True
        
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
        for i in range(ngens):
            perm = [None]*n
            for idx in cosets:
                nbd = neighbors[idx]
                assert nbd[i] is not None, nbd
                src, tgt = lookup[idx], lookup[self.get_label(nbd[i])]
                assert perm[src] is None
                perm[src] = tgt
            assert None not in perm
            perm = Perm(perm)
            gens.append(perm)
        return gens

    def get_group(self):
        gens = self.get_gens()
        G = Group(None, gens)
        return G
    
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
    
        #for d in range(self.ngens):
        #    print("g%d ="%d, cycle(perms[d]))

    
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
    
    ngens = 3
    a, ai, b, bi, c, ci = range(2*ngens)
    rels = [ (ai, a), (bi, b), (c,ci)]
    rels += [ (a,a), (b,b), (c,c), (a,c)*2, (a,b)*3, (b,c)*3 ]
    graph = Schreier(2*ngens, rels)
    graph.build()
    assert len(graph) == 24 # S_4
    
    ngens = 4
    a, ai, b, bi, c, ci, d, di = range(2*ngens)
    rels = [ (ai, a), (bi, b), (c,ci), (d,di)]
    rels += [ (a,a), (b,b), (c,c), (d,d), (a,c)*2, (a,d)*2, (b,d)*2, (a,b)*3, (b,c)*4, (c,d)*3 ]
    graph = Schreier(2*ngens, rels)
    graph.build()
    assert len(graph) == 1152 # F_4
    
    if argv.slow:
        ngens = 4
        a, ai, b, bi, c, ci, d, di = range(2*ngens)
        rels = [ (ai, a), (bi, b), (c,ci), (d,di)]
        rels += [ (a,a), (b,b), (c,c), (d,d), (a,c)*2, (a,d)*2, (b,d)*2, (a,b)*5, (b,c)*3, (c,d)*3 ]
        graph = Schreier(2*ngens, rels)
        graph.build()
        assert len(graph) == 14400 # H_4
    

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
    hgens = []
    graph = Schreier(ngens, rels)
    graph.build()

    assert len(graph) == 240


    
def make_random_modular():
    ngens = 4
    a, ai, b, bi = range(ngens)
    rels = [ (ai, a), (bi, b), (a,)*2, (a,b,)*3, ] # modular group 
    rels += [ (b,)*7, ] # 2,3,7 triangle group

    for _ in range(10000):
        graph = Schreier(ngens, rels)
        hgens = []
        for i in range(2):
            gen = tuple(randint(0, ngens-1) for k in range(20))
            hgens.append(gen)
        if graph.build(hgens, maxsize=2000):
            if len(graph) < 3:
                continue
            #print("hgens:", hgens)
            gens = graph.get_gens()
            G = mulclose(gens, maxsize=10000)
            if len(G) < 10000:
                print(len(graph), end=" ", flush=True)
                print(len(G))


def make_random_55():

    #seed(4)

    ngens = 6
    a, ai, b, bi, c, ci = range(ngens)
    i_rels = [
        (ai, a), (bi, b), (ci, c), 
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*5, (b,c)*5, (a,c)*2,
    ]
    a1 = (b,a)
    b1 = (a,c)
    rel1 = (3*a1+b1)*3 # Bring's curve

    for _ in range(10000):
        rels = list(i_rels)
        for i in range(1):
            #gen = tuple([a, b, c][randint(0, 2)] for k in range(30))
            gen = ()
            for k in range(50):
                gen += choice([(a,b), (a,c), (b,c)])
            #gen = rel1
            rels.append(gen)
        graph = Schreier(ngens, rels)
        if graph.build(maxsize=640000):
            print(".", flush=True, end="")
            n = len(graph)
            if n <= 720:
                continue
            #print("rels:", rels)
            #graph.show()
            print(len(graph.neighbors))
            try:
                gens = graph.get_gens()
            except AssertionError:
                print("XXX")
                continue
            G = mulclose(gens, maxsize=10000)
            if len(G) < 10000:
                print(gen)
                print(len(graph), end=" ", flush=True)
                print(len(G))
            print()



def make_bring():

    # ------------------------------------------------

    # Bring's curve rotation group
    ngens = 2
    a, ai, b, bi = range(2*ngens)
    rels = [ (ai, a), (bi, b), (a,)*5, (b,)*2]
    rels += [ (a,a,a,b)*3 ]
    graph = Schreier(2*ngens, rels)
    #graph.DEBUG = True
    graph.build()
    assert len(graph) == 60 # == 12 * 5

    # --------------- Group ------------------
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
    #print(shortstr(Hz))
    #print()
    #print(shortstr(Hxt))
    #print()

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
    rels += [ (3*a1+b1)*3 ]
    graph = Schreier(2*ngens, rels)
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

    # Here we choose an orientation on each edge:
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

    #print(shortstr(Hz))
    #print()
    #print(shortstr(Hxt))
    #print()

    assert alltrue(dot(Hz, Hxt)==0)


def make_codes_55():
    # start with hyperbolic Coxeter reflection group: a--5--b--5--c
    # Then we add another generator "d" that halves these triangles.
    ngens = 8
    a, ai, b, bi, c, ci, d, di = range(ngens)
    rels_552 = [
        (ai, a), (bi, b), (ci, c), (di, d),
        (a,)*2, (b,)*2, (c,)*2, (d,)*2,
        (a,b)*5, (b,c)*5, (a,c)*2,
        (ci,d,a,d),
        (a,d)*4, (b,d)*2, (d,c)*4,
    ]

    order = argv.get("order")

    # Bring's curve & some random generators found from running make_random_55.
    for rel in [
        (3*(b,a)+(a,c))*3,           # 120  [[30,8,3]]
        (0, 2, 0, 2, 2, 4, 2, 4, 0, 2, 0, 2, 0, 4, 2, 4, 0, 4,
            0, 4, 0, 2, 2, 4, 0, 2, 0, 4, 0, 4, 0, 2, 0, 4, 0, 2,
            0, 4, 0, 2),             # 160  [[40,10,?]]
        (2, 4, 0, 2, 0, 2, 0, 2, 0, 2, 2, 4, 0, 2, 0, 4,
            0, 2, 2, 4, 2, 4, 2, 4, 0, 2, 0, 2, 0, 4, 2, 4, 0, 4,
            0, 4, 0, 2, 0, 4, 0, 2, 2, 4, 0, 2, 0, 2, 0, 4, 0, 2,
            0, 4, 0, 4, 2, 4, 0, 4), # 320  [[80,18,?]]
        (0, 2, 0, 2, 0, 4, 0, 2, 2, 4, 0, 2, 0, 4, 0, 2, 0, 4,
            2, 4, 2, 4, 0, 2, 0, 4, 2, 4, 0, 4, 0, 2, 2, 4, 0, 4,
            2, 4, 0, 4, 2, 4, 0, 2, 0, 4, 0, 2, 0, 4, 0, 4, 0, 4,
            2, 4, 0, 4, 2, 4),       # 600 [[150,32,?]]
        (0, 2, 0, 4, 0, 2, 0, 2, 0, 4, 0, 2, 2, 4, 0, 2, 2, 4,
            0, 2, 2, 4, 0, 4, 0, 2, 0, 2, 0, 4, 0, 2, 0, 2, 0, 2,
            2, 4, 0, 4, 2, 4, 0, 4, 2, 4, 0, 4, 0, 2, 0, 2, 2, 4,
            2, 4, 2, 4, 2, 4),       # 720  [[180,38,?]]
        (2, 4, 0, 4, 0, 4, 0, 4, 2, 4, 0, 2, 2, 4, 2, 4, 0, 2,
            0, 2, 2, 4, 0, 2, 2, 4, 0, 4, 0, 4, 0, 2, 2, 4, 0, 4,
            0, 2, 0, 4, 0, 4, 0, 2, 0, 4, 0, 2, 2, 4, 2, 4, 2, 4,
            0, 2, 0, 2, 0, 2, 2, 4, 2, 4, 0, 2, 0, 2, 0, 2, 2, 4,
            0, 2, 2, 4, 0, 2, 0, 4, 0, 2, 2, 4, 0, 4, 0, 4, 2, 4,
            0, 2, 0, 2, 2, 4, 0, 4, 0, 2), # 1320
    ]:
        graph = Schreier(ngens, rels_552 + [rel])
        graph.build()
        G = graph.get_group()

        if order is not None and len(G) != order*2:
            print(len(G))
            continue

        make_surface_54(G)

        if argv.bring:
            break


class Surface(object):
    def __init__(self, G):
        pass


def make_surface_54(G0):
    from bruhat.solve import shortstr, zeros2, dot2
    from numpy import alltrue, zeros, dot

    print()
    print("|G0| =", len(G0))

    a, ai, b, bi, c, ci, d, di = G0.gens

    # 5,5 coxeter subgroup
    G_55 = Group.generate([a, b, c])
    assert G0.is_subgroup(G_55)

    # orientation subgroup
    L = Group.generate([a*b, a*d])
    #assert len(L) == 60, len(L)
    print("L:", len(L))
    orients = G0.left_cosets(L)
    assert len(orients) == 2
    orients = list(orients)

    H = Group.generate([a, b])
    faces = G_55.left_cosets(H)
    fold_faces = G0.left_cosets(H)
    print("faces:", len(faces))
    print("fold_faces:", len(fold_faces))
    #hom = G.left_action(faces)
    o_faces = [face.intersect(orients[0]) for face in faces] # oriented faces
    r_faces = [face.intersect(orients[1]) for face in faces] # reverse oriented faces

    K = Group.generate([b, c])
    vertices = G_55.left_cosets(K)
    #assert len(vertices) == 12 # unoriented vertices
    o_vertices = [vertex.intersect(orients[0]) for vertex in vertices] # oriented vertices
    print("o_vertices:", len(o_vertices))

    J0 = Group.generate([a, d])
    assert len(J0) == 8

    assert G0.is_subgroup(J0)
    act_edges = G0.action_subgroup(J0)
    act_fold_faces = G0.action_subgroup(H)
    #print("stab:", len(act_edges.get_stabilizer()))
    print("|G0| =", len(G0))
    print("edges:", len(G0) // len(J0))
    edges = G0.left_cosets(J0)
    #for e in edges:
    #    k = choice(e.perms)
    #    e1 = e.left_mul(~k)
    #    assert set(e1.perms) == set(J0.perms)
    for i, g in enumerate(G0):
        assert g in G0
        h_edge = act_edges.tgt[act_edges.send_perms[i]]
        h_fold_faces = act_fold_faces.tgt[act_fold_faces.send_perms[i]]
        n_edge = len(h_edge.fixed())
        n_fold_faces = len(h_fold_faces.fixed())
        count = 0
        perms = [] # extract action on the fixed edges
        for e0 in edges:
            e1 = e0.left_mul(g)
            if e0 != e1:
                continue
            perm = e0.left_mul_perm(g)
            #print(perm, end=" ")
            perms.append(perm)
            count += 1
        assert count == n_edge
        print([p.order() for p in perms] or '', end="")
        print("[|g|=%s,%s.%s.%s.%s]"%(
            g.order(),
            n_edge or ' ',  # num fixed edges
            n_fold_faces or ' ',  # num fixed faces+vertexes
            ["    ", "pres"][int(g in G_55)], # preserves face/vertex distinction
            ["refl", "rot "][int(g in L)], # oriented or un-oriented
        ), end=" ", flush=True)
        print()
    print()

    J = Group.generate([a, c])
    u_edges = G_55.left_cosets(J)
    print("u_edges:", len(u_edges))

    J = Group.generate([c])
    edges = G_55.left_cosets(J)
    print("edges:", len(edges))

    # Here we choose an orientation on each edge:
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

    #print(shortstr(Hz))
    #print()
    #print(shortstr(Hxt))
    #print()

    assert alltrue(dot(Hz, Hxt)==0)

    import qupy.ldpc.solve
    import bruhat.solve
    qupy.ldpc.solve.int_scalar = bruhat.solve.int_scalar
    from qupy.ldpc.css import CSSCode
    Hz = Hz % 2
    Hx = Hxt.transpose() % 2
    code = CSSCode(Hz=Hz, Hx=Hx)
    print(code)



#def make_surface(G0):
#    from bruhat.solve import shortstr, zeros2, dot2
#    from numpy import alltrue, zeros, dot
#
#    print()
#
#    a, ai, b, bi, c, ci, d, di = G0.gens
#    G = Group.generate([a, b, c])
#
#    a1 = b*a
#    b1 = a*c
#
#    bc = b*c
#    ba = b*a
#    L = Group.generate([a1, b1])
#    #assert len(L) == 60, len(L)
#    print("L:", len(L))
#    orients = G.left_cosets(L)
#    #assert len(orients) == 2
#    #L, gL = orients
#    orients = list(orients)
#
#    H = Group.generate([a, b])
#    #assert len([g for g in H if g in L]) == 5
#    faces = G.left_cosets(H)
#    #assert len(faces) == 12 # unoriented faces
#    print("faces:", len(faces))
#    #hom = G.left_action(faces)
#    o_faces = [face.intersect(orients[0]) for face in faces] # oriented faces
#    r_faces = [face.intersect(orients[1]) for face in faces] # reverse oriented faces
#
#    K = Group.generate([b, c])
#    vertices = G.left_cosets(K)
#    #assert len(vertices) == 12 # unoriented vertices
#    o_vertices = [vertex.intersect(orients[0]) for vertex in vertices] # oriented vertices
#    print("o_vertices:", len(o_vertices))
#
#    J0 = Group.generate([a, d])
#    assert len(J0) == 8
#
#    assert G0.is_subgroup(J0)
#    gset = G0.action_subgroup(J0)
#    #print("stab:", len(gset.get_stabilizer()))
#    print("|G0| =", len(G0))
#    print("edges:", len(G0) // len(J0))
#    for i, g in enumerate(gset.src):
#        #if g in L:
#        #    continue
#        j = gset.send_perms[i]
#        h = gset.tgt[j]
#        n = len(h.fixed())
#        print(n or '.', end=" ")
#    print()
#
#    J = Group.generate([a, c])
#    u_edges = G.left_cosets(J)
#    print("u_edges:", len(u_edges))
#
#    J = Group.generate([c])
#    edges = G.left_cosets(J)
#    print("edges:", len(edges))
#
#    # Here we choose an orientation on each edge:
#    pairs = {}
#    for e in u_edges:
#        pairs[e] = []
#        for e1 in edges:
#            if e.intersect(e1):
#                pairs[e].append(e1)
#        assert len(pairs[e]) == 2
#
#    def shortstr(A):
#        s = str(A)
#        s = s.replace("0", '.')
#        return s
#
#    Hz = zeros((len(o_faces), len(u_edges)), dtype=int)
#    for i, l in enumerate(o_faces):
#      for j, e in enumerate(u_edges):
#        le = l.intersect(e)
#        if not le:
#            continue
#        e1, e2 = pairs[e]
#        if e1.intersect(le):
#            Hz[i, j] = 1
#        elif e2.intersect(le):
#            Hz[i, j] = -1
#        else:
#            assert 0
#
#    Hxt = zeros((len(u_edges), len(o_vertices)), dtype=int)
#    for i, e in enumerate(u_edges):
#      for j, v in enumerate(o_vertices):
#        ev = e.intersect(v)
#        if not ev:
#            continue
#        e1, e2 = pairs[e]
#        if e1.intersect(ev):
#            Hxt[i, j] = 1
#        elif e2.intersect(ev):
#            Hxt[i, j] = -1
#        else:
#            assert 0
#
#    #print(shortstr(Hz))
#    #print()
#    #print(shortstr(Hxt))
#    #print()
#
#    assert alltrue(dot(Hz, Hxt)==0)
#
#    import qupy.ldpc.solve
#    import bruhat.solve
#    qupy.ldpc.solve.int_scalar = bruhat.solve.int_scalar
#    from qupy.ldpc.css import CSSCode
#    Hz = Hz % 2
#    Hx = Hxt.transpose() % 2
#    code = CSSCode(Hz=Hz, Hx=Hx)
#    print(code)
#


if __name__ == "__main__":

    profile = argv.profile
    name = argv.next()

    if profile:
        import cProfile as profile
        profile.run("make_bring()")

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        test()
        make_bring()

    print("OK\n")

        
