#!/usr/bin/env python3

# Todd-Coxeter algorithm: compute finite group from generators and relations.
# From: https://math.berkeley.edu/~kmill/notes/todd_coxeter.html

    
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
    DEBUG = False

    def __init__(self, ngens, rels, hgens=[]):
        self.labels = []
        self.neighbors = []
        self.ngens = ngens
        self.rels = list(rels)
        self.hgens = list(hgens)

        #The labels variable is a list of numbers, with the property
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

    def find_label(self, c):
        labels = self.labels
        c2 = labels[c]
        if c == c2:
            return c
        else:
            c2 = self.find_label(c2) # recurse
            labels[c] = c2
            return c2
    
    # The identification procedure is called unify. 
    
    def unify(self, c1, c2):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        c1 = self.find_label(c1)
        c2 = self.find_label(c2)
        if c1 == c2:
            return
        c1, c2 = min(c1, c2), max(c1, c2)
        if self.DEBUG:
            print("unify:", c2, "-->", c1)
        labels[c2] = c1
        for d in range(2*ngens):
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
        neighbors.append((2*ngens)*[None])
        return c

    # For following paths, we have the following two functions: 
    
    def follow_step(self, c, d):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        c = self.find_label(c)
        ns = neighbors[c]
        if ns[d] == None:
            ns[d] = self.add_vertex()
        return self.find_label(ns[d])
    
    def follow_path(self, c, ds):
        c = self.find_label(c)
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
        ngens = self.ngens
        print(labels)
        for idx, row in enumerate(neighbors):
            jdx = self.find_label(idx)
            if idx == jdx:
                print('\t', idx, row)
        print("-"*70)

    def __len__(self):
        return len([k for k in enumerate(self.labels) if k[0]==k[1]])
    
    #    Finally, there is the core algorithm: 
    def build(self):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
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
            c = self.find_label(to_visit)
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
        
        perms = [[cosets.index(self.follow_step(c, 2*d)) for i, c in enumerate(cosets)]
                 for d in range(self.ngens)]
        
        #    to enumerate the cosets (which are vertices which have
        #    the least label in their equivalency classes), and then
        #    taking edges of the Schreier graph for each generator
        #    to construct a permutation.
        #    
    
        for d in range(self.ngens):
            print("g%d ="%d, cycle(perms[d]))

    
def test():
    ngens = 2
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
    graph = Schreier(ngens, rels, hgens)
    graph.build()
    assert len(graph) == 6 # S_3
    
    ngens = 3
    a, ai, b, bi, c, ci = range(2*ngens)
    rels = [ (ai, a), (bi, b), (c,ci)]
    hgens = [ (a,a), (b,b), (c,c), (a,c)*2, (a,b)*3, (b,c)*3 ]
    graph = Schreier(ngens, rels, hgens)
    graph.build()
    assert len(graph) == 24 # S_4
    
    # Bring's curve rotation group
    ngens = 2
    a, ai, b, bi = range(2*ngens)
    rels = [ (ai, a), (bi, b), (a,)*5, (b,b), (b,a)*5]
    hgens = [ (a,a,a,b)*3 ]
    graph = Schreier(ngens, rels, hgens)
    #graph.DEBUG = True
    graph.build()
    assert len(graph) == 60

    ngens = 2
    a, ai, b, bi = range(2*ngens)
    rels = [ (ai, a), (bi, b), (a,)*5, (b,b), (b,a)*5]
    hgens = [ (a,a,a,b)*4 ]
    graph = Schreier(ngens, rels, hgens)
    #graph.DEBUG = True
    graph.build()
    assert len(graph) == 360

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
    graph = Schreier(ngens, rels, hgens)
    graph.build()
    assert len(graph) == 120


if __name__ == "__main__":

    test()

    print("OK\n")

        
