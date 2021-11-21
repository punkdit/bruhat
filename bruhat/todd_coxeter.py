#!/usr/bin/env python3

# From: https://math.berkeley.edu/~kmill/notes/todd_coxeter.html
# Example of Todd-Coxeter to compute S_3/<b>=<a,b|a^3,b^2,abab> from relations

# The permutations can be written in cycle notation fairly
# easily. One way is with
    
def cycle(perm):
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
    def __init__(self, ngens, rels, hgens):
        self.labels = []
        self.neighbors = []
        self.ngens = ngens
        self.rels = rels
        self.hgens = hgens

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

    def find(self, c):
        labels = self.labels
        c2 = labels[c]
        if c == c2:
            return c
        else:
            c2 = self.find(c2) # recurse
            labels[c] = c2
            return c2
    
    # The identification procedure is called unify. 
    
    def unify(self, c1, c2):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        c1 = self.find(c1)
        c2 = self.find(c2)
        if c1 == c2:
            return
        c1, c2 = min(c1, c2), max(c1, c2)
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
    
    # For following paths, we have the following two functions: 
    
    def follow(self, c, d):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        c = self.find(c)
        ns = neighbors[c]
        if ns[d] == None:
            ns[d] = self.new()
        return self.find(ns[d])
    
    def followp(self, c, ds):
        c = self.find(c)
        for d in reversed(ds):
            c = self.follow(c, d)
        return c
    
    #    The first takes a vertex and finds the neighbor in the
    #    d direction, creating a new vertex in that direction
    #    if needed, and the second follows a list of directions
    #    to find the end of a path. The follow function creates
    #    new neighbors with the new() function:
    
    def new(self):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        c = len(labels)
        labels.append(c)
        neighbors.append((2*ngens)*[None])
        return c
    
    #    Finally, there is the core algorithm: 
    def build(self):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        rels = self.rels
        hgens = self.hgens

        start = self.new()
        
        for hgen in hgens:
            self.unify(self.followp(start, hgen), start)
        
        to_visit = 0
        while to_visit < len(labels):
            c = self.find(to_visit)
            if c == to_visit:
                for rel in rels:
                    self.unify(self.followp(c, rel), c)
            to_visit += 1
        
        print("done")

    #    It creates the start vertex, adds all of the relations
    #    for H as relators at this basepoint, and then for each
    #    unmarked vertex, adds all of the relations from rels.
    #    Notice how unify is being used to turn paths into relators.
    #    
    
    def dump(self):
        #    After this, the data structures contain the Schreier
        #    graph for G/H .  This can be interpreted as a permutation
        #    representation, for instance using
        
        cosets = [c for i, c in enumerate(self.labels) if i == c]
        
        perms = [[cosets.index(self.follow(c, 2*d)) for i, c in enumerate(cosets)]
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
        (2,), # b
    ]
    
    graph = Schreier(ngens, rels, hgens)
    graph.build()
    graph.dump()


if __name__ == "__main__":

    test()


        
