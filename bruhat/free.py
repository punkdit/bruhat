#!/usr/bin/env python3

"""
i dunno what i'm doing here...

"""

import sys
from random import randint, seed, choice
import string

from bruhat.gset import Perm, Group, Coset, mulclose
from bruhat.argv import argv

# Todd-Coxeter algorithm: compute finite group from generators and relations.
# Implementation from: 
# https://math.berkeley.edu/~kmill/notes/todd_coxeter.html

# See also:
# https://arxiv.org/pdf/2012.09271.pdf p10.

   
class Word(object):

    def __init__(self, table, path=()):
        self.table = table
        self.path = tuple(path)
        self.idx = table.follow_path(0, path)

    def __str__(self):
        items = []
        for idx in self.path:
            name = string.ascii_lowercase[idx//2]
            if idx%2:
                name = "~"+name # inverse
            items.append(name)
        return "*".join(items)
    __repr__ = __str__

    def __invert__(self):
        path = [idx^1 for idx in reversed(self.path)]
        return Word(self.table, path)

    def __mul__(self, other):
        assert self.table is other.table
        path = self.path + other.path
        return Word(self.table, path)

    def __pow__(g, e):
        if e < 0:
            return (~g).__pow__(-e) # recurse
        op = Word(g.table)
        for i in range(e):
            op = g*op
        return op

    def order(g):
        I = Word(self.table)
        h = g
        count = 1
        while h != I:
            h = g*h
            assert h != g
            count += 1
            assert count < 99999, "???"
        return count

    def __hash__(self):
        table = self.table
        idx = table.follow_path(0, self.path)
        return idx # Warning: table can mutate  XXX

    def __eq__(self, other):
        table = self.table
        assert other.table is table
        left = table.follow_path(0, self.path)
        right = table.follow_path(0, other.path)
        table.dump()
        return left == right



class Table(object):
    """ 
    A Schreier coset graph coming from a group acting on the cosets
    of a subgroup.
    """
    DEBUG = False

    def __init__(self, ngens):
        self.labels = []
        self.neighbors = []
        self.ngens = 2*ngens # include inverses
        start = self.add_label()
        assert start == 0
        row = self.neighbors[start]
        gens = []
        assert ngens < len(string.ascii_lowercase)
        for idx in range(ngens):
            gen = Word(self, (idx,))
            gens.append(gen)
            #self.add_label()
            #self.add_label()
        self.dump()
        self.gens = gens
        self.identity = Word(self)

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
    
#    def unify_recursive(self, c1, c2):
#        labels = self.labels
#        neighbors = self.neighbors
#        ngens = self.ngens
#        c1 = self.get_label(c1)
#        c2 = self.get_label(c2)
#        if c1 == c2:
#            return
#        c1, c2 = min(c1, c2), max(c1, c2)
#        if self.DEBUG:
#            print("unify:", c2, "-->", c1)
#        labels[c2] = c1
#        to_unify = []
#        for d in range(ngens):
#            n1 = neighbors[c1][d]
#            n2 = neighbors[c2][d]
#            if n1 == None:
#                neighbors[c1][d] = n2
#            elif n2 != None:
#                #self.unify(n1, n2) # recurse
#                to_unify.append((n1, n2))
#        for n1, n2 in to_unify:
#            self.unify_recursive(n1, n2) # recurse
    
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

    def add_label(self):
        labels = self.labels
        neighbors = self.neighbors
        ngens = self.ngens
        c = len(labels)
        labels.append(c)
        neighbors.append(ngens*[None])
        return c

    def follow_step(self, c, d):
        labels = self.labels
        neighbors = self.neighbors
        c = self.get_label(c)
        ns = neighbors[c]
        if ns[d] == None:
            ns[d] = self.add_label()
        return self.get_label(ns[d])
    
    def follow_path(self, c, ds):
        c = self.get_label(c)
        for d in reversed(ds):
            c = self.follow_step(c, d)
        return c
    
    def dump(self):
        labels = self.labels
        neighbors = self.neighbors
        print(list(enumerate(labels)))
        for idx, row in enumerate(neighbors):
            jdx = self.get_label(idx)
            if idx == jdx:
                print('\t', idx, row)
        print("-"*70)

    def __len__(self):
        return len([k for k in enumerate(self.labels) if k[0]==k[1]])
    
    def build(self, rels, hgens=[], maxsize=None):
        labels = self.labels
        neighbors = self.neighbors

        for rel in rels+hgens:
            for gen in rel:
                assert 0 <= gen < self.ngens, repr(rel)

        start = 0
    
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
    

    
def test():
    ngens = 4 # include the inverses
    rels = [
        (1, 0), # a^-1a
        (3, 2), # b^-1b
        (0, 0, 0), #a^3
        (2, 2), # b^2
        (0, 2, 0, 2) # abab
    ]
    table = Table(2)
    a, b = table.gens
    I = table.identity
    print(a, b, a**3)

    assert I == a*~a
    assert I == (a*b)*~(a*b)

    #rels = [a**3, 
    #graph.build(rels)
    #print(len(graph))
    #graph.show()
    #assert len(graph) == 6 # S_3








if __name__ == "__main__":

    profile = argv.profile
    name = argv.next() or "test"

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    else:
        fn = eval(name)
        fn()

    print("OK\n")

        

