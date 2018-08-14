#!/usr/bin/env python3

"""
Calculate the genus of hecke modular curves.
"""

import sys, os

from argv import argv


class Field(object):

    def __init__(self, p):
        self.p = p
    
    def neg(self, i):
        p = self.p
        return (p-i)%p
    
    def recip(self, i):
        p = self.p
        assert 0<i<p
        for j in range(1, p):
            if (j*i)%p == 1:
                return j
        assert 0
    
    def negr(self, i):
        return self.neg(self.recip(i))


#def Vertex(object):
#    def __init__(self


class Edge(object):
    def __init__(self, name, v0, v1):
        self.name = name
        self.v0 = v0
        self.v1 = v1

    def __str__(self):
        return "(%s)[%s--%s]"%(self.name, self.v0, self.v1)
    __repr__ = __str__

    #def __eq__(self, other):


class Curve(object):
    def __init__(self, edges):
        self.edges = list(edges)
        self.lookup = {}
        for e in edges:
            self.lookup[e.name] = e

    def get(self, name):
        return self.lookup[name]

    def remove(self, e):
        self.edges.remove(e)
        del self.lookup[e.name]

    def __str__(self):
        #return "Curve(%d, %d)"%(len(self.edges), len(self.vertices))
        return "Curve(%s, %s)"%(self.edges, self.vertices)

    def send(self, vmap):
        for e in self.edges:
            e.v0 = vmap.get(e.v0, e.v0)
            e.v1 = vmap.get(e.v1, e.v1)

    def collapse(self, i):
        e = self.lookup[i]
        #print("collapse", e)
        self.remove(e)
        v0, v1 = e.v0, e.v1
        self.send({v0:v1})

    def identify(self, i, j):
        if i==j:
            return self.collapse(i) # <-- return
        e0 = self.lookup[i]
        e1 = self.lookup[j]
        #print("identify", e0, e1)
        edges = self.edges
        assert e0 in edges
        assert e1 in edges
        if e0.v0 == e0.v1:
            assert 0 # TODO
        if e1.v0 == e1.v1:
            assert 0 # TODO
        if e0.v0 == e1.v1:
            e0, e1 = e1, e0
        if e0.v1 == e1.v0:
            self.remove(e0)
            #self.remove(e1)
            #assert e0.v0 != e1.v1, (e0, e1)
            self.send({e0.v0:e1.v1})
            return  # <---- return
            
        self.remove(e1)
        vmap = {}
        assert e0.v0 != e0.v1
        assert e1.v0 != e1.v1
        vmap[e1.v1] = e0.v0
        vmap[e1.v0] = e0.v1
        # rewrite vertices
        self.send(vmap)

    @property
    def vertices(self):
        vs = set()
        for e in self.edges:
            vs.add(e.v0)
            vs.add(e.v1)
        v = list(vs)
        v.sort()
        return v


def test():
    f = Field(5)
    assert f.negr(1)==4
    assert f.negr(2)==2
    assert f.negr(3)==3
    assert f.negr(4)==1


def get_genus(p, verbose=False):
    assert p>2

    f = Field(p)
    
    pairs = []
    for i in range(1, p):
        pair = [i, f.negr(i)]
        pair.sort()
        if pair in pairs:
            continue
        pairs.append(pair)
        #print("%d <-> %d," % tuple(pair), end=' ')

    edges = []
    verts = range(p)
    for i in range(p):
        v0 = verts[(i-1)%p]
        v1 = verts[i]
        edge = Edge(i, v0, v1)
        edges.append(edge)

    curve = Curve(edges)

    if verbose:
        print(curve)

    curve.collapse(0)

    if verbose:
        print(curve)

    for (i, j) in pairs:
        if verbose:
            print(i, j)
        curve.identify(i, j)
        if verbose:
            print(curve)

    # Euler characteristic
    chi = len(curve.vertices) - len(curve.edges) + 1
    #print("chi:", chi)
    assert chi%2==0

    genus = (2-chi)//2
    #print("genus:", genus)

    return genus


from theta import all_primes

def main():

    p = argv.get("p")

    if p:
        genus = get_genus(p)
    
        print("genus:", genus)
    
    else:
        for p in all_primes(104):
            if p<5:
                continue
            genus = get_genus(p)
            print("p=%d, genus=%d"%(p, genus))
            


if __name__ == "__main__":

    main()


