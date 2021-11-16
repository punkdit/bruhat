#!/usr/bin/env python3

"""
Construct the five platonic solids, their incidence geometry and
symmetry group.
"""

import sys, os

from bruhat.argv import argv
from bruhat.util import write, uniqtuples
from bruhat import isomorph
from bruhat.action import Group, Perm



EPSILON = 1e-8


class Vector(object):
    def __init__(self, *cs):
        self.cs = cs
        self.n = len(cs)
        self._hash = hash(self.__str__())

    def __str__(self):
        cs = ["%.4f"%c for c in self.cs]
        return "(%s)"%(','.join(cs))

    def __hash__(self):
        return self._hash
        
    def __eq__(self, other):
        return self.__str__() == other.__str__()
        #for i in range(self.n):
        #    if abs(self.cs[i] - other.cs[i]) > EPSILON:
        #        return False
        #return True

    def __getitem__(self, idx):
        return self.cs[idx]
        
    def __ne__(self, other):
        return self.__str__() != other.__str__()
        return not (self==other)

    def __add__(self, other):
        cs = [self[0]+other[0],self[1]+other[1],self[2]+other[2]]
        return Vector(*cs)

    def __sub__(self, other):
        cs = [self[0]-other[0],self[1]-other[1],self[2]-other[2]]
        return Vector(*cs)

    def reflect(self, i):
        cs = list(self.cs)
        cs[i] = -cs[i]
        return Vector(*cs)

    def cross(self, other):
        cs = [
            self[1]*other[2] - self[2]*other[1],
            self[2]*other[0] - self[0]*other[2],
            self[0]*other[1] - self[1]*other[0]]
        return Vector(*cs)

    def dot(self, other):
        r = self[0]*other[0]+self[1]*other[1]+self[2]*other[2]
        return r



class Item(object):

    items = None
    def __init__(self, verts):
        self.verts = list(verts)
        self.s_verts = set(verts)
        self._hash = hash(tuple(self.verts))

        if self.__class__.items is not None:
            self.__class__.items.append(self)

    @classmethod
    def push(cls):
        cls.items = []

    @classmethod
    def pop(cls):
        items = cls.items
        cls.items = None
        return items

    def __eq__(self, other):
        return self.verts == other.verts

    def __ne__(self, other):
        return self.verts != other.verts

    def __hash__(self):
        #return hash(tuple(self.verts))
        return self._hash

    def __getitem__(self, idx):
        return self.verts[idx]

    def reversed(self):
        verts = list(reversed(self.verts))
        return Face(verts)

    def reflect_x(self):
        verts = [v.reflect(0) for v in self.verts]
        return self.__class__(verts)

    def reflect_y(self):
        verts = [v.reflect(1) for v in self.verts]
        return self.__class__(verts)

    def reflect_z(self):
        verts = [v.reflect(2) for v in self.verts]
        return self.__class__(verts)

    def intersection(self, other):
        return (self.s_verts).intersection((other.s_verts))

    def contains(self, other):
        return (self.s_verts).issuperset((other.s_verts))


class Edge(Item):
    pass

class Face(Item):
    pass

class Vertex(Item):
    pass


def fix_orientation(face):
    v0, v1, v2 = face[:3]
    a = v1-v0
    b = v2-v0
    n = a.cross(b)
    items = [v.dot(n) for v in face]
    reverse = items[0] < 0
    if reverse:
        face = face.reversed()
    return face


def make_cube():
    Face.push()

    vs = []
    
    for i in [-1, 1]:
        vs = [
            Vector(i, -1, -1),
            Vector(i, -1, +1),
            Vector(i, +1, +1),
            Vector(i, +1, -1)]
        Face(vs)
    
    for j in [-1, 1]:
        vs = [
            Vector(-1, j, -1),
            Vector(-1, j, +1),
            Vector(+1, j, +1),
            Vector(+1, j, -1)]
        Face(vs)
    
    for k in [-1, 1]:
        vs = [
            Vector(-1, -1, k),
            Vector(-1, +1, k),
            Vector(+1, +1, k),
            Vector(+1, -1, k)]
        Face(vs)

    return [fix_orientation(face) for face in Face.pop()]


def make_dodecahedron():
    Face.push()
    phi = 0.5*(1 + 5**0.5) # golden ratio

    f0 = Face([
        Vector(+1/phi, +phi, 0),
        Vector(-1/phi, +phi, 0),
        Vector(-1, 1, 1),
        Vector(0, 1/phi, phi),
        Vector(+1, 1, 1)])
    f1 = f0.reflect_z()
    f2 = f0.reflect_y()
    f3 = f1.reflect_y()

    f0 = Face([
        Vector(+1/phi, +phi, 0),
        Vector(+1, 1, 1),
        Vector(phi, 0, 1/phi),
        Vector(phi, 0, -1/phi),
        Vector(1, 1, -1)])
    f1 = f0.reflect_x()
    f2 = f0.reflect_y()
    f3 = f1.reflect_y()

    f0 = Face([
        Vector(+1, 1, 1),
        Vector(phi, 0, 1/phi),
        Vector(+1, -1, 1),
        Vector(0, -1/phi, phi),
        Vector(0, +1/phi, phi)])
    f1 = f0.reflect_x()
    f2 = f0.reflect_z()
    f3 = f1.reflect_z()

    return [fix_orientation(face) for face in Face.pop()]


def make_icosahedron():
    Face.push()
    phi = 0.5*(1 + 5**0.5) # golden ratio

    f0 = Face([
        Vector(0,    1, phi),
        Vector(0,   -1, phi),
        Vector(phi,  0, 1)])
    f1 = f0.reflect_x()
    f2 = f0.reflect_z()
    f3 = f1.reflect_z()

    f0 = Face([
        Vector( 1,  phi, 0),
        Vector(-1,  phi, 0),
        Vector( 0,  1,   phi)])
    f1 = f0.reflect_z()
    f2 = f0.reflect_y()
    f3 = f1.reflect_y()

    f0 = Face([
        Vector(phi, 0, 1),
        Vector(phi, 0, -1),
        Vector(1,   phi, 0)])

    f1 = f0.reflect_y()
    f2 = f0.reflect_x()
    f3 = f1.reflect_x()

    f0 = Face([
        Vector(0,    1, phi),
        Vector(1,   phi,  0),
        Vector(phi,  0, 1)])
    f1 = f0.reflect_x()
    f2 = f0.reflect_y()
    f3 = f0.reflect_z()

    f0 = Face([
        Vector(0,   -1,-phi),
        Vector(-1,  -phi,  0),
        Vector(-phi,  0,-1)])
    f1 = f0.reflect_x()
    f2 = f0.reflect_y()
    f3 = f0.reflect_z()

    return [fix_orientation(face) for face in Face.pop()]



def make_octahedron():
    Face.push()

    f0 = Face([
        Vector(1, 0, 0),
        Vector(0, 1, 0),
        Vector(0, 0, 1)])

    f1 = f0.reflect_x()
    f2 = f1.reflect_y()
    f3 = f2.reflect_x()

    f0 = f0.reflect_z()
    f1 = f0.reflect_x()
    f2 = f1.reflect_y()
    f3 = f2.reflect_x()

    return [fix_orientation(face) for face in Face.pop()]


def make_simplex():
    Face.push()

    v0 = Vector(1, 1, 1)
    v1 = Vector(1, -1, -1)
    v2 = Vector(-1, -1, 1)
    v3 = Vector(-1, 1, -1)

    Face([v0, v1, v2])
    Face([v0, v1, v3])
    Face([v0, v2, v3])
    Face([v1, v2, v3])

    return [fix_orientation(face) for face in Face.pop()]
make_tetrahedron = make_simplex


def make_geometry(name):
    #print "make_geometry", name
    faces = eval("make_%s"%name)()
    #print "faces:", len(faces)

    edges = []
    n = len(faces)
    for i in range(n):
      for j in range(i+1, n):
        #items = set(faces[i].verts).intersection(faces[j].verts)
        items = faces[i].intersection(faces[j])
        assert len(items)<=2
        if len(items) == 2:
            edge = Edge(items)
            edges.append(edge)
    #print "edges:", len(edges)

    verts = set()
    n = len(edges)
    for i in range(n):
      for j in range(i+1, n):
        #items = set(edges[i].verts).intersection(edges[j].verts)
        items = edges[i].intersection(edges[j])
        assert len(items)<=1
        if len(items) == 1:
            vert = Vertex(items)
            verts.add(vert)
    #print "verts:", len(verts)

    items = faces + edges + list(verts)
    tpmap = {}
    for face in faces:
        tpmap[face] = "face"
    for edge in edges:
        tpmap[edge] = "edge"
    for vert in verts:
        tpmap[vert] = "vert"
    incidence = []
    for item in items:
      for jtem in items:
        if tpmap[item] != tpmap[jtem]:
            if item.contains(jtem):
                incidence.append((item, jtem))
    for item in items:
        incidence.append((item, item))
    from bruhat.geometry import Geometry
    geometry = Geometry(incidence, tpmap)

    return geometry
    


def test():

    for name in "simplex cube octahedron dodecahedron icosahedron".split():
        make_geometry(name)
        print()
    print("OK")



def orbiplex(G, k=4):

    #import numpy
    #from gelim import zeros, dotx, rank, nullity

    print("orbiplex: |G|=%d" % len(G))
    items = G.items

    k = argv.get("k", k)
    for n in range(k):

        write("|C_%d| ="%n)
        #if n > len(items):
        #    break

        tpls = list(uniqtuples(items, n))

        perms = []
        for perm in G:
            _perm = {}
            for key in tpls:
                value = tuple(perm[i] for i in key)
                _perm[key] = value
            _perm = Perm(_perm, tpls)
            perms.append(_perm)

        G2 = Group(perms, tpls)

        orbits = list(G2.orbits())
        d = len(orbits)
        write("%d,"%d)


def main():
    name = argv.next() or "icosahedron"
    geometry = make_geometry(name)

    graph = geometry.get_bag()
    graph1 = geometry.get_bag()
    #print graph
    n = len(graph)
    items = [i for i in range(n) if graph[i].desc=="vert"]
    #print items

    perms = []
    for f in isomorph.search(graph, graph1):
        #print f
        f = dict((i, f[i]) for i in items)
        perm = Perm(f, items)
        perms.append(perm)
        write('.')
    print("symmetry:", len(perms))
    assert len(set(perms)) == len(perms)

    G = Group(perms, items)

    orbiplex(G)



if __name__ == "__main__":

    if argv.profile:
       import cProfile as profile
       profile.run("test()")

    elif argv.test:
        test()

    else:
       main()


    
