#!/usr/bin/env python3

"""

"""


import sys, os
from time import time
start_time = time()
import random
from random import randint, choice
from functools import reduce
from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)
from operator import add, mul
from math import prod

import numpy

#scalar = numpy.int64
scalar = numpy.int8 # CAREFUL !!

from bruhat.action import mulclose, mulclose_hom
from bruhat.spec import isprime
from bruhat.argv import argv
from bruhat.solve import parse, enum2, row_reduce, span, shortstr, rank, shortstrx, pseudo_inverse, intersect
from bruhat.solve import zeros2, identity2
from bruhat import solve 
from bruhat.dev import geometry
from bruhat.util import cross, allperms, choose
from bruhat.smap import SMap
from bruhat.qcode import Geometry, get_adj

EPSILON = 1e-8

DEFAULT_P = argv.get("p", 2)



class Matrix(object):
    def __init__(self, A, p=DEFAULT_P, shape=None, name="?"):
        if type(A) == list or type(A) == tuple:
            A = numpy.array(A, dtype=scalar)
        else:
            A = A.astype(scalar) # makes a copy
        if shape is not None:
            A.shape = shape
        self.A = A
        #n = A.shape[0]
        #assert A.shape == (n, n)
        assert int(p) == p
        assert p>=0
        self.p = p
        #self.n = n
        if p>0:
            self.A %= p
        self.key = (self.p, self.A.tobytes())
        self._hash = hash(self.key)
        self.shape = A.shape
        self.name = name

    @classmethod
    def perm(cls, items, p=DEFAULT_P, name="?"):
        n = len(items)
        A = numpy.zeros((n, n), dtype=scalar)
        for i, ii in enumerate(items):
            A[ii, i] = 1
        return Matrix(A, p, name=name)

    @classmethod
    def identity(cls, n, p=DEFAULT_P):
        A = numpy.identity(n, dtype=scalar)
        return Matrix(A, p, name="I")

    def get_bitkey(self):
        assert self.p == 2
        A = numpy.packbits(self.A)
        return A.tobytes()

    def __str__(self):
        return str(self.A)

    def __repr__(self):
        return "Matrix(%s)"%str(self.A)

    def shortstr(self):
        return shortstr(self.A)

    def __hash__(self):
        return self._hash

    def is_zero(self):
        return self.A.sum() == 0

    def __len__(self):
        return len(self.A)

    def __eq__(self, other):
        assert self.p == other.p
        return self.key == other.key

    def __ne__(self, other):
        assert self.p == other.p
        return self.key != other.key

    def __lt__(self, other):
        assert self.p == other.p
        return self.key < other.key

    def __add__(self, other):
        A = self.A + other.A
        return Matrix(A, self.p)

    def __sub__(self, other):
        A = self.A - other.A
        return Matrix(A, self.p)

    def __neg__(self):
        A = -self.A
        return Matrix(A, self.p)

    def __mul__(self, other):
        if isinstance(other, Matrix):
            A = numpy.dot(self.A, other.A)
            return Matrix(A, self.p, name=self.name+other.name)
        else:
            return NotImplemented

    def __rmul__(self, r):
        A = r*self.A
        return Matrix(A, self.p)

    def __getitem__(self, idx):
        A = self.A[idx]
        #print("__getitem__", idx, type(A))
        if type(A) is scalar:
            return A
        return Matrix(A, self.p)

    def transpose(self):
        A = self.A
        return Matrix(A.transpose(), self.p)

    def sum(self):
        return self.A.sum()


def test():
    M = Matrix([[1,1,0],[0,1,0]])
    M1 = Matrix([[1,0,0],[0,1,0]])


def make_genons():

    #solve.int_scalar = numpy.int64
    #print(solve.int_scalar)
    import qupy.ldpc.solve
    solve.int_scalar = qupy.ldpc.solve.int_scalar
    from qupy.ldpc.css import CSSCode

    key = (3, 6)
    idx = argv.get("idx", 12)
    geometry = Geometry(key, idx, True)
    #graph = geometry.build_graph(desc)
    G = geometry.G
    print("|G| = %d, idx = %d" % (len(G), idx))

    faces = geometry.get_cosets([0,1,1])
    edges = geometry.get_cosets([1,0,1])
    verts = geometry.get_cosets([1,1,0])
    print("faces=%d, edges=%d, verts=%d"%(len(faces), len(edges), len(verts)))

    print(G)
    h_face = G.left_action(faces)
    h_edge = G.left_action(edges)
    h_vert = G.left_action(verts)

    lookup = {v:idx for (idx,v) in enumerate(edges)}
    for g in G:
        perm = h_edge[g]
        idxs = [lookup[perm[edge]] for edge in edges]
        print(idxs)
        M = Matrix.perm(idxs)
        print(M.shortstr())
        break

    return

    A = get_adj(faces, verts)
    Az = get_adj(faces, edges)
    Ax = get_adj(verts, edges)
    print(A.shape, Az.shape, Ax.shape)
    print(rank(Az))
    print(rank(Ax))

    code = CSSCode(Hx=A, Hz=A)
    print(code)

    code = CSSCode(Hx=Ax, Hz=Az)
    print(code)
    print(shortstrx(Ax, Az.transpose()))
    #print(shortstr(code.Lx))
    #print()
    #print(shortstr(code.Lz))





if __name__ == "__main__":
    fn = argv.next() or "test"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))





