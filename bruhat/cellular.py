#!/usr/bin/env python

"""
cellular complexes

"""


from time import time
from random import choice, randint

import numpy

from bruhat.argv import argv


class Cell(object):
    def __init__(self, dim=0, children=[]):
        self.dim = dim
        if type(children) in [list, tuple, set]:
            children = {child:1 for child in children}
        self.children = dict(children)

    def __str__(self):
        return "Cell(%d, %d)"%(self.dim, len(self.children))
    __repr__ = __str__


class Complex(object):
    def __init__(self, grades=None):
        if grades is None:
            grades = [[], [], []]
        self.grades = [list(g) for g in grades]
        self.lookup = {}

    def __str__(self):
        return "Complex(%s)"%(self.grades,)

    def bdy(self, dim):
        tgt, src = self.grades[dim:dim+2]
        m, n = len(tgt), len(src)
        A = numpy.zeros((m, n), dtype=int)
        #print("bdy", dim, src)
        for j,cell in enumerate(src):
            #print(j, cell)
            for dell,r in cell.children.items():
                i = tgt.index(dell)
                #print('\t', (i, j), "<--", r)
                A[i, j] = r
        return A

    def cell(self, dim=0, children=[]):
        cell = Cell(dim, children)
        self.grades[dim].append(cell)
        return cell

    def face(self, children):
        cell = self.cell(2, children)
        return cell

    def edge(self, v0, v1):
        cell = self.cell(1, {v0:-1, v1:+1})
        return cell

    def vertex(self):
        cell = self.cell(0)
        return cell




def make_ball():
    cx = Complex()
    vertex, edge, face = cx.vertex, cx.edge, cx.face
    v0 = vertex()
    v1 = vertex()
    e0 = edge(v0, v1)
    e1 = edge(v1, v0)
    f0 = face({e0, e1})
    f1 = face({e0:-1, e1:-1})
    return cx


def main():

    cx = make_ball()
    print(cx)
    H0 = cx.bdy(0)
    print(H0)

    H1 = cx.bdy(1)
    print(H1)

    HH = numpy.dot(H0, H1)
    print(HH)


if __name__ == "__main__":

    start_time = time()


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
        main()

    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")

