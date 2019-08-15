#!/usr/bin/env python3

"""
Here we turn string _diagrams into polynomials.

See:
"FINITE DIMENSIONAL HILBERT SPACES ARE COMPLETE FOR
DAGGER COMPACT CLOSED CATEGORIES" PETER SELINGER
https://arxiv.org/pdf/1207.6972.pdf
"""

from functools import reduce
from operator import add
from string import ascii_lowercase as letters

import numpy

from bruhat import element
#from bruhat.vec import Space, Hom, Map
from bruhat.poly import Poly


class Vertex(object):
    r""" The name and type of a _vertex.
        f : B --> A\otimes A
    """
    def __init__(self, name, src, tgt):
        self.name = name
        self.src = src
        self.tgt = tgt

    def __str__(self):
        return "(%s:%s --> %s)"%(self.name, self.src, self.tgt)

    def __repr__(self):
        return "Vertex(%s, %s, %s)"%(self.name, self.src, self.tgt)

    def clone(self):
        return Vertex(self.name, self.src, self.tgt)

    def get_inputs(self):
        return [(self, i) for i in range(len(self.src))]

    def get_outputs(self):
        return [(self, i) for i in range(len(self.tgt))]

def is_uniq(items):
    return len(set(items)) == len(items)


class Diagram(object):
    def __init__(self, verts, links=[]):
        assert is_uniq(verts)
        assert is_uniq(links)
        sig = {}
        for v in verts:
            assert isinstance(v, Vertex)
            tp = sig.get(v.name)
            assert tp is None or tp == (v.src, v.tgt), \
                "%s: %s already has type %s --> %s" % (v, v.name, tp[0], tp[1])
            sig[v.name] = v.src, v.tgt
        self.sig = sig
        self.verts = list(verts)
        self.lookup = dict((v, i) for (i, v) in enumerate(verts))
        self.links = list(links)

    def __str__(self):
        return "Diagram(%s, %s)"%(self.verts, self.links)

    def link(self, f, g, i, j):
        "link i-th output of f to j-th input of g"
        link = (f, g, i, j)
        assert link not in self.links
        assert f in self.lookup
        assert g in self.lookup
        assert f.tgt[i] == g.src[j]
        self.links.append(link)

    def get_inputs(self):
        items = reduce(add, (v.get_inputs() for v in self.verts))
        assert is_uniq(items)
        items = set(items)
        for (f, g, i, j) in self.links:
            assert (g, j) in items
            items.remove((g, j))
        return items

    def get_outputs(self):
        items = reduce(add, (v.get_outputs() for v in self.verts))
        assert is_uniq(items)
        items = set(items)
        for (f, g, i, j) in self.links:
            assert (f, i) in items
            items.remove((f, i))
        return items

    def is_closed(self):
        return not self.get_inputs() and not self.get_outputs()

    def interpret(self):
        links = self.links
        wires = {}
        for link in self.links:
            f, g, i, j = link
            label = f.tgt[i]
            assert label == g.src[j]
            occurs = wires.setdefault(label, [])
            occurs.append(link)

        # dimension of each vector space
        shapes = list(wires.items())
        shapes.sort()
        shapes = [(k, len(v)) for (k, v) in shapes]
        print("shapes:", shapes)
        shapes = dict(shapes)

        ring = element.Z
        polys = [Poly(letters[i], ring) for i in range(len(self.verts))]

        ops = {} # interpret each vertex name
        for vi, v in enumerate(self.verts):
            print("interpret", v)
            # shape = (out[0], out[1], ..., in[0], in[1], ...)
            shape = []
            for label in v.tgt + v.src: # out + in
                shape.append(shapes[label])
            F = numpy.zeros(shape=shape, dtype=object)
            F[:] = 0
            idxs = []
            for i, label in enumerate(v.tgt):
                occurs = wires[label] # all the links where this label occurs 
                idx = None
                for _idx, link in enumerate(occurs):
                    if link[0] == v and link[2] == i:
                        assert idx is None
                        idx = _idx
                assert idx is not None, "missing link"
                idxs.append(idx)
            for i, label in enumerate(v.src):
                occurs = wires[label] # all the links where this label occurs 
                idx = None
                for _idx, link in enumerate(occurs):
                    if link[1] == v and link[3] == i:
                        assert idx is None
                        idx = _idx
                assert idx is not None, "missing link"
                idxs.append(idx)
            #F[tuple(idxs)] = polys[vi]
            F[tuple(idxs)] = 1
            print(F)
            op = ops.get(v.name)
            if op is None:
                op = F
            else:
                op = op + F
            ops[v.name] = op

        #contract = [] # for einsum
        args = []
        for vi, v in enumerate(self.verts):
            op = ops[v.name].astype(numpy.int)
            args.append(op)
            idxs = []
            for i, label in enumerate(v.tgt):
              for idx, link in enumerate(self.links):
                if link[0] == v and link[2] == i:
                    idxs.append(idx)
            for i, label in enumerate(v.src):
              for idx, link in enumerate(self.links):
                if link[1] == v and link[3] == i:
                    idxs.append(idx)
            assert len(idxs) == len(op.shape)
            args.append(idxs)
        print(args)
        value = numpy.einsum(*args)
        return value
        


def test():

    ring = element.Z

    x = Poly("x", ring)
    y = Poly("y", ring)

    #print(x*x*y)


def main():

    A = "A"
    B = "B"
    f = Vertex("f", A, A)

    d = Diagram([f])
    assert d.get_inputs() == set([(f, 0)])
    assert d.get_outputs() == set([(f, 0)])

    d.link(f, f, 0, 0)

    assert not d.get_inputs()
    assert not d.get_outputs()
    s = str(d)

    print(d)
    print(d.interpret())

    # -------------------------

    f0, f1 = f, f.clone()
    d = Diagram([f0, f1])
    d.link(f0, f1, 0, 0)
    d.link(f1, f0, 0, 0)
    assert d.is_closed()

    print(d)
    print(d.interpret())

    # -------------------------

    A = "A"
    B = "B"
    f = Vertex("f", A, B)
    g = Vertex("g", B, A)

    d = Diagram([f, g])
    d.link(f, g, 0, 0)
    d.link(g, f, 0, 0)
    assert d.is_closed()

    print(d)
    print(d.interpret())



if __name__ == "__main__":

    main()
    
    



