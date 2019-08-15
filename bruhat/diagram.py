#!/usr/bin/env python3

"""
Here we turn string _diagrams into polynomials.

See:
"FINITE DIMENSIONAL HILBERT SPACES ARE COMPLETE FOR
DAGGER COMPACT CLOSED CATEGORIES" PETER SELINGER
https://arxiv.org/pdf/1207.6972.pdf
"""

from random import choice
from functools import reduce
from operator import add
#from string import ascii_lowercase
letters = "xyzuvwabcdefghijklmnopqrst"

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

    def interpret(self, verbose=False):
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
        #print("shapes:", shapes)
        shapes = dict(shapes)

        ring = element.Z
        polys = [Poly(letters[i], ring) for i in range(len(self.verts))]

        ops = {} # interpret each vertex name
        for vi, v in enumerate(self.verts):
            #print("interpret", v)
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
            F[tuple(idxs)] = polys[vi]
            #F[tuple(idxs)] = 1 # numpy.einsum can't do object dtype...
            #print(F)
            op = ops.get(v.name)
            if op is None:
                op = F
            else:
                op = op + F
            ops[v.name] = op

        args = []
        for vi, v in enumerate(self.verts):
            op = ops[v.name]
            #op = op.astype(numpy.int)
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
        ops = [ops[v.name] for v in self.verts]
        #value = numpy.einsum(*args)
        value = einsum(ops, args)
        return value
        

def einsum(ops, addrs):
    n = len(ops)
    assert len(addrs) == n

    all_addrs = reduce(add, addrs)
    all_addrs = list(set(all_addrs))
    all_addrs.sort()
    #print("addrs:", addrs)
    #print("all_addrs:", all_addrs)

    # shape lookup for each addr:
    lookup = dict((i, None) for i in all_addrs)
    for i, op in enumerate(ops):
        op = ops[i]
        for j, addr in enumerate(addrs[i]):
            k = lookup[addr]
            if k is None:
                lookup[addr] = op.shape[j]
            else:
                assert k == op.shape[j]

    #print(lookup)
    shape = tuple(lookup[i] for i in all_addrs)
    #print("einsum shape:", shape)

    value = 0
    for idx in numpy.ndindex(shape):
        val = 1
        #print("idx:", idx)
        for i, op in enumerate(ops):
            j = tuple(idx[all_addrs.index(addr)] for addr in addrs[i])
            #print("j:", j)
            _val = op[j]
            if _val == 0:
                break
            val = val * _val
        else:
            value += val

    return value


def test():

    ring = element.Z

    x = Poly("x", ring)
    y = Poly("y", ring)

    #print(x*x*y)


def main():

    A, B, C = "ABC"
    f = Vertex("f", A, A)

    d = Diagram([f])
    assert d.get_inputs() == set([(f, 0)])
    assert d.get_outputs() == set([(f, 0)])

    d.link(f, f, 0, 0)

    assert not d.get_inputs()
    assert not d.get_outputs()
    s = str(d)

    #print(d)
    assert str(d.interpret()) == "x"

    # -------------------------

    f0, f1 = f, f.clone()
    d = Diagram([f0, f1])
    d.link(f0, f1, 0, 0)
    d.link(f1, f0, 0, 0)
    assert d.is_closed()

    #print(d)
    assert str(d.interpret()) == "2*x*y"

    # -------------------------

    f0, f1, f2 = f, f.clone(), f.clone()
    d = Diagram([f0, f1, f2])
    d.link(f0, f1, 0, 0)
    d.link(f1, f2, 0, 0)
    d.link(f2, f0, 0, 0)
    assert d.is_closed()

    #print(d)
    assert str(d.interpret()) == "3*x*y*z"

    # -------------------------

    f = Vertex("f", A, B)
    g = Vertex("g", B, A)

    d = Diagram([f, g])
    d.link(f, g, 0, 0)
    d.link(g, f, 0, 0)
    assert d.is_closed()

    #print(d)
    assert str(d.interpret()) == "x*y"

    # -------------------------

    f = Vertex("f", B, A+A)
    g = Vertex("g", A+A, B)
    h = Vertex("h", B+A, A+B)

    d = Diagram([f, g, h])
    d.link(f, g, 0, 0)
    d.link(f, g, 1, 1)
    d.link(g, h, 0, 0)
    d.link(h, h, 0, 1)
    d.link(h, f, 1, 0)
    assert d.is_closed()

    #print(d)
    assert str(d.interpret())  == "x*y*z"

    # -------------------------

    f = Vertex("f", A, A+A)
    g = Vertex("g", A+A, A)

    for trial in range(10):

        items = [f.clone() for i in range(2)]
        items += [g.clone() for i in range(2)]
        d = Diagram(items)
        while not d.is_closed():
            ilinks = list(d.get_inputs())
            olinks = list(d.get_outputs())
            #print(ilinks, olinks)
            if len(ilinks)*len(olinks) == 0:
                break
            #srcs = [f.tgt[i] for (f, i) in olinks]
            #tgts = [f.src[i] for (f, i) in ilinks]
            f, i = choice(olinks)
            g, j = choice(ilinks)
            d.link(f, g, i, j)
        else:
            #print(d)
            #print("einsum:")
            s = d.interpret()
            print(s)





if __name__ == "__main__":

    main()
    
    



