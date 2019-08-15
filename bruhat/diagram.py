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

    # -------------------------

    f0, f1 = f, f.clone()
    d = Diagram([f0, f1])
    d.link(f0, f1, 0, 0)
    d.link(f1, f0, 0, 0)
    assert d.is_closed()

    # -------------------------

    A = "A"
    B = "B"
    f = Vertex("f", A, B)
    g = Vertex("f", B, A)

    d = Diagram([f, g])
    d.link(f, g, 0, 0)
    d.link(g, f, 0, 0)
    assert d.is_closed()



if __name__ == "__main__":

    main()
    
    



