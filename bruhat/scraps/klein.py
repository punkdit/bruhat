#!/usr/bin/env python3

"""
try to build the klein quartic from heptagons in the
hyperbolic plane. fail.
"""

from bruhat.argv import argv



class Item(object):
    "Elements of a graded poset, or cellular complex."
    def __init__(self, dim, desc=None):
        self.dim = dim
        self.desc = desc
        self.bdy = []
        self.cobdy = []

    def add_bdy(parent, child):
        assert child.dim == parent.dim-1
        assert child not in parent.bdy
        parent.bdy.append(child)
        assert parent not in child.cobdy
        child.cobdy.append(parent)

    def remove_bdy(parent, child):
        assert child in parent.bdy
        for item in list(child.bdy):
            child.bdy.remove_bdy(item)
        parent.bdy.remove(child)

    def intersect(self, other):
        if other.dim > self.dim:
            self, other = other, self # swap
        result = None
        if other.dim == self.dim:
            if other is self:
                result = self
            for item in self.bdy:
                if item in other.bdy:
                    assert result is None
                    result = item
        elif other.dim == self.dim-1:
            if other in self.bdy:
                result = other
        elif other.dim < self.dim-1:
            for child in self.bdy:
                item = child.intersect(other)
                if item:
                    result = item
        else:
            assert 0
        return result


def join(parent0, child0, parent1, child1):
    assert child0 in parent0.bdy
    assert child1 in parent1.bdy
    assert parent0.dim == parent1.dim



def polyhedron(n):
    face = Item(2)
    verts = [Item(0) for i in range(n)]
    edges = [Item(1) for i in range(n)]
    for i, edge in enumerate(edges):
        v0 = verts[i]
        v1 = verts[(i+1)%n]
        edge.add_bdy(v0)
        edge.add_bdy(v1)
        face.add_bdy(edge)
    return face


class Figure(object):
    def __init__(self, items=[]):
        self.items = list(items)



def main():

    p = argv.get(p, "5")

    face = Item(2, 


if __name__ == "__main__":

    main()



