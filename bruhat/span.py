#!/usr/bin/env python3

"""

attempted implementation of:
http://math.ucr.edu/home/baez/week188.html

Schubert cell decomposition of Grassmanians...
"""

import numpy


from bruhat.argv import argv
from bruhat.util import cross
from bruhat.action import mulclose

scalar = numpy.int64



class Matrix(object):
    def __init__(self, v, p):
        self.v = v
        self.p = int(p)
        self.shape = v.shape
        _str = str(v)
        _str = _str.replace("\n", ',')
        self._str = _str

    def __eq__(self, other):
        return self._str == other._str

    def __ne__(self, other):
        return self._str != other._str

    def __hash__(self):
        return hash(self._str)

    def __mul__(self, other):
        if not isinstance(other, Matrix):
            return NotImplemented
        v = numpy.dot(self.v, other.v)
        v %= self.p
        return Matrix(v, self.p)

    def __str__(self):
        return self._str
    __repr__ = __str__

    def transpose(self):
        v = self.v.transpose()
        return Matrix(v, self.p)


def SL(n, p): 
    "special linear group"
    assert int(n)==n
    assert int(p)==p
    assert n>0 

    I = numpy.identity(n, scalar)
    gen = []
    for i in range(n):
        for j in range(n):
            if i==j:
                continue
            A = I.copy()
            A[i, j] = 1 
            gen.append(Matrix(A, p))
    G = mulclose(gen)
    return G



def parse(spec):
    rows = spec.split('\n')
    rows = [row.strip() for row in rows]
    rows = [row for row in rows if row]
    rows = [row.replace(" ", "") for row in rows]
    rows = [','.join(list(row)) for row in rows]
    rows = ["[%s]"%row for row in rows]
    m = ','.join(rows)
    m = "[%s]"%m
    return m


def subs(spec, p):
    spec = parse(spec)
    idxs = []
    for i, c in enumerate(spec):
        if c=="*":
            idxs.append(i)
    if not idxs:
        yield spec
        return
    spec = list(spec)
    items = list(range(p))
    for sub in cross([items]*len(idxs)):
        for i, val in enumerate(sub):
            spec[idxs[i]] = str(val)
        yield ''.join(spec)


def generate(spec, p, transpose=False):
    ms = set()
    for item in subs(spec, p):
        value = eval(item)
        value = numpy.array(value)
        if transpose:
            value = value.transpose()
        value = Matrix(value, p)
        ms.add(value)
    assert len(ms) == (p**spec.count("*"))
    return Span(ms)


class Span(object):
    "a subspace spanned by a set of Matrix objects"
    def __init__(self, items):
        self.items = set(items) # a set of Matrix objects

    def __str__(self):
        items = [str(v) for v in self.items]
        items.sort()
        return "Span([%s])"%(', '.join(items))

    def incident(self, other):
        a, b = self.items, other.items
        return a.issubset(b) or b.issubset(a)

    def __eq__(self, other):
        a, b = self.items, other.items
        return a==b

    def __ne__(self, other):
        a, b = self.items, other.items
        return a!=b

    def __len__(self):
        return len(self.items)

    def __iter__(self):
        return iter(self.items)

    def __mul__(self, other):
        items = set()
        for a in self.items:
          for b in other.items:
            items.add(a*b)
        return Span(items)

    def __rmul__(self, other):
        if not isinstance(other, Matrix):
            return NotImplemented
        # other * self
        items = set()
        for b in self.items:
            items.add(other*b)
        return Span(items)

    def get_spec(self):
        items = [g.v for g in self.items]
        g = items[0]
        for h in items[1:]:
            g = g + h
        a = numpy.ones(g.shape, scalar)
        a = numpy.minimum(a, g)
        s = str(a)
        s = s.replace("1", "*")
        return s



def test():

    p = 2

    s = parse(
    """
    * * *
    """)
    assert s == "[[*,*,*]]"

    spec = """
    1 * 1
    1 * *
    """

    assert len(generate(spec, p)) == 8

test()


def main():

    n = 4
    p = 2

    point = generate("*000", p, transpose=True)
    line  = generate("**00", p, transpose=True)
    plane = generate("***0", p, transpose=True)

    assert point.incident(line)
    assert point.incident(plane)
    assert line.incident(plane)

    M = generate(
    """
        ****
        ****
        ****
        ****
    """, p)

    #G = SL(n, p)
    #print("|G| =", len(G))

    print("|M| =", len(M))

    G = [g for g in M if abs(numpy.linalg.det(g.v)%p - 1.) < 1e-10]
    print("|SL| =", len(G))

    items = []
    for g in G:
        if g*point == point and g*line == line and g*plane == plane:
            items.append(g)
    B = Span(items)

    print("|B| =", len(B))

    print(B.get_spec())

    P = []
    for g in G:
        if g*point == point:
            items.append(g)
    P = Span(items)
    print("|P| =", len(P))
    print(P.get_spec())

    


if __name__ == "__main__":

    main()



