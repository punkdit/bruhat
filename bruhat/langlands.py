#!/usr/bin/env python3

"""
Try to construct the Langlands correspondence
between conjugacy classes of elements of G = GL(n, q)
with complex irreducible representations of G.
"""

from bruhat.element import FiniteField, PolynomialRing, GaloisField, Linear
from bruhat.action import mulclose
from bruhat.util import cross
from bruhat.argv import argv


def GL2(field):

    n = 2
    GL = Linear(n, field)

    els = field.elements
    zero = field.zero

    pairs = []
    for a in els:
      for b in els:
        if a==zero and b==zero:
            continue
        pairs.append((a, b))
    G = []
    for (a, b) in pairs:
      for (c, d) in pairs:
        if a*d != b*c:
            G.append(GL.get([[a, b], [c, d]]))

    print("GL(%d, %d), order = %d" % (n, len(field), len(G)))

    I = GL.get([[1, 0], [0, 1]])
    remain = set(G)
    remain.remove(I)
    inv = {I:I}
    while remain:
        g = iter(remain).__next__()
        for h in remain:
            if g*h == I:
                break
        else:
            assert 0
        assert h*g == I
        inv[g] = h
        remain.remove(g)
        if g!=h:
            inv[h] = g
            remain.remove(h)

    print("conjugacy classes =", end=" ", flush=True)
    remain = set(G)
    cgys = []
    while remain:
        g = iter(remain).__next__()
        remain.remove(g)
        cls = set([g])
        for h in G:
            g1 = inv[h] * g * h
            if g1 in remain:
                cls.add(g1)
                remain.remove(g1)
        cgys.append(list(cls))
    print(len(cgys))


    if 0:
        for cls in cgys:
            for g in cls:
                print(str(g).replace("\n", ""))
            print()
    


def main():

    p = argv.get("p", 3)
    deg = argv.get("deg", 1)

    field = FiniteField(p)

    if deg == 1:
        G = GL2(field)
        return

    base, field = field, field.extend(deg)

    assert len(field.elements) == p**deg
    print("GF(%d)"%len(field.elements))

    print(field.mod)
    print()

    ring = PolynomialRing(field)
    poly = ring.evaluate(field.mod)
    assert str(poly) == str(field.mod)
    assert poly(field.x) == 0
    subs = []
    for a in field.elements:
        if poly(a) == field.zero:
            subs.append(a)
            print(a)

#    subs = [ring.evaluate(a) for a in subs]
#    for a in subs:
#      for b in subs:
#        print(ring.evaluate(a, b))

    return

    lin = Linear(deg, base)

    zero = field.zero
    one = field.one
    x = field.x
    a = one
    basis = []
    for i in range(deg):
        basis.append(a)
        a = x*a

    for a in field.elements:
        if a == zero:
            continue
        op = [[0 for j in range(deg)] for i in range(deg)]
        for j, b in enumerate(basis):
            c = a*b
            poly = c.value # unwrap
            for i in range(deg):
                op[i][j] = poly[i]
        op = lin.get(op)
        print(str(op).replace("\n", ""), a)
        print()

    if argv.check:
        zero = field.zero
        div = {}
        for a in field.elements:
            for b in field.elements:
                c = a*b
                div[c, a] = b
                div[c, b] = a
    
        for a in field.elements:
            for b in field.elements:
                if b != zero:
                    assert (a, b) in div




if __name__ == "__main__":

    main()


