#!/usr/bin/env python3

from bruhat import element
from bruhat.vec import Space, Hom, Map
from bruhat.action import mulclose, mulclose_hom, Perm, Group
from bruhat.argv import argv
from bruhat.util import factorial


def get_perms(items, parts):
    perms = []
    for part in parts:
        for i in range(len(part)-1):
            perm = dict((item, item) for item in items)
            perm[part[i]] = part[i+1]
            perm[part[i+1]] = part[i]
            perm = Perm(perm, items)
            perms.append(perm)
    if not perms:
        # identity
        perms = [Perm(dict((item, item) for item in items), items)]
    G = Group.generate(perms)
    return G


class Young(object):

    def __init__(self, part, labels=None):
        assert part
        n = sum(part)
        i = part[0]
        for j in part[1:]:
            assert j <= i
            i = j
        if labels is None:
            labels = list(range(n))
        assert len(labels)==n
        rows = []
        idx = 0
        for i in part:
            row = []
            for j in range(i):
                row.append(labels[idx])
                idx += 1
            rows.append(row)
        cols = []
        for i in range(len(rows[0])):
            col = []
            for row in rows:
                if i < len(row):
                    col.append(row[i])
            cols.append(col)
        self.rows = rows
        self.cols = cols
        self.labels = labels
        self.part = tuple(part)
        self.n = n

    def get_rowperms(self):
        return get_perms(self.labels, self.rows)

    def get_colperms(self):
        return get_perms(self.labels, self.cols)

    def __str__(self):
        lines = []
        for row in self.rows:
            pre = ' ' if lines else '['
            line = pre + '[%s]' %(' '.join("%d"%i for i in row))
            lines.append(line)
        return '\n'.join(lines) + ']'


def test():

    t = Young((4, 3, 1, 1))
    assert t.rows == [[0, 1, 2, 3], [4, 5, 6], [7], [8]]
    assert t.cols == [[0, 4, 7, 8], [1, 5], [2, 6], [3]]

    G = t.get_rowperms()
    assert len(G) == factorial(4) * factorial(3)



def main():

    #ring = element.Z
    ring = element.Q

    qubit = Space(2, ring)

    q2 = qubit @ qubit

    hom = Hom(qubit, qubit)

    I = Map.from_array([[1, 0], [0, 1]], hom)
    X = Map.from_array([[0, 1], [1, 0]], hom)
    Z = Map.from_array([[1, 0], [0, -1]], hom)
    H = X+Z
    E = Map.from_array([[0, 1], [0, 0]], hom) # raising 
    F = Map.from_array([[0, 0], [1, 0]], hom) # lowering
    II = I@I
    XI = X@I
    IX = I@X
    XX = X@X
    assert XI*IX == XX
    
    CNOT = Map.from_array(
        [[1, 0, 0, 0],
         [0, 1, 0, 0],
         [0, 0, 0, 1],
         [0, 0, 1, 0]], hom@hom)

    SWAP = Map.from_array(
        [[1, 0, 0, 0],
         [0, 0, 1, 0],
         [0, 1, 0, 0],
         [0, 0, 0, 1]], hom@hom)

    assert SWAP * SWAP == II
    assert CNOT * CNOT == II
    
    G = mulclose([XI, IX, CNOT]) # order 8
    G = mulclose([XI, IX, CNOT, SWAP]) # order 24.. must be S_4
    G = mulclose([CNOT, SWAP])
#    print(len(G))
#    for g in G:
#        print(g)

    A = SWAP @ I
    B = I @ SWAP
    S_3 = list(mulclose([A, B]))
    assert len(S_3) == 6
#    for g in G:
#        print(g)
#    print(g.hom)

    hom = A.hom
    space = hom.src
    N = 2**3
    basis = space.basis()
    orbits = set()
    for v in basis:
        v1 = space.zero_vector()
        for g in S_3:
            u = g*v
            v1 = v1 + u
        orbits.add(v1)
    orbits = list(orbits)
#    for v in orbits:
#        print(v)

    HHH = H@H@H
    v = basis[7]
    #print(v)
    #print(HHH * v)

    # ----------------------------------------------------------

    items = [0, 1, 2]
    g1 = Perm({0:1, 1:0, 2:2}, items)
    g2 = Perm({0:0, 1:2, 2:1}, items)

    s1 = SWAP @ I
    s2 = I @ SWAP

    G = mulclose([g1, g2])
    hom = mulclose_hom([g1, g2], [s1, s2])
    for g in G:
      for h in G:
        assert hom[g*h] == hom[g]*hom[h] # check it's a group hom

    # ----------------------------------------------------------

    items = [0, 1, 2, 3]
    g1 = Perm({0:1, 1:0, 2:2, 3:3}, items)
    g2 = Perm({0:0, 1:2, 2:1, 3:3}, items)
    g3 = Perm({0:0, 1:1, 2:3, 3:2}, items)
    gen1 = [g1, g2, g3]

    s1 = SWAP @ I @ I
    s2 = I @ SWAP @ I
    s3 = I @ I @ SWAP
    gen2 = [s1, s2, s3]

    G = mulclose(gen1)
    hom = mulclose_hom(gen1, gen2)
    for g in G:
      for h in G:
        assert hom[g*h] == hom[g]*hom[h] # check it's a group hom

    projs = []
    for part in [(4,), (3,1), (2,2)]:
        t = Young(part)

        rowG = t.get_rowperms()
        colG = t.get_colperms()
        horiz = None
        for g in rowG:
            P = hom[g]
            horiz = P if horiz is None else (horiz + P)
    
        vert = None
        for g in colG:
            P = hom[g]
            s = g.sign()
            P = s*P
            vert = P if vert is None else (vert + P)
        A = vert * horiz

        assert vert*vert == len(colG) * vert
        assert horiz*horiz == len(rowG) * horiz
        projs.append(A)

    for A in projs:
      for B in projs:
        assert A*B == B*A

    for A in projs:
        print("proj:")
        im = A.image()
        src, tgt = im.hom
        for i in src:
          for j in tgt:
            v = im[j, i]
            if v != ring.zero:
                print("%s*%s"%(v, j), end=" ")
          print()


if __name__ == "__main__":

    name = argv.next()
    fn = eval(name)
    fn()


