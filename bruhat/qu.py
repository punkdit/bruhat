#!/usr/bin/env python3

from functools import reduce

from bruhat import element
from bruhat.vec import Space, Hom, Map
from bruhat.action import mulclose, mulclose_hom, Perm, Group
from bruhat.argv import argv
from bruhat.util import factorial, partitions
from bruhat.rep import get_perms, Young


class Specht(object):

    def __init__(self, n, space):
        assert n>=2
        d = len(space)
        ring = space.ring

        items = list(range(n))
        gen1 = []
        for i in range(n-1):
            perm = dict((item, item) for item in items)
            perm[items[i]] = items[i+1]
            perm[items[i+1]] = items[i]
            perm = Perm(perm, items)
            gen1.append(perm)

        perms = mulclose(gen1)
        G = Group(perms, items)

        #I = space.ident

        # tensor power of the space
        tensor = lambda x, y : x@y
        tspace = reduce(tensor, [space]*n)

        # build action of symmetric group on the tensor power
        thom = Hom(tspace, tspace)
        gen2 = []
        for g in gen1:
            items = []
            for idxs in tspace.gen:
                jdxs = tuple(idxs[g[i]] for i in range(len(idxs)))
                items.append(((idxs, jdxs), ring.one))
                #print(idxs, "->", jdxs)
            #print()
            swap = Map(items, thom)
            gen2.append(swap)

        action = mulclose_hom(gen1, gen2)
        for g in G:
          for h in G:
            assert action[g*h] == action[g]*action[h] # check it's a group hom


        # Build the young symmetrizers
        projs = []
        parts = []
        for part in partitions(n):
            if len(part) > d:
                continue
            parts.append(part)
            t = Young(G, part)
    
            rowG = t.get_rowperms()
            colG = t.get_colperms()
            horiz = None
            for g in rowG:
                P = action[g]
                horiz = P if horiz is None else (horiz + P)
        
            vert = None
            for g in colG:
                P = action[g]
                s = g.sign()
                P = s*P
                vert = P if vert is None else (vert + P)
            A = vert * horiz + horiz*vert
            #A = horiz * vert
    
            assert vert*vert == len(colG) * vert
            assert horiz*horiz == len(rowG) * horiz
            #A = A.transpose()
            projs.append(A)

            #print(part)
            #print(t)
            #print(A)

        self.projs = projs
        self.parts = parts



def test():

    n = 9
    G = Group.symmetric(n)
    t = Young(G, (4, 3, 1, 1))
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
    basis = space.get_basis()
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

    specht = Specht(4, qubit)

    for i, A in enumerate(specht.projs):
        print("part:", specht.parts[i])
        print("proj:")
        im = A.image()
        #im = A
        src, tgt = im.hom
        for i in src:
          for j in tgt:
            v = im[j, i]
            if v != ring.zero:
                print("%s*%s"%(v, j), end=" ")
          print()

    return

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
        #A = horiz * vert

        assert vert*vert == len(colG) * vert
        assert horiz*horiz == len(rowG) * horiz
        #print(t)
        #print(A)
        projs.append(A)

    for A in projs:
      for B in projs:
        assert A*B == B*A

#    print()
#    for A in projs:
#        for g in G:
#            P = hom[g]
#            if P*A != A*P:
#                print(P*A)
#                print(A*P)
#                return

    for A in projs:
        print("proj:")
        im = A.image()
        #im = A
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


