#!/usr/bin/env python3

"""
binary linear codes from quadratic residues.
eg. p=23 gives the Golay code.

See:
Sphere Packings, Lattices and Groups
J. H. ConwayN. J. A. Sloane
p84.

The Theory of Error-Correcting Codes
F.J. MacWilliams
N.J.A. Sloane
p45.
"""


import sys, os

import numpy

from argv import argv
from solve import zeros2, shortstr, find_kernel, dot2, array2, span, eq2, rank
from solve import linear_independent
from action import Perm, Group


def even_code():

    p = argv.get("p", 7)

    #assert (p%8) in [1, 7]
    # equivalent to "2 (binary!) is a quadratic _residue mod p"
    # eg. 7 17 23 31 41 47 71 73 79 89 97

    def neginv(i):
        if i==0:
            return p
        if i==p:
            return 0
        for j in range(1, p):
            if (i*j)%p == 1:
                break
        else:
            assert 0
        return (-j)%p

    nonresidues = set(range(1, p))
    residues = set()
    for i in range(1, p):
        j = (i*i)%p
        if j not in residues:
            residues.add(j)
            nonresidues.remove(j)
    print("residues:", residues)
    print("non-residues:", nonresidues)

    # extended binary quadratic _residue code
    N = p+1
    G = zeros2(N, N)

    G[p, :] = 1
    for u in range(p):
        G[u, p] = 0 if (p-1)%8 == 0 else 1
        for v in range(p):
            if u==v:
                i = 0
            elif (v-u)%p in residues:
                i = 1
            else:
                i = 0
            G[u, v] = i

    G = linear_independent(G)
    print("G =")
    print(shortstr(G))

    from qupy.ldpc.css import CSSCode
    from qupy.ldpc.gallagher import classical_distance
    G = G.astype(numpy.int32)
    code = CSSCode(Hx=G, Hz=G)
    print(code)

    from bruhat.codes import strong_morthogonal
    for genus in range(1, 4):
        print("genus:", genus, "strong_morthogonal:", strong_morthogonal(G, genus))

    def double(G):
        M, N = G.shape
        DG = zeros2(M+1, 2*N)
        DG[1:, 0:N] = G
        DG[1:, N:2*N] = G
        DG[0, 0:N] = 1
        DG = DG.astype(numpy.int32)
        return DG

    DG = G
    DG = DG.astype(numpy.int32)
    print("distance:", classical_distance(DG))
    for _ in range(2):

        DG = double(DG)
        DG = linear_independent(DG)
        print(shortstr(DG))
        
        for genus in range(1, 5):
            print("genus:", genus, "strong_morthogonal:", strong_morthogonal(DG, genus))
    
        code = CSSCode(Hx=DG, Hz=DG)
        print(code)
        #print("distance:", classical_distance(DG))
    


def main():

    p = argv.get("p", 7)

    #assert (p%8) in [1, 7]
    # equivalent to "2 (binary!) is a quadratic _residue mod p"
    # eg. 7 17 23 31 41 47 71 73 79 89 97

    def neginv(i):
        if i==0:
            return p
        if i==p:
            return 0
        for j in range(1, p):
            if (i*j)%p == 1:
                break
        else:
            assert 0
        return (-j)%p

    nonresidues = set(range(1, p))
    residues = set()
    for i in range(1, p):
        j = (i*i)%p
        if j not in residues:
            residues.add(j)
            nonresidues.remove(j)
    print("residues:", residues)
    print("non-residues:", nonresidues)

    # extended binary quadratic _residue code
    N = p+1
    G = zeros2(N, N)

    G[p, :] = 1
    for u in range(p):
        G[u, p] = 0 if (p-1)%8 == 0 else 1
        for v in range(p):
            if u==v:
                i = 0
            elif (v-u)%p in residues:
                i = 1
            else:
                i = 0
            G[u, v] = i

    print("G =")
    print(shortstr(G))

    GG = dot2(G, G.transpose())
    if GG.sum()==0:
        print("self-dual code, p=%d mod 4" % (p%4))
    else:
        print("not a self-dual code, p=%d mod 4" % (p%4))

    m = rank(G)
    #assert m == N/2
    print("F_2 rank =", m)
    #print("det:", numpy.linalg.det(G))

    H = find_kernel(G)
    H = array2(H)
    print()
    print("H =")
    print(shortstr(H))

    # -----------------------------------

    print()
    print("non extended:")
    print("G =")
    G1 = G[:-1, :-1]
    print(shortstr(G1))

    #print("det:", numpy.linalg.det(G1.astype(numpy.float)))

    m = rank(G1)
    print("rank =", m)
    #assert m == N/2

    H1 = find_kernel(G1)
    H1 = array2(H1)
    print()
    print("H =")
    print(shortstr(H1))

    GG = dot2(G, G.transpose())
    if GG.sum()==0:
        print("self-dual code")
    else:
        print("not a self-dual code")

    i = iter(nonresidues).__next__()
    idxs = [(j*i)%p for j in range(p)]
    assert len(set(idxs)) == p
    NG1 = G1[:, idxs]
    print("G~ =")
    print(shortstr(NG1))
    print("G . G~ =")
    print(dot2(G1, NG1.transpose()).sum())

    NH1 = find_kernel(NG1)
    NH1 = array2(NH1)
    print()
    print("NH =")
    print(shortstr(NH1))
    #print(dot2(H1, NH1.transpose()))

    print("linear_independent(G):")
    print(shortstr(linear_independent(G1)))
    print("linear_independent(G~):")
    print(shortstr(linear_independent(NG1)))

    # -----------------------------------
    # code should be fixed by PSL(2, p).

    G1 = zeros2(N, N)
    for u in range(N):
      for v in range(N):
        G1[u, v] = G[u, neginv(v)]

    #print()
    #print(shortstr(G1))
    #print()

    # still in the codespace:
    assert dot2(G1, H.transpose()).sum() == 0

    # -----------------------------------
    # build PSL(2, p) via action on P(F_p)

    basis = list(range(N))
    perm = dict((i, (i+1)%p) for i in range(p))
    perm[p] = p # infty
    A = Perm(perm, basis)

    perm = {}
    for u in range(N):
        v = neginv(u)
        perm[u] = v
    B = Perm(perm, basis)
    PSL = Group.generate([A, B])
    print("|PSL(2,%d)| = %d"%(p, len(PSL)))

    print(residues)
    count = 0
    for g in PSL:
        r1 = set(g[r] for r in residues)
        if r1==residues:
            count += 1
    print("count =", count)

#    g = PSL[5]
#    print(g.orbits())

#    for u in span(G):
#        v = array2([u[A[i]] for i in basis])
#        if not eq2(u, v):
#            continue
#        v = array2([u[B[i]] for i in basis])
#        if not eq2(u, v):
#            continue
#        print(u)

    if p==23:
        # Conway & Sloane, p274
        cycles = [(p,), list(range(23))]
        alpha = Perm.fromcycles(cycles, basis)
        cycles = [(p,), (15, 7, 14, 5, 10, 20, 17, 11, 22, 21, 19), 
            (0,), (3, 6, 12, 1, 2, 4, 8, 16, 9, 18, 13)]
        beta = Perm.fromcycles(cycles, basis)
        cycles = [(p, 0), (15, 3), (7, 13), (14, 18), (5, 9), (10, 16),
            (20, 8), (17, 4), (11, 2), (22, 1), (21, 12), (19, 6)]
        gamma = Perm.fromcycles(cycles, basis)
        cycles = [(p,), (14, 17, 11, 19, 22), (15,), (20, 10, 7, 5, 21,),
            (0,), (18, 4, 2, 6, 1), (3,), (8, 16, 13, 9, 12)]
        delta = Perm.fromcycles(cycles, basis)

        G = Group.generate([alpha, beta, gamma]) # PSL(2,23)
        assert len(G) == 6072

        # way too big to generate:
        #G = Group.generate([alpha, beta, gamma, delta]) # M_24 
        #print(len(G))


if __name__ == "__main__":

    even_code()
    #main()


