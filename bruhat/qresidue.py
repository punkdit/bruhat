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
from action import Perm, Group


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

    residues = set()
    for i in range(1, p):
        residues.add((i*i)%p)
    print("residues:", residues)

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
    print("rank =", m)
    print("det:", numpy.linalg.det(G))

    H = find_kernel(G)
    H = array2(H)
    print()
    print(shortstr(H))

    return

    # -----------------------------------

    print()
    print("non extended:")
    G1 = G[:-1, :-1]
    print(shortstr(G1))

    print("det:", numpy.linalg.det(G1.astype(numpy.float)))

    m = rank(G1)
    print("rank =", m)
    #assert m == N/2

    H1 = find_kernel(G1)
    H1 = array2(H1)
    print()
    print(shortstr(H1))

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

    main()


