#!/usr/bin/env python


import sys, os
import random
from random import randint, choice
from functools import reduce
from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)
from operator import add, mul, lshift
from math import prod

import numpy

from sage.all_cmdline import PolynomialRing, ZZ, factor
from sage import all_cmdline as sage


from bruhat.matrix_sage import Matrix as SMatrix
from bruhat.algebraic import Algebraic, Matrix
from bruhat.gset import mulclose, Perm, Group, cayley
from bruhat.argv import argv
from bruhat import solve as lin
from bruhat.util import cross, choose

p = argv.get("p", 2)


@cache
def get_bits(n, arity=2):
    bits = list(numpy.ndindex((arity,)*n))
    assert len(bits) == arity**n
    bits.sort(key = sum)
    bits = tuple(bits)
    return bits


@cache
def get_idxs(n):
    idxss = []
    for bits in get_bits(n):
        idxs = tuple(i for (i,ii) in enumerate(bits) if ii==1)
        idxss.append(idxs)
    return tuple(idxss)



def get_lower(H):
    #print("get_lower")
    #print(H.A)

    ring = sage.GF(p)
    H = SMatrix(ring, H)
    #print(H)
    m, nn = H.shape
    n = nn//2

    counts = []
    for idxs in get_idxs(n):
        A = numpy.zeros((2*len(idxs), nn))
        I = numpy.array([[1,0],[0,1]])
        for i,ii in enumerate(idxs):
            A[2*i:2*i+2, 2*ii:2*ii+2] = I
        A = SMatrix(ring, A)
        A = A.intersect_rowspace(H)
        A = A.to_numpy().astype(int)
        #print(A, idxs, A.sum(1), end=' ')
        for rank, i in enumerate(A.sum(1)):
            if i==0:
                break
        else:
            rank += 1
        #print(rank)
        #print(A, A.shape, A.sum(1))
        counts.append(rank)
        #print()
    return counts




def sp_orbits():

    n = argv.get("n", 2)
    m = argv.get("m", n)

    Cliff = Cliff = Algebraic.Sp(2*n, p)
    F = Cliff.invariant_form
    U = Cliff.get_zip_uturn()

    Cliff1 = Algebraic.Sp(2, p)
    U1 = Cliff1.get_zip_uturn()
    print("Cliff1", len(Cliff1))

    if p < 5 and n <= 4:

        count = 0
        orbit = []
        for H in Cliff.qchoose(m):
            assert H == H.normal_form()
            #print(H)
            count += 1
            orbit.append(H)

    #elif p == 5 and n==4:
    else:
        orbit = set()
        for H in Cliff.qchoose(m):
            assert H == H.normal_form()
            break
        orbit.add(H)
        bdy = list(orbit)
        while bdy:
            _bdy = []
            for H in bdy:
                for g in Cliff.gen:
                    J = (H*g.t).normal_form()
                    if J in orbit:
                        continue
                    _bdy.append(J)
                    orbit.add(J)
            bdy = _bdy
            print("[%d:%d]"%(len(orbit),len(bdy)), end='', flush=True)
        print()

#    else:
#        orbit = set()
#        for H0 in Cliff.qchoose(n):
#            assert H0 == H0.normal_form()
#            break
#        orbit.add(H0)
#        while len(orbit) < 1000:
#            H = H0
#            assert Cliff.is_isotropic(H)
#            for _ in range(2+len(orbit)):
#                g = choice(Cliff.gen)
#                assert g.t*F*g == F
#                H = H*g
#                assert Cliff.is_isotropic(H), g
#            H = H.normal_form()
#            orbit.add(H)


    I = Cliff1.I
    gen = []
    for a in Cliff1.gen:
        for i in range(n):
            ops = [I]*n
            ops[i] = a
            op = reduce(lshift, ops)
            op = U*op*U.t
            #assert op in Cliff
            gen.append(op)
    print("LCliff:", len(Cliff1) ** n )

    limit = argv.get("limit", None)
    verbose = argv.verbose

    print("total:", len(orbit))

    remain = set(orbit) 
    del orbit
    count = 0
    while remain:
        if verbose:
            print("\t(%d)" % len(remain), end="")
        H = remain.pop()
        assert Cliff.is_isotropic(H)
        #orbit = {(H*g.t).normal_form() for g in LCliff}
        bdy = [H]
        orbit = set(bdy)
        while bdy:
            _bdy = []
            for H in bdy:
                for g in gen:
                    J = (H*g.t).normal_form()
                    if J in orbit:
                        continue
                    _bdy.append(J)
                    orbit.add(J)
                    if J in remain:
                        remain.remove(J)
            bdy = _bdy
            if verbose:
                print("[%d:%d]"%(len(orbit),len(bdy)), end='', flush=True)
            if limit and len(orbit) > limit:
                print("FAIL")
                break
        else:
            if verbose:
                print()
            H = iter(orbit).__next__()
            sig = get_lower(H*U)
            print("orbit:", len(orbit), sig, flush=True)
            count += 1
        #remain.difference_update(orbit)
    print()
    print("total orbits =", count)
    
test_sp = sp_orbits


def sp_stab():

    n = argv.get("n", 2)
    m = argv.get("m", n)

    Cliff = Cliff = Algebraic.Sp(2*n, p)
    F = Cliff.invariant_form
    U = Cliff.get_zip_uturn()

    Cliff1 = Algebraic.Sp(2, p)
    U1 = Cliff1.get_zip_uturn()
    print("Cliff1", len(Cliff1))

    if p < 5 and n <= 4:

        count = 0
        orbit = []
        for H in Cliff.qchoose(m):
            assert H == H.normal_form()
            #print(H)
            count += 1
            orbit.append(H)

    #elif p == 5 and n==4:
    else:
        orbit = set()
        for H in Cliff.qchoose(m):
            assert H == H.normal_form()
            break
        orbit.add(H)
        bdy = list(orbit)
        while bdy:
            _bdy = []
            for H in bdy:
                for g in Cliff.gen:
                    J = (H*g.t).normal_form()
                    if J in orbit:
                        continue
                    _bdy.append(J)
                    orbit.add(J)
            bdy = _bdy
            print("[%d:%d]"%(len(orbit),len(bdy)), end='', flush=True)
        print()

    I = Cliff1.I
    gen = []
    for a in Cliff1.gen:
        for i in range(n):
            ops = [I]*n
            ops[i] = a
            op = reduce(lshift, ops)
            op = U*op*U.t
            #assert op in Cliff
            gen.append(op)
    print("LC:", len(Cliff1) ** n )

    N = len(Cliff1) ** n
    LC = mulclose(gen, verbose=True)
    assert len(LC) == N

    print("LC:", N)

    print(len(orbit))
    for H in orbit:
        stab = []
        for g in LC:
            J = (H*g.t).normal_form()
            if H==J:
                stab.append(g)
        G = cayley(stab)
        s = G.structure_description().replace(" ", "")
        print("%s:%d"%(s,N//len(stab)), end=' ', flush=True)

    print()
        


def sample_sp():

    p = argv.get("p", 2)
    n = argv.get("n", 2)

    Cliff = Cliff = Algebraic.Sp(2*n, p)
    F = Cliff.invariant_form
    U = Cliff.get_zip_uturn()

    Cliff1 = Algebraic.Sp(2, p)
    U1 = Cliff1.get_zip_uturn()
    print("Cliff1", len(Cliff1))

    def sample():
        for H in Cliff.qchoose(n):
            assert H == H.normal_form()
            break
        assert Cliff.is_isotropic(H)
        for _ in range(argv.get("size", 10)):
            g = choice(Cliff.gen)
            #assert g.t*F*g == F
            H = H*g.t
        H = H.normal_form()
        assert Cliff.is_isotropic(H), g
        return H

    I = Cliff1.I
    gen = []
    for a in Cliff1.gen:
        for i in range(n):
            ops = [I]*n
            ops[i] = a
            op = reduce(lshift, ops)
            op = U*op*U.t
            #assert op in Cliff
            gen.append(op)
    #print("LCliff:", len(Cliff1) ** n )

    limit = argv.get("limit", None)
    verbose = argv.verbose

    #while 1:
    found = set()
    for _ in range(argv.get("trials", 10)):
        H = sample()
        bdy = [H]
        orbit = set(bdy)
        while bdy:
            _bdy = []
            for H in bdy:
                for g in gen:
                    J = (H*g.t).normal_form()
                    if J in orbit:
                        continue
                    _bdy.append(J)
                    orbit.add(J)
            bdy = _bdy
            if verbose:
                print("[%d:%d]"%(len(orbit),len(bdy)), end='', flush=True)
            if limit and len(orbit) > limit:
                print("FAIL")
                break
        else:
            if verbose:
                print()
            c = len(orbit)
            if c in found:
                print("/", end='', flush=True)
            else:
                print("\n(%d)" % c, end=' ', flush=True)
                found.add(c)
    print()


def test_qpoly():

    R = PolynomialRing(ZZ, "q")
    q = R.gens()[0]

    cs = [int(i) for i in "11122222111"]
    f = sum(j*q**i for (i,j) in enumerate(cs))
    print("f =", f)
    for p in [2,3,5,7]:
        print("f(%d) = %s" % (p, f.subs(q=p)))
    print()

    v = 0
    items = [
        (1, (q+1)**4),
        (6, (q+1)**3 * q * (q-1)),
        (3, ((q-1)*(q)*(q+1))**2 ),
        (4, (q-1)**2 * q * (q+1)**4),
        (1, (q+1)**4 * q * (q-1)**3),
        ((q - 2), (q - 1)**3 * q**3 * (q + 1)**3),
        (3, (q+1)**4 * q**2 * (q-1)**3),
    ]
    for m,item in items:
        v = v + m*item
    assert v == f

    for p in [2,3,5,7]:
        print("p = %s"%p)
        for m,item in items:
            print("%10d"%item.subs(q=p), "\t%s" % factor(item))
        print()

    #print(v)
    #print(factor(f-v))



def lc(M, i):
    #print("lc", i)
    n = len(M)
    nbd = [j for j in range(n) if M[i,j]]
    #print("nbd:", nbd)

    B = M[:, nbd][nbd, :]
    #print(B)
    B = Matrix(B.A + 1)
    #print(B)

    A = M.A.copy()
    for i in nbd:
      for j in nbd:
        if i!=j:
            A[i,j] = 1-A[i,j]
    M = Matrix(A)
    return M
    


def graph_lc():

    n = argv.get("n", 4)
    p = n*(n-1)//2
    count = 2**p

    items = set()
    idxs = [(i,j) for i in range(n) for j in range(i+1,n)]
    #print(idxs)
    for bits in numpy.ndindex((2,)*p):
        A = numpy.zeros((n,n), dtype=int)
        #print(bits)
        for (i,j),bit in zip(idxs, bits):
            A[i,j] = bit
            A[j,i] = bit
        M = Matrix(A)
        items.add(M)
        #print(M)
    #print()

    assert len(items) == count

    items = list(items)
    items.sort(key=str)

    remain = set(items)
    for M in remain:
        for i in range(n):
            M1 = lc(M,i)
            #print(M1)
            assert M1 in remain

    orbits = []
    while remain:
        M = remain.pop()
        found = {M}
        bdy = [M]
        while bdy:
            _bdy = []
            for M in bdy:
              for i in range(n):
                M1 = lc(M, i)
                assert lc(M1, i) == M
                if M1 in remain:
                    found.add(M1)
                    _bdy.append(M1)
                    remain.remove(M1)
            bdy = _bdy
        #print(len(found))
        orbits.append(list(found))
    orbits.sort(key = len)
    print([len(o) for o in orbits])

    gens = []
    lookup = {M:i for (i,M) in enumerate(items)}
    for i in range(n):
        perm = [lookup[lc(M,i)] for M in items]
        perm = Perm(perm)
        gens.append(perm)
        #print(perm)

    #G = mulclose(gens, verbose=True)
    #print(len(G))

    G = Group(gens = gens, build=False)
    #print(G.structure_description())

    #f = open("lc_graphs.gap", "w")
    #print("G := %s;"%G.gapstr(), file=f)
    #f.close()




def search_poly():

    R = sage.PolynomialRing(sage.ZZ, 'q')
    q = R.gens()[0]

    vals = argv.vals
    vals = tuple(2*v for v in vals)
    print("vals =", vals)
    qs = [3, 5, 7, 11][:len(vals)]

    lookup = {}
    for degree in range(8):
        #print("degree", degree)
        #items = range(-5, 6)
        items = [-1, 0, 1]
        for cs in cross( [items] * (1+degree) ):
            if cs[-1] == 0:
                continue
            p = sum( c*q**i for (i,c) in enumerate(cs) )
            key = tuple(p.subs(q=val) for val in qs)
            if key == vals:
                print(p, "=", str(sage.factor(p)).replace(" ", ""))
                return
            #lookup[key] = p
            #print(key, end=" ")
        #print()
    
    #if vals in lookup:
    #    print("found:", lookup[vals])
    #else:
    #    print("not found")


def interpolate():

    R = sage.PolynomialRing(sage.QQ, 'q')
    q = R.gens()[0]

    vals = argv.vals
    print("vals =", vals)

    qs = argv.qs or [3, 5, 7, 11, 13]
    qs = qs[:len(vals)]

    points = list(zip(qs, vals))
    print(points)
    p = R.lagrange_polynomial(points)
    print(p, "=", sage.factor(p))


def find_poly():

    vals = argv.vals
    print("vals =", vals)

    qs = [2,3,5,7,11,13][:len(vals)]

    exps = range(5)
    for i in exps:
     for j in exps:
      for k in exps:
        items = tuple(((q-1)**i) * (q**j) * ((q+1)**k) for q in qs)
        assert len(items) == len(vals)
        if items == vals:
            print(i,j,k)

def tripartite():
    q = argv.get("q", 2)
    val = argv.get("val", 168)
    exp = argv.get("exp")
    
    p = lambda q : ((q-1)**j) * (q**i) * ((q+1)**k)
    if exp is not None:
        i,j,k = exp
        print(exp, p(q))
        return
    exps = range(15)
    for i in exps:
     for j in exps:
      for k in exps:
        #if i>j or j>k:
        #    continue
        if p(q)==val:
            print((i,j,k), p(q), p(2))




def quadratic_matrix(): # quadratic_form is much faster

    n = argv.get("n", 2)
    q = argv.get("q", 3)
    p = n*(n+1)//2

    print(q**p)

    R = sage.GF(q)

    items = set()
    idxs = [(i,j) for i in range(n) for j in range(i,n)]
    assert len(idxs) == p
    #print(idxs)
    for bits in numpy.ndindex((q,)*p):
        A = numpy.zeros((n,n), dtype=int)
        #print(bits)
        for (i,j),bit in zip(idxs, bits):
            A[i,j] = bit
            A[j,i] = bit
        Q = SMatrix(R, A)
        items.add(Q)
        #print(Q)
    #print()

    print(len(items))

    G = sage.GL(n, R)
    gens = []
    for g in G.gens():
        M = SMatrix(R, g)
        gens.append(M)

    #G = mulclose(gens)
    #print(len(G))

    orbits = []
    remain = set(items)
    while remain:
        Q = remain.pop()
        orbit = [Q]
        bdy = [Q]
        while bdy:
            _bdy = []
            for Q in bdy:
                for M in gens:
                    R = M*Q*M.t
                    if R not in remain:
                        continue
                    _bdy.append(R)
                    remain.remove(R)
            bdy = _bdy
            orbit += bdy
            #print(len(orbit), len(remain))
        orbits.append(orbit)

    orbits.sort(key=len)
    counts = [len(orbit) for orbit in orbits]
    print(counts, sum(counts))



def quadratic_form():

    n = argv.get("n", 2)
    q = argv.get("q", 5)
    p = n*(n+1)//2

    K = sage.GF(q)
    u, = K.gens()

    if u==1:
        elements = [K(i) for i in range(q)]
    else:
        elements = [K(0)] + [u**i for i in range(q-1)]
    #print(elements, len(elements))

    vs = ["x%d"%i for i in range(n)]
    R = sage.PolynomialRing(K, vs)
    vs = R.gens()

    items = set()
    idxs = [(i,j) for i in range(n) for j in range(i,n)]
    assert len(idxs) == p
    #print(idxs)

    pairs = [(x,x) for x in vs]
    for (x,y) in choose(vs, 2):
        pairs.append( (x,y) )
    terms = [a*b for (a,b) in pairs]

    items = []
    for idxs in numpy.ndindex((q,)*p):
        F = sum( elements[idx]*term for (idx,term) in zip(idxs, terms) )
        items.append(F)
    print(len(items))

    V = SMatrix(R, vs).t

    G = sage.GL(n, K)
    gens = []
    for g in G.gens():
        M = SMatrix(R, g)
        MV = M*V
        subs = {v:MV.M[i,0] for (i,v) in enumerate(vs)}
        gens.append(subs)

    orbits = []
    remain = set(items)
    while remain:
        Q = remain.pop()
        orbit = [Q]
        bdy = [Q]
        while bdy:
            _bdy = []
            for Q in bdy:
                for replace in gens:
                    R = Q.subs(replace)
                    if R not in remain:
                        continue
                    _bdy.append(R)
                    remain.remove(R)
            bdy = _bdy
            orbit += bdy
            #print(len(orbit), len(remain))
        orbits.append(orbit)

    orbits.sort(key=len)
    counts = [len(orbit) for orbit in orbits]
    print(counts, sum(counts))


def geometry():

    n = argv.get("n", 2)
    q = argv.get("q", 3)
    p = n*(n+1)//2

    K = sage.GF(q)
    u, = K.gens()

    if u==1:
        elements = [K(i) for i in range(q)]
    else:
        elements = [K(0)] + [u**i for i in range(q-1)]
    #print(elements, len(elements))

    vs = ["x%d"%i for i in range(n)]
    R = sage.PolynomialRing(K, vs)
    vs = R.gens()

    items = set()
    idxs = [(i,j) for i in range(n) for j in range(i,n)]
    assert len(idxs) == p
    #print(idxs)

    pairs = [(x,x) for x in vs]
    for (x,y) in choose(vs, 2):
        pairs.append( (x,y) )
    terms = [a*b for (a,b) in pairs]

    items = []
    for idxs in numpy.ndindex((q,)*p):
        F = sum( elements[idx]*term for (idx,term) in zip(idxs, terms) )
        items.append(F)
    assert len(items) == q**p
    print("%d^%d = %d"%(q,p,q**p))

    nn = 2*n
    G = sage.Sp(nn,K)
    #print(G)
    print("|Sp(%d)| ="%nn, len(G))
    #print(dir(G))
    F = G.invariant_form()
    F = SMatrix(K, F)

    # convert block order to u-turn order
    U = numpy.empty((nn,nn), dtype=object)
    U[:] = K(0)
    cols = [i for i in range(n)] + list(reversed([n+i for i in range(n)]))
    for i in range(nn):
        U[i, cols[i]] = K(1)
    U = SMatrix(K, U)

    B = SMatrix.identity(K, n)
    found = set()
    for item in items:
        diffs = [sage.derivative(item, v) for v in vs]
        A = numpy.empty((n,n), dtype=object)
        for i in range(n):
            diff = sage.derivative(item, vs[i])
            for j in range(n):
                subs = {vs[k]:K(0) for k in range(n)}
                subs[j] = K(1)
                A[i,j] = diff.subs(subs)
        A = SMatrix(K, A)
        found.add(A)
        #print(item, A*B.t == B*A.t )
        assert A*B.t == B*A.t
        A = B.augment(A)
        assert A == A.row_reduce()
        #print(A)
        A = A*U # convert block order to u-turn order
        AFA = (A*F*A.t)
        assert AFA.is_zero()

    print(len(items), len(found))


if __name__ == "__main__":
    from time import time
    start_time = time()
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))

    

    




