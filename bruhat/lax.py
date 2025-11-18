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
from bruhat.algebraic import Algebraic
from bruhat.action import mulclose
from bruhat.argv import argv

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




def test_sp():

    n = argv.get("n", 2)

    Cliff = Cliff = Algebraic.Sp(2*n, p)
    F = Cliff.invariant_form
    U = Cliff.get_zip_uturn()

    Cliff1 = Algebraic.Sp(2, p)
    U1 = Cliff1.get_zip_uturn()
    print("Cliff1", len(Cliff1))

    if p < 5 and n <= 4:

        count = 0
        orbit = []
        for H in Cliff.qchoose(n):
            assert H == H.normal_form()
            #print(H)
            count += 1
            orbit.append(H)

    #elif p == 5 and n==4:
    else:
        orbit = set()
        for H in Cliff.qchoose(n):
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

    orbits = []
    remain = set(orbit) 
    del orbit
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
            print("orbit:", len(orbit), sig)
        #remain.difference_update(orbit)
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
    for _ in range(10):
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
            print("orbit:", len(orbit))


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
        #(1, (q-2)*(q+1)**4 * q**3 * (q-1)**1),
        (1, (q - 2) * (q - 1)**3 * q**3 * (q + 1)**3),
        (3, (q+1)**4 * q**2 * (q-1)**3),
    ]
    for m,item in items:
        v = v + m*item

    for p in [2,3,5,7]:
        print("p = %s"%p)
        for m,item in items:
            print("%10d"%item.subs(q=p), "\t%s" % factor(item))
        print()

    #print(v)
    #print(factor(f-v))

    poly = (q - 2) * (q - 1)**3 * q**3 * (q + 1)**3
    #print( poly )
    #print( poly / (q-2) )

    poly = (q - 1)**3 * q**3 * (q + 1)**3
    for p in [2,3,5,7]:
        print("p=%d"%p, poly.subs(q=p))






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

    

    




