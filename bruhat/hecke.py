#!/usr/bin/env python

"""
build dual containing codes (self-dual CSS) from hecke operators
"""


from time import time
start_time = time()
from random import choice
from functools import reduce
from operator import mul
from random import random

import numpy

from bruhat.gset import Group, Perm, mulclose
from bruhat.argv import argv
from bruhat.solve import zeros2, dot2, span, shortstr, linear_independent, parse, enum2
from bruhat.util import choose
from bruhat.orthogonal import get_logops, row_reduce
from bruhat.oeqc import is_triorthogonal
from bruhat.todd_coxeter import Schreier


def colour(H):
    m, n = H.shape
    #print(m, n)
    A = numpy.dot(H, H.transpose())
    #print(A.shape)
    #print(shortstr(A!=0))
    nbds = {i:[] for i in range(m)} # map idx -> nbd list
    for i in range(m):
     for j in range(m):
        if A[j,i] and i!=j:
            nbds[i].append(j)
    #print(nbds)

    found = {} # map idx -> colour
    remain = list(range(m))
    found[0] = 0
    j = nbds[0][0]
    found[j] = 1
    #print(found)
    bdy = nbds[0][1:]
#    maxk = 0
#    while bdy:
#        _bdy = []
#        for i in bdy:
#            if i in found:
#                continue
#            ks = [found.get(j) for j in nbds[i]]
#            k = 0
#            while k in ks:
#                k += 1
#            found[i] = k
#            maxk = max(k, maxk)
#            for j in nbds[i]:
#                _bdy.append(j)
#        bdy = _bdy
#        #print(found, bdy)

    while len(found) < m:
        N = len(found)
        for i in range(m):
            if i in found:
                continue
            cs = [found.get(j) for j in nbds[i] if found.get(j) is not None]
            if len(cs) < 2:
                continue
            for k in range(m):
                if k not in cs:
                    found[i] = k
                    break
            else:
                assert 0
        if len(found) == N:
            break
    maxk = max(found.values())

    if len(found) < m:
        print(len(found), " <", m)
    #ks = list(set((found.values())))
    #ks.sort()
    ks = {i:[] for i in range(maxk+1)}
    for i in found:
        ks[found[i]].append(i)
    print("colour:", [len(k) for k in ks.values()])
    print("%d-colourable checks" % len(ks))


def test_colour():
    H = parse("""
    """)
    colour(H)
    print(H.shape)


def make_hecke(X, Y):
    #print("make")
    G = X.src
    assert Y.src is G
    N = len(G)
    m, n = X.tgt.rank, Y.tgt.rank
    #print(X.send_perms)
    #print(X.tgt)
    remain = set((i,j) for i in range(m) for j in range(n))
    while remain:
        A = zeros2(m, n)
        ij = iter(remain).__next__()
        i, j = ij
        for idx in range(N):
            x = X.tgt[X.send_perms[idx]]
            y = Y.tgt[Y.send_perms[idx]]
            ij = x[i], y[j]
            if A[ij] == 0:
                A[ij] = 1
                remain.remove(ij)
        yield A


def get_lower_distance(H, L, max_d=4):
    m, n = H.shape
    idxs = list(range(n))
    for w in range(2, max_d):
      for ii in choose(idxs, w):
        v = zeros2(n, 1)
        for i in ii:
            v[i, 0] = 1
        if dot2(H, v).sum() == 0 and dot2(L, v).sum():
            return w
    return None


def find_lower_distance(H, L, d):
    # search all weight d vectors
    m, n = H.shape
    idxs = list(range(1, n))
    for ii in choose(idxs, d-1):
        v = zeros2(n, 1)
        v[0,0] = 1
        for i in ii:
            v[i, 0] = 1
        if dot2(H, v).sum() == 0 and dot2(L, v).sum():
            return True # distance <= d
    return False # distance > d


def get_distance_xz(H, L, min_d=2):

    m, n = H.shape
    if len(L) > 24:
        d = min(L.sum(1))
        return None, d

    if len(L) + len(H) > 24:
        d = n
        for u in span(L):
            d1 = u.sum()
            if d1 == 0:
                continue
            d = min(d, d1)
            if d<=min_d:
                return None, d
        return None, d

    vs = list(span(H))

    d = n
    for u in span(L):
        if u.sum() == 0:
            continue
        for v in vs:
            uv = (u+v)%2
            d = min(d, uv.sum())
        if d<=min_d:
            return d,d
    return d,d


def find_distance_random(A, L):
    m, n = A.shape
    w = n
    for v in L:
        d = v.sum()
        w = min(w, d)
        for trial in range(100):
          for a in A:
            v1 = (v+a)%2
            w1 = v1.sum()
            if w1 < w:
                v = v1
                w = w1
            elif w1 == w and random() < 0.5:
                v = v1
    return w


def minimize(A, v, min_d=2):
    bdy = {v.tobytes() : v}
    found = set(bdy.keys())
    d = v.sum()
    while bdy and d>=min_d:
        #print(d, end=" ", flush=True)
        _bdy = []
        for v in bdy.values():
            for u in A:
                v1 = (v+u)%2
                if v1.tobytes() in found:
                    continue
                found.add(v1.tobytes())
                d1 = v1.sum()
                if d1 < min_d:
                    print("!",end="",flush=True)
                    _bdy = [v1]
                    d = d1
                    break
                elif d1 == min_d:
                    _bdy.append(v1)
                    #print("+",end="",flush=True)
                    # keep going
                #else:
                #    print("[%d]"%d1, end="")
            else:
                continue
            #print("!!")
            break
        bdy = {v.tobytes(): v for v in _bdy}
        #print("(%s)"%len(bdy), end="\n")
    return d


def find_distance_dynamic(A, L, min_d=2):
    #print(A.shape)
    k, n = L.shape
    assert k<=28, "really??"
    count = 0
    d = n
    for v in enum2(k):
        count += 1
        #if count>100:
        #    break
        u = dot2(v, L)
        if u.sum() == 0:
            continue
        d1 = minimize(A, u, min_d)
        if d1 < d:
            d = d1
            print(d)
    return d

    

def get_distance(H, L=None, min_d=2):
    # return lower & upper bound on distance
    m, n = H.shape
    if L is None:
        L = get_logops(H)

    max_d = 3 if n>100 else 4
    ld = get_lower_distance(H, L, max_d)
    if ld is not None:
        return ld, ld

    return get_distance_xz(H, L, min_d)


def get_distance_css(code):
    Hx, Hz, Lx, Lz = code.Hx, code.Hz, code.Lx, code.Lz
    dx = get_lower_distance(Hx, Lz)
    dz = get_lower_distance(Hz, Lx)
    if dx is None:
        _, dx = get_distance_xz(Hx, Lx)
    if dz is None:
        _, dz = get_distance_xz(Hz, Lz)

    return dx, dz

    
def get_group(name):
    if "*" in name:
        Gs = [get_group(name) for name in name.split("*")]
        G = reduce(mul, Gs)
        return G

    graph = None
    stem = name[0]
    n = int(name[1:]) if name[1:] else None
    if stem == "A":
        G = Group.alternating(n)
    elif stem == "S":
        G = Group.symmetric(n)
    elif stem == "I":
        G = Group.dihedral(n)
    elif stem == "C":
        G = Group.cyclic(n)
    elif stem == "K":
        ngens = 2
        a, ai, b, bi = range(2*ngens)
    
        # Klein Quartic
        rels = [ (ai, a), (bi, b), (a,)*3, (b,)*7, (a,b)*2 ]
        rels += [ (a,bi,bi)*4 ]
        graph = Schreier(2*ngens, rels)
        graph.build()
        assert len(graph) == 168, len(graph)

    elif stem == "D":
        graph = Schreier.make_D(n) # D4 has order 192

    elif stem == "B":
        #graph = Schreier.make_B(3) # nada
        graph = Schreier.make_B(n) # B4 has order 384

    elif stem == "Z":
        ngens = 3
        r, g, b = (0,), (1,), (2,)
        rels = [r*2, g*2, b*2, (r+b)*3, (r+g)*3, (b+g)*3]
        rels += [(b+r+g)*4]
        graph = Schreier(ngens, rels)
        graph.build(maxsize=100000)

    elif stem == "G":
        from bruhat.qcode import Geometry
        key = argv.get("key", (4, 4))
        idx = argv.get("idx", 24)
        geometry = Geometry(key, idx)
        graph = geometry.build_graph()

    else:
        assert 0, name

    if graph is not None:
        print("graph:", len(graph))
        G = graph.get_group()
        N = len(G)
        perms = []
        for g in G:
            perm = [g[i] for i in range(N)]
            perm = Perm(perm)
            perms.append(perm)
        G = Group(perms)

    return G


def get_actions(G, Hs):
    Xs = [G.action_subgroup(H) for H in Hs]
    print("|Xs| =", len(Xs))

    i = 0
    while i < len(Xs):
        j = i+1
        Xi = Xs[i]
        while j < len(Xs):
            Xj = Xs[j]
            print("/", end="", flush=True)
            if Xi.rank==Xj.rank and Xs[i].isomorphic(Xs[j]):
                Xs.pop(j)
                print(".", end="", flush=True)
            else:
                j += 1
        i += 1
    print()
    return Xs


def main_rand():
    name = argv.next()
    G = get_group(name)
    print("|G| =", len(G))
    gen = [g for g in G if g*g!=g]

    n = argv.n
    Hs = []
    while len(Hs) < argv.get("subgroups", 5):
        H = Group.generate([choice(gen) for i in range(1)])
        nn = len(G) // len(H)
        if n is not None and nn>n:
            continue
        if 4 < nn:
            Hs.append(H)
    #while len(Hs) < 40:
    #    H = Group.generate([choice(gen) for i in range(2)])
    #    if 4*len(H) < len(G):
    #        Hs.append(H)
    Hs.sort(key = len, reverse=True)
    print([len(H) for H in Hs])

    #if n is not None:
    #    Hs = [H for H in Hs if len(G)//len(H)<=n]
    #    print([len(H) for H in Hs])
    
    #Xs = get_actions(G, Hs)
    Xs = [G.action_subgroup(H) for H in Hs]
    print("|Xs| =", len(Xs))


    print("make_hecke...")
    ops = {}
    count = 0
    N = len(Xs)
    for i in range(N):
      for j in range(i, N):
        ops[i,j] = []
        X, Y = Xs[i], Xs[j]
        for H in make_hecke(X, Y):
            rw, cw = H.sum(1).max(), H.sum(0).max()
            if rw == 1 or cw == 1 or rw>20 or cw>20:
                continue
            if dot2(H, H.transpose()).max() == 0:
                ops[i,j].append(H)
                print(H.shape, end=' ', flush=True)
                count += 1
        #print(H.shape, end='', flush=True)
    
    print("\nops:", count)

    min_d = argv.get("min_d", 3)
    found = set()
    for key in ops:
      for A in ops[key]:
        #print(A, A.shape)
        H = row_reduce(A)
        m, n = H.shape
        k = n-2*m
        if k<=0:
            continue
        #ws = H.sum(0)
        #assert ws.min()
        #HHt = dot2(H, H.transpose())
        #assert HHt.max() == 0
        #if HHt.max():
        #    continue

        dmin, dmax = get_distance(H)
        if dmax<min_d:
            print(choice("/\\"), end='', flush=True)
            continue
        w = A.sum(1)[0]
        desc = "[[%d, %d, %s<=d<=%s]]_%d"%(n, k, dmin, dmax, w)
        if desc in found:
            continue
        found.add(desc)
        print()
        print("A.shape =", A.shape, "rank =", len(H))
        #A = linear_independent(A)
        print(shortstr(A))
        print(desc)
        print("is_triorthogonal:", is_triorthogonal(A))
        colour(H)
        print()




def main():

    name = argv.next()
    G = get_group(name)
    print("|G| =", len(G))
    Hs = list(G.subgroups())
    Hs = [H for H in Hs if len(H)<len(G)]

    print("proper subgroups:", len(Hs))

    Xs = get_actions(G, Hs)
    print("|Xs| =", len(Xs))

    print("make_hecke...")
    ops = {}
    N = len(Xs)
    for i in range(N):
      for j in range(N):
        ops[i,j] = []
    for i, X in enumerate(Xs):
      for j, Y in enumerate(Xs):
        for A in make_hecke(X, Y):
            A = linear_independent(A)
            m, n = A.shape
            if m<3 or n<5:
                continue
            ops[i,j].append(A)
    print("ops:", len(ops))

    min_d = argv.get("min_d", 3)

    if argv.css:
        print("CSSCode's")
        from qupy.ldpc.css import CSSCode
        for i in range(N):
         for j in range(N):
          for k in range(j, N):
            Hxs = ops[j,i]
            Hzs = ops[k,i]
            for Hx in Hxs:
             for Hz in Hzs:
                assert Hx.shape[1] == Hz.shape[1]
                mx, n = Hx.shape
                mz, n = Hz.shape
                k = n-mx-mz
                if k==0 or k>=n//2:
                    continue
                if Hx.tobytes() == Hz.tobytes(): # self-dual
                    continue
                if dot2(Hx, Hz.transpose()).sum() == 0:
                    code = CSSCode(Hx=Hx, Hz=Hz)
                    #dx, dz = code.distance(2)
                    dx, dz = get_distance_css(code)
                    if (dx is not None and dx<min_d) or (dz is not None and dz < min_d):
                        print(".", end="", flush=True)
                        continue
                    print()
                    print(type(dx) is int, dx<min_d)
                    print("Hx=")
                    print(shortstr(Hx))
                    print("Hz=")
                    print(shortstr(Hz))
                    print("[[%d, %d, %s %s]"%(n, k, dx, dz), min_d)
        return

    found = set()
    for key in ops:
      for A in ops[key]:
        #print(A, A.shape)
        H = row_reduce(A)
        m, n = H.shape
        k = n-2*m
        if k<=0:
            continue
        ws = H.sum(0)
        if ws.min() == 0 or ws.max() == 1:
            continue
        HHt = dot2(H, H.transpose())
        if HHt.max():
            continue

        dmin, dmax = get_distance(H)
        if dmax<min_d:
            print(choice("/\\"), end='', flush=True)
            continue
        w = A.sum(1)[0]
        desc = "[[%d, %d, %s<=d<=%s]]_%d"%(n, k, dmin, dmax, w)
        if desc in found:
            continue
        found.add(desc)
        print()
        print("A.shape =", A.shape)
        A = linear_independent(A)
        print(shortstr(A))
        print(desc)
        print()




if __name__ == "__main__":
    fn = "main"
    if argv.rand:
        fn = "main_rand"
    if argv.test_colour:
        fn = "test_colour"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds.\n"%(time() - start_time))



