#!/usr/bin/env python

from time import time
start_time = time()

import numpy

from bruhat.solve import parse, span, shortstr, array2, rank, solve, intersect, row_reduce, zeros2
from bruhat.isomorph import Tanner, search
from bruhat.equ import Equ, quotient
from bruhat.action import Perm, Group, mulclose
from bruhat.util import factorize
from bruhat.qcode import QCode
from bruhat.argv import argv


def preproc():
    data = open(argv.next()).read()
    
    data = data.replace(".", "")
    data = data.replace("3", "")
    data = data.replace(" ", "")
    data = data.replace("\n\n", "\n")
    data = data.replace("[[", "[")
    data = data.replace("]]", "]\n")
    data = data.replace("[", "")
    data = data.replace("]", "")
    
    #data = data.replace("\n\n", ".")
    #data = data.replace("\n", "")
    #data = data.replace(".", "\n")
    
    f = open('oeqc.txt', 'w')
    f.write(data)
    f.close()


def find_isos(H1, H2):
    src = Tanner.build2(H1, H1)
    tgt = Tanner.build2(H2, H2)
    G = [f for f in search(src, tgt)]
    return G

def is_iso(H1, H2):
    src = Tanner.build2(H1, H1)
    tgt = Tanner.build2(H2, H2)
    for f in search(src, tgt):
        return True
    return False


def read_codes():
    print("read_codes")
    data = open("oeqc.txt").read()
    items = data.split("\n\n")
    #print(len(items))
    
    found = set()
    for item in items:
        item = item.strip()
        if not item or item[0] not in "01":
            #print("skipping:", repr(item))
            continue
        try:
            A = parse(item)
        except ValueError:
            print("parse failed:", item)
            raise
        yield A

def get_dist():
    for A in read_codes():
        n = A.shape[1]
        dist = [0 for i in range(n)]
        for v in span(A):
            dist[v.sum()] += 1
        dist = tuple(dist)
        if dist not in found:
            print(dist)
            found.add(dist)


def get_autos(H):
    rows = [v for v in span(H) if v.sum()]
    V = array2(rows)
    print(shortstr(V))
    count = 0
    for f in find_isos(V, V):
        #print(f)
        print(".", end="", flush=True)
        count += 1
    print(count)


            
def gap_fmt(perm):
    cs = perm.cycles()
    ss = []
    for c in cs:
        if len(c)==1:
            continue
        c = [i+1 for i in c]
        ss.append(tuple(c))
    return ''.join(str(s) for s in ss)

def gap_code(perms):
    s = [gap_fmt(perm) for perm in perms]
    s = "Group(%s);"%(', '.join(s))
    return s

def get_autos_nauty(H):
    print("get_autos_nauty", H.shape)
    #dist = {i:0 for i in range(H.shape[1]+1)}
    #for v in span(H):
    #    dist[v.sum()] += 1
    #print(dist)
    #return
    rows = [v for v in span(H) if v.sum()]
    V = array2(rows)
    #print(shortstr(V))
    m, n = V.shape
    from pynauty import Graph, autgrp
    g = Graph(m+n) # checks + bits
    for check in range(m):
        bits = [m + bit for bit in range(n) if V[check, bit]]
        g.connect_vertex(check, bits)
    for bit in range(n):
        checks = [check for check in range(m) if V[check, bit]]
        g.connect_vertex(m+bit, checks)
    g.set_vertex_coloring([set(range(m)), set(range(m, m+n))])
    #print(g)
    aut = autgrp(g)
    #print(aut)

    gen = aut[0]
    items = list(range(m+n))
    perms = []
    for perm in gen:
        #print(perm)
        perm = Perm(perm, items)
        perms.append(perm)
    #print(gap_code(perms))
    return perms

    #G = Group.generate(perms)
    G = mulclose(perms, verbose=True)
    N = len(G) # 322560, (C2 x C2 x C2 x C2) : A8
    print("autos:", N, factorize(N))



def main_rm():
    H = parse("""
    11111111........
    ....11111111....
    ........11111111
    11..11..11..11..
    .11..11..11..11.
    """)
    assert rank(H) == len(H)
    print(H.shape)
    
    perms = get_autos_nauty(H)

    #G = mulclose(perms, verbose=True)
    #N = len(G) # 322560, (C2 x C2 x C2 x C2) : A8
    #print("autos:", N, factorize(N))

    v = '1111............'
    n = len(v)
    m = len(perms[0].items)
    items = list(range(m-n, m))
    print(items)
    perms = [g.restrict(items) for g in perms]
    items = list(range(n))
    #perms = [Perm(perm.perm, items) for perm in perms]
    #perms = [Perm(dict((k-m+n,v-m+n) for (k,v) in g.perm.items()), items) for g in perms]
    perms = [dict((k-m+n,v-m+n) for (k,v) in g.perm.items()) for g in perms]
    perms = [[g[i] for i in items] for g in perms]
    print(perms)
    orbit = {v}
    bdy = {v}
    while bdy:
        _bdy = set()
        for g in perms:
            for v in bdy:
                v1 = ''.join(v[i] for i in g)
                if v1 not in orbit:
                    _bdy.add(v1)
                    orbit.add(v1)
        bdy = _bdy
    print(len(orbit))
    orbit = list(orbit)
    orbit.sort()
    orbit = "\n".join(orbit)
    J0 = parse(orbit)
    J0 = row_reduce(J0)
    print(shortstr(J0), J0.shape)
    JX = zeros2(*J0.shape, 2)
    JX[:,:,0] = J0
    JZ = zeros2(*J0.shape, 2)
    JZ[:,:,1] = J0
    J = numpy.concatenate((JX, JZ))
    print(J.shape)

    code = QCode.build_gauge(J)
    print(code)
    L = code.get_logops()
    print(L.shape)



def main_golay():
    """
    012345678901234567890123456789
    1...........11...111.1.1
    .1...........11...111.11
    ..1.........1111.11.1...
    ...1.........1111.11.1..
    ....1.........1111.11.1.
    .....1......11.11..11..1
    ......1......11.11..11.1
    .......1......11.11..111
    ........1...11.111...11.
    .........1..1.1.1..1.111
    ..........1.1..1..11111.
    ...........11...111.1.11
    """

    # Golay code
    H0 = parse("""
    1...........11...111.1.1
    .1...........11...111.11
    ..1.........1111.11.1...
    ...1.........1111.11.1..
    ....1.........1111.11.1.
    .....1......11.11..11..1
    ......1......11.11..11.1
    .......1......11.11..111
    ........1...11.111...11.
    .........1..1.1.1..1.111
    ..........1.1..1..11111.
    ...........11...111.1.11
    """)

    print(H0.shape)
    n = H0.shape[1]
    words = {i:[] for i in range(n+1)}
    dist = [0 for i in range(n+1)]
    for v in span(H0):
        d = v.sum()
        dist[d] += 1
        words[d].append(v)
    dist = tuple(dist)
    print(dist)

    Hs = []
    for word in words[16]:
        rows = []
        for u in words[8]:
            if numpy.alltrue(u*word == u):
                rows.append(u)
        H = array2(rows)
        H = row_reduce(H)
        Hs.append(H)
        #print(shortstr(H))
        #print()
    print(len(Hs))

    for H1 in Hs:
        w = intersect(Hs[0], H1)
        print(len(w), end=' ')
    print()

    return

    Hs = []
    for i in range(12):
        rows = [j for j in range(12) if H0[j, i+12]==0]
        H1 = H0[rows]
        Hs.append(H1)
        #print(shortstr(H1), H1.shape)
        #print()

    for Ha in Hs:
      for Hb in Hs:
        Hab = intersect(Ha, Hb)
        print(Hab.shape[0], end=' ')
      print()

    #get_autos_nauty(H)

        
def find_iso():
    # broken: must take span of H
    codes = list(read_codes())
    print(len(codes))

    #codes = codes[:100]
    #codes = [H.tobytes() for H in codes]
    equs = quotient(codes, is_iso, verbose=True)
    found = set()
    for equ in equs:
        found.add(equ.top)
    print()
    print("distinct codes:", len(found))
    for equ in found:
        print(shortstr(equ.items[0]))
        print()

if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "test_all"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds"%(time() - start_time))
    print()


