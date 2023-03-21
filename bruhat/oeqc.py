#!/usr/bin/env python

from time import time
start_time = time()

from bruhat.argv import argv

from bruhat.solve import parse, span, shortstr, array2, rank, solve, intersect
from bruhat.isomorph import Tanner, search
from bruhat.equ import Equ, quotient
from bruhat.action import Perm, Group, mulclose
from bruhat.util import factorize


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
    dist = {i:0 for i in range(H.shape[1]+1)}
    for v in span(H):
        dist[v.sum()] += 1
    print(dist)
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
    print(gap_code(perms))

    return
    #G = Group.generate(perms)
    G = mulclose(perms, verbose=True)
    N = len(G) # 322560, (C2 x C2 x C2 x C2) : A8
    print("autos:", N, factorize(N))



def main_autos():
    H = parse("""
    101001010011
    011010010011
    000111010011
    000000110101
    000000001111
    """)

    H = parse("""
    11111111........
    ....11111111....
    ........11111111
    11..11..11..11..
    .11..11..11..11.
    """)
    assert rank(H) == len(H)

    from bruhat.triply_even import codes24
    H_golay = codes24.get("g_{24}")
    print(shortstr(H_golay))

    W = intersect(H, H_golay)
    print(W.shape)

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


