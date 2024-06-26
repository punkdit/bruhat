#!/usr/bin/env python

"""
Use z3 to look for logical operators (for CSS codes),
thereby bounding/finding the distance.
"""


from time import time
start_time = time()
import os
from random import choice
from functools import reduce
from operator import mul, add
from random import random, randint

import numpy

from bruhat.gset import Group, Perm, mulclose
from bruhat.argv import argv
from bruhat.solve import zeros2, dot2, span, shortstr, linear_independent, parse, enum2, rank
from bruhat.util import choose
from bruhat.orthogonal import get_logops, row_reduce
from bruhat.oeqc import is_triorthogonal
from bruhat.todd_coxeter import Schreier



def make_checksum(H):
    m, n = H.shape

    lines = ["def check(v):"]
    for i in range(m):
        line = []
        for j in range(n):
            if H[i,j]==0:
                continue
            line.append("v[%d]"%(j,))
        line = "  if (%s)%%2:" % ("+".join(line),)
        lines.append(line)
        lines.append("    return False")
    lines.append("  return True")
    co = ("\n".join(lines))
    ns = {}
    fn = exec(co, ns)
    return ns["check"]


def make_check(H):
    m, n = H.shape

    lines = ["def check(v):"]
    lines.append("  count = 0")
    for i in range(m):
        line = []
        for j in range(n):
            if H[i,j]==0:
                continue
            line.append("v[%d]"%(j,))
        line = "  if (%s)%%2:" % ("+".join(line),)
        lines.append(line)
        lines.append("    count += 1")
    lines.append("  return count")
    co = ("\n".join(lines))
    ns = {}
    fn = exec(co, ns)
    return ns["check"]


def perf():
    name = argv.load
    H = parse(open(name).read())
    print(H.shape)
    m, n = H.shape
    check = make_check(H)

    import numba
    check = numba.njit(check)

    print("warmup...")
    v = zeros2(n)
    for _ in range(10):
        assert check(v) == dot2(H, v).sum()
        i = randint(0,n-1)
        v[i] = (v[i]+1)%2

    print("go...")
    N=10000
    t0 = time()
    v = zeros2(n)
    count = 0
    for _ in range(N):
        count += check(v)
        i = randint(0,n-1)
        v[i] = (v[i]+1)%2
    print("time:", time() - t0)

    t0 = time()
    v = zeros2(n)
    count = 0
    for _ in range(N):
        count += dot2(H,v).sum()
        i = randint(0,n-1)
        v[i] = (v[i]+1)%2
    print("time:", time() - t0)


def find_lower_distance(H, L, d):
    # search all weight d vectors
    m, n = H.shape
    idxs = list(range(1, n))
    v = zeros2(n, 1)
    for ii in choose(idxs, d-1):
        v[:] = 0
        v[0,0] = 1
        for i in ii:
            v[i, 0] = 1
        if dot2(H, v).sum() == 0 and dot2(L, v).sum():
            return True # distance <= d
    return False # distance > d


def search_distance(H, L, d, homogeneous=False):
    import z3
    from z3 import Bool, And, Or, Xor, Not, Implies, Sum, If, Solver

    print("search_distance: d=%d"%d)

    m, n = H.shape
    k, n1 = L.shape
    assert n==n1

    solver = Solver()
    add = solver.add
    v = [Bool("v%d"%i) for i in range(n)]

    #add(Or(*v)) # non-zero
    term = Sum([If(vi,1,0) for vi in v]) == d
    #print(term)
    add(term)

    if homogeneous:
        add(v[0])

    # parity checks
    for i in range(m):
        terms = [v[j] for j in range(n) if H[i, j] == 1]
        assert len(terms)>1
        add(Not(reduce(Xor, terms)))

    # non-trivial logical 
    terms = []
    for i in range(k):
        items = [v[j] for j in range(n) if L[i, j] == 1]
        terms.append( reduce(Xor, items) )
    term = reduce(Or, terms)
    add(term)

    result = solver.check()
    print(result)
    if result != z3.sat:
        return

    model = solver.model()
    v = [model.evaluate(v[i]) for i in range(n)]
    v = [int(eval(str(vi))) for vi in v]
    print(v, sum(v))

    print(dot2(H, v))
    assert dot2(H, v).sum() == 0

    print(dot2(L, v))
    assert sum(v) == d, "z3 bug found"

    return True


def _sparse_distance(H, L, d, v, remain, max_check):

    m, n = H.shape
    k, _ = L.shape

    w0 = len(v)

    remain = set(remain)
    for i in list(remain):
        assert i not in v
        v.append(i)

        weight = H[:, v]
        weight = weight.sum(1) % 2
        weight = weight.sum()
        Lv = L[:, v]
        Lv = Lv.sum(1)
        Lv = Lv%2
        if weight==0 and Lv.sum():
            return True

        #if weight==0: we made a stabilizer

        if weight <= max_check and w0 < d:
            remain.remove(i)
            if _sparse_distance(H, L, d, v, remain, max_check): # <----- recurse
                return True
            remain.add(i)

        j = v.pop()
        assert j==i


def sparse_distance(H, L, d, homogeneous=True, verbose=False):
    # faster version of tree_distance
    if verbose:
        print("sparse_distance: d=%d"%d)
    m, n = H.shape
    k, n1 = L.shape
    assert n==n1

    assert homogeneous
    v = [0]
    remain = set(range(1, n))

    max_check = H[:, 0].sum()
    if argv.slow:
        max_check += 1 # hmmm... adding 1 here, see test()
    if verbose:
        print("max_check:", max_check)

    if _sparse_distance(H, L, d-1, v, remain, max_check):
        if verbose:
            print("found!")
            print(v)
        return True

    if verbose:
        print("not found")


def _tree_distance(H, L, d, v, remain, max_check):

    m, n = H.shape
    k, _ = L.shape

    w0 = v.sum()

    #print(" "*w0 + "_tree_distance", max_check)

    remain = set(remain)
    for i in list(remain):
        assert v[i] == 0
        v[i] = 1 # <------- set

        weight = dot2(H, v).sum()
        if weight==0 and dot2(L, v).sum():
            return True

        #if weight==0: we made a stabilizer

        if weight <= max_check and w0 < d:
            remain.remove(i)
            if _tree_distance(H, L, d, v, remain, max_check): # <----- recurse
                return True
            remain.add(i)

        v[i] = 0 # <------- reset


def tree_distance(H, L, d, homogeneous=False):
    # this works for codes where we never violate more checks
    # than the first bit, ie. codes with string-like logical operators,
    # eg. 2d topological codes.
    #print("tree_distance: d=%d"%d)
    m, n = H.shape
    k, n1 = L.shape
    assert n==n1

    check = lambda v : dot(H, v).sum() == 0

    v = zeros2(n)
    v[0] = 1
    remain = set(range(1, n))

    max_check = dot2(H, v).sum()
    if argv.slow:
        max_check += 1 # hmmm... adding 1 here, see test()
    #print("max_check:", max_check)

    if _tree_distance(H, L, d-1, v, remain, max_check):
        #print("found!")
        #print(v, v.sum())
        #print(dot2(H, v))
        assert dot2(H, v).sum() == 0
    
        #print(dot2(L, v))
        assert sum(v) == d

        return True

    #print("not found")


def dynamic_distance(H, L, verbose=False):
    if verbose:
        print("dynamic_distance")

    m, n = H.shape
    k, n1 = L.shape
    assert n==n1

    check = make_check(H)
    #check = lambda v : dot2(H, v).sum()
    if argv.fastnumba:
        import numba
        check = numba.njit(check)
    #else:

    #print(shortstr(H))
    support = [list(numpy.where(h)[0]) for h in H]
    #print(support)

    v = zeros2(n)
    v[0] = 1
    remain = set(range(1, n))

    max_check = dot2(H, v).sum()
    assert max_check == check(v)
    if argv.slow:
        max_check += 1 # hmmm... adding 1 here, see test()
    if verbose:
        print("max_check:", max_check)

    found = {(0,)}
    bdy = list(found)
    d = 2
    while bdy:
        #print()
        #print("="*79)
        if verbose:
            print("[d = %d]"%d, end="", flush=True)
            print("[bdy: %d]"%len(bdy), end="", flush=True)
        _bdy = []
        for path in bdy:
            #print("path:", path)
            v = zeros2(n)
            v[list(path)] = 1
            syndrome = dot2(H, v)
            #print("syndrome:", syndrome)
            idxs = numpy.where(syndrome)[0]
            #print("idxs:", idxs)
            nbd = set(reduce(add, [support[i] for i in idxs]))
            nbd = [i for i in nbd if v[i]==0]
            #print("nbd:", nbd)
            for i in nbd:
                v[i] = 1 # <----- set
                #Hv = dot2(H, v)
                #w = Hv.sum()
                w = check(v)
                if w == 0:
                    k = dot2(L, v).sum()
                    if k:
                        if verbose:
                            print()
                        return v.sum() # <------- return
                elif w <= max_check:
                    p = list(path)
                    p.append(i)
                    p.sort()
                    p = tuple(p)
                    if p not in found:
                        found.add(p)
                        _bdy.append(p)
                v[i] = 0 # <------ unset
        bdy = _bdy
        d += 1


def test():
    H = parse("""
......1.11.........1.....1........1..1.1................
.1................................11..111........1...1..
..........11...........1.........1.......1..1.....1.1...
.1..........11.......1.....1....1..............1.1......
.....1......1.....1.1...................1.1......1.....1
1...1.1........1....1....1.....1..........1.............
..1........1..11.1....1........1.........1..............
...11..1....1................1..1.........1..1..........
.......1.............1.....11.......1......1.11.........
.....1............1.....1....1.............111....1.....
.1.......1......1....1........1........1......1....1....
........1..1.....1.1...1....1.......1.................1.
1.................1.....111..........1..1............1..
....1.11.1....1..1..........1.1.........................
..11..........1...............1.1..............11..1....
.....1.......1.............1.....1..1.......1.........11
1..1......................1..1.1.........1......1.1.....
..........1.....1.........1...........1.........1..111..
........1.1.....1......11............1.....1..1.........
...............1...11.1...........11..................11
..1..........1........1..........1.1..1........1....1...
    """)
    L = parse("""
...................1...111.11.111.1.1..1.1...1..........
......................1..11....1..1...1.................
.......................11........11111.....1............
........................1........111.1......1...........
.........................1.11.1..1..1..11.1.11.1........
..........................11.1.11111...1111..1.....1....
...........................1.1.1.111.11.1111...11.......
............................1.1111111..1.11......1......
.............................1.......1111....11111......
...............................111....1..11.1...111.....
................................111..1..1....1111.11...1
.................................111..11....1...1.111...
..................................1.111.1.....11.1.1..1.
...................................11111......1......11.
    """)
    idxs = [0, 18, 29, 42]
    m, n = H.shape
    for w in range(1, 5):
        for jdxs in choose(idxs, w):
          v = zeros2(n)
          for i in jdxs:
            v[i] = 1
          print(jdxs, dot2(H, v).sum() )
    assert dot2(H, v).sum() == 0
    assert dot2(L, v).sum()
    print("OK")
    return


def update():
    names = os.listdir("codes")
    names.sort(key = lambda name:int(name.split("_")[1]))
    #print(names)

    idx = argv.get("idx", 10)

    for name in names:
        stem = name[:-4]
        n, k, _idx = [int(s) for s in stem.split("_")[1:4]]
        if idx is not None and _idx < idx:
            continue
        s = open("codes/"+name).read()
        print(name, end=" ", flush=True)
        H = parse(s)
        L = get_logops(H)
        d = argv.get("d", 2)
        while 1:
            print("[%d]"%d, end="", flush=True)
            result = sparse_distance(H, L, d)
            if result is not None:
                print(" d=%s"%d)
                break
            d += 2


def main_css():
    key = argv.get("key", (5,4))
    idx = 4
    while 1:
        css_geometry(key, idx)
        idx += 1


def css_geometry(key, idx):
    from bruhat.qcode import Geometry, get_adj
    geometry = Geometry(key, idx)

    G = geometry.G
    print("|G| = %4d, idx = %d" % (len(G), idx), end=" ", flush=True)

    faces = geometry.get_cosets([0,1,1])
    edges = geometry.get_cosets([1,0,1])
    verts = geometry.get_cosets([1,1,0])
    print("faces=%d, edges=%d, verts=%d"%(len(faces), len(edges), len(verts)), end=" ", flush=True)

    if argv.homology:
        Hz = get_adj(faces, edges)
        Hx = get_adj(verts, edges)
    elif argv.bicolour:
        assert 0
    else:
        assert 0

    _, n = Hz.shape
    if n < 10:
        print()
        return
    #print(shortstr(Hz), Hz.shape)
    #print()
    #print(shortstr(Hx), Hx.shape)

    from qumba.csscode import CSSCode
    code = CSSCode(Hx=Hx, Hz=Hz)
    print(code, end=" ", flush=True)

    dist = []
    for (H, L) in [(code.Hx, code.Lz), (code.Hz, code.Lx)]:
        d = dynamic_distance(H, L)
#        d = 1
#        while 1:
#            result = tree_distance(H, L, d, homogeneous=True)
#            if result:
#                break
#            if d == 20:
#                break
#    
#            d += 1
        dist.append(d)
    print("dist =", dist)


def main():

    H = parse("""
    ..11...1.111..11
    .....11..11.1111
    ....1..11..11111
    .1.1..1.1.11.1.1
    1.1.11.1.1..1.1.
    """)
    #L = get_logops(H)
    #print(shortstr(L))
    L = parse("""
    ...1..11....1...
    ......1.1.11....
    .......1.111....
    ........1..11.1.
    .........11.11..
    ............1111
    """)


    name = argv.next() or argv.load
    if name is not None:
        print(name)
        H = parse(open(name).read())
        L = get_logops(H)
    print(H.shape)
    m = rank(H)
    n = H.shape[1]
    k = len(L)
    print("[[%d, %d, --]]"%(n, k))
    assert n == 2*m + k

    if argv.dynamic_distance:
        result = dynamic_distance(H, L)
        print("d =", result)
        return

    if 1:
      for d in range(1, 4):
        result = find_lower_distance(H, L, d)
        if result:
            print("d =", d)
            return
        print("d >", d)

    d = argv.get("d", 4)
    max_d = argv.get("max_d", H.shape[1])
    homogeneous = argv.get("homogeneous", True)
    while 1:

        if argv.tree:
            result = tree_distance(H, L, d, homogeneous)
        elif argv.sparse:
            result = sparse_distance(H, L, d, homogeneous)
        else:
            result = search_distance(H, L, d, homogeneous)
        if result:
            break
        if d == max_d:
            break

        d += 1



if __name__ == "__main__":
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds.\n"%(time() - start_time))


