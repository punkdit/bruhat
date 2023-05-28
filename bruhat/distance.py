#!/usr/bin/env python

"""
Use z3 to look for logical operators (for CSS codes),
thereby bounding/finding the distance.
"""


from time import time
start_time = time()
from random import choice
from functools import reduce
from operator import mul
from random import random, randint

import numpy

from bruhat.gset import Group, Perm, mulclose
from bruhat.argv import argv
from bruhat.solve import zeros2, dot2, span, shortstr, linear_independent, parse, enum2
from bruhat.util import choose
from bruhat.orthogonal import get_logops, row_reduce
from bruhat.oeqc import is_triorthogonal
from bruhat.todd_coxeter import Schreier



def make_check(H):
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
        assert check(v) == (dot2(H, v).sum() == 0)
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
    print("tree_distance: d=%d"%d)
    m, n = H.shape
    k, n1 = L.shape
    assert n==n1

    #print(shortstr(H))
    #print(H.sum(0))
    #print(H.sum(1))

    check = lambda v : dot(H, v).sum() == 0

    v = zeros2(n)
    v[0] = 1
    remain = set(range(1, n))

    max_check = dot2(H, v).sum()
    if argv.slow:
        max_check += 1 # hmmm... adding 1 here, see test()
    print("max_check:", max_check)

    if _tree_distance(H, L, d-1, v, remain, max_check):
        print("found!")
        print(v, v.sum())
        print(dot2(H, v))
        assert dot2(H, v).sum() == 0
    
        print(dot2(L, v))
        assert sum(v) == d

        return True

    print("not found")


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


    name = argv.load
    if name:
        H = parse(open(name).read())
        print(H.shape)
        L = get_logops(H)

    if 0:
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


