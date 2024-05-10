#!/usr/bin/env python

"""

"""

import string

from sage.all_cmdline import GF, NumberField
from sage import all_cmdline 

from bruhat.argv import argv
from bruhat.smap import SMap



def slow_count_points(p, l):
    q = p**l

    F = GF(q)

    e = F.gen()
    if l==1:
        # e is additive gen
        items = [i*e for i in range(q)]
    else:
        # e is multiplicative gen
        items = [e**i for i in range(1, q)] + [0*e]
    #print(p, l, len(set(items)), len(items))
    assert len(set(items)) == len(items) == q, [e, e**2]

    f = lambda x,y : y**2 == (x-1)*x*(x+1)
    sols = [f(x,y) for x in items for y in items]
    n = sols.count(True)

    lhs = lambda y : y**2
    rhs = lambda x : (x-1)*x*(x+1)

    lsend = {a : 0 for a in items}
    for y in items:
        lsend[lhs(y)] += 1

    rsend = {a : 0 for a in items}
    for x in items:
        rsend[rhs(x)] += 1

    count = 0
    for a in items:
        count += lsend[a]*rsend[a]
    assert count == n

    return n+1 # add projective point


def count_points(p, l):
    q = p**l

    F = GF(q)

    e = F.gen()
    if l==1:
        # e is additive gen
        items = [i*e for i in range(q)]
    else:
        # e is multiplicative gen
        items = [e**i for i in range(1, q)] + [0*e]
    #print(p, l, len(set(items)), len(items))
    #assert len(set(items)) == len(items) == q, [e, e**2]
    assert len(items) == q, [e, e**2]

    #f = lambda x,y : y**2 == (x-1)*x*(x+1)
    #sols = [f(x,y) for x in items for y in items]
    #n = sols.count(True)
    print('_', end='', flush=True)

    lhs = lambda y : y**2
    rhs = lambda x : (x-1)*x*(x+1)

    lsend = {a : 0 for a in items}
    rsend = {a : 0 for a in items}
    for a in items:
        lsend[lhs(a)] += 1
        rsend[rhs(a)] += 1
        #lsend[a**2] += 1
        #rsend[(a-1)*a*(a+1)] += 1
    print('_', end='', flush=True)

    count = 0
    for a in items:
        count += lsend[a]*rsend[a]
    #assert count == n

    return count+1 # add projective point


def spec():

    #F = GF(5)
    #print(F, F.gen())
    #return

    p = argv.get("p")
    l = argv.get("l", 1)
    if p is None:
        ps = [2, 3, 5]
        ls = list(range(1,9))
    else:
        ps = [p]
        ls = [l]
    for p in ps:
      for l in ls:
        print(p, l, p**l, end=" ", flush=True)
        n = count_points(p, l)
        print(n, n-(p**l+1))


def main():
    n = argv.get("n")

    K0 = NumberField(x**2 + 1, "x")
    print(K0)

    K1 = K0.extension(x**2+x+1, "w")
    print(K1)

    print(K1.absolute_field("a").galois_group())



if __name__ == "__main__":

    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()


