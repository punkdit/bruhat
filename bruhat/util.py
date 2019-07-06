#!/usr/bin/env python

import sys


def all_subsets(n):

    if n==0:
        yield []
        return

    if n==1:
        yield []
        yield [0]
        return

    for subset in all_subsets(n-1):
        yield subset
        yield subset + [n-1] # sorted !!

assert len(list(all_subsets(5))) == 2**5




def factorial(n):
    r = 1
    for i in range(1, n+1):
        r *= i
    return r

assert factorial(0) == 1
assert factorial(1) == 1
assert factorial(2) == 2
assert factorial(3) == 2*3
assert factorial(4) == 2*3*4


def choose(items, n):
    if n > len(items):
        return
    if n == 0:
        yield ()
        return
    if n == 1:
        for item in items:
            yield (item,)
        return
    for i, item in enumerate(items):
        for rest in choose(items[i+1:], n-1):
            yield (item,)+rest

assert len(list(choose(range(4), 1))) == 4
assert len(list(choose(range(4), 2))) == 6
assert len(list(choose(range(4), 3))) == 4



def allperms(items):
    items = tuple(items)
    if len(items)<=1:
        yield items
        return
    n = len(items)
    for i in range(n):
        for rest in allperms(items[:i] + items[i+1:]):
            yield (items[i],) + rest

assert list(allperms("abc")) == [
    ('a', 'b', 'c'), 
    ('a', 'c', 'b'), 
    ('b', 'a', 'c'), 
    ('b', 'c', 'a'), 
    ('c', 'a', 'b'), 
    ('c', 'b', 'a')]

all_perms = allperms

def allders(items):
    "all derangements of items"
    # hack this
    for perm in allperms(items):
        for a, b in zip(perm, items):
            if a==b:
                break
        else:
            yield perm

assert list(allders("abc")) == [('b', 'c', 'a'), ('c', 'a', 'b')]


all_ders = allders



def allsignedperms(items):
    items = tuple(items)
    if len(items)<=1:
        yield +1, items
        return
    n = len(items)
    sign = 1
    for i in range(n):
        for _sign, rest in allsignedperms(items[:i] + items[i+1:]):
            yield sign*_sign, (items[i],) + rest
        sign *= -1

assert list(allsignedperms("ab")) == [(1, ('a', 'b')), (-1, ('b', 'a'))]


def write(s):
    sys.stdout.write(str(s)+' ')
    sys.stdout.flush()




def cross(itemss):
    if len(itemss)==0:
        yield ()
    else:
        for head in itemss[0]:
            for tail in cross(itemss[1:]):
                yield (head,)+tail


def all_parts2(items):
    items = list(items)
    n = len(items)
    if n==0:
        yield items, items
        return

    if n==1:
        yield [], items
        yield items, []
        return

    if n==2:
        yield [], items
        yield [items[0]], [items[1]]
        yield [items[1]], [items[0]]
        yield items, []
        return

    bits = [(0, 1)]*n
    for idxs in cross(bits):
        left = [items[i] for i in range(n) if idxs[i]==0] 
        right = [items[i] for i in range(n) if idxs[i]==1] 
        yield left, right

assert list(all_parts2([0, 1, 2])) == [
    ([0,1,2], []), ([0,1], [2]), ([0,2], [1]), ([0], [1, 2]), 
    ([1, 2], [0]), ([1], [0, 2]), ([2], [0, 1]), ([], [0, 1, 2])]


def uniqtuples(items, n):
    if n==0:
        yield ()
        return # <-- return
    assert n>0
    if n > len(items):
        return # <-- return
    if len(items)==1:
        assert n==1
        yield (items[0],)
        return # <-- return

    m = len(items)
    for i in range(m):
        item = items[i]
        for tail in uniqtuples(items[:i] + items[i+1:], n-1):
            yield (item,)+tail


def alltuples(items, n):
    if n==0:
        yield ()
        return # <-- return
    assert n>0
    if n > len(items):
        return # <-- return
    if len(items)==1:
        assert n==1
        yield (items[0],)
        return # <-- return

    m = len(items)
    for i in range(m):
        item = items[i]
        for tail in uniqtuples(items, n-1):
            yield (item,)+tail


def partitions(n, maxn=None):
    if maxn is None:
        maxn = n
    if n==0:
        yield ()
    else:
        for i in range(1, min(n, maxn)+1):
            for rest in partitions(n-i, i):
                yield (i,) + rest


assert list(partitions(0)) == [()]
assert list(partitions(1)) == [(1,)]
assert list(partitions(2)) == [(1, 1), (2,)]
assert list(partitions(3)) == [(1, 1, 1), (2, 1), (3,)]
assert list(partitions(4)) == [(1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,)]


def all_primes(n, ps=None):
    "list of primes < n"

    items = [0]*n
    p = 2 

    while p**2 < n:
        i = 2 
        while p*i < n:
            items[p*i] = 1 
            i += 1

        p += 1
        while p < n and items[p]:
            p += 1

    ps = [i for i in range(2, n) if items[i]==0]
    return ps



def all_binary_trees(els):
    if len(els)==0:
        return

    els = list(els)
    if len(els)==1:
        #yield (els[0],) # a leaf
        yield els[0] # a leaf
        return

    els.sort()
    n = len(els)
    # iterate over size of left sub-tree
    for left in range(1, n):
        right = n-left # size of right sub-tree
        rtrees = list(all_binary_trees(els[left:]))
        for ltree in all_binary_trees(els[:left]):
            for rtree in rtrees:
                yield (ltree, rtree)

# https://en.wikipedia.org/wiki/Catalan_number
# https://oeis.org/A000108
# https://en.wikipedia.org/wiki/Associahedron
assert list(all_binary_trees([1])) == [1]
assert list(all_binary_trees([1,2])) == [(1, 2),]
assert list(all_binary_trees([1,2,3])) == [(1, (2, 3)), ((1, 2), 3)]
assert list(all_binary_trees([1,2,3,4])) == [
    (1, (2, (3, 4))), (1, ((2, 3), 4)), ((1, 2), (3, 4)), ((1, (2, 3)), 4), (((1, 2), 3), 4)]
assert len(list(all_binary_trees([1,2,3,4,5]))) == 14




def set_seed(i=0):
    import numpy
    from random import seed
    seed(i)
    numpy.random.seed(i)

