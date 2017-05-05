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


