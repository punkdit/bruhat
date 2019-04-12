#!/usr/bin/env python3

"""
See https://golem.ph.utexas.edu/category/2019/03/entropy_mod_p.html

"""

from bruhat.action import Perm, Group
from bruhat.argv import argv


p = argv.get("p", 5)


if 0:
    items = [i for i in range(1, p**2) if i%p]
    print(items)
    
    perms = []
    for i in items:
        perm = dict((j, (j*i)%p**2) for j in items)
        perms.append(Perm(perm, items))
    
    G = Group(perms, items)
    
    print(len(G))
    print(G.is_cyclic())


def div(a, b):
    # a/b = c mod p
    # a = c*b mod p
    assert 0<=a<p
    assert 0<b<p
    for c in range(p):
        if a == (c*b)%p:
            return c
    assert 0


def h(x):
    total = 0
    for r in range(1, p):
        v = (x**r)%p
        v = div(v, r)
        total += v
    return total%p


for x in range(1, p):
    print("h(%d)=%d" % (x, h(x)), end=" ")
print()


