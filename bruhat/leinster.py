#!/usr/bin/env python3

"""
See https://golem.ph.utexas.edu/category/2019/03/entropy_mod_p.html

"""

from bruhat.action import Perm, Group
from bruhat.util import factorial
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


def h2(x1, x2):
    x1 %= p
    x2 %= p
    total = 0
    for r1 in range(p):
      for r2 in range(p):
        if r1+r2 != p:
            continue
        total += - div( (x1**r1) * (x2**r2) % p, factorial(r1)*factorial(r2)%p )
    return total%p


for x in range(1, p):
    print("h(%d)=%d" % (x, h(x)), end=" ")
print()

for x in range(1, p):
    print("h2(%d)=%d" % (x, h2(x, 1-x)), end=" ")
print()

print("h(1,1)", h2(1, 1))


