#!/usr/bin/env python3

"""
First attempt to do Galois theory with cyclotomic numbers...
"""

from bruhat.action import Perm, Group
from bruhat.element import CyclotomicRing
from bruhat.argv import argv

n = argv.get("n", 16)

items = list(range(1, n))

mul = set(items)
for i in items:
    for j in items:
        if (i*j)%n:
            continue
        if i in mul:
            mul.remove(i)

items = list(mul)
items.sort()

perms = []
for i in items:
    perm = dict((j, (j*i)%n) for j in items)
    perm = Perm(perm, items)
    perms.append(perm)

G = Group(perms, items)
print("   ", end="")
for jdx, qerm in enumerate(perms):
    j = items[jdx]
    print("%3d"%j, end="")
print()
print("  ", "---"*len(items))
for idx, perm in enumerate(perms):
    i = items[idx]
    print("%3d"%i, end="")
    for jdx, qerm in enumerate(perms):
        j = items[jdx]
        kdx = perms.index(perm*qerm)
        k = items[kdx]
        print("%3d"%k, end="")
    print()

print()

#for H in G.subgroups():
#    print(H)

ring = CyclotomicRing(n)
x = ring.x

r2 = x**2 - x**6
assert r2**2 == 2
            
