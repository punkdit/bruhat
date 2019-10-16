#!/usr/bin/env python3

"""
Verifying that the normalizer of a subgroup is
the automorphism group of the corresponding G-set.
"""

from bruhat.action import Group, Perm
from bruhat.argv import argv


G = Group.symmetric(4)

print(len(G))

Hs = G.subgroups()

print(len(Hs))


def setpromote(items):
    if isinstance(items, Perm):
        items = {items}
    else:
        items = set(items)
    for item in items:
        assert isinstance(item, Perm)
    return items


def setmul(*itemss):
    items = itemss[0]
    items = setpromote(items)

    idx = 1
    while idx < len(itemss):
        jtems = itemss[idx]
        jtems = setpromote(jtems)
        items = set(a*b for a in items for b in jtems)
        idx += 1

    return items


for H in Hs:

    for n in G:

        # nH = HnH

        lhs = setmul(n, H)
        rhs = setmul(H, n, H)

        i = int(lhs == rhs)

        rhs = setmul(H, n)
        j = int(lhs == rhs)

        #print("%d==%d"%(i, j), end=" ")
        assert i==j

        if not j:
            continue

        

    print()



