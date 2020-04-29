#!/usr/bin/env python3

# From: https://math.berkeley.edu/~kmill/notes/todd_coxeter.html


# Example of Todd-Coxeter to compute S_3 from relations

idents = []
neighbors = []
to_visit = 0

ngens = 2
rels = [
    (1, 0), # a^-1a
    (3, 2), # b^-1b
    (0, 0, 0), #a^3
    (2, 2), # b^2
    (0, 2, 0, 2) # abab
]
hgens = [
    (2,), # b
]

def find(c):
    c2 = idents[c]
    if c == c2:
        return c
    else:
        c2 = find(c2)
        idents[c] = c2
        return c2

def new():
    c = len(idents)
    idents.append(c)
    neighbors.append((2*ngens)*[None])
    return c

def unify(c1, c2):
    c1 = find(c1)
    c2 = find(c2)
    if c1 == c2:
        return
    c1, c2 = min(c1, c2), max(c1, c2)
    idents[c2] = c1
    for d in range(2*ngens):
        n1 = neighbors[c1][d]
        n2 = neighbors[c2][d]
        if n1 == None:
            neighbors[c1][d] = n2
        elif n2 != None:
            unify(n1, n2)

def follow(c, d):
    c = find(c)
    ns = neighbors[c]
    if ns[d] == None:
        ns[d] = new()
    return find(ns[d])

def followp(c, ds):
    c = find(c)
    for d in reversed(ds):
        c = follow(c, d)
    return c

start = new()

for hgen in hgens:
    unify(followp(start, hgen), start)

while to_visit < len(idents):
    c = find(to_visit)
    if c == to_visit:
        for rel in rels:
            unify(followp(c, rel), c)
    to_visit += 1

print("done")

cosets = [c for i, c in enumerate(idents) if i == c]

perms = [[cosets.index(follow(c, 2*d)) for i, c in enumerate(cosets)]
         for d in range(ngens)]

def cycle(perm):
    parts = []
    for i in range(len(perm)):
        part = [str(i+1)]
        k = perm[i]
        while k != i:
            if k < i: break
            part.append(str(k+1))
            k = perm[k]
        else:
            parts.append(" ".join(part))
    return "("+")(".join(parts)+")"

for d in range(ngens):
    print("g%d ="%d, cycle(perms[d]))


