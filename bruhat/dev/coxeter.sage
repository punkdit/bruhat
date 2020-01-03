#!/usr/bin/env sage

from sage.monoids.hecke_monoid import HeckeMonoid

<<<<<<< HEAD
G = CoxeterGroup(["B", 3])
=======
#from bruhat.action import mulclose

def mulclose(gen, verbose=False, maxsize=None):
    els = set(gen)
    bdy = list(els)
    changed = True
    while bdy:
        _bdy = []
        for A in gen:
            for B in bdy:
              for C in [A*B, B*A]:
                if C not in els:
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
        bdy = _bdy
    return els


def mulclose_names(gen, names, verbose=False, maxsize=None):
    bdy = list(set(gen))
    assert len(names) == len(gen)
    names = dict((gen[i], (names[i],)) for i in range(len(gen)))
    changed = True
    while bdy:
        _bdy = []
        for A in gen:
            for B in bdy:
                C = A*B
                if C not in names:
                    #els.add(C)
                    names[C] = names[B] + names[A] # <-------- reverse the order
                    _bdy.append(C)
                    if maxsize and len(names)>=maxsize:
                        return list(names)
        bdy = _bdy
    return names


G = CoxeterGroup(["A", 3])
>>>>>>> a5dfec032fcba9c554e5271a8154e5d1167b8ce3

M = HeckeMonoid(G)
gen = M.monoid_generators()

a, b, c = gen
assert a*b*a == b*a*b
assert a*c == c*a
assert b*c*b == c*b*c

gen = list(gen[i+1] for i in range(len(gen)))
names = mulclose_names(gen, "GBR")
values = list(names.values())
values.sort()
print(values)
print(len(values))

print(len(G))

idem = []
for a in M:
    if a*a == a:
        idem.append(a)
#idem = list(reversed(idem))
print(len(idem))
#print(idem)

def getname(g):
    if g not in names:
        names[g] = "I"
    name = names[g]
    return ''.join(name)

for g in idem:
    print(getname(g))

#idem = idem[:4]

for g in idem:
  for h in idem:
    print("%s --> %s:" % (getname(g), getname(h))),
    count = 0
    n = 0
    #H = list(mulclose([g, h]))
    gen = [g, h]
    for x in M:
        if x == g*x*h:
            count += 1
            print getname(x),
            gen.append(x)
    #n = len(mulclose([g, h]))
    n = len(mulclose(gen))
    #print "%2d(%2d)"%(count, n),
    print
  #print




