#!/usr/bin/env sage


from sage.monoids.hecke_monoid import HeckeMonoid

G = CoxeterGroup(["B", 3])

M = HeckeMonoid(G)
gen = M.monoid_generators()

print(len(G))

idem = []
for a in M:
    if a*a == a:
        idem.append(a)
#idem = list(reversed(idem))
print(len(idem))
#print(idem)

for g in idem:
  for h in idem:
    count = 0
    for x in M:
        if x == g*x*h:
            count += 1
    print "%2d"%count,
  print




