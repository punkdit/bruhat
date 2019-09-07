#!/usr/bin/env sage

#from __future__ import print_function

import sys
print(sys.argv)


K.<e8> = NumberField(x^4+1)

r2 = e8+e8^7
assert r2**2 == 2

G = K.galois_group()

assert len(G) == 4

for g in G:
    H = G.subgroup([g])
    J = H.fixed_field()
    e = J.gen()



print("OK", 1)


