#!/usr/bin/env sage

#from __future__ import print_function

#import sys
#print(sys.argv)


if 0:
    K.<e8> = NumberField(x^4+1)
    
    r2 = e8+e8^7
    assert r2**2 == 2
    
    G = K.galois_group()
    assert len(G) == 4
    
    #for g in G:
    #    H = G.subgroup([g])
    #    result = H.fixed_field()
    #    J = result[0]
    #    print(J)
    #    e = J.gen()

else:
    K.<e16> = NumberField(x^8+1)
    assert K.is_galois()
    G = K.galois_group()
    assert len(G) == 8
    print(G)
    print(list(G))
    #print(dir(G))

for g in G:
    H = G.subgroup([g])
    J = H.fixed_field()
    print(J)

for g in G:
  for h in G:
    if h == g:
        continue
    H = G.subgroup([g, h])
    print("H:", len(H))
    J = H.fixed_field()
    print(J[0])
    

#e8 = e16**2
#r2 = e8+e8^7
#assert r2**2 == 2
#assert e16**16 == 1


#K.<e32> = NumberField(x^16+1)

#for item in K.subfields():
#    J = item[0]
#    print(J)


print("OK", 1)


