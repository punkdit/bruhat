#!/usr/bin/env python3
"""
Here we look at universal quantum gate sets over a finite field.

See also: clifford.py
"""

from random import shuffle

from bruhat.argv import argv
from bruhat.element import FiniteField, PolynomialRing
from bruhat.vec import Space, Hom, Map
from bruhat.action import mulclose
from bruhat.util import all_primes

# Finite field notation in gap
# https://www.gap-system.org/Manuals/doc/ref/chap59.html
# [ [ Z(3)^0, Z(3)^0,   Z(3) ], [   Z(3), 0*Z(3),   Z(3) ], [ 0*Z(3),   Z(3), 0*Z(3) ] ]


def get_gen(p):
    for i in range(1, p):
        j = (i*i)%p
        count = 2
        while j != i:
            count += 1
            j = (i*j) % p
        assert count <= p, (count,)
        if count == p:
            return i
    assert 0

def get_order(p, x):
    j = i = get_gen(p)
    order = 1
    while j != x:
        j = i*j
        order += 1
    return order


def element_str(p, x):
    if x==0:
        return "0*Z(%s)"%p
    return "Z(%s)^%s" % (p, get_order(p, x))


def gapstr(A):
    p = A.ring.p
    s = A.str(element_str=(lambda x,p=p:element_str(p, x)), sep=', ')
    s = s.replace("\n", " ")
    return s


def build_gates(ring, i4, root2, i8):
    qubit = Space(2, ring)
    hom = Hom(qubit, qubit)

    I = Map.from_array([[1, 0], [0, 1]], hom)
    X = Map.from_array([[0, 1], [1, 0]], hom)
    Z = Map.from_array([[1, 0], [0, -1]], hom)
    H = (1/root2)*(X+Z)

    assert H*H == I

    S = Map.from_array([[1, 0], [0, i4]], hom)
    assert S*S == Z

    T = Map.from_array([[1, 0], [0, i8]], hom)
    assert S*S == Z

    #gen = [X, Z, S] # generates 32 elements
    gen = [X, Z, S, H] # generates 192 elements
    #gen = [X, Z, S, H, T]

    if argv.gap:
        print("G := Group(")
        ms = []
        for g in gen:
            m = gapstr(g)
            ms.append(m)
        print(",\n".join(ms)+");;")

    if ring.p <= 41:
        G = mulclose(gen)
        print("|G| =", len(G))




def test():

    p = argv.get("p")
    ps = [p] if p else all_primes(200)

    for p in ps:

        if p==2:
            continue
    
        field = FiniteField(p)
        #ring = PolynomialRing(field)
        #x = ring.x
    
        items = [field.promote(i) for i in range(p)]
        has_imag = [i for i in items if i*i == -1]
        
        #print(p, has_imag)
        if not has_imag:
            continue

        i4 = has_imag[0]
        assert -i4 == has_imag[1]

        has_root2 = [i for i in items if i*i == 2]
        if not has_root2:
            continue
        r2 = has_root2[0]
        assert -r2 == has_root2[1]

        has_i8 = [i for i in items if i*i == i4]
        if not has_i8:
            continue
        i8 = has_i8[0]
        assert -i8 == has_i8[1]

        print("p =", p)
        #print(i4, r2, i8)
        build_gates(field, i4, r2, i8)


def clifford():

    p = argv.get("p", 17)
    
    field = FiniteField(p)
    #ring = PolynomialRing(field)
    #x = ring.x

    for fgen in range(1, p):
        items = set()
        j = fgen
        while j not in items:
            items.add(j)
            j = (j*fgen)%p
        if len(items) == p-1:
            break
    else:
        assert 0
    finv = {}
    for i in range(1, p):
      for j in range(1, p):
        if (i*j)%p == 1:
          finv[i] = j
          finv[j] = i
    assert len(finv) == p-1

    #print("fgen:", fgen)

    items = [field.promote(i) for i in range(p)]
    has_imag = [i for i in items if i*i == -1]
    
    #print(p, has_imag)
    if not has_imag:
        assert 0

    i4 = has_imag[0]
    assert -i4 == has_imag[1]

    has_root2 = [i for i in items if i*i == 2]
    if not has_root2:
        assert 0
    r2 = has_root2[0]
    assert -r2 == has_root2[1]

    has_i8 = [i for i in items if i*i == i4]
    if not has_i8:
        assert 0
    i8 = has_i8[0]
    assert -i8 == has_i8[1]

    print("p =", p)
    #print(i4, r2, i8)

    qubit = Space(2, field)
    hom = Hom(qubit, qubit)

    I = Map.from_array([[1, 0], [0, 1]], hom)
    X = Map.from_array([[0, 1], [1, 0]], hom)
    Z = Map.from_array([[1, 0], [0, -1]], hom)
    H = (1/r2)*(X+Z)

    assert H*H == I

    S = Map.from_array([[1, 0], [0, i4]], hom)
    assert S*S == Z

    T = Map.from_array([[1, 0], [0, i8]], hom)
    assert S*S == Z

    C1 = mulclose([X, Z])
    assert len(C1) == 8

    # C1 is Pauli group + phases
    P = fgen*I # phase
    C1 = mulclose([X, Z, P])  # add phases
    assert len(C1) == 64, len(C1)
    C1_lookup = set(C1)

    #gen = [X, Z, S, H]
    #C2 = mulclose(gen)
    #assert len(C2) == 192

    G = []
    for a in range(p):
     for b in range(p):
      for c in range(p):
       for d in range(p):
        if (a*d - b*c)%p:
            G.append(Map.from_array([[a, b], [c, d]], hom))
    G_lookup = set(G)
    print("|GL(%d, 2)|=%d" % (p, len(G)))

    gen = [X, Z, S, H, P]

    # Clifford group + phases
    C2 = mulclose(gen)
    assert len(C2) == 384
    C2_lookup = set(C2)
    print("|C2| =", len(C2))

    for g in C2:
        assert g in G_lookup

    inv = {I:I}
    for a in G:
      if a in inv:
        continue
      for b in G:
        ab = a*b
        ci = inv.get(ab)
        if ci is None:
            continue
        inv[a] = b*ci
        inv[b] = ci*a
      print(len(inv), end=" ", flush=True)
    print()

    if 0:
        for g2 in C2:
            for g in C1:
                assert g2 * g * inv[g2] in C1
    
        C2 = []
        for g2 in G:
            for g in C1:
                if g2*g*inv[g2] not in C1_lookup:
                    break
            else:
                C2.append(g2)
        assert len(C2) == 384 # same as above
        C2_lookup = set(C2)

    C3 = []
    for g3 in G:
        for g in C1:
            if g3*g*inv[g3] not in C2_lookup:
                break
        else:
            C3.append(g3)
    print("|C3| =", len(C3))
    C3_lookup = set(C3)

    shuffle(C3)

#    count = 0
#    for a in C3:
#      for b in C3:
#        if a*b in C3_lookup:
#            count += 1
#    print(count)

    def src(a):
        return set([b for b in C3 if a*b in C3_lookup])

    def tgt(a):
        return set([b for b in C3 if b*a in C3_lookup])

    if 0:
        items = iter(C3)
        a = items.__next__()
    
        src_a = src(a)
        print("|src_a| =", len(src_a))

    srcs = []
    for b in C3:
        src_b = src(b)
        if src_b not in srcs:
            print("|src_b| = ", len(src_b))
            srcs.append(src_b)
        if len(srcs)==4: # there is only 4 of these to be found
            break

    obs = list(srcs)
    tgts = []
    for b in C3:
        tgt_b = tgt(b)
        if tgt_b not in obs:
            obs.append(tgt_b)
        if tgt_b not in tgts:
            print("|tgt_b| = ", len(tgt_b))
            tgts.append(tgt_b)
        if len(tgts)==4: # there is only 4 of these to be found
            break

    done = False
    while not done:
        done = True
        print("obs:", len(obs))
    
        obs1 = list(obs)
        for s in obs:
          for t in obs:
            st = s.intersection(t)
            a = ' '
            if st not in obs1:
                a = '*'
                obs1.append(st)
                done = False
            print("%4d"%(len(st)), end=a)
          print()
    
        obs = obs1

    print("obs:", len(obs))
    for ob in obs:
        print(len(ob), end=" ")



if __name__ == "__main__":

    fn = argv.next()
    if fn is None:
        test()
    else:
        fn = eval(fn)
        fn()



