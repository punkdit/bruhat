#!/usr/bin/env python3
"""
Here we look at universal quantum gate sets over a finite field.

Derived from gates.py 
"""

from random import shuffle

import numpy

from bruhat.argv import argv
from bruhat.util import all_primes

print(all_primes(100))

# Finite field notation in gap
# https://www.gap-system.org/Manuals/doc/ref/chap59.html
# [ [ Z(3)^0, Z(3)^0,   Z(3) ], [   Z(3), 0*Z(3),   Z(3) ], [ 0*Z(3),   Z(3), 0*Z(3) ] ]


def test():

    p = argv.get("p", 17)
    print("p =", p)
    
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

    # look for a square root of -1
    has_i4 = [i for i in items if i*i%p == p-1]
    
    if not has_i4:
        assert 0

    i4 = has_i4[0]
    assert p-i4 == has_i4[1]

    # look for a square root of 2
    has_r2 = [i for i in items if i*i%p == 2]
    if not has_r2:
        assert 0
    r2 = has_r2[0]
    assert p-r2 == has_r2[1]

    # look for a square root of i4
    has_i8 = [i for i in items if i*i%p == i4]
    if not has_i8:
        assert 0
    i8 = has_i8[0]
    assert p-i8 == has_i8[1]

    # look for a square root of i8... ?
    has_i16 = [i for i in items if i*i%p == i8]
    if has_i16:
        i16 = has_i16[0]
        assert p-i16 == has_i16[1]
        #print(i16**15%p)
        print(i4, r2, i8, i16)


    def array(items):
        A = numpy.array(items, dtype=int) % p
        return A

    def mul(A, B):
        AB = numpy.dot(A, B) % p
        return AB

    def inv(A):
        a, b = A[0]
        c, d = A[1]
        r = finv[(a*d - b*c)%p]
        B = array([[d, -b], [-c, a]])*r
        return B

    def equ(A, B):
        A = (A-B)%p
        return numpy.abs(A).sum() == 0

    def fromstring(s):
        A = numpy.fromstring(s, int)
        A.shape = 2,2
        return A

    def mulclose(gen, verbose=False, maxsize=None):
        els = set(g.tostring() for g in gen)
        bdy = list(gen)
        changed = True
        while bdy:
            #if verbose:
            #    print "mulclose:", len(els)
            _bdy = []
            for A in gen:
                for B in bdy:
                    C = mul(A, B)
                    s = C.tostring()
                    if s not in els:
                        els.add(s)
                        _bdy.append(C)
            bdy = _bdy
        return [fromstring(s) for s in els]

    I = array([[1, 0], [0, 1]])
    X = array([[0, 1], [1, 0]])
    Z = array([[1, 0], [0, -1]])
    ir2 = finv[r2]
    H = (X+Z)*ir2%p

    assert equ(mul(H, H), I)

    S = array([[1, 0], [0, i4]])
    assert equ(mul(S, S), Z)

    T = array([[1, 0], [0, i8]])
    assert equ(mul(T, T), S)

    C1 = mulclose([X, Z])
    assert len(C1) == 8

    # C1 is Pauli group + phases
    P = mul(fgen, I) # phase
    C1 = mulclose([X, Z, P])  # add phases
    assert len(C1)/(p-1) == 4, len(C1)
    C1_lookup = set(g.tostring() for g in C1)

    #gen = [X, Z, S, H]
    #C2 = mulclose(gen)
    #assert len(C2) == 192

    G = []
    for a in range(p):
     for b in range(p):
      for c in range(p):
       for d in range(p):
        if (a*d - b*c)%p:
            G.append(array([[a, b], [c, d]]))
    G_lookup = set(g.tostring() for g in G)
    print("|GL(%d, 2)|=%d" % (p, len(G)))

    gen = [X, Z, S, H, P]
    for g in gen:
        assert numpy.allclose(g%p, g)

    # Clifford group + phases
    C2 = mulclose(gen)
    for g in C2:
        assert numpy.allclose(g%p, g)
    assert len(C2)/(p-1) == 24, len(C2)
    C2_lookup = set(g.tostring() for g in C2)
    print("|C2| =", len(C2))
    print("|C2|/(p-1) =", len(C2)/(p-1))

    for g in C2:
        assert g.tostring() in G_lookup

    C3 = []
    for g3 in G:
        for g in C1:
            if mul(mul(g3, g), inv(g3)).tostring() not in C2_lookup:
                break
        else:
            C3.append(g3)
    print("|C3| =", len(C3))
    print("|C3|/(p-1) =", len(C3)/(p-1))
    C3_lookup = set(g.tostring() for g in C3)

    shuffle(C3)

    if 0:
        a = ir2*i16**15%p
        b = ir2*i16**3%p
        H = array([[a, b], [b, a]])
        print(H)
        print(H.tostring() in C3_lookup) # True
        return
    

    def src(a):
        return set([b.tostring() for b in C3 if mul(a, b).tostring() in C3_lookup])

    def tgt(a):
        return set([b.tostring() for b in C3 if mul(b, a).tostring() in C3_lookup])

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
            
    print()

    for ob in obs:
        if len(ob) == len(C2):
            assert set(ob) == set(C2_lookup)


if __name__ == "__main__":

    fn = argv.next()
    if fn is None:
        test()
    else:
        fn = eval(fn)
        fn()



