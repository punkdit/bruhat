#!/usr/bin/env python3
"""
Here we look at universal quantum gate sets over a finite field.
"""

from bruhat.argv import argv
from bruhat.element import FiniteField, PolynomialRing
from bruhat.vec import Space, Hom, Map
from bruhat.action import mulclose
from bruhat.util import all_primes


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
    #gen = [X, Z, S, H] # generates 192 elements
    gen = [X, Z, S, H, T]

    print("G := Group(")
    ms = []
    for g in gen:
        m = gapstr(g)
        ms.append(m)
    print(",\n".join(ms)+");;")

    if ring.p <= 41:
        G = mulclose(gen)
        print("|G| =", len(G))



def main():

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


if __name__ == "__main__":

    main()


