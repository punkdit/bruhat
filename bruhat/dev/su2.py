#!/usr/bin/env python3
"""
Finite subgroups of su2... FAIL
"""

from bruhat.argv import argv
from bruhat.element import FiniteField, PolynomialRing, CyclotomicField, Linear
from bruhat.vec import Space, Hom, Map
from bruhat.action import mulclose
from bruhat.util import all_primes

# Finite field notation in gap
# https://www.gap-system.org/Manuals/doc/ref/chap59.html
# [ [ Z(3)^0, Z(3)^0,   Z(3) ], [   Z(3), 0*Z(3),   Z(3) ], [ 0*Z(3),   Z(3), 0*Z(3) ] ]



def main():

    # ----------------------------------
    # binary tetrahedral group ... maybe?

    field = CyclotomicField(4)
    i = field.x

    print(i)
    M = Linear(2, field)
    gen = [
        M.get([[(-1+i)/2, (-1+i)/2], [(1+i)/2, (-1-i)/2]]),
        M.get([[0,i], [-i,0]])
    ]

    G = mulclose(gen)
    assert len(G)==48 # ???


    # ----------------------------------
    # binary octahedral group

    field = CyclotomicField(8)
    x = field.x

    a = x-x**3 #sqrt(2)
    i = x**2

    print((1+i)/a) # FAIL

    M = Linear(2, field)
    gen = [
        M.get([[(-1+i)/2,(-1+i)/2], [(1+i)/2,(-1-i)/2]]), 
        M.get([[(1+i)/a,0], [0,(1-i)/a]])
    ]
    print(gen[0])
    print(gen[1])


    return

    G = mulclose(gen)
    print(len(G))






if __name__ == "__main__":

    fn = argv.next()
    if fn is None:
        main()
    else:
        fn = eval(fn)
        fn()



