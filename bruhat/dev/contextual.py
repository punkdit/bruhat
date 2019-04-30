#!/usr/bin/env python3

from bruhat.util import cross

"""

The Mermin-Peres (MP) square

    XI    IX    XX   -->  I

    IZ    ZI    ZZ   -->  I

    XZ    ZX    YY   -->  I
    |     |     |     
    v     v     v     
    I     I     -I    

Theorem: It is not possible to assign values phi in {-1, 1} to 
every member of the 2-qubit pauli group such that
    phi(XI)*phi(IX)*phi(XX) = 1
    phi(IZ)*phi(ZI)*phi(ZZ) = 1
    phi(XZ)*phi(ZX)*phi(YY) = 1
    phi(XI)*phi(IZ)*phi(XZ) = 1
    phi(IX)*phi(ZI)*phi(ZX) = 1
    phi(XX)*phi(ZZ)*phi(YY) = -1
"""


items = "XI IX XX IZ ZI ZZ XZ ZX YY".split()
XI, IX, XX, IZ, ZI, ZZ, XZ, ZX, YY = items

n = len(items)

for vals in cross([(-1, 1)]*n):
    phi = {}
    for idx, item in enumerate(items):
        phi[items[idx]] = vals[idx]

    if (
        phi[XI]*phi[IX]*phi[XX] == 1 and
        phi[IZ]*phi[ZI]*phi[ZZ] == 1 and
        phi[XZ]*phi[ZX]*phi[YY] == 1 and
        phi[XI]*phi[IZ]*phi[XZ] == 1 and
        phi[IX]*phi[ZI]*phi[ZX] == 1 and
        phi[XX]*phi[ZZ]*phi[YY] == -1 ):

        for idx, item in enumerate(items):
            print("phi(%s)=%s "%(item, phi[item]), end=" ")
            if (idx+1)%3==0:
                print()
        break



