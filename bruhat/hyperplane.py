#!/usr/bin/env python

from sage.all_cmdline import *


def make_parent(base_ring, dimension, names=None):
    if names is None:
        names = tuple('t'+str(i) for i in range(dimension))
    else:
        names = tuple(map(str, names))
        if len(names) != dimension:
            raise ValueError('number of variable names does not match dimension')
    return HyperplaneArrangements(base_ring, names=names)




def Coxeter(data, K=QQ, names=None):
    from sage.combinat.root_system.cartan_type import CartanType
    from sage.combinat.root_system.root_system import RootSystem
    from sage.combinat.root_system.weyl_group import WeylGroup
    
    if data in NN:
        cartan_type = CartanType(["A", data - 1])
    else:
        cartan_type = CartanType(data)
    if not cartan_type.is_crystallographic():
        raise NotImplementedError(
        "Coxeter arrangements are not implemented for non crystallographic Cartan types")
    W = WeylGroup(cartan_type)
    Ra = RootSystem(cartan_type).ambient_space()
    PR = Ra.positive_roots()
    d = Ra.dimension()
    H = make_parent(K, d, names)
    x = H.gens()
    hyperplanes = []
    
    for a in PR:
        hyperplanes.append(sum(a[j] * x[j] for j in range(d)))
    A = H(*hyperplanes)
    x = polygen(QQ, 'x')
    charpoly = prod(x - d + 1 for d in W.degrees())
    A.characteristic_polynomial.set_cache(charpoly)
    return A


h = Coxeter("B4")

print(h)

print(' '.join(dir(h)))

""" B3
t0 + + + 0
+ t1 + + 0
+ + t2 + 0
t0 - t1 + + 0
t0 + t1 + + 0
t0 + - t2 + 0
t0 + + t2 + 0
+ t1 - t2 + 0
+ t1 + t2 + 0
"""



