#!/usr/bin/env python3

"""
Find power series of rational functions
"""

from bruhat.util import factorial
from bruhat.element import Q
from bruhat.poly import Poly


def choose(m, n):
    return factorial(m) // (factorial(n) * factorial(m-n))




def leibniz(base, p, q, x):
    # f(x) == p(x) / q(x)
    # q(x)*f(x) == p(x)
    p = Poly.promote(p, base)
    q = Poly.promote(q, base)

    print("leibniz(%s, %s)" % (p, q))
    
    #for n in range(

    zero, one = base.zero, base.one

    fs = [] # coeff's

    p0, q0 = p(x=zero), q(x=zero)
    f0 = p0/q0
    fs.append(f0)

    print(fs)

    for n in range(1, 5):
        rhs = p.diff(x, n)(x=zero)
        for r in range(n):
            q0 = q.diff(x, n-r)(x=zero)
            rhs -= choose(n, r) * fs[r] * q0
        f0 = rhs / q(x=zero)
        print(f0)
        fs.append(f0)






def test():

    ring = Q
    zero = Poly({}, ring)
    one = Poly({():1}, ring)
    x = Poly("x", ring)
    y = Poly("y", ring)
    z = Poly("z", ring)
    xx = Poly({(("x", 2),) : 1}, ring)
    xy = Poly({(("x", 1), ("y", 1)) : 1}, ring)


    leibniz(ring, one, (1-x)**2, x)



if __name__ == "__main__":

    test()

