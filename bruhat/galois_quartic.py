#!/usr/bin/env python

"""
Galois theory of polynomial equations ...

"""

import string

from sage import all_cmdline as sage

from bruhat.argv import argv
from bruhat.smap import SMap


def render():

    from huygens.namespace import RGB, Canvas, path
    
    yellow = RGB(1.0, 0.8, 0.2)
    orange = RGB(252/255, 132/255, 3/255) * (1./0.8)


    
def main_galois():

    K = sage.PolynomialRing(sage.QQ, "x")
    x, = K.gens()

    #p = x**4 - x + 1 # Galois group 4T5 (S4) with order 24
    p = x**3 - x + 1 # Galois group 3T2 (S3) with order 6 
    #p = x**6 - 6*x**4 + 9*x**2 + 23 # Galois group 6T2 ([3]2) with order 6 

    E = sage.NumberField(p, "x")
    print(E)
    G = E.galois_group()
    print(G)
    print()
    Hs = G.subgroups()
    for H in Hs:
        print("|H| =", len(H))
        print(H)
        F, hom = (H.fixed_field())
        print(F)
        print(hom)
        print()

main = main_galois

def main_cubic():

    B = sage.CyclotomicField(3)
    w, = B.gens()

    K = sage.PolynomialRing(B, list("rstx"))
    vs = K.gens()
    r, s, t, x = vs

    y = (x-r)*(x-s)*(x-t)

    #print(y)
    #for item in y:
    #    print(item)

    a = 1
    b = -(r+s+t)
    c = r*s + r*t + s*t
    d = -r*s*t
    assert y == a*x**3 + b*x**2 + c*x + d

    p = r**2*s + w*s**2*t + w**2*t**2*r
    p = p**3
    for coeff,monic in p:
        print("%14s"%coeff, monic)

    """
       1 r^6*s^3
       1 s^6*t^3
       1 r^3*t^6
     3*w r^4*s^4*t
     3*w r^4*s*t^4
     3*w r*s^4*t^4
  3*w**2 r^5*s^2*t^2 +
  3*w**2 r^2*s^5*t^2 +
  3*w**2 r^2*s^2*t^5   == 
       6 r^3*s^3*t^3   == -6*d**3 

    """
    #print(p)
    #p += 6*d**3
    #print(p)
    #print(c**3)

    vals = [-2,-1,0,1,2]
    S3 = True
    for vr in vals:
      for vs in vals:
        for vt in vals:
            #print((r,s,t,v), end=' ')
            v = p(r=vr,s=vs,t=vt)
            assert v == p(r=vs,s=vt,t=vr)
            assert v == p(r=vt,s=vr,t=vs)
            w = p(r=vr,s=vt,t=vs) # != v in general
            assert w == p(r=vs,s=vr,t=vt)
            assert w == p(r=vt,s=vs,t=vr)
            if w != v:
                S3 = False
    print("S3:", S3)
    

def main_quartic():

    B = sage.CyclotomicField(3)
    w, = B.gens()

    K = sage.PolynomialRing(B, list("rstx"))
    vs = K.gens()
    r, s, t, x = vs

    y = (x-r)*(x-s)*(x-t)

    #print(y)
    #for item in y:
    #    print(item)

    a = 1
    b = -(r+s+t)
    c = r*s + r*t + s*t
    d = -r*s*t
    assert y == a*x**3 + b*x**2 + c*x + d

    X, Y, Z = r, s, t
    blue = X*Y**2*Z**3
    red = Y*Z**2*X**3
    green = Z*X**2*Y**3

    p = (blue + w*red + w**2*green)**3
    for coeff,monic in p:
        print("%14s"%coeff, monic)
    #print(len(p))
    """
    p = 
    1 * (
        r^3*s^6*t^9 + 
        r^6*s^9*t^3 + 
        r^9*s^3*t^6 ) == d**3 * (s^3*t^6 + r^3 ... ) == d**3 * ( ? * b**3 ? - ... )

    -3*w**2 * (
        r^4*s^7*t^7 + 
        r^7*s^4*t^7 + 
        r^7*s^7*t^4 ) == d**4 * (r**3*s**3 + ... ) == d**4 * (c**3 - ... )

    3*w * (
        r^5*s^5*t^8 + 
        r^5*s^8*t^5 + 
        r^8*s^5*t^5 ) == d**5 * (r**3 + s**3 + t**3) == d**6 * (b**3 - ... )

    6 * 
        r^6*s^6*t^6 == d**6
    """

    print()
    q = p - 6 * d**6 
    for coeff,monic in q:
        print("%14s"%coeff, monic)

    return

    vals = [-2,-1,0,1,2]
    for vr in vals:
      for vs in vals:
        for vt in vals:
            #print((r,s,t,v), end=' ')
            v = p(r=vr,s=vs,t=vt)
            assert v == p(r=vs,s=vt,t=vr)
            assert v == p(r=vt,s=vr,t=vs)
            w = p(r=vr,s=vt,t=vs)
            assert w == p(r=vs,s=vr,t=vt)
            assert w == p(r=vt,s=vs,t=vr)
    
    for vr in range(10):
        for vs in range(10):
            v = p(r=vr, s=vs, t=0)
            assert v==0
            
        #print()
    



if __name__ == "__main__":

    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()


