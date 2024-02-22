#!/usr/bin/env python
"""
build modular curves X_0(N) 
and find their genus.
"""

from bruhat.argv import argv
from bruhat.action import Group, Perm
from bruhat.smap import SMap
from bruhat.solve import shortstr, zeros2

def hstr(n, H):
    s = SMap()
    for i in range(n):
     for j in range(n):
        s[i,j] = '.'
    for g in H:
        perm = g.perm
        (i,j) = perm[0,0]
        s[i,j] = '1'
    return str(s)
    

def mkop(n, H):
    A = zeros2(n,n)
    for g in H:
        A[g.perm[0,0]] = 1
    return A


def send_blue(A):
    n = len(A)
    B = zeros2(n, n)
    for j in range(n):
      for i in range(n):
        i1,j1 = j, -i
        B[i1%n, j1%n] = A[i, j]
    return B


def send_green(A):
    n = len(A)
    B = zeros2(n, n)
    for j in range(n):
      for i in range(n):
        B[(i+j)%n, j] = A[i, j]
    return B


def get_genus(n):
    G = Group.cyclic(n)
    GG = G.direct_product(G)
    assert len(GG) == n**2

    found = []
    for g in GG:
        H = Group.generate([g])
        if len(H) == n and H not in found:
            found.append(H)

    #print(len(found))
    lookup = {}
    ops = []
    for H in found:
        A = mkop(n, H)
        ops.append(A)
        lookup[shortstr(A)] = A
        #print(shortstr(A))
        #print("send_blue:")
        #print(shortstr(send_blue(A)))
        #print("send_green:")
        #print(shortstr(send_green(A)))
        #print()

    items = [shortstr(A) for A in ops]
    b = blue = Perm({shortstr(A) : shortstr(send_blue(A)) for A in ops}, items)
    g = green = Perm({shortstr(A) : shortstr(send_green(A)) for A in ops}, items)
    r = red = g*b*g*b

    #print("blue:", len(blue.orbits()))
    #print("green:", len(green.orbits()))
    #print("red:", len(red.orbits()))
    #print("lozenges:", len(ops))

    nblue = len(blue.orbits())
    ngreen = len(green.orbits())
    nred = len(red.orbits())
    nops = len(ops)
    chi = nred+nblue+ngreen - nops
    #print("chi:", chi)
    assert chi%2 == 0

    #G = Group.generate([blue,green], items)
    #print(len(G))

    return (2-chi)//2


    

def main():
    print(" N=   genus=")
    for n in range(1,33):
        genus = get_genus(n)
        print("%4s       "%n, genus)


if __name__ == "__main__":

    from time import sleep, time
    start_time = time()
    profile = argv.profile
    name = argv.next()
    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%(_seed))
        seed(_seed)

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)

    elif name is not None:
        fn = eval(name)
        fn()

    else:
        main()

    print("OK: finished in %.3f seconds\n"%(time() - start_time))



