#!/usr/bin/env python
"""
build modular curves X_0(N) 
and find their genus.
"""

import os

from bruhat.argv import argv
from bruhat.action import Group, Perm
from bruhat.smap import SMap
from bruhat.solve import shortstr, zeros2
from bruhat import elim, element


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


BODY = """
digraph {

    node [
        shape = circle
        color = "black"
        width = 0.2
        height = 0.2
        label = ""
    ]
    edge [ penwidth = 2.0 ]

%s

}
"""

def todot(N, red, green, blue, items):
    print("blue:", len(blue.orbits()))
    print("green:", len(green.orbits()))
    print("red:", len(red.orbits()))
    print("lozenges:", len(items))
    lines = []
    lookup = dict((v,k) for (k,v) in enumerate(items))
    n = len(items)
    for i,item in enumerate(items):
        j = lookup[blue[item]]
        lines.append('  %d -> %d [color="#0000b2"];'%(i,j))
        j = lookup[green[item]]
        lines.append('  %d -> %d [color="#009200"];'%(i,j))
        j = lookup[red[item]]
        lines.append('  %d -> %d [color="#d20000"];'%(i,j))
    stem = "X0_%d"%N
    body = BODY % ('\n'.join(lines))
    f = open("%s.dot"%stem, "w")
    print(body, file=f)
    f.close()
    rval = os.system("dot %s.dot -Tpdf > %s.pdf"%(stem, stem))
    assert rval == 0, rval
    print("wrote %s.pdf"%stem)


def get_genus(N, render=False):
    G = Group.cyclic(N)
    GG = G.direct_product(G)
    assert len(GG) == N**2

    # order N cyclic subgroups of GG
    found = []
    for g in GG:
        H = Group.generate([g])
        if len(H) == N and H not in found:
            found.append(H)

    lookup = {}
    ops = []
    for H in found:
        A = mkop(N, H)
        ops.append(A)
        lookup[shortstr(A)] = A
        #print(shortstr(A))
        #print("send_blue:")
        #print(shortstr(send_blue(A)))
        #print("send_green:")
        #print(shortstr(send_green(A)))
        #print()

    # blue  : z -> -1/z
    # green : z -> z+1
    # red   : z -> 1/(1-z)

    items = [shortstr(A) for A in ops]
    b = blue = Perm({shortstr(A) : shortstr(send_blue(A)) for A in ops}, items)
    g = green = Perm({shortstr(A) : shortstr(send_green(A)) for A in ops}, items)
    r = red = g*b*g*b

    if render:
        todot(N, red, green, blue, items)

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

def render():
    N = argv.get("N", 11)
    get_genus(N, True)


def test_atkin_lehner():
    from elim import eq, pseudo_inverse, shortstr, dot, array, identity
    # See: Harada 2010, p13
    # https://www.lmfdb.org/knowledge/show/cmf.atkin-lehner
    # generators of gamma_0(11):
    # https://swc-math.github.io/notes/files/01Weston1.pdf
    ring = element.Q

    N = 11

    p = ring.promote
    promote = lambda A : array([[p(A[0][0]), p(A[0][1])], [p(A[1][0]), p(A[1][1])]])
    def find():
        n = 7
        for x in range(-n, n+1):
         for y in range(-n, n+1):
          for z in range(-n, n+1):
            #W = [[N*x, y], [N*z, N]]
            d = N*x - y*z
            if d == 1:
                W = promote([[N*x, y], [N*z, N]])
                yield W
                #print(x, y, z)

    I = identity(ring, 2)

    # T,U,V generate gamma_0(11)
    T = promote([[1,1],[0,1]])
    U = promote([[7,-2],[11,-3]])
    V = promote([[8,-3],[11,-4]])
    for W in find():
        #print(W)
        Wi = pseudo_inverse(ring, W)
        W2 = dot(ring, W, W)
        assert eq(I, dot(ring, Wi, W))
        assert W2[1,0]%N == 0
        assert W2[0,0]%N == 0
        assert W2[1,1]%N == 0
        W2i = pseudo_inverse(ring, W2)
        assert eq(W2i, dot(ring, Wi, Wi))

        conj = lambda A : dot(ring, Wi, dot(ring, A, W))
        assert conj(T)[1,0] % N == 0
        #print(shortstr(conj(U)))
        assert conj(U)[1,0] % N == 0
        assert conj(V)[1,0] % N == 0

        conj2 = lambda A : dot(ring, W2i, dot(ring, A, W2))
        print("W:")
        print(shortstr(W))
        print("W2:")
        print(shortstr(W2))
        #print(shortstr(conj2(U)))
        print()


def main():
    print(" N=   genus=")
    for N in range(1,50):
        genus = get_genus(N)
        print("%4s       "%N, genus)


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



