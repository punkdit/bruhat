#!/usr/bin/env python

from operator import mul, matmul, add
from functools import reduce

from sage.all_cmdline import ZZ, QQ

#from bruhat.matrix_sage import Matrix
from bruhat import matrix_sage
from bruhat.argv import argv
from bruhat.gset import mulclose, Perm, Group


ring = QQ
one = ring.one()

def Matrix(v):
    return matrix_sage.Matrix(ring, v)

def identity(n):
    return matrix_sage.Matrix.identity(ring, n)

def make_free(n):
    if not n:
        return []
    r = one / n
    vs = []
    for i in range(n):
        v = []
        for j in range(n):
            if j==i:
                v.append((n-1)*one/n)
            else:
                v.append((-1)*one/n)
        assert sum(v) == 0
        v = Matrix(v).t
        vs.append(v)
    assert reduce(add, vs).sum() == 0
    return vs


def conv(a, b, r=one/2):
    return r*a + (1-r)*b


#def stack(vs):
#    u = None
#    for v in vs:


class Conv:
    def __init__(self, gens):
        self.gens = list(gens)
        # put gens in the columns of u:
        u = None
        for v in gens:
            u = v if u is None else u.augment(v)
        self.u = u
        shape = None
        for v in gens:
            assert v.sum() == 0
            assert shape is None or v.shape == shape
            shape = v.shape

    def __getitem__(self, idx):
        return self.gens[idx]

    def __len__(self):
        return len(self.gens)

    def __str__(self):
        #s = str(self.gens)
        s = str(self.u)
        #s = s.replace("\n", "")
        #s = s.replace(" ", "")
        return "Conv(\n%s, dim=%d)"%(s, self.dim)

    @property
    def dim(self):
        return self.u.t.rank()

    @classmethod
    def free(cls, n):
        vs = make_free(n)
        return Conv(vs)

    def _quotient1(self, l, r, rels):
        u = l-r
        uu = u@u.d       
        K = uu.cokernel()
        #print("_quotient1:")
        #print(K)
        #print(self)
        gens = [K*v for v in self.gens]
        rels = [(K*l, K*r) for (l,r) in rels]
        return Conv(gens), rels

    def quotient(self, *rels):
        #print("\nquotient:")
        other = self
        rels = list(rels)
        while rels:
            rel = rels.pop(0)
            other, rels = other._quotient1(*rel, rels)
            #print(other, len(rels))
        return other

    def _fail_quotient(self, *rels):
        # argh, how to do this in one go?
        u = None
        for (l,r) in rels:
            lr = l-r
            op = lr@lr.d 
            u = op if u is None else u.stack(op)
        K = u.cokernel()
        print("quotient:")
        print(u)
        print(K)
        gens = [K*v for v in self.gens]
        return Conv(gens)


def test():

    v = Matrix([1,0,0]).t

    a, b, c = make_free(3)
    assert (a+b+c).sum() == 0
    
    assert conv(a, b) == (one/2)*(a+b)
    assert conv(a, b, one/3) == (one/3)*(a+2*b)


    vs = make_free(4)
    a, b, c, d = vs
    u = conv(a,b,one/3) - conv(c, d)

    uu = ( u@u.d )
    assert uu.rank() == 1

    K = uu.cokernel()

    ws = [K*v for v in vs]
    a1, b1, c1, d1 = ws

    assert conv(a1,b1,one/3) == conv(c1, d1)

    for w in ws:
        assert w.sum() == 0

    # -------------------------------------------------------

    space = Conv.free(4)
    assert space.dim == 3

    # now construct the square by imposing one relator
    a, b, c, d = space
    other = space.quotient( (conv(a,b), conv(c,d) ) )
    a, b, c, d = other
    assert conv(a,b) == conv(c,d)
    assert other.dim == 2

    a, b, c, d = space
    other = space.quotient( (conv(a,b,one/3), conv(c,d) ) )
    a, b, c, d = other
    assert conv(a,b,one/3) == conv(c,d)
    assert other.dim == 2

    a, b, c, d = space
    other = space.quotient( (conv(a,b), conv(c,d)), (a, d) )
    assert other.dim == 1
    a, b, c, d = other
    assert a==d

    #return

    # now we construct the tensor product of the square with
    # itself, giving an 8 dimensional conv space
    space = Conv.free(16)
    lookup = {(i,j):space[i+4*j] for i in range(4) for j in range(4)}

    rels = []
    for i in range(4):
        a, b, c, d = [lookup[i,j] for j in range(4)]
        rels.append( (conv(a,b), conv(c,d)) )

        a, b, c, d = [lookup[j,i] for j in range(4)]
        rels.append( (conv(a,b), conv(c,d)) )

    other = space.quotient(*rels)
    assert other.dim == 8

    # -------------------------------------------------------

    space = Conv.free(6)
    a, b, c, d, e, f = space

    # construct an Octahedron == Bloch Octahedron
    octa = space.quotient( 
        (conv(a,b), conv(c,d)),
        (conv(a,b), conv(e,f)),
    )
    assert octa.dim == 3

    # tensor square of this is the separable states
    N = 6
    space = Conv.free(N*N)
    lookup = {(i,j):space[i+N*j] for i in range(N) for j in range(N)}

    rels = []
    for i in range(N):
        a, b, c, d, e, f = [lookup[i,j] for j in range(N)]
        rels.append( (conv(a,b), conv(c,d)) )
        rels.append( (conv(a,b), conv(e,f)) )

        a, b, c, d, e, f = [lookup[j,i] for j in range(N)]
        rels.append( (conv(a,b), conv(c,d)) )
        rels.append( (conv(a,b), conv(e,f)) )

    assert len(rels) == 24
    other = space.quotient(*rels)
    assert other.dim == 15

    print(len(space.gens))
    print(len(set(space.gens)))

    # -------------------------------------------------

    names = [(pm, a, b) for pm in (1,-1) for a in "XYZ" for b in "XYZ"]
    print(names)
    N = len(names)

    space = Conv.free(N)
    lookup = {name:space[i] for (i,name) in enumerate(names)}

    mixed = []
    for a in "XYZ":
      for b in "XYZ":
        l = lookup[(+1, a, b)]
        r = lookup[(-1, a, b)]
        mixed.append( conv(l,r) )
    l = mixed[0]
    rels = [(l,r) for r in mixed[1:]]

    other = space.quotient(*rels)
    print(other.dim)


def test_clifford():
    from bruhat.clifford_sage import Clifford, K

    n = argv.get("n", 1)

    c = Clifford(n)
    S, H, CX, CZ = c.S, c.H, c.CX, c.CZ
    X, Y, Z = c.X, c.Y, c.Z
    wI = c.w()
    I = c.I

    # use generators in SU(N)?
    gen = []
    for i in range(n):
        #gen.append(wI*S(i))
        gen.append(S(i))
        #gen[-1].name = ("S",)
        gen.append(H(i))
        #gen[-1].name = ("H",)
        for j in range(i+1, n):
            gen.append(CZ(i,j))
    #gen = [wI*S(0), wI*S(1), wI*wI*H(0), wI*wI*H(1), wI*CZ()]
    for g in gen:
        #print(g)
        #assert g.determinant() == 1
        assert g*g.d == I

#    # for n=2 we get 46080 or 92160 depending on generators...
#    G = mulclose(gen, verbose=True)
#    print("|G| =", len(G))

    N = 2**n

    def canonical(v):
        for i in range(N):
            r = v.M[i,0]
            if r:
                break
        else:
            assert 0
        v = (1/r)*v
        return v


    v = [0]*N
    v[0] = 1
    v = matrix_sage.Matrix(K, [v]).t
#    stab = [g for g in G if canonical(g*v)==v]
#    print("stab =", len(stab))
#
#    for g in stab:
#        print(g)
#    return

    bdy = [v]
    found = set(bdy)
    while bdy:
        _bdy = []
        for v in bdy:
          for g in gen:
            gv = g*v
            #gv = canonical(gv)
            if gv in found:
                continue
            found.add(gv)
            _bdy.append(gv)
        bdy = _bdy
        print(len(found), end=" ", flush=True)
    print()

    found = list(found)
    items = set(canonical(v) for v in found)
    items = list(items)

    lookup = {v:i for (i,v) in enumerate(items)}

    def getop(op):
        idxs = []
        for (i,v) in enumerate(items):
            u = op*v
            j = lookup[canonical(u)]
            idxs.append(j)
        return Perm(idxs)

    pauli = []
    for i in range(n):
        pauli.append( getop(X(i)) )
        pauli.append( getop(Y(i)) )
        pauli.append( getop(Z(i)) )
    G = Group.generate(pauli)
    print("Pauli:", len(G))

    M = len(items)
    space = Conv.free(M)

    if 0:
        mixed = []
        for i in range(M):
            orbit = set(g[i] for g in G)
            assert len(orbit) == N
            v = reduce(add, [space[j] for j in orbit])
            v = (one/N)*v
            mixed.append(v)
        l = mixed[0]
        rels = [(l,r) for r in mixed[1:]]
    
        #space = space.quotient(*rels)
        #print(space.dim)

    op = getop(S(0))
    perms = [getop(g) for g in gen]
    #G = mulclose(perms, verbose=True)
    #orbit = { g*op*~g for g in G }

    bdy = [op]
    orbit = set(bdy)
    while bdy:
        _bdy = []
        for op in bdy:
          for g in perms:
            gop = g*op*~g
            if gop not in orbit:
                orbit.add(gop)
                _bdy.append(gop)
        bdy = _bdy

    print("S ops:", len(orbit))

    squares = set()
    found = set()
    for s in orbit:
        assert s.order() == 4
        for item in s.get_orbits():
            item = tuple(item)
            if len(item)!=4:
                continue
            i = item[0]
            item = (i, s[i], s[s[i]], s[s[s[i]]])
            #if item in found:
            #    continue
            squares.add(item)
    print("squares:", len(squares))
    rels = []
    for item in squares:
        if 0 in item:
            print(item)

        i, j, k, l = item
        rels.append((conv(space[i], space[k]), conv(space[j], space[l])))

    print("quotient:")
    other = space.quotient(*rels)
    print(other.dim)

    #for g in G:
    #    f = g.fixed()
    #    print("fixed:", len(f), "order:", g.order())

    return

    ops = []
    for i in range(n):
        ops.append(X(i))
        ops.append(Y(i))
        ops.append(Z(i))

    found = list(found)
    for v in found:
        for op in ops:
            if op*v==v:
                print("+1", end=" ")
            elif op*v==-v:
                print("-1", end=" ")
            else:
                print(" .", end=" ")
        #print([int(op*v==v) for op in ops], end=" ") 
        #print([int(op*v==-v) for op in ops]) 
        #print(v.t)
        print()
    
#    for v in found:
#        print(v.t, [int(v==u) for u in found])
#    for v in found:
#        print(hash(v), end=" ")
#    print()




if __name__ == "__main__":

    from time import time
    start_time = time()

    profile = argv.profile
    name = argv.next() or "test"
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
        test()


    t = time() - start_time
    print("OK! finished in %.3f seconds\n"%t)




