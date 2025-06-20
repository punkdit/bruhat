#!/usr/bin/env python

from operator import mul, matmul, add
from functools import reduce

from sage.all import ZZ, QQ

import numpy

#from bruhat.matrix_sage import Matrix
from bruhat import matrix_sage
from bruhat.argv import argv
from bruhat.gset import mulclose, Perm, Group
from bruhat.solve import shortstr


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
            #assert v.sum() == 0 
            assert shape is None or v.shape == shape
            shape = v.shape

    def __getitem__(self, idx):
        return self.gens[idx]

    def __len__(self):
        return len(self.gens)

    def longstr(self):
        s = str(self.u)
        return "Conv(\n%s, dim=%d)"%(s, self.dim)

    def __str__(self):
        return "Conv(dim=%d)"%(self.dim)

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
        #print("cokernel", K.shape)
        #print("_quotient1:")
        #print(K)
        #print(self)
        gens = [K*v for v in self.gens]
        rels = [(K*l, K*r) for (l,r) in rels]
        return Conv(gens), rels

    def slow_quotient(self, *rels):
        #print("\nquotient:")
        other = self
        rels = list(rels)
        while rels:
            rel = rels.pop(0)
            other, rels = other._quotient1(*rel, rels)
            print("quotient1:", other.dim)
            #print(other, len(rels))
        return other

    def fast_quotient(self, *rels):
        #print("fast_quotient:", self)
        if not rels:
            return self
        rows = []
        for (l,r) in rels:
            row = (l-r).t.M[0]
            rows.append(row)
        u = Matrix(rows)
        #print("kernel:")
        K = u.kernel().t
        #print(K.shape)
        gens = [K*v for v in self.gens]
        return Conv(gens)

    quotient = fast_quotient

    def get_polyhedron(self):
        from sage.all import Polyhedron
        vs = [v.t.M[0] for v in self.gens]
        p = Polyhedron(vertices=vs)
        return p


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


def canonical(v):
    N = v.shape[0]
    for i in range(N):
        r = v.M[i,0]
        if r:
            break
    else:
        assert 0
    v = (1/r)*v
    return v



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

    other = space.quotient(*rels)
    assert other.dim == N**2-1


def get_clifford_hull(n, verbose=False):
    "find convex stabilizer polytope using density matrices"

    from bruhat.clifford_sage import Clifford, K, w4

    c = Clifford(n)
    S, H, CX, CZ = c.S, c.H, c.CX, c.CZ
    X, Y, Z = c.X, c.Y, c.Z
    wI = c.w()
    I = c.I

    gen = []
    for i in range(n):
        gen.append(S(i))
        gen.append(H(i))
        for j in range(i+1, n):
            gen.append(CZ(i,j))
    for g in gen:
        assert g*g.d == I

    #G = mulclose(gen)
    #print(len(G))

    N = 2**n
    v = [0]*N
    v[0] = 1
    v = matrix_sage.Matrix(K, [v]).t

#    orbit = set()
#    for g in G:
#        u = g*v
#        u = canonical(u)
#        orbit.add(u)
#    print(len(orbit))

    rho = v@v.d
    #orbit = set([g*rho*~g for g in G])
    #print(len(orbit))

    pairs = [(g,~g) for g in gen]

    bdy = [rho]
    orbit = set(bdy)
    while bdy:
        _bdy = []
        for rho in bdy:
          #for g in gen:
          for (g,ig) in pairs:
            grho = g*rho*ig
            if grho not in orbit:
                orbit.add(grho)
                _bdy.append(grho)
        bdy = _bdy
        if verbose:
            print(len(orbit), end=" ", flush=True)
    if verbose:
        print()

    M = (2**n)**2
    zero = Matrix([0]*M)
    verts = []
    for rho in orbit:
        #print(rho)
        A = rho.real().coerce(QQ)
        #print(A.shape, M)
        A = A.reshape(1,M)
        B = ((-w4)*rho.imag()).coerce(QQ).reshape(1,M)
        A = A.stack(zero).reshape(1,2*M)
        B = zero.stack(B).reshape(1,2*M)
        v = A+B
        verts.append(v.t)

    u = reduce(add, verts)
    u = (QQ.one()/len(verts)) * u
    verts = [v-u for v in verts]

    space = Conv(verts)
    return space


def test_orbit():

    print("test_orbit")
    n = argv.get("n", 1)
    space = get_clifford_hull(n)
    print(space)
    p = space.get_polyhedron()

    print(p)
    #print(" ".join(dir(p)))
    #print(p.Hrepresentation())
    #print(p.face_lattice())
    i = 0
    while 1:
        fs = p.faces(i)
        print(i, len(fs))
        if len(fs)==0:
            break
        i += 1


def get_autos(A, B):
    m, n = A.shape
    k, l = B.shape
    assert l==n

    from pynauty import Graph, autgrp
    g = Graph(n+m+k) # bits + checks

    for bit in range(n):
        checks = [n+check for check in range(m) if A[check, bit]]
        g.connect_vertex(bit, checks)

    for check in range(m):
        bits = [bit for bit in range(n) if A[check, bit]]
        g.connect_vertex(n+check, bits)

    for face in range(k):
        idxs = [i for i in range(n) if B[face, i]]
        g.connect_vertex(n+m+face, idxs)

    g.set_vertex_coloring([set(range(n)), set(range(n, m+n)), set(range(m+n, m+n+k))])
    aut = autgrp(g)

    gen = aut[0]
    print(aut[1])
    assert aut[1] == int(aut[1]), str(aut[1])
    N = int(aut[1])
    items = list(range(n))
    perms = []
    for perm in gen:
        perm = perm[:n]
        #print(perm)
        perm = Perm(perm)
        perms.append(perm)
    #print(gap_code(perms))
    return N, perms




def test_autos():
    
    print("test_autos")
    n = argv.get("n", 1)
    space = get_clifford_hull(n)
    print(space)

    p = space.get_polyhedron()
    print(p)

    verts = p.faces(0)
    edges = p.faces(1)
    faces = p.faces(2)

    print(len(verts))

    #for v in verts:
    v = verts[0]
    e = edges[0]
    #for name in dir(e):
    #    print(name)

    #print(e.ambient_Vrepresentation())

    #return

    verts = [v.vertices()[0] for v in verts]
    edges = [set(e.vertices()) for e in edges]
    faces = [set(f.vertices()) for f in faces]

    lookup = {v:i for (i,v) in enumerate(verts)}

    A = numpy.zeros((len(edges), len(verts)), dtype=int)
    B = numpy.zeros((len(faces), len(verts)), dtype=int)

    for idx,(v0,v1) in enumerate(edges):
        A[idx, lookup[v0]] = 1
        A[idx, lookup[v1]] = 1

    for idx,f in enumerate(faces):
        for v in f:
            B[idx, lookup[v]] = 1

    print(A.shape)
    print(shortstr(A))
    print()
    print(shortstr(B))

    #A = A.transpose()
    N, perms = get_autos(A, B)

    from bruhat.oeqc import gap_code
    print(gap_code(perms))




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




