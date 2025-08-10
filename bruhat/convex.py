#!/usr/bin/env python

from operator import mul, matmul, add
from functools import reduce

from sage.all import ZZ, QQ

import numpy

#from bruhat.matrix_sage import Matrix
from bruhat import matrix_sage
from bruhat.argv import argv
from bruhat.gset import mulclose, Perm, Group
from bruhat.action import mulclose_names
from bruhat.solve import shortstr

from bruhat.clifford_sage import Clifford, K, w4

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


class Convex:
    def __init__(self, gens, rhos=None):
        if rhos is None:
            rhos = gens
        assert len(gens) == len(rhos)
        gens = list(gens)
        # put gens in the columns of u:
        lookup = {}
        u = None
        for v in gens:
            lookup[v] = len(lookup)
            u = v if u is None else u.augment(v)
        shape = None
        for v in gens:
            #assert v.sum() == 0 
            assert shape is None or v.shape == shape
            shape = v.shape
        self.u = u
        self.gens = gens
        self.rhos = list(rhos)
        self.lookup = lookup

    def __getitem__(self, idx):
        return self.gens[idx]

    def __len__(self):
        return len(self.gens)

    def longstr(self):
        s = str(self.u)
        return "Convex(\n%s, dim=%d)"%(s, self.dim)

    def __str__(self):
        return "Convex(dim=%d)"%(self.dim)

    @property
    def dim(self):
        return self.u.t.rank()

    @classmethod
    def free(cls, n):
        vs = make_free(n)
        return Convex(vs)

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
        return Convex(gens), rels

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
        return Convex(gens)

    quotient = fast_quotient

    def get_polyhedron(self):
        from sage.all import Polyhedron
        vs = [v.t.M[0] for v in self.gens]
        p = Polyhedron(vertices=vs)
        return p

    def get_verts(self, face):
        verts = face.ambient_Vrepresentation()
        verts = [Matrix(v).t for v in verts]
        return verts

    def get_idxs(self, face):
        verts = self.get_verts(face)
        idxs = [self.lookup[v] for v in verts]
        return idxs


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

    space = Convex.free(4)
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
    space = Convex.free(16)
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

    space = Convex.free(6)
    a, b, c, d, e, f = space

    # construct an Octahedron == Bloch Octahedron
    octa = space.quotient( 
        (conv(a,b), conv(c,d)),
        (conv(a,b), conv(e,f)),
    )
    assert octa.dim == 3

    # tensor square of this is the separable states
    N = 6
    space = Convex.free(N*N)
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
    #from bruhat.clifford_sage import Clifford, K

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
    space = Convex.free(M)

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


def get_clifford_gens(n, local=False):

    c = Clifford(n)
    S, H, CX, CZ = c.S, c.H, c.CX, c.CZ
    X, Y, Z = c.X, c.Y, c.Z
    wI = c.w()
    I = c.I

    gen = []
    for i in range(n):
        gen.append(S(i))
        gen.append(H(i))
        if not local:
            for j in range(i+1, n):
                gen.append(CZ(i,j))
    for g in gen:
        assert g*g.d == I

    #G = mulclose(gen)
    #print(len(G))

    return gen


def get_pauli_gens(n):

    c = Clifford(n)
    X, Y, Z = c.X, c.Y, c.Z
    I = c.I

    gen = []
    for i in range(n):
        gen.append(X(i))
        gen.append(Y(i))
        gen.append(Z(i))
    for g in gen:
        assert g*g.d == I

    #G = mulclose(gen)
    #print(len(G))

    return gen




def get_clifford_states(n, local=False, verbose=False):
    cliff_gen = get_clifford_gens(n, local)
    pauli_gen = get_pauli_gens(n)

    N = 2**n
    v = [0]*N
    v[0] = 1
    v = matrix_sage.Matrix(K, [v]).t

    rho = v@v.d

    pairs = [(g,~g) for g in cliff_gen]

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

    orbit = list(orbit)
    orbit.sort(key = str)
    lookup = {rho:i for (i,rho) in enumerate(orbit)}

    perms = []
    for g,ig in pairs:
        perm = Perm([lookup[g*rho*ig] for rho in orbit])
        perms.append(perm)

    pauli = []
    for g in pauli_gen:
        perm = Perm([lookup[g*rho*~g] for rho in orbit])
        pauli.append(perm)

    #G = Group.generate(perms)
    #print(len(G))

    return orbit, perms, pauli


def get_clifford_hull(n, local=False, verbose=False):
    "find convex stabilizer polytope using density matrices"

    orbit, perms, pauli = get_clifford_states(n, local, verbose)

    orbit = list(orbit)

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

    space = Convex(verts, orbit)
    space.perms = perms # Clifford gens
    space.pauli = pauli # Pauli gens
    return space


def test_CZ_state():

    n = 2
    gens = get_clifford_gens(n)

    G = mulclose(gens, verbose=True, maxsize=None)
    print(len(G))

    v = Matrix([1,1,1,0]).t
    for g in G:
        if g*v != v:
            continue
        #print('.', end='', flush=True)
        print(g, g==~g)
    print()



def test_fixed():
    n = 2

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


    N = 2**n
    def get(idx):
        v = [0]*N
        v[idx] = 1
        v = matrix_sage.Matrix(K, [v]).t
        return v

    v = get(0) + get(1) + get(2)
    rho = v@v.d

    g = CZ(0,1)
    assert g*rho*~g == rho

    G = mulclose(gen, verbose=True)
    count = 0
    for g in G:
        grho = g*rho*~g
        if grho == rho:
            count += 1

    print("stab:", count)



def test_stabilizer_states():
    c = Clifford(1)
    I, X, Y, Z = c.I, c.X(), c.Y(), c.Z()

    n = 2
    
    pauli = [I, X, Y, Z]
    lookup = {}
    gens = []

    II = I@I

    for idx in numpy.ndindex((4,)*n):
        #lr = pauli[i]@pauli[j]
        #lookup[lr] = "IXYZ"[i] + "IXYZ"[j]
        #if i==j==0:
            #continue
        ops = [pauli[i] for i in idx]
        op = reduce(matmul, ops)
        name = ''.join("IXYZ"[i] for i in idx)
        lookup[op] = name
        lookup[-op] = "-"+name
        if sum(idx):
            gens.append(op)
            print(name, end=' ')
    print()
    print("gens:", len(gens))

    found = set()
    for a in gens:
      for b in gens:
        if a*b != b*a:
            continue
        if a==b:
            continue
        S = mulclose([a,b])
        assert len(S) == 4
        assert -II not in S
        S = list(S)
        S.sort(key = str)
        S = tuple(S)
        found.add(S)

    print(len(found))

    for S in found:
        print(" ".join(lookup[op] for op in S))



def test_bruhat():

    n = 2
    
    orbit, perms, pauli = get_clifford_states(n)
    print("orbit:", len(orbit))

    orbit = list(orbit)
    values = set()
    send = {}
    N = len(orbit)
    for i in range(N):
      #for j in range(i, N):
      for j in range(N):
        v = (orbit[i]*orbit[j]).trace()
        values.add(v)
        send.setdefault(v, []).append((i,j))
    print(len(values))
    values = list(values)
    values.sort()
    print(values)
    for v in values:
        j = len(send[v])
        print(v, j, j/60)


def get_clifford_orbit(n, k):

    assert 0<=k<=n
    N = 2**n
    gen = get_clifford_gens(n)

    v = [0]*N
    v[0] = 1
    v = Matrix(v).t
    rho = v.d @ v

    rho = numpy.zeros((N,N), dtype=int)
    for i in range(2**k):
        rho[i,i] = 1
    rho = matrix_sage.Matrix(K, rho)
    rho = (one / (2**k)) * rho
    assert rho.trace() == 1

    print("orbit:", end=' ')
    bdy = [rho]
    orbit = set(bdy)
    while bdy:
        _bdy = []
        for sho in bdy:
          for g in gen:
            tho = (~g)*sho*g
            if tho not in orbit:
                orbit.add(tho)
                _bdy.append(tho)
        bdy = _bdy
        print(len(orbit), end=" ", flush=True)
    print()
    return orbit



def test_density():

    n = argv.get("n", 3)
    k = argv.get("k", 0)

    orbit = get_clifford_orbit(n, k)

    print("total orbit size:", len(orbit))
    
    orbit = list(orbit)
    values = set()
    send = {}
    N = len(orbit)
    i = 0
    for j in range(N):
        v = (orbit[i]*orbit[j]).trace()
        values.add(v)
        send.setdefault(v, []).append((i,j))
    values = list(values)
    values.sort(key = lambda v:-eval(str(v)))
    print("overlap:", values)
    for v in values:
        j = len(send[v])
        #assert j%N == 0
        print(v, "\t", j)


def test_density_bundle():

    c = Clifford(1)
    I, X, Y, Z = c.I, c.X(), c.Y(), c.Z()
    pauli = [I, X, Y, Z]

    assert( (X@X) * (Z@Z) == -(Y@Y) )
    assert( (X@Y) * (Z@X) == -(Y@Z) )
    assert( (X@Z) * (Z@Y) == -(Y@X) )

    k = argv.get("k", 0)
    n = argv.get("n", 2)
    N = 2**n
    
    lookup = {}
    gens = []
    for idx in numpy.ndindex((4,)*n):
        ops = [pauli[i] for i in idx]
        op = reduce(matmul, ops)
        name = ''.join("IXYZ"[i] for i in idx)
        lookup[op] = name
        lookup[-op] = "-"+name
        #lookup[w4*op] = "i"+name
        #lookup[-w4*op] = "-i"+name
        if sum(idx):
            gens.append(op)

    print("gens:", len(gens))

    orbit = get_clifford_orbit(n, k)
    print("orbit:", len(orbit))

    G = mulclose(gens)
    In = reduce(matmul, [I]*n)

    def mkpoint(ops):
        point = [lookup[op] for op in ops]
        point.sort(key = lambda s:(-s.count("I"),s.count("-"),s))
        point = tuple(point)
        return point

    total_space = {} # base point to fiber
    for i, rho in enumerate(orbit):
        ops = [op for op in G if rho*op==rho]
        ops.remove(In)

        assert len(ops) == N-1

        base = [lookup[op].replace("-", "") for op in ops]
        base.sort()
        key = tuple(base)

        point = mkpoint(ops)
        total_space.setdefault(key, []).append(ops)

#        name = ','.join(point)
#        print("|%s>"%name, end=" ")
#    print()

    print(len(total_space))

    bases = list(total_space.keys())
    bases.sort(key = lambda base : (-str(base).count("I"), str(base)))
    for row,base in enumerate(bases):
        fiber = total_space[base]
        fiber = [mkpoint(ops) for ops in fiber]
        fiber.sort(key = lambda p:(str(p).count("-"), p))
        for i,point in enumerate(fiber):
            name = ','.join(point)
            #print(r"\ket{%s}"%name, end=" & " if i<N-1 else r" \\")
            print("%4s"%row, "|%s>"%name)
            break

            


def test_bundle():

    c = Clifford(1)
    I, X, Y, Z = c.I, c.X(), c.Y(), c.Z()
    pauli = [I, X, Y, Z]

    assert( (X@X) * (Z@Z) == -(Y@Y) )
    assert( (X@Y) * (Z@X) == -(Y@Z) )
    assert( (X@Z) * (Z@Y) == -(Y@X) )

    n = argv.get("n", 2)
    assert n<=2
    
    lookup = {}
    gens = []
    for idx in numpy.ndindex((4,)*n):
        ops = [pauli[i] for i in idx]
        op = reduce(matmul, ops)
        name = ''.join("IXYZ"[i] for i in idx)
        lookup[op] = name
        lookup[-op] = "-"+name
        #lookup[w4*op] = "i"+name
        #lookup[-w4*op] = "-i"+name
        if sum(idx):
            gens.append(op)

    print("gens:", len(gens))

    orbit, perms, pauli = get_clifford_states(n)
    print("orbit:", len(orbit))

    G = mulclose(gens)
    II = I@I

    def mkpoint(ops):
        assert len(ops)==3
        point = [lookup[op] for op in ops]
        point.sort(key = lambda s:(-s.count("I"),s.count("-"),s))
        point = tuple(point)
        return point

    total_space = set()
    for i, rho in enumerate(orbit):
        ops = [op for op in G if rho*op==rho]
        ops.remove(II)

        a, b, c = ops
        assert a*b == c
        assert a*c == b # etc

        fiber = []
        for i in [-1, 1]:
          for j in [-1, 1]:
            other = [i*a, j*b, i*j*c]
            point = mkpoint(other)
            fiber.append(point)
        assert len(fiber) == 4
        fiber.sort(key = lambda p:(str(p).count("-"), p))
        fiber = tuple(fiber)
        total_space.add(fiber)

        point = mkpoint(ops)
        name = ','.join(point)
        print("|%s>"%name, end=" ")
    print()

    found = list(set(total_space))
    print(len(found))
    for fiber in total_space:
        #print(fiber)
        for i,point in enumerate(fiber):
            name = ','.join(point)
            if argv.latex:
                print(r"\ket{%s}"%name, end=" & " if i<3 else r" \\")
            else:
                print(("|%s>"%name).ljust(13), end=" ")
        print()
            

        
def test_hull():

    print("test_hull")
    n = argv.get("n", 1)
    local = argv.get("local", False)
    space = get_clifford_hull(n, local=local)
    print(space)
    
    assert n < 3, "too big.."
    p = space.get_polyhedron()

    print(p)
    #print(" ".join(dir(p)))

    for dim in range(16):
        faces = p.faces(dim)
        N = len(faces)
        print("dim %d, N=%d" % (dim, N))


def render_hull(): # FAIL

    print("render_hull")
    n = argv.get("n", 2)
    local = argv.get("local", False)
    space = get_clifford_hull(n, local=local)
    print(space)

    #p = space.get_polyhedron()
    #faces = p.faces(0)
    #print(len(faces))

    N = 2**n
    shape = (N, N, 2)
    points = []
    for rho in space.rhos:
        #print(rho)
        A = numpy.zeros(shape, dtype=float)
        #print(type(rho.M), rho.M[0,0])
        for i in range(N):
         for j in range(N):
            a = rho.M[i,j]
            z = complex(a)
            A[i,j,0] = z.real
            A[i,j,1] = z.imag
        A = A.reshape(N*N*2)
        #print(A)
        points.append(A)
    M = N*N*2
    print(M)

    from math import sin, pi

#    u = numpy.zeros((M,))
#    #u[:M//2] = 1
#    #u[M//2:] = -1/4
#    for i in range(M):
#        u[i] = sin(2*pi*i/M)
#
#    v = numpy.zeros((M,))
#    #v[0::3] = 1
#    #v[1::3] = -1
#    #v[2::3] = 1/3
#    for i in range(M):
#        v[i] = sin(3*2*pi*i/M)

    # subtract center
    u = points[0]
    for p in points[1:]:
        u = u+p
    u = (1/len(points))*u
    points = [p-u for p in points]

    #u = points[7]
    #v = points[12] + points[44]

    # fail

    from huygens.namespace import Canvas, path

    cvs = Canvas()

    J = len(points)
    u = points[0]+points[1]
    for i in range(2, J):
      for j in range(i+1, J):
        v = points[i]  + points[j]

        fg = Canvas()
    
        r = 3
        for p in points:
            x = r*p.dot(u)
            y = r*p.dot(v)
            fg.fill(path.circle(x, y, 0.1))

        bb = fg.get_bound_box()
        fg.stroke(path.rect(bb.llx, bb.lly, bb.width, bb.height))

        cvs.insert(2*j*r, 2*i*r, fg)
        #cvs.translate(2*bb.width, 0)

    #cvs.writePDFfile("stabilizer_states")
    


def test_faces():
    print("test_faces")
    n = argv.get("n", 2)
    local = argv.get("local", False)
    space = get_clifford_hull(n, local=local)
    print(space)
    print("perms:", len(space.perms))

    #G = mulclose(space.perms) # 11520
    #print(len(G))

    names = []
    for i in range(n):
        names.append("X%d"%i)
        names.append("Y%d"%i)
        names.append("Z%d"%i)
    pauli = mulclose_names(space.pauli, names)
    identity = Perm(list(range(space.perms[0].rank)))
    assert identity in pauli
    pauli[identity] = ()
    print(len(pauli)) # projective pauli group

    #return
    
    assert n < 3, "too big.."
    p = space.get_polyhedron()

    print(p)
    #print(" ".join(dir(p)))

    found = []
    #for dim in [0,1,2,3]:
    for dim in [3]:
        faces = p.faces(dim)
        N = len(faces)
        print("dim %d, N=%d" % (dim, N))
        items = []
        for face in faces:
            idxs = space.get_idxs(face)
            if len(idxs) != 6:
                continue
            found.append(idxs)

    print("found:", len(found))

    lookup = {
        "" : "II",
        "I" : "II",
        "X0" : "XI",
        "Y0" : "YI",
        "Z0" : "ZI",
        "X1" : "IX",
        "Y1" : "IY",
        "Z1" : "IZ",
        "X0*X1" : "XX",
        "Y0*X1" : "YX",
        "Z0*X1" : "ZX",
        "X0*Y1" : "XY",
        "Y0*Y1" : "YY",
        "Z0*Y1" : "ZY",
        "X0*Z1" : "XZ",
        "Y0*Z1" : "YZ",
        "Z0*Z1" : "ZZ",
    }

    stabs = set()
    for idxs in found:
        #print(idxs)
        #idxs = set(idxs)
        ops = []
        for op in pauli:
            #jdxs = set(op[i] for i in idxs)
            jdxs = [op[i] for i in idxs]
            if jdxs == idxs:
                name = '*'.join(pauli[op]) or "I"
                ops.append(lookup.get(name, name))
        ops.sort()
        #assert len(ops) == 8
        name = ",".join(ops)
        stabs.add(name)
        #for g in space.pauli:
        #    print("\t", [g[i] for i in idxs])
    stabs = list(stabs)
    stabs.sort()
    for stab in stabs:
        print(stab)
    



def test_orbit():

    print("test_orbit")
    n = argv.get("n", 1)
    local = argv.get("local", False)
    space = get_clifford_hull(n, local=local)
    print(space)
    
    assert n < 3, "too big.."
    p = space.get_polyhedron()

    print(p)
    #print(" ".join(dir(p)))

    for dim in [0,1,2,3]:
        faces = p.faces(dim)
        N = len(faces)
        print("dim %d, N=%d" % (dim, N))
        items = []
        for face in faces:
            idxs = space.get_idxs(face)
            #idxs = set(idxs)
            idxs.sort()
            idxs = tuple(idxs)
            items.append(idxs)
        print("6's:", len([item for item in items if len(item)==6]))
        lookup = {item:idx for (idx,item) in enumerate(items)}
        perms = []
        for g in space.perms:
            idxs = []
            for item in items:
                jtem = [g[i] for i in item]
                jtem.sort()
                jtem = tuple(jtem)
                idxs.append(lookup[jtem])
            perms.append( Perm(idxs) )

        orbits = []
        remain = set(range(N))
        while remain:
            i = remain.pop()
            orbit = {i}
            bdy = list(orbit)
            while bdy:
                _bdy = []
                for g in perms:
                    for i in bdy:
                        j = g[i]
                        if j not in orbit:
                            orbit.add(j)
                            _bdy.append(j)
                            remain.remove(j)
                bdy = _bdy
            orbits.append(orbit)
        orbits.sort(key = len)
        print("orbits:", [len(orbit) for orbit in orbits])
        assert sum([len(orbit) for orbit in orbits]) == N

    return

    for i in range(N):
        v = vs[i]
        print(v.t)
        print(space.gens[idx].t)
        print()

    m = argv.get("m")
    if m is not None:
        fs = p.faces(m)
        print(len(fs))
        print(fs[0])
        return

    return

    #p = p*p

    #print(p.Hrepresentation())
    #print(p.face_lattice())
    counts = []
    total = -1
    i = 0
    while 1:
        c = len(p.faces(i))
        print(i, c)
        total += ((-1)**i)*c
        counts.append(c)
        if c==0:
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


def test_real_cliff():

    for n in [1,2]:
        gen = get_clifford_gens(n)

        G = mulclose(gen, verbose=True)

        G = set( g@g.d for g in G )

        print("|PCliff| = ", len(G))

        RG = [g for g in G if g==g.conjugate()]
        print("|RPCliff| = ", len(RG))

        print("factor =", len(G)//len(RG))

        #return
        

def test_double_cosets():
    from bruhat.matrix_sage import CyclotomicField, Matrix

    n = argv.get("n", 2)
    k = argv.get("k", 0)

    # -----------------------------------
    # build Pauli group

    ring = CyclotomicField(8)
    w8 = ring.gens()[0]
    w = w8**2
    r2 = w8 + w8.conjugate()
    ir2 = 1/r2
    assert r2**2 == 2
    assert 2*ir2**2 == 1

    I = Matrix(ring, [[1,0],[0,1]])
    wI = w*I
    X = Matrix(ring, [[0,1],[1,0]])
    Z = Matrix(ring, [[1,0],[0,-1]])
    S = Matrix(ring, [[1,0],[0,w]])
    H = Matrix(ring, [[ir2,ir2],[ir2,-ir2]])
    CZ = Matrix(ring, [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])

    gen = []
    for i in range(n):
        for op in [X,Z,wI]:
            g = [I]*n
            g[i] = op
            gen.append(reduce(matmul, g))
    
    Pauli = mulclose(gen)
    Pauli = list(Pauli)
    lookup = dict((g,i) for (i,g) in enumerate(Pauli))
    print("Pauli:", len(Pauli))
    assert len(Pauli) == 4**(n+1)

    #gen = [HI, IH, SI, IS, CZ]
    gen = get_clifford_gens(n)

    perms = []
    for g in gen:
        ig = ~g
        idxs = [lookup[ig*h*g] for h in Pauli]
        perm = Perm(idxs)
        perms.append(perm)
    gen = perms

    stab = []
    for i in range(n-k):
        op = [I]*n
        op[i] = X
        stab.append(reduce(matmul, op))
    stab = mulclose(stab)
    stab = [lookup[g] for g in stab]

    del lookup

    def mkpoint(stab):
        stab = list(stab)
        stab.sort()
        stab = tuple(stab)
        return stab

    stab = mkpoint(stab)
    codes = {stab}
    bdy = list(codes)
    while bdy:
        _bdy = []
        for g in gen:
            for p in bdy:
                q = mkpoint([g[i] for i in p])
                if q in codes:
                    continue
                codes.add(q)
                _bdy.append(q)
        bdy = _bdy
    if n==2:
        assert len(codes) == 60
    print("codes:", len(codes))
    codes = list(codes)
    codes.sort()
    lookup = {c:idx for (idx,c) in enumerate(codes)}
    N = len(codes)

    #pairs = {(p,q) for p in codes for q in codes}
    #print("pairs:", len(pairs))
    mask = numpy.zeros((N,N), dtype=numpy.uint8)
    mask[:] = 1

    #orbits = []
    row = col = 0
    #while pairs:
    counts = []
    while 1:
        #pair = iter(pairs).__next__()
        while row<N and mask[row,col] == 0:
            col += 1
            if col==N:
                row += 1
                col = 0
        if row==N:
            break
        pair = codes[row], codes[col]
        orbit = {pair}
        bdy = list(orbit)
        while bdy:
            _bdy = []
            for g in gen:
                for (p,q) in bdy:
                    p1 = mkpoint([g[i] for i in p])
                    q1 = mkpoint([g[i] for i in q])
                    if (p1,q1) in orbit:
                        continue
                    orbit.add((p1,q1))
                    _bdy.append((p1,q1))
            bdy = _bdy
        #orbits.append(orbit)
        #assert orbit.issubset(pairs)
        #pairs.difference_update(orbit)
        for (p,q) in orbit:
            i,j = lookup[p], lookup[q]
            assert mask[i,j]
            mask[i,j] = 0
        k = len(orbit)
        assert k%N == 0
        print("orbit:", k//N)
        counts.append(k//N)

    #orbits.sort(key=len)
    #print([len(o)//N for o in orbits], len(orbits))
    counts.sort()
    print(counts, len(counts))


def test_hecke(): # another version of test_double_cosets
    from bruhat.matrix_sage import CyclotomicField, Matrix

    n = argv.get("n", 2)
    k = argv.get("k", 0)

    # -----------------------------------
    # build Pauli group

    ring = CyclotomicField(8)
    w8 = ring.gens()[0]
    w = w8**2
    r2 = w8 + w8.conjugate()
    ir2 = 1/r2
    assert r2**2 == 2
    assert 2*ir2**2 == 1

    I = Matrix(ring, [[1,0],[0,1]])
    wI = w*I
    X = Matrix(ring, [[0,1],[1,0]])
    Z = Matrix(ring, [[1,0],[0,-1]])
    S = Matrix(ring, [[1,0],[0,w]])
    H = Matrix(ring, [[ir2,ir2],[ir2,-ir2]])
    CZ = Matrix(ring, [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]])

    gen = []
    for i in range(n):
        for op in [X,Z,wI]:
            g = [I]*n
            g[i] = op
            gen.append(reduce(matmul, g))
    
    Pauli = mulclose(gen)
    Pauli = list(Pauli)
    lookup = dict((g,i) for (i,g) in enumerate(Pauli))
    print("Pauli:", len(Pauli))
    assert len(Pauli) == 4**(n+1)

    #gen = [HI, IH, SI, IS, CZ]
    gen = get_clifford_gens(n)

    perms = []
    for g in gen:
        ig = ~g
        idxs = [lookup[ig*h*g] for h in Pauli]
        perm = Perm(idxs)
        perms.append(perm)
    gen = perms

    stab = []
    for i in range(n-k):
        op = [I]*n
        op[i] = X
        stab.append(reduce(matmul, op))
    stab = mulclose(stab)
    stab = [lookup[g] for g in stab]

    del lookup

    def mkpoint(stab):
        stab = list(stab)
        stab.sort()
        stab = tuple(stab)
        return stab

    stab = mkpoint(stab)
    codes = {stab}
    bdy = list(codes)
    while bdy:
        _bdy = []
        for g in gen:
            for p in bdy:
                q = mkpoint([g[i] for i in p])
                if q in codes:
                    continue
                codes.add(q)
                _bdy.append(q)
        bdy = _bdy
    if n==2 and k==0:
        assert len(codes) == 60
    print("codes:", len(codes))
    codes = list(codes)
    codes.sort()
    lookup = {c:idx for (idx,c) in enumerate(codes)}
    N = len(codes)

    mask = numpy.zeros((N,N), dtype=numpy.uint8)
    mask[:] = 1

    row = col = 0
    counts = []
    while 1:
        while row<N and mask[row,col] == 0:
            col += 1
            if col==N:
                row += 1
                col = 0
        if row==N:
            break
        pair = codes[row], codes[col]
        bdy = [(row,col)]
        count = 1
        while bdy:
            if argv.verbose:
                print("%d:%d"%(count,len(bdy)), end=" ", flush=True)
            _bdy = []
            while bdy:
                i,j = bdy.pop()
                p, q = codes[i], codes[j]
                for g in gen:
                    p1 = mkpoint([g[i] for i in p])
                    q1 = mkpoint([g[i] for i in q])
                    i,j = lookup[p1], lookup[q1]
                    if mask[i,j]:
                        _bdy.append((i,j))
                        mask[i,j] = 0
                        count += 1
            bdy = _bdy
        if argv.verbose:
            print()
        print("orbit:", count//N)
        counts.append(count//N)

    #orbits.sort(key=len)
    #print([len(o)//N for o in orbits], len(orbits))
    counts.sort()
    print(counts, len(counts))





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

