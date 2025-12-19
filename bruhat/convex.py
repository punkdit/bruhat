#!/usr/bin/env python

from operator import mul, matmul, add
from functools import reduce
from math import sin, cos, pi
from os import popen

from sage.all import ZZ, QQ

import numpy

#from bruhat.matrix_sage import Matrix
from bruhat import matrix_sage
from bruhat.argv import argv
from bruhat.gset import mulclose, Perm, Group, gap_code
from bruhat.action import mulclose_names, mulclose_hom
from bruhat.solve import shortstr
from bruhat.smap import SMap

from bruhat.clifford_sage import (
    Clifford, K, w4, r2, w8, get_clifford_gens, get_pauli_gens)

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



def get_clifford_states(n, local=False, verbose=False, name=False):
    cliff_gen = get_clifford_gens(n, local, name)
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

    return orbit, perms, pauli, cliff_gen # arggghhh


def get_clifford_hull(n, local=False, verbose=False):
    "find convex stabilizer polytope using density matrices"

    orbit, perms, pauli, cliff_gen = get_clifford_states(n, local, verbose)

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


def test_clifford_states():
    n = 1

    local = False
    verbose = False
    orbit, perms, pauli, cliff_gen = get_clifford_states(n, local, verbose)

    orbit = list(orbit)

    M = (2**n)**2
    zero = Matrix([0]*M)
    verts = []
    for rho in orbit:
        print(rho)


def test_coxeter():

    print("test_orbit")
    n = argv.get("n", 2)
    local = argv.get("local", False)
    space = get_clifford_hull(n, local=local)
    print(space)
    
    assert n < 3, "too big.."
    p = space.get_polyhedron()

    print(p)
    #print(" ".join(dir(p)))

    D = 2*((2**n)**2)
    I = identity(D)

    faces = p.faces(1)
    N = len(faces)
    print("dim %d, N=%d" % (1, N))
    gen = set()
    items = []
    for face in faces:
        idxs = space.get_idxs(face)
        vs = face.ambient_Vrepresentation()
        ms = [Matrix(v._vector) for v in vs]
        d = ms[1] - ms[0]
        r = (d*d.t).M[0,0]
        #print(str(d).replace(" ",""), r)
        #if r != 1:
        #    continue
        #assert r == 1
        #print(d, d*d.t)
        P = I - (2/r)*d.t*d
        assert P*P == I
        gen.add(P)
    print("gen:", len(gen))

    maxsize = argv.get("maxsize", 100000)
    G = mulclose(gen, maxsize=maxsize, verbose=True)

    print("\n|G| =", len(G))



def test_gap():
    n = 2

    orbit, gens, pauli, cliff_gen = get_clifford_states(n)
    N = len(orbit)
    assert N == 60
    #for gen in gens:
    #    print(gen)

    #PCliff = mulclose(gens, verbose=True)
    PCliff = Group.generate(gens, verbose=True)
    assert len(PCliff) == 11520
    print(PCliff)

    items = []
    for g in gens:
        item = g.gapstr()
        items.append(item)

    f = open("/tmp/bruhat_convex.gap", "w")
    print('SizeScreen([1000,1000]);;', file=f);
    gapstr = "Group(%s)"%(','.join(items))
    print(gapstr)
    return

    #from bruhat.tom import load_tom
    #tom = load_tom(gapstr=gapstr)
    #print(tom)
    #return

    print('G := %s;;'%gapstr, file=f)
    print('Hs := MaximalSubgroupClassReps(G);;', file=f)
    print('for H in Hs do Print(H, ";\\n"); od;', file=f)
    f.close()

    Xs = []
    X = PCliff.left_action(PCliff)
    Xs.append(X)

    data = popen("gap --norepl /tmp/bruhat_convex.gap").read()
    #open("gapdump.out", "w").write(data)
    lines = data.split("\n")
    for line in lines:
        line = line.strip()
        if not line.startswith("Group"):
            continue
        #print(line)
        line = line.replace(" ", "")
        END = ")]);"
        START = "Group([("
        assert line.endswith(END)
        assert line.startswith(START)
        line = line[len(START):-len(END)]
        cycs = line.split("),(")
        cycs = [cyc.split(")(") for cyc in cycs]
        cycs = [[tuple(int(i)-1 for i in c.split(","))
            for c in cyc] for cyc in cycs]
        hens = [Perm.from_cycles(N, cyc) for cyc in cycs]
        #print(cycs[0])
        #print(hens[0])
        H = Group.generate(hens)

        X = PCliff.left_action(H)
        print(H, X, len(X.gens))

        Xs.append(X)

    from bruhat.hecke import Builder
    builder = Builder(len(PCliff.gens))
    for X in Xs:
        builder.add(X)

    builder.get_hecke()

    #builder.build()
    builder.get_tom()
    #builder.dump()


def test_outer():

    from bruhat.algebraic import Algebraic, get_permrep
    G = Algebraic.Sp(4)
    print(len(G))

    G = get_permrep(G)
    print(G)
    N = G.rank

    from bruhat.gap import Gap
    gap = Gap()
    #_G = gap.define(G)
    auto = gap.AutomorphismGroup(G)
    auto = gap.to_group(N, auto) # FAIL XXX
    print(len(auto))

    tom = gap.tom(G)

    print("tom:", len(tom))


def test_tom():
    n = 2

    orbit, gens, pauli, cliff_gen = get_clifford_states(n)
    N = len(orbit)
    assert N == 60
    #for gen in gens:
    #    print(gen)
    #print(orbit[0])
    #print(pauli[0])
    #print(len(pauli))

    #PCliff = mulclose(gens, verbose=True)
    PCliff = Group.generate(gens, verbose=True)
    assert len(PCliff) == 11520
    print(PCliff)

    #PCliff = Group(gens=gens)

    from bruhat.gap import Gap
    gap = Gap()
    #_PCliff = gap.define(PCliff)
    auto = gap.AutomorphismGroup(PCliff)
    auto = gap.to_group(N, auto)
    print(len(auto))

    tom = gap.tom(PCliff)

    print("tom:", len(tom))

    n = len(tom)
    names = tom.names

    #Gap.DEBUG = True
    He = tom.get_stabilizer(N, "C")
    print(He)

    Hf = tom.get_stabilizer(N, "D")
    print(Hf)

    assert He != Hf
    assert not gap.IsConjugate(PCliff, He, Hf, get=True)

#    return
#
#    for g in auto:
#        if g in PCliff:
#            continue
#        #if g.order() != 2:
#        #    continue
#
#        He1 = He.conjugate(g)
#        if He1 == He:
#            continue
#        #_He1 = gap.define(He1)
#        #_Hf = gap.define(Hf)
#        #_PCliff = gap.define(PCliff)
#        if gap.IsConjugate(PCliff, He1, Hf, get=True) == "true":
#            print("T",end='',flush=True)
#        else:
#            print("_", end='', flush=True)
#
#    return

    def dump(count):
        print("%d:"%count, end=' ')
        total = 0
        for i in range(n):
            if tom[i, n-1] == count:
                print(names[i], end='')
                total += 1
        print("", total)

    dump(6) # C D  <--- colour / mutant-colour
    dump(15) # H I <--- fiber / mutant-fiber
    dump(30) # M N O P Q R S T U V
    dump(60) # c d e f g h i j k l m n 
    #dump(180) 
    #dump(11520 // 12)

    # colours(6) fibers(15) octahedra(30) vertices(60)
    #    D <------> H  <-----> M  <--------> e, f
    #    C <------> I  <-----> O  <<<-- mutant

    # e(60) <---> U(30)
    # f(60) <---> T(30)

    # Argh, the outer automorpism of PCliff(2) doesn't seem to swap the 
    # colours ?!?!? just swaps e & f structure ???!? check this.. XXX

    """

    * | M                        O
    ------------------------------------------------------
    I | C_1(90)+J_4(360)         O(30)+a_2(180)+z_2(240)
    H | M(30)+X_2(180)+y_2(240)  C_1(90)+J_4(360)
    
    * | e                        f
    ------------------------------------------------------
    I | q_1(180)+U_5(720)        i_1(180)+T_5(720)
    H | e(60)+Y_3(360)+n_4(480)  f(60)+m_3(360)+n_4(480)
    
    * | e                  f
    ------------------------------------------
    C | S_3(360)           o_3(360)
    D | R_1(120)+f_2(240)  T_1(120)+f_2(240)
    
    * | I            H
    ------------------------------
    C | P(30)+c(60)  x(90)
    D | D_1(90)      N(30)+d(60)


    How do the 30's relate to the 15's ?
    * | H                          I
    ----------------------------------------------------------
    N | N(30)+2*a_1(120)+L_2(180)  D_1(90)+b_4(360)
    P | x(90)+a_4(360)             P(30)+2*b_1(120)+M_2(180)
    M | M(30)+X_2(180)+y_2(240)    C_1(90)+J_4(360)
    O | C_1(90)+J_4(360)           O(30)+a_2(180)+z_2(240)
    S | S(30)+w_1(180)+y_2(240)    E_1(90)+a_4(360)
    V | F_1(90)+b_4(360)           V(30)+w_1(180)+z_2(240)
    Q | s(90)+e_3(360)             Q(30)+u_1(180)+r_2(240)
    R | r(90)+j_3(360)             R(30)+s_1(180)+u_2(240)
    T | y(90)+v_3(360)             T(30)+y_1(180)+r_2(240)
    U | B_1(90)+t_3(360)           U(30)+A_2(180)+u_2(240)

    # these are the 30's
    N*N = 2*N(30)+4*a_1(120)+b_4(360)
    P*P = 2*P(30)+4*b_1(120)+a_4(360)
    M*M = 2*M(30)+2*X_2(180)+2*y_2(240) <--- Octahedron ?
    O*O = 2*O(30)+2*a_2(180)+2*z_2(240) <--- Octahedron ?
    S*S = 2*S(30)+2*y_2(240)+a_4(360)
    V*V = 2*V(30)+2*z_2(240)+b_4(360)
    Q*Q = 2*Q(30)+Q_4(360)+x_4(480)
    R*R = 2*R(30)+N_4(360)+z_4(480)
    T*T = 2*T(30)+d_4(360)+x_4(480)
    U*U = 2*U(30)+c_4(360)+z_4(480)

    Octahedron * Octahedron = 2*Octahedron + 180 + ...

    # these are the 60's
    c*c = c(60)+b_1(120)+3*M_2(180)+6*a_4(360)+k_6(720)
    d*d = d(60)+a_1(120)+3*L_2(180)+6*b_4(360)+k_6(720)
    e*e = e(60)+q_1(180)+I_6(720)+e_6(720)+G_8(1920) # <------- vertex !
    f*f = f(60)+i_1(180)+l_5(720)+A_6(720)+G_8(1920) # <------- vertex !
    g*g = 4*g(60)+4*z_4(480)+x_7(1440)
    h*h = 4*h(60)+8*z_2(240)+2*k_6(720)
    i*i = 4*i(60)+8*y_2(240)+2*k_6(720)
    l*l = 4*l(60)+4*x_4(480)+w_7(1440)
    j*j = 4*j(60)+2*A_7(960)+A_8(1440)
    k*k = 4*k(60)+2*A_7(960)+t_7(1440)
    m*m = 4*m(60)+2*i_6(720)+2*A_7(960)
    n*n = 4*n(60)+2*g_6(720)+2*A_7(960)

    two vertices make an edge (or a vertex)
    there are 90 internal edges and 2*360+960=1680 external edges (1770 total)
    so the oriented edges: 180 internal, 2*720 + 1920 external: structure "e" or "f" above

    How do the vertices relate to the 6's ?
    * | C                        D
    ------------------------------------------------------
    e | S_3(360)                 R_1(120)+f_2(240)
    f | o_3(360)                 T_1(120)+f_2(240)


    * | H                        I
    ------------------------------------------------------
    e | e(60)+Y_3(360)+n_4(480)  q_1(180)+U_5(720)
    f | f(60)+m_3(360)+n_4(480)  i_1(180)+T_5(720)

    * | H                       I
    ----------------------------------------------------
    H | H(15)+F_1(90)+a_1(120)  b(45)+w_1(180)
    I | b(45)+w_1(180)          I(15)+E_1(90)+b_1(120)

    """

    for s in "MNOPQRSTUV":
        print("%s*%s ="%(s,s), tom.get_desc(tom[s]*tom[s], True))
    for s in "cdefghijklmn":
        print("%s*%s ="%(s,s), tom.get_desc(tom[s]*tom[s], True))
    #return

    #for l in "CD":
    #  for r in "HI":
    #    print("%s*%s ="%(l,r), tom.get_desc(tom[l] * tom[r]) )

    tabulate = tom.tabulate

    print(tabulate("MNOPQRSTUV", "HI", True))
    print()
    
    print(tabulate("cdefghijklmn", "CD", True))
    print()
    
    #print(tabulate("MNOPQRSTUV", "cdefghijklmn"))
    #print()

    for i in [tom.names.index(c) for c in "HI"]:
        assert tom[i,i] == 1 # um..

    print(tabulate("HI", "HI", True))
    print()

    return tom
    

def test_pcliff():
    n = 2

    orbit, perms, pauli, cliff_gen = get_clifford_states(n)
    N = len(orbit)
    assert N == 60
    #for gen in perms:
    #    print(gen)
    PCliff = mulclose(perms, verbose=True)
    assert len(PCliff) == 11520

    Pauli = mulclose(pauli)

    fibers = []
    remain = list(range(N))
    exclude = set()
    while remain:
        i = remain[0]
        fiber = {g[i] for g in Pauli}
        fiber = list(fiber)
        fiber.sort()
        fibers.append(fiber)
        #print(fiber)
        for i in fiber:
            remain.remove(i)
        for i in fiber:
          for j in fiber:
            exclude.add((i,j))
    print(fibers, len(fibers))

    twelves = []
    found = set()
    for g in PCliff:
        i = g.order()
        if i==12:
            twelves.append(g)
        if i not in found:
            #print(i)
            found.add(i)

    print("twelves:", len(twelves))
    #for g in twelves:
    #    for h in [g**i for i in range(12)]:
    #        print(h.order() == 12) # 4 of these
    #    break

    found = set()
    for g in twelves:
        orbits = g.get_orbits()
        if len(orbits) == 5:
            print("found!")
        counts = [len(o) for o in orbits]
        counts.sort()
        counts = tuple(counts)
        found.add(counts)
        #for o in orbits:
        #    print("\t", o, len(o))
    print(found)

    g = twelves[0]
    orbits = g.get_orbits()
    orbits.sort(key = len)

    _orbits = []
    idx = 0
    for o in orbits:
        i = o[idx%len(o)]
        #idx += 7
        p = [i]
        while g[i] not in p:
            p.append(g[i])
            i = g[i]
        print(o, p)
        #for j in o[1:]:
        #    assert g[i] == j
        #    i = j
        _orbits.append(p)
    orbits = _orbits

    edges = []
    for i in range(N):
        for j in range(i+1,N):
            if (i,j) in exclude:
                edges.append((i,j))

    print(len(edges))

    from huygens.namespace import Canvas, path, grey, st_thick

    cvs = Canvas()

    Rs = [0.5,1,1.5,2.0,3,5,7,9]
    #Rs = [1,2,3,4,5,6,7,8]
    dtheta = 0.
    coords = {}
    for i,R in enumerate(Rs):
        orbit = orbits[i]
        n = len(orbit)
        for j,idx in enumerate(orbit):
            theta = 2*pi*j/n + dtheta
            x, y = R*sin(theta), R*cos(theta)
            coords[idx] = (x,y)
        dtheta += pi/n

    for (i,j) in edges:
        x0, y0 = coords[i]
        x1, y1 = coords[j]
        cvs.stroke(path.line(x0, y0, x1, y1), st_thick+[grey])

    for idx in range(N):
        x, y = coords[idx]
        cvs.fill(path.circle(x, y, 0.05), [])
        

    cvs.writePDFfile("stabilizer_states")



def test_CZ_state():

    n = 2
    #gens = get_clifford_gens(n, name=True)
    orbit, perms, pauli, cliff_gen = get_clifford_states(n, name=True)

    mat = lambda v : matrix_sage.Matrix(K,v)

    v = mat([1,1,1,0]) # |CZ>
    #v = mat([1,1-r2]) @ mat([1,+w4])  # unentangled
    #v = mat([1,-w4,0,2*w4]) "D12" stabilizer
    v = v.t

    #return

    #G = mulclose(cliff_gen, verbose=True, maxsize=None)
    #print(len(G))

    hom = mulclose_hom(cliff_gen, perms, verbose=True)
    G = set(hom.keys())
    print(len(G))

    stab = []
    for g in G:
        if g*v != v:
            continue
        #print('.', end='', flush=True)
        #print(g, g==~g)
        print(g.name)
        stab.append(g)
    print()
    print(len(stab))

    perms = [hom[s] for s in stab]
    #return perms

    G = Group(perms)
    #print(G.gapstr())
    from bruhat.gap import Gap
    gap = Gap()
    print( gap.StructureDescription(G, get=True) )

    R = stab[0]
    for g in stab[1:]:
        R = R+g
    R = (one/len(stab))*R
    print(R*R == R)
    spaces = R.eigenvectors()
    for (val, vecs, dim) in spaces:
        if val != 1:
            continue
        print(vecs.t, dim)


def test_edges():
    print("test_edges")
    n = argv.get("n", 2)

    rhos, gens, pauli, cliff_gen = get_clifford_states(n)

    N = gens[0].rank
    verts = list(range(N))
    pairs = {(p,q) for p in verts for q in verts}
    print("pairs:", len(pairs))

    orbits = []
    while pairs:
        pair = pairs.pop()
        orbit = {pair}
        bdy = list(orbit)
        while bdy:
            _bdy = []
            for g in gens:
              for (p,q) in bdy:
                pair = (g[p], g[q])
                if pair in orbit:
                    continue
                orbit.add(pair)
                pairs.remove(pair)
                _bdy.append(pair)
            bdy = _bdy
        orbit = list(orbit)
        orbit.sort()
        print("orbit:", orbit[0])
        orbits.append(orbit)

    G = mulclose(gens)
    print(len(G))

    evals = lambda rho : [s[0] for s in rho.eigenvectors()]

    def purify(rho):
        spaces = rho.eigenvectors()
        for (val, vecs, dim) in spaces:
            if val==0:
                continue
            #v = vecs[0, :]
            for i in range(dim):
                yield vecs[:, i].t

    v = Matrix([1,1,1,0])
    magic = v.t@v

    i, j = orbits[2][0]
    s, t = rhos[i], rhos[j]
    for v in purify(s+t):
        print("magic:", v)
        break
    magic = v.t@v
    #print(magic)
    #return v

    rs = [(magic*rho).trace() for rho in rhos]
    top = max(rs)
    magic = [i for i in range(N) if rs[i]==top]
    magic = tuple(magic)
    print("magic idxs:", magic)


    #print(len(orbits))
    for o in orbits:
        i, j = o[0]
        s, t = rhos[i], rhos[j]
        A = (s*t).trace()
        print()
        print("evals:", evals((one/2)*(s+t)))
        print(p, A, len(o), len(o)//N)
    
        stab = []
        for g in G:
            gi, gj = g[i], g[j]
            if gi==i and gj==j or gi==j and gj==i:
                stab.append(g)
        print("stab:", len(stab), len(G)//len(stab))
        print("found magic idxs:", magic in o)

#    for o in orbits:
#        print(len(o))
#        for (i,j) in o:
#            if i!=0:
#                continue
#            print( "\t", evals(rhos[i]+rhos[j]) )


    
    

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
    
    orbit, perms, pauli, cliff_gen = get_clifford_states(n)
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

    orbit, perms, pauli, cliff_gen = get_clifford_states(n)
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
    n = argv.get("n", 2)
    local = argv.get("local", False)
    space = get_clifford_hull(n, local=local)
    print(space)
    
    assert n < 3, "too big.."
    p = space.get_polyhedron()

    print(p)
    #print(" ".join(dir(p)))

    dims = [0,1,2,3]
    dim = argv.get("dim")
    if dim is not None:
        dims = [dim]

    print("dims =", dims)

    for dim in dims:
        faces = p.faces(dim)
        N = len(faces)
        print("dim %d, N=%d" % (dim, N))
        items = []
        for face in faces:
            idxs = space.get_idxs(face)
            #print(face.ambient_Vrepresentation())
            #idxs = set(idxs)
            idxs.sort()
            idxs = tuple(idxs)
            items.append(idxs)
        #print("6's:", len([item for item in items if len(item)==6]))
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

def test_qutrit():
    vs = [
        (3,0,0,1,1,1,1,1,1,1,1,1),
        (0,3,0,1,1,1,1,1,1,1,1,1),
        (0,0,3,1,1,1,1,1,1,1,1,1),
        (1,1,1,3,0,0,1,1,1,1,1,1),
        (1,1,1,0,3,0,1,1,1,1,1,1),
        (1,1,1,0,0,3,1,1,1,1,1,1),
        (1,1,1,1,1,1,3,0,0,1,1,1),
        (1,1,1,1,1,1,0,3,0,1,1,1),
        (1,1,1,1,1,1,0,0,3,1,1,1),
        (1,1,1,1,1,1,1,1,1,3,0,0),
        (1,1,1,1,1,1,1,1,1,0,3,0),
        (1,1,1,1,1,1,1,1,1,0,0,3),
    ]

    vs = [Matrix(v).t for v in vs]
    v0 = (one/len(vs))*reduce(add, vs) 
    vs = [v-v0 for v in vs]
    space = Convex(vs)
    assert space.dim == 8

    p = space.get_polyhedron()

    print(p)

    found = []
    for dim in range(9):
        faces = p.faces(dim)
        N = len(faces)
        print("dim %d, N=%d" % (dim, N))


        

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

