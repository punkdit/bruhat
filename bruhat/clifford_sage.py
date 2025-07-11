#!/usr/bin/env python

"""
build clifford/pauli groups in sage
"""

from random import choice
from operator import mul, matmul
from functools import reduce
#from functools import cache
from functools import lru_cache
cache = lru_cache(maxsize=None)

from sage.all_cmdline import FiniteField, CyclotomicField, latex, block_diagonal_matrix
from sage import all_cmdline 

from bruhat.solve import zeros2, identity2
from bruhat.action import mulclose, mulclose_find, mulclose_hom, mulclose_names
from bruhat.matrix_sage import Matrix
from bruhat.gset import Group, Perm, gap_code
from bruhat.util import cross
from bruhat.smap import SMap
from bruhat.argv import argv

K = CyclotomicField(8)
w8 = K.gen()
w4 = w8**2
r2 = w8 + w8**7
assert r2**2 == 2

def simplify_latex(self):
    M = self.M
    m, n = self.shape
    idxs = [(i,j) for i in range(m) for j in range(n)]
    for idx in idxs:
        if M[idx] != 0:
            break
    else:
        assert 0
    scale = M[idx]
    if scale != 1 and scale != -1:
        M = (1/scale) * M
        s = {
            r2 : r"\sqrt{2}",
            1/r2 : r"\frac{1}{\sqrt{2}}",
            #2/r2 : r"\frac{2}{\sqrt{2}}",
            2/r2 : r"\sqrt{2}",
            #r2/2 : r"\sqrt{2}/2",
        }.get(scale, latex(scale))
        if "+" in s:
            s = "("+s+")"
        s = "%s %s"%(s, latex(M))
    else:
        s = latex(M)
    s = s.replace(r"\zeta_{8}^{2}", "i")
    return s


class Coset(object):
    def __init__(self, group, items):
        self.group = group
        self.items = list(items)
        self.g = items[0] # pick one

    def __mul__(self, other):
        assert isinstance(other, Coset)
        gh = self.g * other.g
        result = self.group.lookup[gh]
        return result

    def __hash__(self):
        return id(self)


class FactorGroup(object):
    def __init__(self, G, H):
        remain = set(G)
        cosets = []
        while remain:
            g = iter(remain).__next__()
            coset = []
            for h in H:
                gh = g*h
                remain.remove(gh)
                coset.append(gh)
            coset = Coset(self, coset) # promote
            cosets.append(coset)
            #if len(remain)%1000 == 0:
            #    print(len(remain), end=" ", flush=True)
        #print()
        self.G = G
        self.H = H
        lookup = {g:coset for coset in cosets for g in coset.items}
        self.lookup = lookup
        #cosets = [Coset(self, coset) for coset in cosets] # promote
        self.cosets = cosets

    def __getitem__(self, idx):
        return self.cosets[idx]

    def __len__(self):
        return len(self.cosets)


class Clifford(object):
    "clifford group on n qubits"
    def __init__(self, n):
        self.n = n
        K = CyclotomicField(8)
        self.K = K
        self.I = Matrix.identity(K, 2**n)

    def w(self):
        K = self.K
        w = K.gen()
        return w*self.I

    def mkop(self, i, g):
        n = self.n
        K = self.K
        assert 0<=i<n
        I = Matrix.identity(K, 2)
        items = [I]*n
        items[i] = g
        gi = reduce(matmul, items)
        return gi
        #while len(items)>1:
        #    #items[-2:] = [items[-2] @ items[-1]]
        #    items[:2] = [items[0] @ items[1]]
        #return items[0]
        
    @cache
    def Z(self, i=0):
        K = self.K
        Z = Matrix(K, [[1, 0], [0, -1]])
        Zi = self.mkop(i, Z)
        return Zi
        
    @cache
    def S(self, i=0):
        K = self.K
        w = K.gen()
        S = Matrix(K, [[1, 0], [0, w*w]])
        Si = self.mkop(i, S)
        return Si
        
    @cache
    def X(self, i=0):
        K = self.K
        X = Matrix(K, [[0, 1], [1,  0]])
        Xi = self.mkop(i, X)
        return Xi

    @cache
    def Y(self, i=0):
        K = self.K
        Y = Matrix(K, [[0, w4], [-w4,  0]])
        Yi = self.mkop(i, Y)
        return Yi
        
    @cache
    def H(self, i=0):
        K = self.K
        w = K.gen()
        r2 = w+w.conjugate()
        ir2 = r2 / 2
        H = Matrix(K, [[ir2, ir2], [ir2, -ir2]])
        Hi = self.mkop(i, H)
        return Hi

    @cache
    def CZ(self, idx=0, jdx=1):
        n = self.n
        K = self.K
        assert 0<=idx<n
        assert 0<=jdx<n
        assert idx!=jdx
        N = 2**n
        A = zeros2(N, N)
        ii, jj = 2**(n-idx-1), 2**(n-jdx-1)
        for i in range(N):
            if i & ii and i & jj:
                A[i, i] = -1
            else:
                A[i, i] = 1
        return Matrix(K, A)

    @cache
    def CX(self, idx=0, jdx=1):
        assert idx != jdx
        CZ = self.CZ(idx, jdx)
        H = self.H(jdx)
        return H*CZ*H

    @cache
    def SWAP(self, idx=0, jdx=1):
        assert idx != jdx
        HH = self.H(idx) * self.H(jdx)
        CZ = self.CZ(idx, jdx)
        return HH*CZ*HH*CZ*HH*CZ


def test_clifford():

    c2 = Clifford(2)
    II = c2.I
    XI = c2.X(0)
    IX = c2.X(1)
    ZI = c2.Z(0)
    IZ = c2.Z(1)
    wI = c2.w()

    Pauli = mulclose([wI*wI, XI, IX, ZI, IZ])
    assert len(Pauli) == 64

    SI = c2.S(0)
    IS = c2.S(1)
    HI = c2.H(0)
    IH = c2.H(1)
    CZ = c2.CZ(0, 1)

    #C2 = mulclose([SI, IS, HI, IH, CZ], verbose=True) # slow
    #assert len(C2) == 92160

    CX = c2.CX(0, 1)
    SWAP = c2.SWAP(0, 1)

    assert SWAP * CZ * SWAP == CZ

    c2 = Clifford(3)
    A = c2.CX(0, 1)
    B = c2.CX(1, 2)
    assert A*B != B*A

    A = c2.CZ(0, 1)
    B = c2.CZ(1, 2)
    assert A*B == B*A

    CX01 = c2.CX(0, 1)
    CX10 = c2.CX(1, 0)
    CZ01 = c2.CZ(0, 1)
    CZ10 = c2.CZ(1, 0)
    SWAP01 = c2.SWAP(0, 1)

    assert HI*HI == II
    assert HI*XI*HI == ZI
    assert HI*ZI*HI == XI
    assert HI*IX*HI == IX
    assert SWAP01 == CX01*CX10*CX01



def test_clifford3():
    c3 = Clifford(3)
    wI = c3.w()
    III = c3.I
    XII = c3.X(0)
    IXI = c3.X(1)
    IIX = c3.X(2)
    ZII = c3.Z(0)
    IZI = c3.Z(1)
    IIZ = c3.Z(2)
    SII = c3.S(0)
    ISI = c3.S(1)
    IIS = c3.S(2)
    HII = c3.H(0)
    IHI = c3.H(1)
    IIH = c3.H(2)
    CX01 = c3.CX(0, 1)
    CX10 = c3.CX(1, 0)
    CX02 = c3.CX(0, 2)
    CX20 = c3.CX(2, 0)
    CX12 = c3.CX(1, 2)
    CX21 = c3.CX(2, 1)
    CX12 = c3.CX(1, 2)
    CZ01 = c3.CZ(0, 1)
    CZ02 = c3.CZ(0, 2)
    CZ12 = c3.CZ(1, 2)
    SWAP01 = c3.SWAP(0, 1)
    SWAP02 = c3.SWAP(0, 2)
    SWAP12 = c3.SWAP(1, 2)

    assert HII*HII == III
    assert HII*XII*HII == ZII
    assert HII*ZII*HII == XII
    assert HII*IXI*HII == IXI

    I = Clifford(1).I
    c2 = Clifford(2)
    assert CZ01 == c2.CZ(0,1)@I
    assert SWAP01 == CX01*CX10*CX01

    assert CX01*CX02 == CX02*CX01
    assert CX12*CX01 == CX02*CX01*CX12


def test_bruhat():
    K = CyclotomicField(8)
    w = K.gen()
    one = K.one()
    w2 = w*w
    r2 = w+w.conjugate()
    ir2 = r2 / 2
    
    I = Matrix(K, [[1, 0], [0, 1]])
    w2I = Matrix(K, [[w2, 0], [0, w2]])
    S = Matrix(K, [[1, 0], [0, w2]])
    X = Matrix(K, [[0, 1], [1, 0]])
    Z = Matrix(K, [[1, 0], [0, -1]])
    H = Matrix(K, [[ir2, ir2], [ir2, -ir2]])
    
    Pauli1 = mulclose([w2I, X, Z])
    Cliff1 = mulclose([w2I, S, H])
    
    assert len(Pauli1) == 16
    assert len(Cliff1) == 192
    
    CX = Matrix(K,
       [[1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0]])
    
    CZ = Matrix(K,
       [[1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0,-1]])
    
    SWAP = Matrix(K,
       [[1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1]])
    
    II = I @ I
    w2II = w2*II
    XI = X @ I
    IX = I @ X
    ZI = Z @ I
    IZ = I @ Z
    SI = S @ I
    IS = I @ S
    HI = H @ I
    IH = I @ H

    Pauli2 = mulclose([XI, IX, ZI, IZ])

    lhs = CZ * CX
    rhs = CX * CZ

    g = rhs * lhs.inverse()
    print(g)
    print(g == ZI)
    print(g == IZ)

    return

    Cliff2 = mulclose([SI, IS, HI, IH, CZ], verbose=True)
    #assert len(Cliff2) == 92160

    torus = []
    for g in Cliff2:
        if g.is_diagonal():
            #print(g)
            torus.append(g)
    print("torus:", len(torus))

    while 1:
        gen = [choice(torus) for i in range(4)]
        T = mulclose(gen)
        if len(T) == len(torus):
            break
    print("gen", len(gen))

    # ------------------------------------------------------------------
    # See:
    # https://arxiv.org/abs/2003.09412
    # Hadamard-free circuits expose the structure of the Clifford group
    # Sergey Bravyi, Dmitri Maslov
    # although this is for the projective Clifford group only 

    n = 2
    W = mulclose([HI, IH, SWAP]) # Weyl group
    assert len(W) == 8

    #B = mulclose([XI, IX, CX, CZ, SI, IS])
    B = mulclose([XI, IX, CX, CZ, SI, IS]+gen, verbose=True)
    print("Borel:", len(B))

    #T = [g for g in torus if g in B]
    #print(len(T))
    #return
    #B.extend(torus)

    # build the double cosets:
    dcs = {w:set() for w in W}
    total = 0
    for w in W:
        dc = dcs[w]
        size = 0
        for l in B:
          lw = l*w
          for r in B:
            lwr = lw*r
            dc.add(lwr)
          if len(dc) > size:
            size = len(dc)
            print(size, end=" ", flush=True)
        print("//")
        total += size
    print("total:", total) # == 92160 ?

    for w in W:
        for u in W:
            if u==w:
                continue
            a = dcs[w].intersection(dcs[u])
            #assert len(a) == 0
            print(len(a), end=" ")
    print()

    return dcs


def test_cocycle():
    K = CyclotomicField(8)
    w = K.gen()
    one = K.one()
    w2 = w*w
    r2 = w+w.conjugate()
    ir2 = r2 / 2
    
    I = Matrix(K, [[1, 0], [0, 1]])
    w2I = Matrix(K, [[w2, 0], [0, w2]])
    S = Matrix(K, [[1, 0], [0, w2]])
    X = Matrix(K, [[0, 1], [1, 0]])
    Z = Matrix(K, [[1, 0], [0, -1]])
    H = Matrix(K, [[ir2, ir2], [ir2, -ir2]])
    
    Pauli1 = mulclose([w2I, X, Z])
    Cliff1 = mulclose([w2I, S, H])
    
    assert len(Pauli1) == 16
    assert len(Cliff1) == 192
    
    CZ = Matrix(K,
       [[1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0,-1]])
    
    II = I @ I
    w2II = w2*II
    XI = X @ I
    IX = I @ X
    ZI = Z @ I
    IZ = I @ Z
    SI = S @ I
    IS = I @ S
    HI = H @ I
    IH = I @ H

    pauli_gen = [XI, IX, ZI, IZ]
    Pauli2 = mulclose([w2II] + pauli_gen)
    assert len(Pauli2) == 64

    Phase2 = mulclose([w2II])
    assert len(Phase2) == 4

    F4 = FactorGroup(Pauli2, Phase2)
    assert len(F4) == 16

    sy_gen = [ZI, IZ, XI, IX]
    pauli_lin = lambda g: [int(g*h != h*g) for h in sy_gen]
    
    # no phase generator needed here
    Cliff2 = mulclose([SI, IS, HI, IH, CZ], verbose=True)
    assert len(Cliff2) == 92160

    F2 = FiniteField(2)
    cliff_lin = lambda g:Matrix(F2, [pauli_lin(g*h*g.inverse()) for h in pauli_gen]).transpose()
    Sp4 = mulclose([cliff_lin(g) for g in [SI,IS,HI,IH,CZ]])
    assert len(Sp4) == 720

    Mp4 = FactorGroup(Cliff2, Pauli2)
    assert len(Mp4) == 1440

    hom = {} # Mp4 --> Sp4
    for coset in Mp4:
        g = coset.g # Cliff2
        assert g in Cliff2
        h = cliff_lin(g)
        hom[coset] = h
    
    Sp4_i = Matrix(F2, [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    kernel = [g for g in Mp4 if hom[g]==Sp4_i]
    assert len(kernel) == 2
    
    homi = {g:[] for g in Sp4}
    for g in Mp4:
        homi[hom[g]].append(g)
    #print([len(v) for v in homi.values()])
    def cocyc(g, h): # Sp4xSp4 --> Z/2
        gh = g*h
        gi, hi, ghi = homi[g], homi[h], homi[gh]
        lhs = gi[0]*hi[0]
        assert lhs in ghi
        return ghi.index(lhs)
    items = list(Sp4)
    for _ in range(64):
        g, h, k = [choice(items) for _ in range(3)]
        #print(cocyc(g, h), end=" ")
        lhs = cocyc(g, h*k) + cocyc(h, k)
        rhs = cocyc(g*h, k) + cocyc(g, h)
        assert lhs%2 == rhs%2
        print("%s=%s"%(lhs%2,rhs%2), end=" ")


def test_extension_1():

    n = 1
    c = Clifford(n)
    I = c.I

    gen = [w4*I, c.X(), c.Z()]
    Pauli = mulclose(gen)
    Pauli = list(Pauli)
    assert len(Pauli) == 16

    gen = [c.I, c.S(), c.H()]
    names = "ISH"
    Cliff = mulclose(gen)
    assert len(Cliff) == 192

    lookup = {g:i for (i,g) in enumerate(Pauli)}
    perms = []
    for g in gen:
        ig = g.inverse()
        idxs = [lookup[ig*h*g] for h in Pauli]
        perm = Perm(idxs)
        perms.append(perm)

    #s = gap_code(perms)
    #print(s)
    # gap> G := Group((2, 5, 3, 15)(8, 10, 12, 13), 
    # (2, 3)(4, 13)(5, 16)(6, 15)(8, 12)(10, 11));
    # gap> StructureDescription(G);
    # "S4"

    #hom = mulclose_hom(gen, perms)
    APs = Group.generate(perms)
    assert len(APs) == 24

    lookup = mulclose_names(perms, names)
    lookup = {k:''.join(v) for (k,v) in lookup.items()}

    S4 = Group.symmetric(4)
    iso = S4.isomorphism(APs)
    for g in iso:
        h = iso[g]
        name = lookup[h]
        name = name.replace("SS", "Z")
        name = name.replace("HZH", "X")
        #print(g, "-->", name)

    for trial in range(10):
        iso = S4.isomorphism(APs)
        for g in iso:
          for h in iso:
            assert iso[g*h] == iso[g]*iso[h]
        assert len(set(iso.values())) == 24

    a = Perm([1,0,3,2])
    b = Perm([3,2,1,0])
    assert a in S4
    assert b in S4

    H = Group.generate([a,b])
    assert len(H) == 4
    assert S4.is_subgroup(H)
    assert S4.is_normal(H)

    e = S4.identity
    pairs = [(a,b) for a in S4 for b in S4]
    for (a,b) in pairs:
        if a*a!=e or b*b!=e or a*b == b*a:
            continue
        K = Group.generate([a,b])
        if len(K) != 6:
            continue
        HK = set(h*k for h in H for k in K)
        if len(HK) != 24:
            continue

        print(a)
        print(b)

        HK = [(h,k) for h in H for k in K]
        for (h,k) in HK:
          for (h1,k1) in HK:
            assert h*k*h1*k1 == h*k*h1*~k*k*k1
        break
    else:
        assert 0
        
    return

    I, S, H = perms
    Z = S*S
    X = H*Z*H

    assert X*Z == Z*X # life without phases
    A = Group.generate([X, Z])
    assert len(A) == 2**2

    X = APs.action_subgroup(A)

    # Z2xZ2 >--> APs==S_4 -->> Sp(2,2) == S_3
    #   A   >-->   APs    -->>   tgt

    src, tgt = X.src, X.tgt
    send = X.send_perms # argh
    proj = {}
    for i,g in enumerate(src):
        h = tgt[send[i]]
        proj[g] = h

    #print("section:")
    section = {}
    for g in src:
        pg = proj[g]
        if pg in section:
            h = section[pg]
            if len(lookup[g]) < len(lookup[h]):
                section[pg] = g
        else:
            section[pg] = g
        #print(proj[g], "-->", lookup[g])
        # the proj is a 4-to-1 mapping

    #   A   >-->   APs    -->>   tgt
    #   4   >-->    24    -->>     6

    #print("section:")
    for k,v in section.items():
        assert v in APs
        #print(k, "-->", lookup[v])

    assert len(X.tgt) == 6 # Sp(2,2) == GL(2,2) == S_3
    #print(X.tgt, A, X.tgt is A)
    
    cocycle = lambda g,h : section[g]*section[h]*~section[g*h]

    alookup = {"I":"I", "SS":"Z", "HSSH":"X", "SSHSSH":"XZ"}

    for i,g in enumerate(tgt):
      for j,h in enumerate(tgt):
        u = cocycle(g, h)
        #print("c(%s, %s) = %s" % (
        #    lookup[section[g]], 
        #    lookup[section[h]], 
        #    alookup[lookup[u]]))

    smap = SMap()
    items = list(tgt)
    items.sort()
    for i,g in enumerate(items):
      smap[0, 4+4*i] = lookup[section[g]]
      smap[1, 4+4*i] = '----'
      smap[2+i, 0] = lookup[section[g]]
      smap[2+i, 3] = "|"
      for j,h in enumerate(items):
        u = cocycle(g, h)
        s = alookup[lookup[u]]
        smap[2+i,4+4*j] = s

    print(smap)
    
    for g in tgt:
      for h in tgt:
        u = cocycle(g, h)
        assert proj[u].is_identity()
        #print(u)
        k = choice(tgt)
        a, b = cocycle(g, h*k), cocycle(h,k)
        b = section[g]*b*~section[g]
        assert a*b == b*a
        c, d = cocycle(g*h, k), cocycle(g, h)
        assert c*d == d*c
        #print(int(a*b == c*d), end="")
        assert a*b == c*d

    print()

    # dict'ify this
    pairs = [(g,h)  for g in tgt for h in tgt]
    cocycle = {(g,h) : cocycle(g,h) for (g,h) in pairs}

    # phi: tgt --> A
    N = len(tgt)
    M = len(A)
    print( "phi:", tgt, "-->", A )
    count = 0
    for orig in cross( (A,) * N ):
        phi = {tgt[i]:orig[i] for i in range(N) }
        dphi = {}
        for (g,h) in pairs:
            dphi[g,h] = section[g]*phi[h]*~section[g] * ~phi[g*h] * phi[g]
        count += 1
        if dphi!=cocycle:
            continue
        print("found phi:")
        for i in range(N):
            o = orig[i]
            print("\t", lookup[section[tgt[i]]], "-->", alookup[lookup[o]])
    assert count == M**N



def test_extension():

    c2 = Clifford(2)
    II = c2.I
    XI = c2.X(0)
    IX = c2.X(1)
    ZI = c2.Z(0)
    IZ = c2.Z(1)
    SI = c2.S(0)
    IS = c2.S(1)
    HI = c2.H(0)
    IH = c2.H(1)
    CX01 = c2.CX(0, 1)
    CX10 = c2.CX(1, 0)
    CZ = c2.CZ()
    SWAP = c2.SWAP(0, 1)

    gen = [XI, IX, ZI, IZ, w4*II]
    Pauli = mulclose(gen)
    Pauli = list(Pauli)
    assert len(Pauli) == 64

    gen = [SI, IS, HI, IH, CZ]
    lookup = {g:i for (i,g) in enumerate(Pauli)}
    perms = []
    for g in gen:
        ig = g.inverse()
        idxs = [lookup[ig*h*g] for h in Pauli]
        perm = Perm(idxs)
        perms.append(perm)

    s = gap_code(perms)
    print(s)

    return

    #hom = mulclose_hom(gen, perms)
    APs = Group.generate(perms)
    assert len(APs) == 11520

    #return

    SI, IS, HI, IH, CZ = perms
    ZI = SI*SI
    IZ = IS*IS
    XI = HI*ZI*HI
    IX = IH*IZ*IH

    assert XI*ZI == ZI*XI # life without phases
    A = Group.generate([XI, IX, ZI, IZ])
    assert len(A) == 2**4

    X = APs.action_subgroup(A)

    src, tgt = X.src, X.tgt
    send = X.send_perms # argh
    proj = {}
    for i,g in enumerate(src):
        h = tgt[send[i]]
        proj[g] = h
    section = {}
    for g in src:
        section[proj[g]] = g
    
    cocycle = lambda g,h : section[g]*section[h]*~section[g*h]
    
    for trial in range(100):
        g = choice(tgt)
        h = choice(tgt)
        u = cocycle(g, h)
        assert proj[u].is_identity()
        #print(u)
        k = choice(tgt)
        a, b = cocycle(g, h*k), cocycle(h,k)
        b = section[g]*b*~section[g]
        assert a*b == b*a
        c, d = cocycle(g*h, k), cocycle(g, h)
        assert c*d == d*c
        print(int(a*b == c*d), end="")
        assert a*b == c*d

    print()
    return X


def test_CCZ():

    I = Clifford(1).I

    c2 = Clifford(2)
    wI = c2.w()
    II = c2.I
    XI = c2.X(0)
    IX = c2.X(1)
    ZI = c2.Z(0)
    IZ = c2.Z(1)
    SI = c2.S(0)
    IS = c2.S(1)
    HI = c2.H(0)
    IH = c2.H(1)
    CX01 = c2.CX(0, 1)
    CX10 = c2.CX(1, 0)
    CZ01 = c2.CZ(0, 1)
    CZ10 = c2.CZ(1, 0)
    SWAP01 = c2.SWAP(0, 1)

    assert HI*HI == II
    assert HI*XI*HI == ZI
    assert HI*ZI*HI == XI
    assert HI*IX*HI == IX
    assert SWAP01 == CX01*CX10*CX01

    #print( CX01 )
    #print( CX10 )
    #print( CZ01 )
    #print( IH * CZ01 * IH )
    #print( HI * CZ01 * HI )

    c3 = Clifford(3)
    wI = c3.w()
    III = c3.I
    XII = c3.X(0)
    IXI = c3.X(1)
    IIX = c3.X(2)
    ZII = c3.Z(0)
    IZI = c3.Z(1)
    IIZ = c3.Z(2)
    SII = c3.S(0)
    ISI = c3.S(1)
    IIS = c3.S(2)
    HII = c3.H(0)
    IHI = c3.H(1)
    IIH = c3.H(2)
    CX01 = c3.CX(0, 1)
    CX10 = c3.CX(1, 0)
    CX02 = c3.CX(0, 2)
    CX20 = c3.CX(2, 0)
    CX12 = c3.CX(1, 2)
    CX21 = c3.CX(2, 1)
    CX12 = c3.CX(1, 2)
    CZ01 = c3.CZ(0, 1)
    CZ02 = c3.CZ(0, 2)
    CZ12 = c3.CZ(1, 2)
    SWAP01 = c3.SWAP(0, 1)
    SWAP02 = c3.SWAP(0, 2)
    SWAP12 = c3.SWAP(1, 2)

    assert HII*HII == III
    assert HII*XII*HII == ZII
    assert HII*ZII*HII == XII
    assert HII*IXI*HII == IXI

    assert CZ01 == c2.CZ(0,1)@I
    assert SWAP01 == CX01*CX10*CX01

    assert CX01*CX02 == CX02*CX01
    assert CX12*CX01 == CX02*CX01*CX12

    N = 2**c3.n
    rows = []
    for i in range(N):
        row = [0]*N
        row[i] = -1 if i==N-1 else 1
        rows.append(row)
    CCZ = Matrix(c3.K, rows)

    g = CCZ * XII*IXI*IIX * CCZ.inverse()
    print("g =")
    print(g)

    half = c3.K.one() / 2
    op = half*(IIX + ZII*IIX + IZI*IIX - ZII*IZI*IIX)

    ns =locals()

    #names = "wI SII ISI IIS HII IHI IIH CZ01 CZ02 CZ12".split()
    names = "XII IXI IIX CZ01 CZ02 CZ12".split()
    gen = [ns[name] for name in names]
    name = mulclose_find(gen, names, g, verbose=True)
    print(name)


def test_su2():
    K = CyclotomicField(8)
    w = K.gen()
    one = K.one()
    w2 = w*w
    r2 = w+w.conjugate()
    ir2 = r2 / 2
    
    I = Matrix(K, [[1, 0], [0, 1]])
    w2I = Matrix(K, [[w2, 0], [0, w2]])
    S = Matrix(K, [[1, 0], [0, w2]])
    X = Matrix(K, [[0, 1], [1, 0]])
    Z = Matrix(K, [[1, 0], [0, -1]])
    H = Matrix(K, [[ir2, ir2], [ir2, -ir2]])
    
    Pauli1 = mulclose([w2I, X, Z])
    Cliff1 = mulclose([w2I, S, H])
    
    assert len(Pauli1) == 16
    assert len(Cliff1) == 192

    G = [g for g in Cliff1 if g.determinant() == 1]
    print(len(G))


    
def test():
    test_clifford()
    test_clifford3()
    test_bruhat()
    #test_cocycle() # sloooow
    test_CCZ()



if __name__ == "__main__":

    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "test"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()


