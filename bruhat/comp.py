#!/usr/bin/env python

"""
Chain complexes.

rewrite of chain.py, reversing order of Chain lins

"""


from time import time
start_time = time()

from functools import reduce
from operator import matmul, mul

import numpy
#where = lambda u : numpy.where(u)[0]
def where(u):
    if isinstance(u, Lin):
        u = u.A
    return numpy.where(u==1)[0]

def uwhere(u): # unique where
    idxs = where(u)
    assert len(idxs)==1
    return idxs[0]

from bruhat import elim
from bruhat.lin import Lin, Space, AddSpace, MulSpace, element
from bruhat.smap import SMap
from bruhat.argv import argv
from bruhat.util import cross, allperms, distinct
from bruhat.action import Perm, Group


def shortstr(A):
    s = elim.shortstr(A)
    s = str(s)
    s = s.replace(" 0 ", " . ")
    if len(A.shape)>1 and A.shape[1] == 1:
        s = s.replace("\n", "")
    return s


def xstr(A):
    s = shortstr(A)
    for c in " [],":
        s = s.replace(c, "")
    s = s.replace("1", "X")
    return s

def zstr(A):
    s = shortstr(A)
    for c in " [],":
        s = s.replace(c, "")
    s = s.replace("1", "Z")
    return s


# does not need hashable operators
def mulclose(gen, verbose=False, maxsize=None):
    if verbose:
        print("mulclose")
    ops = list(gen)
    bdy = gen
    while bdy:
        _bdy = []
        for g in bdy:
            for h in gen:
                k = g*h
                if k not in ops:
                    ops.append(k)
                    _bdy.append(k)
                    if maxsize and len(ops) >= maxsize:
                        return ops # <---------- return
        bdy = _bdy
        if verbose:
            print("mulclose:", len(ops))
    return ops


class Seq(object):
    def __init__(self, lins):
        assert lins
        self.ring = lins[0].ring
        self.lins = list(lins)

    def __len__(self):
        return len(self.lins)

    def __getitem__(self, idx):
        return self.lins[idx]





class Chain(Seq):

    def __init__(self, lins):
        Seq.__init__(self, lins)
        prev = None
        grade = None
        for idx, lin in enumerate(lins):
            tgt, src = lin.hom
            if grade is None:
                grade = tgt.grade
            assert grade == tgt.grade
            grade += 1
            assert grade == src.grade
            assert prev is None or prev.src == lin.tgt
            prev = lin

    def get(self, idx):
        if idx < len(self):
            return self[idx].tgt
        assert idx == len(self)
        return self[idx-1].src

    def identity(chain):
        n = len(chain)
        lins = [chain.get(i).identity() for i in range(n+1)]
        return ChainMap(chain, chain, lins)

    def __str__(self):
        #return str([lin.homstr() for lin in self])
        spaces = [lin.tgt for lin in self.lins] + [self.lins[-1].src]
        return "%s(%s)"%(self.__class__.__name__,
            "<--".join(str(space) for space in spaces))

    def __eq__(self, other):
        if len(self.lins) != len(other.lins):
            return False
        for lhs, rhs in zip(self.lins, other.lins):
            if lhs != rhs:
                return False
        return True

    def sub(self, start, stop):
        lins = self.lins[start : stop]
        return Chain(lins)

    @classmethod
    def promote(cls, item):
        if isinstance(item, Chain):
            return item
        if isinstance(item, Lin):
            return Chain([item])
        assert 0, "what's this: %s ?"%item

    def longstr(self):
        lins = self.lins
        ss = [elim.shortstr(lin.A) for lin in lins]
        rows = max([s.rows for s in ss])
        #print("rows", rows, [s.rows for s in ss])
        smap = SMap()
        col = 2
        for s, lin in zip(ss, lins):
            tgt, src = lin.hom
            smap.paste(s, rows-s.rows, col)
            col1 = col + s.cols + 4
            for i in range(col, col1-2):
                smap[rows+1, i] = '-'
            smap[rows+1, i] = '>'
            smap[rows+1, col-2] = str(src.n)
            col = col1 + 1
        smap[rows+1, col-2] = str(tgt.n)
        result = str(smap)
        result = result.replace(" 0 ", " . ")
        return result

    cache = {} # TODO
    @classmethod
    def product(cls, ring, chains):
        key = tuple(id(f) for f in chains)
        if key in cls.cache:
            return cls.cache[key]

        chains = [Chain.promote(f) for f in chains]
        for chain in chains:
            assert len(chain) == 1, "not implemented"
        fs = [chain[0] for chain in chains]

        N = len(chains)
        grades = [[] for i in range(N+1)]
        #slookup = {}
        #tlookup = {}
        for items in cross([f.hom for f in fs]):
            space = MulSpace(ring, *items)
            grades[space.grade].append(space)
        grades = [AddSpace(ring, *spaces) for spaces in grades]
        #print([space.name for space in grades])
    
        srcs = [[] for i in range(N+1)]
        tgts = [[] for i in range(N+1)]
        items = [(f.src.identity(), f.tgt.identity()) for f in fs]
        for idx in range(N):
            _items = list(items)
            _items[idx] = (fs[idx],)
            for lins in cross(_items):
                lin = reduce(matmul, lins)
                srcs[lin.src.grade].append(lin)
                tgts[lin.tgt.grade].append(lin)
                #print([lin.homstr() for lin in lins])
                #print(lin.homstr())
                #print()
        
        #print([len(grade) for grade in srcs])
    
        chain = []
        for i in range(N):
            tgt = grades[i]
            src = grades[i+1]
            #print([space.name for space in tgt.get_add_items()])
            #print([space.name for space in src.get_add_items()])
            A = elim.zeros(ring, tgt.n, src.n)
            for lin in tgts[i]:
                #print("\t", lin.homstr())
                cols = src.get_slice(lin.src)
                rows = tgt.get_slice(lin.tgt)
                assert numpy.alltrue(A[rows, cols] == ring.zero)
                A[rows, cols] = lin.A
            lin = Lin(tgt, src, A)
            #print(shortstr(A))
            #print()
            chain.append(lin)
    
        chain = Chain(chain)
        cls.cache[key] = chain # doesn't help ...??? XXX profile
        return chain
    


class ChainMap(Seq):
    CHECK = True
    def __init__(self, tgt, src, lins):
        tgt = Chain.promote(tgt)
        src = Chain.promote(src)
        assert isinstance(tgt, Chain)
        assert isinstance(src, Chain)
        n = len(src)
        assert len(tgt) == n
        assert len(lins) == n+1
        self.ring = tgt.ring
        self.tgt = tgt
        self.src = src
        self.hom = (tgt, src) # yes it's backwards, just like shape is.
        Seq.__init__(self, lins)
        if self.CHECK:
            self.check()

    def check(self):
        tgt, src = self.hom
        lins = self.lins
        n = len(src)
        for i in range(n):
            lhs = lins[i] * src[i]
            rhs = tgt[i] * lins[i+1]
            assert lhs == rhs, "not a chain map"

    def __eq__(self, other):
        assert self.hom == other.hom
        assert len(self) == len(other)
        for lhs, rhs in zip(self, other):
            if lhs != rhs:
                return False
        return True

    def sub(self, start, stop):
        tgt = self.tgt.sub(start, stop)
        src = self.src.sub(start, stop)
        lins = self.lins[start : stop+1]
        return ChainMap(tgt, src, lins)

    def __mul__(self, other):
        assert isinstance(other, ChainMap)
        assert other.tgt == self.src
        lins = [l*r for (l,r) in zip(self, other)]
        return ChainMap(self.tgt, other.src, lins)

    def add(self, other):
        assert isinstance(other, ChainMap)
        assert other.hom == self.hom
        lins = [l.add(r) for (l,r) in zip(self, other)]
        return ChainMap(self.tgt, self.src, lins)

    @classmethod
    def zero(cls, tgt, src):
        n = len(src)
        assert len(tgt) == n
        # zero map
        lins = [Lin(tgt.get(i), src.get(i)) for i in range(n+1)]
        return ChainMap(tgt, src, lins)

    @classmethod
    def product(cls, ring, cmaps): # tensor product
        for cmap in cmaps:
            assert isinstance(cmap, ChainMap)
            assert len(cmap) == 2
            f0, f1 = cmap
            #print(f0.homstr(), f1.homstr())
        N = len(cmaps)
    
        linss = [[] for i in range(N+1)]
        for items in cross(cmaps):
            #print([f.homstr() for f in items])
            f = reduce(matmul, items)
            #print(f.homstr())
            assert f.src.grade == f.tgt.grade
            linss[f.src.grade].append(f)
        cmap = [reduce(Lin.direct_sum, lins) for lins in linss]
        #print([lin.homstr() for lin in chain])

        tgt = Chain.product(ring, [item.tgt for item in cmaps])
        src = Chain.product(ring, [item.src for item in cmaps])
        cmap = ChainMap(tgt, src, cmap)

        return cmap

    @classmethod
    def perm(cls, ring, fs, perm): # perm'ute tensor factors
        #print("ChainMap.perm")

        fs = [(f[0] if isinstance(f, Chain) else f) for f in fs]
                
        N = len(fs)
        assert N == len(perm)
        gs = [fs[i] for i in perm]
    
        def product(fs):
            grades = [[] for i in range(N+1)]
            for items in cross([f.hom for f in fs]):
                space = MulSpace(ring, *items)
                grades[space.grade].append(space)
            grades = [AddSpace(ring, *spaces) for spaces in grades]
            return grades
    
        sgrades = product(fs)
        tgrades = product(gs)
    
        cmap = []
        for i in range(N+1):
            tgt, src = tgrades[i], sgrades[i]
            #print(tgt.name, "<---", src.name)
            tgts = tgt.get_add_items()
            srcs = src.get_add_items()
            lins = [src.get_mul_swap(perm) for src in srcs]
    
            f = reduce(Lin.direct_sum, lins)
            #print(f.homstr())
            assert f.src == src
    
            swap = [None]*len(lins)
            for i, lin in enumerate(lins):
                assert tgts.count(lin.tgt) == 1 # careful !
                idx = tgts.index(lin.tgt)
                swap[idx] = i
            #print("swap:", swap)
            swap = f.tgt.get_add_swap(swap)
            f = swap * f
    
            #print(f.homstr())
    
            assert f.hom == (tgt, src)
            cmap.append(f)
            #print()

        src = Chain.product(ring, fs)
        tgt = Chain.product(ring, [fs[i] for i in perm])
        cmap = cls(tgt, src, cmap)
    
        return cmap


def test_chain_symmetry(ring, lins):
    N = len(lins)
    src = Chain.product(ring, lins)
    #print([lin.homstr() for lin in src])
    for perm in allperms(list(range(N))):
        #print("perm:", perm)
        tgt = Chain.product(ring, [lins[i] for i in perm])
        #print([lin.homstr() for lin in tgt])
        cmap = ChainMap.perm(ring, lins, perm)


def test_chain_product():

    p = 2
    ring = element.FiniteField(p)

    C1 = Space(ring, 1, 1, "C_1")
    C0 = Space(ring, 1, 0, "C_0")
    f = Lin(C0, C1, [[1]])

    chain = Chain.product(ring, [f, f])

    D1 = Space(ring, 1, 1, "D_1")
    D0 = Space(ring, 1, 0, "D_0")
    g = Lin(D0, D1, [[1]])

    E1 = Space(ring, 1, 1, "E_1")
    E0 = Space(ring, 1, 0, "E_0")
    h = Lin(E0, E1, [[1]])

    test_chain_symmetry(ring, [f, g])
    test_chain_symmetry(ring, [f, g, h])
    test_chain_symmetry(ring, [f, f])
    test_chain_symmetry(ring, [f, f, h])
    test_chain_symmetry(ring, [f, f, f])


    C1 = Space(ring, 3, 1, "C_1")
    C0 = Space(ring, 2, 0, "C_0")
    f = Lin(C0, C1, [[1,1,0],[0,1,1]])

    chain = Chain.product(ring, [f, f, f])
    #test_chain_symmetry(ring, [f, f, f])

    C1 = Space(ring, 3, 1, "C_1")
    C0 = Space(ring, 3, 0, "C_0")
    f = Lin(C0, C1, [[1,1,0],[0,1,1],[1,0,1]])

    D1 = Space(ring, 3, 1, "D_1")
    D0 = Space(ring, 3, 0, "D_0")
    g = Lin(D0, D1, [[1,1,0],[0,1,1],[1,0,1]])

    cmap = [
        Lin(D0, C0, [[0,1,0],[0,0,1],[1,0,0]]),
        Lin(D1, C1, [[0,1,0],[0,0,1],[1,0,0]]),
    ]

    cmap = ChainMap(g, f, cmap)

    ident = [
        Lin(D0, C0, [[1,0,0],[0,1,0],[0,0,1]]),
        Lin(D1, C1, [[1,0,0],[0,1,0],[0,0,1]]),
    ]
    ident = ChainMap(g, f, ident)

    src = Chain.product(ring, [f, f, f])
    tgt = Chain.product(ring, [g, g, g])

    c3 = ChainMap.product(ring, [cmap, cmap, cmap])
    c3 = ChainMap.product(ring, [cmap, cmap, ident])


def get_point(tgt, idx=0):
    ring = tgt.ring
    src = Space(ring, 1, tgt.grade, "src")
    A = elim.zeros(ring, tgt.n, src.n)
    A[idx,0] = ring.one
    return Lin(tgt, src, A)


def test_kagome():
    print("test_kagome")

    if 1:
        p = 2
        ring = element.FiniteField(p)
    else:
        ring = element.Z

    dims = argv.get("dims", 3)
    L = argv.get("L", 3)


    C1 = Space(ring, L, 1, "C_1")
    C0 = Space(ring, L, 0, "C_0")
    f = elim.zeros(ring, L, L)
    for i in range(L):
        f[i,i] = ring.one
        f[i,(i+1)%L] = ring.one
    f = Lin(C0, C1, f)
    f = Chain([f])

    A = elim.zeros(ring, L, L)
    for i in range(L):
        A[i,(i+1)%L] = ring.one

    send = ChainMap(f, f, [ Lin(C0, C0, A), Lin(C1, C1, A), ])
    inv = ChainMap(f, f, [lin.transpose() for lin in send])

    lins = [f]*dims
    toric = Chain.product(ring, lins)
    perm = tuple((i+1)%dims for i in range(dims))
    rotate = ChainMap.perm(ring, lins, perm)
    assert rotate.src == rotate.tgt

    ChainMap.CHECK = argv.get("CHECK", False)

    print(toric)
    #grade = argv.get("grade", dims//2)
    grade = argv.get("grade", 1)
    print("grade:", grade)
    Hx, Hzt = toric[grade-1 : grade+1]
    #Hx, Hzt = Hx.A, Hzt.A
    Hz = Hzt.transpose()
    print("Hx:")
    print(Hx.shape) #, shortstr(Hx.sum(0)), shortstr(Hx.sum(1)))
    print("Hz:")
    print(Hz.shape) #, shortstr(Hz.sum(0)), shortstr(Hz.sum(1)))

    i = f.identity()
    assert i*i == i
    assert send != i
    assert send*send != i
    assert reduce(mul, [send]*L) == i

    assert send != inv
    assert inv != i
    assert inv*inv != i
    assert reduce(mul, [inv]*L) == i

    assert send*inv == i
    assert inv*send == i

    I = ChainMap.product(ring, [i]*dims)
    gens = []
    for j in range(dims):
        g = [i]*dims
        g[j] = send
        g[(j+1)%dims] = inv
        g = ChainMap.product(ring, g)
        gens.append(g)

        assert g.src == g.tgt

        #for count in range(1, L+1):
        #    op = reduce(mul, [g]*count)
        #    assert (op == I) == (count==L)

    if dims == 3:
        G = []
        a, b = gens[:2]
        for ia in range(L):
         for ib in range(L):
            ma = reduce(mul, [a]*ia, I)
            mb = reduce(mul, [b]*ib, I)
            G.append(ma*mb)
        rr = rotate*rotate
        G += [rotate*g for g in G] + [rr*g for g in G]

    elif dims == 4:
        G = []
        a, b, c = gens[:3]
        for ia in range(L):
         for ib in range(L):
          for ic in range(L):
            ma = reduce(mul, [a]*ia, I)
            mb = reduce(mul, [b]*ib, I)
            mc = reduce(mul, [c]*ic, I)
            G.append(ma*mb*mc)
        rr = rotate*rotate
        rrr = rr*rotate
        G += [rotate*g for g in G] + [rr*g for g in G] + [rrr*g for g in G]

    else:
        assert 0

#    if dims==3 and L==3:
#        a, b, c = gens[1:]
#        G = [
#            I, a, a*a, 
#            b, a*b, a*a*b, 
#            b*b, a*b*b, a*a*b*b, 
#        ]
#        rr = rotate*rotate
#        G += [rotate*g for g in G] + [rr*g for g in G]
#    else:
#        G = mulclose(gens, verbose=True)

    print("|G| =", len(G))

    if argv.check:
        assert distinct(G)
    assert( len(G) == dims * (L**(dims-1)) ) # assert 

    # Hx: weight 6 checks
    # Hz: weight 4 checks
    Hx = Hx.A
    Hz = Hz.A

    tgt = toric
    Cn = toric.get(grade)
    point = get_point(Cn)

    #points = [g[1]*point for g in gens]
    points = []
    for g in G:
        A = (g[grade]*point).A
        A.shape = len(A),
        points.append(A)
    points = elim.array(points)
    print("points:")
    print(points.shape)
    zstabs = elim.intersect(ring, points, Hz)
    print("zstabs:")
    print(zstabs.shape)
    print(shortstr(zstabs))
    Hz1t = elim.dot(ring, points, zstabs.transpose())
    Hz1 = Hz1t.transpose()
    print("Hz1t:")
    print(Hz1t.shape)
    print(shortstr(Hz1t))

    cmap1 = points.transpose()

    inv = elim.pseudo_inverse(ring, Hz.transpose())
    cmap2 = elim.dot(ring, inv, zstabs.transpose())
    
    lhs = elim.dot(ring, cmap1, Hz1t)
    rhs = elim.dot(ring, Hz.transpose(), cmap2)
    assert elim.eq(lhs, rhs)

    Hx1 = elim.dot(ring, Hx, cmap1)
    print("Hx1:")
    print(Hx1.shape)
    print(shortstr(Hx1))
    H = (Hx1==1).astype(int)
    idxs = []
    cmap0 = []
    for idx, row in enumerate(H):
        if row.sum() == 0:
            continue
        idxs.append(idx)
        v = [ring.zero]*H.shape[0]
        v[idx] = ring.one
        cmap0.append(v)

    cmap0 = elim.array(cmap0).transpose()
    print("cmap0:")
    print(cmap0.shape)
    print(shortstr(cmap0))
            
    Hx1 = Hx1[idxs, :]
    print("Hx1:")
    print(Hx1.shape)
    print(shortstr(Hx1))

    lhs = elim.dot(ring, cmap0, Hx1)
    rhs = elim.dot(ring, Hx, cmap1)
    assert elim.eq(lhs, rhs)

    lhs = elim.dot(ring, cmap1, Hx1.transpose())
    rhs = elim.dot(ring, Hx.transpose(), cmap0)
    #print(elim.eq(lhs, rhs)) # False

    lhs = elim.dot(ring, cmap2, Hz1)
    rhs = elim.dot(ring, Hz, cmap1)
    #print(elim.eq(lhs, rhs)) # False

    return locals()

    return

    N0 = Space(ring, 0, 0, "N0")
    N1 = Space(ring, 0, 1, "N1")
    N2 = Space(ring, 0, 2, "N2")
    N3 = Space(ring, 0, 3, "N3")

    #toric = toric.sub(0, 2)
    src = Chain([Lin.zero(N0, point.src), Lin.zero(point.src, N2)])
    point = ChainMap(toric, src, [
        #Lin.zero(toric.get(0), src.get(0)),
        toric[0] * point,
        point,
        #Lin.zero(toric.get(1), src.get(1)),
        Lin.zero(toric.get(2), src.get(2)),
    ])

    return

    if 0: # SLOW
        G = mulclose([a, b], verbose=True, maxsize=9)
        print("|G| =", len(G))
        assert len(G) == 9
    
        for g in G:
          for h in G:
            assert g*h == h*g
    

def test_kagome_floquet():
    print("test_kagome_floquet")

    if 1:
        p = 2
        ring = element.FiniteField(p)
    else:
        ring = element.Z

    dims = argv.get("dims", 3)
    L = argv.get("L", 3)

    dot = lambda *args : elim.dot(ring, *args)

    C1 = Space(ring, L, 1, "C_1")
    C0 = Space(ring, L, 0, "C_0")
    f = elim.zeros(ring, L, L)
    for i in range(L):
        f[i,i] = ring.one
        f[i,(i+1)%L] = ring.one
    f = Lin(C0, C1, f)
    f = Chain([f])

    A = elim.zeros(ring, L, L)
    for i in range(L):
        A[i,(i+1)%L] = ring.one

    send = ChainMap(f, f, [ Lin(C0, C0, A), Lin(C1, C1, A), ])
    inv = ChainMap(f, f, [lin.transpose() for lin in send])

    lins = [f]*dims
    toric = Chain.product(ring, lins)
    perm = tuple((i+1)%dims for i in range(dims))
    rotate = ChainMap.perm(ring, lins, perm)
    assert rotate.src == rotate.tgt

    ChainMap.CHECK = argv.get("CHECK", False)

    print(toric)
    #grade = argv.get("grade", dims//2)
    grade = argv.get("grade", 1)
    print("grade:", grade)

    i = f.identity()
    assert i*i == i
    assert send != i
    assert send*send != i
    assert reduce(mul, [send]*L) == i

    assert send != inv
    assert inv != i
    assert inv*inv != i
    assert reduce(mul, [inv]*L) == i

    assert send*inv == i
    assert inv*send == i

    Cn = toric.get(grade)
    n = Cn.n

    I = ChainMap.product(ring, [i]*dims)

    if 0:
        # this group is _transitive on the cubes
        gens = []
        for j in range(dims):
            g = [i]*dims
            g[j] = send
            g = ChainMap.product(ring, g)
            gens.append(g)
        G = mulclose(gens, verbose=True)
        assert len(G) == Hz.shape[1] // 3

    # translate in Z direction
    TZ = ChainMap.product(ring, [i, i, send])

    Hx, Hzt = toric[grade-1 : grade+1]
    Hz = Hzt.transpose()

    if 0:
        for i in range(n):
            p0 = get_point(Cn, i)
            assert uwhere(p0) == i
            p1 = TZ[grade] * p0
            #p1 = rotate[grade] * p0
            j = uwhere(p1)
            #hz = Hx[:, [i,j]]
            hz = Hz[:, [i,j]]
            for row in hz:
                if (row==1).sum()==2:
                    print(shortstr(row), "%d --> %d"%(i,j))
    
        return
    
    # this group is _transitive on the qubits (edges)
    gens = [rotate, rotate*rotate]
    for j in range(dims):
        g = [i]*dims
        g[j] = send
        g = ChainMap.product(ring, g)
        gens.append(g)
    gens = [g[grade] for g in gens]
    G0 = mulclose(gens, maxsize=n, verbose=True)
    assert len(G0) == n
    print("|G0| =", len(G0))

    # build a lookup
    lookup = {} # idx -> g such that g[idx] == 0
    for g in G0:
        #assert g * g.transpose() == Cn.identity()
        gi = g.transpose()
        point = get_point(Cn)
        point = gi*point
        idx = uwhere(point.A==ring.one)
        assert idx not in lookup
        lookup[idx] = g
        assert 0 == uwhere((g*get_point(Cn, idx)).A==1)
        #print(idx, end=" ")
    #print()

    perm = {}
    found = set()
    for idx in range(n):
        g = lookup[idx]
        g = g.transpose() * TZ[grade] * rotate[grade] * g
        point = get_point(Cn, idx)
        point = g*point
        jdx = uwhere(point)
        #print(jdx, end=" ")
        perm[idx] = jdx
        assert jdx not in found
        found.add(jdx)
    assert perm[0] == 29, perm
    assert perm[27] == 60, perm
    assert perm[54] == 18, perm
    assert perm == {0: 29, 1: 27, 2: 28, 3: 32, 4: 30, 5:
        31, 6: 35, 7: 33, 8: 34, 9: 38, 10: 36, 11: 37, 12: 41,
        13: 39, 14: 40, 15: 44, 16: 42, 17: 43, 18: 47, 19: 45,
        20: 46, 21: 50, 22: 48, 23: 49, 24: 53, 25: 51, 26: 52,
        27: 60, 28: 61, 29: 62, 30: 54, 31: 55, 32: 56, 33: 57,
        34: 58, 35: 59, 36: 69, 37: 70, 38: 71, 39: 63, 40: 64,
        41: 65, 42: 66, 43: 67, 44: 68, 45: 78, 46: 79, 47: 80,
        48: 72, 49: 73, 50: 74, 51: 75, 52: 76, 53: 77, 54: 18,
        55: 19, 56: 20, 57: 21, 58: 22, 59: 23, 60: 24, 61: 25,
        62: 26, 63: 0, 64: 1, 65: 2, 66: 3, 67: 4, 68: 5, 69:
        6, 70: 7, 71: 8, 72: 9, 73: 10, 74: 11, 75: 12, 76: 13,
        77: 14, 78: 15, 79: 16, 80: 17}

    T = Perm(perm, list(range(n)))
    print("T.order() ==", T.order())

    gens = []
    for j in range(dims):
        g = [i]*dims
        g[j] = send
        g[(j+1)%dims] = inv
        g = ChainMap.product(ring, g)
        gens.append(g)

        assert g.src == g.tgt

        #for count in range(1, L+1):
        #    op = reduce(mul, [g]*count)
        #    assert (op == I) == (count==L)

    assert dims==3, "todo"
    #h = [send, send, send] # automorphism of the slice...
    h = [send, send, inv]
    h = ChainMap.product(ring, h)
    assert (h != I)
    assert (h*h != I)
    assert (h*h*h == I)

    if dims == 3:
        G = []
        a, b = gens[:2]
        for ia in range(L):
         for ib in range(L):
            ma = reduce(mul, [a]*ia, I)
            mb = reduce(mul, [b]*ib, I)
            G.append(ma*mb)
        rr = rotate*rotate
        G += [rotate*g for g in G] + [rr*g for g in G]

    elif dims == 4:
        G = []
        a, b, c = gens[:3]
        for ia in range(L):
         for ib in range(L):
          for ic in range(L):
            ma = reduce(mul, [a]*ia, I)
            mb = reduce(mul, [b]*ib, I)
            mc = reduce(mul, [c]*ic, I)
            G.append(ma*mb*mc)
        rr = rotate*rotate
        rrr = rr*rotate
        G += [rotate*g for g in G] + [rr*g for g in G] + [rrr*g for g in G]

    else:
        assert 0

    print("|G| =", len(G))

    if argv.check:
        assert distinct(G)
    assert( len(G) == dims * (L**(dims-1)) ) # assert 

    perm = [perm[i] for i in range(n)]
    T = elim.zeros(ring, n, n)
    for i, j in enumerate(perm):
        T[j, i] = ring.one
    T = Lin(Cn, Cn, T)

    cmaps = []
    codes = []
    T1s = []
    for order in range(3):

        Hx, Hzt = toric[grade-1 : grade+1]
        #Hx, Hzt = Hx.A, Hzt.A
        Hz = Hzt.transpose()

        # Hx: weight 6 checks
        # Hz: weight 4 checks
        Hx = Hx.A.copy()
        Hz = Hz.A.copy()

        T1 = T**order
        #point = get_point(Cn, (T**order)[0])
        #point = (h[grade]**order)*point
        point = T1*get_point(Cn)
    
        #points = [g[1]*point for g in gens]
        points = []
        for g in G:
            A = (g[grade]*point).A
            A.shape = len(A),
            points.append(A)
        points = elim.array(points)
        print("points:", points.shape)
        print(shortstr(points.sum(0)))
        zstabs = elim.intersect(ring, points, Hz)
        print("zstabs:", zstabs.shape)
    
        Hz1t = dot(points, zstabs.transpose())
        Hz1 = Hz1t.transpose()
    
        cmap1 = points.transpose()
        print("cmap1", cmap1.shape)
    
        inv = elim.pseudo_inverse(ring, Hz.transpose())
        cmap2 = dot(inv, zstabs.transpose())
        
        lhs = dot(cmap1, Hz1t)
        rhs = dot(Hz.transpose(), cmap2)
        assert elim.eq(lhs, rhs)
    
        Hx1 = dot(Hx, cmap1)
        print("Hx1:", Hx1.shape)
        H = (Hx1==1).astype(int)
        idxs = []
        cmap0 = []
        for idx, row in enumerate(H):
            if row.sum() == 0:
                continue
            idxs.append(idx)
            v = [ring.zero]*H.shape[0]
            v[idx] = ring.one
            cmap0.append(v)
    
        cmap0 = elim.array(cmap0).transpose()
        print("cmap0:", cmap0.shape)
                
        Hx1 = Hx1[idxs, :]
        print("Hx1:", Hx1.shape)

        lhs = dot(cmap0, Hx1)
        rhs = dot(Hx, cmap1)
        assert elim.eq(lhs, rhs)

        print("Hx1:")
        print(xstr(Hx1))
        print("Hz1:")
        print(zstr(Hz1))
        codes.append((Hx1, Hz1))
        cmaps.append((cmap0, cmap1, cmap2))
        T1s.append(T1)

    return codes, cmaps, toric

    return

    Hzt = toric[grade]
    Hz = Hzt.transpose()
    #print([list(row).count(1) for row in Hz]) # 4
    #print([list(row).count(1) for row in Hzt]) # 4

    if 0:
        for i in range(27):
            print("i =", i)
            idxs = []
            for cmap1 in cmap1s:
                jdx = uwhere(cmap1[:, i])
                idxs.append(jdx)
            hz = Hz[:, idxs]
            for row in hz:
                if (row==1).sum() == 2:
                    print(shortstr(row))
            #print(shortstr(Hz[:, idxs]))
    
        return

    for idx in range(3):
      for jdx in range(idx+1,3):
        l, r = codes[idx][0], codes[jdx][0]
        H = elim.intersect(ring, l, r)
        print(len(H), end=" ")
        l, r = codes[idx][1], codes[jdx][1]
        H = elim.intersect(ring, l, r)
        print(len(H))

    a, b, c = [codes[i][0] for i in range(3)]
    H = elim.intersect(ring, a, b)
    H = elim.intersect(ring, H, c)
    print(len(H))

#    for idx in range(3):
#      for jdx in range(3):
#        H = dot(codes[idx][0], codes[jdx][1].transpose())
#        #print(numpy.alltrue(H==0))
#        if idx != jdx:
#            print(idx, jdx)
#            print(shortstr(H))

def test():

    #test_chain_product()
    test_kagome()



if __name__ == "__main__":


    fn = argv.next() or "test"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK: finished in %.3f seconds\n"%(time() - start_time))



