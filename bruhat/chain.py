#!/usr/bin/env python

"""
see also: kvchain.py, comp.py

"""


from random import choice
from functools import reduce

import numpy

from bruhat import elim
from bruhat.lin import Lin, Space, AddSpace, element
from bruhat.smap import SMap
from bruhat.argv import argv
from bruhat.vect2 import Rig, Cell0, Cell1, Cell2

#shortstr = elim.shortstr

def shortstr(A):
    s = str(elim.shortstr(A))
    s = s.replace('0', '.')
    s = s.replace(' ', '')
    return s

def shortstr(A, deco=False, zero='.'):
    if isinstance(A, Lin):
        A = A.A
    #A = elim.array(A)
    A = A.view()
    assert len(A.shape), A.shape
    assert len(A.shape)<=2
    if 1 in A.shape:
        A.shape = (1, A.shape[0]*A.shape[1])
    if len(A.shape)==1:
        A.shape = (1, A.shape[0])
    m, n = A.shape
    rows = []
    for row in range(m):
        idx = row
        row = list(A[row, :])
        row = ''.join(str(i) for i in row)
        row = row.replace('0', zero)
        if deco:
            row = "%3d "%idx + row
        rows.append(row)
    if deco:
        row = ''.join(str(i%10) for i in range(n))
        rows.append('    '+row)
    return '\n'.join(rows)


def int_array(A):
    if isinstance(A, Lin):
        A = A.A
    B = numpy.zeros(A.shape, dtype=int)
    for idx in numpy.ndindex(A.shape):
        B[idx] = int(str(A[idx]))
    return B


# ------------------------------------------------------------

# ------------------------ graded stuff ----------------------


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
        for lin in lins:
            tgt, src = lin.hom
            if prev is not None:
                assert prev.tgt == src
                assert (lin*prev).is_zero()
            prev = lin

    def get(self, idx):
        lins = self.lins
        n = len(lins)
        if idx < n:
            return lins[idx].src
        elif idx == n:
            return lins[n-1].tgt
        raise IndexError

    def get_grades(self):
        lins = self.lins
        if not lins:
            return [] # ?
        return [lin.src for lin in lins] + [lins[-1].tgt]

    def identity(chain):
        n = len(chain)
        lins = [chain.get(i).identity() for i in range(n+1)]
        return ChainMap(chain, chain, lins)

    def __str__(self):
        spaces = [lin.src for lin in self.lins] + [self.lins[-1].tgt]
        return "%s(%s)"%(self.__class__.__name__,
            "-->".join(str(space) for space in spaces))

    def __eq__(self, other):
        if len(self.lins) != len(other.lins):
            return False
        for lhs, rhs in zip(self.lins, other.lins):
            if lhs != rhs:
                return False
        return True

    def __add__(self, other):
        assert isinstance(other, Chain)
        return AddChain(self, other)

    def __matmul__(self, other):
        assert isinstance(other, Chain)
        return MulChain(self, other) # cache these ?

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

    def encode(self, rig):
        # encode self as a "circulant" 2-vector space matrix
        grades = self.get_grades()
        grades.sort(key = lambda space:space.grade)
        #print("grades:", grades)
        bdys = list(reversed(self))
        n = len(grades)+1
        # f : B <--- A
        A = numpy.empty((n, n), dtype=object)
        B = numpy.empty((n, n), dtype=object)
        f = numpy.empty((n, n), dtype=object)
        for i in range(n):
          for j in range(n):
            A[i, j] = rig.zero
            B[i, j] = rig.zero
        for i in range(n):
          for j in range(i, n):
            idx = j-i
            if idx < len(grades):
                A[i, j] = grades[idx]
            if idx>0:
                B[i, j] = grades[idx-1]
        for i in range(n):
          for j in range(n):
            if 0<=j-i-1<len(bdys):
                f[i,j] = bdys[j-i-1]
                #print("f[i,j]", i, j)
            else:
                f[i,j] = Lin.zero(B[i,j], A[i,j])
        cell = Cell0(rig, n, "n")
        A = Cell1(cell, cell, A)
        #print(A)
        B = Cell1(cell, cell, B)
        #print(B)
        #for i in range(n):
        #  for j in range(n):
        #    print("%12s"%f[i,j].homstr(), end=" ")
        #  print()
        f = Cell2(B, A, f)
        return f
            

class AddChain(Chain):
    "direct sum of Chain complexes"
    def __init__(self, lhs, rhs):
        assert lhs.lins
        assert rhs.lins
        assert lhs.ring is rhs.ring
        ring = lhs.ring

        lins = []
        for l, r in zip(lhs.lins, rhs.lins):
            lin = l.direct_sum(r)
            lins.append(lin)
        Chain.__init__(self, lins)


class MulChain(Chain):
    "tensor product of Chain complexes"
    def _init(self):
        self.spaces = set()  # found spaces
        self.srcs = {}       # Space -> [Lin]
        self.tgts = {}       # Space -> [Lin]
        self.grades = {}     # int -> [Space]

    def _addspace(self, space):
        spaces = self.spaces
        grades = self.grades
        if space not in spaces:
            spaces.add(space)
            grades.setdefault(space.grade, []).append(space)

    def _addlin(self, lin):
        #print("_addlin", lin.homstr())
        tgt, src = lin.hom
        self._addspace(src)
        self.srcs.setdefault(lin.src, []).append(lin)
        self._addspace(tgt)
        self.tgts.setdefault(lin.tgt, []).append(lin)

    def __init__(self, lhs, rhs):
        assert lhs.lins
        assert rhs.lins
        assert lhs.ring is rhs.ring
        ring = lhs.ring
        self._init()
        for g in rhs.lins: # cols
            for f in lhs.lins: # rows
                self._addlin( f @ g.src.identity() ) # vertical arrow
                sign = -1 if f.src.grade % 2 else 1
                self._addlin( sign * f.src.identity() @ g ) # _horizontal arrow
            sign = -1 if f.tgt.grade % 2 else 1
            self._addlin( sign * f.tgt.identity() @ g ) # _horizontal arrow
  
        for f in lhs.lins: # rows
            self._addlin( f @ g.tgt.identity() ) # vertical arrow

        keys = list(self.grades.keys())
        keys.sort(reverse=True)
        #print("grades:", keys)

        N = len(keys)
        lins = []
        for idx in range(N-1):
            i = keys[idx]
            assert keys[idx+1] == i-1, keys
            tgt = AddSpace(ring, *self.grades[i-1])
            src = AddSpace(ring, *self.grades[i])
            #print(tgt, "<------", src)
            A = elim.zeros(ring, tgt.n, src.n)
            #print(shortstr(A))
            items = src.get_add_items()
            for s in items:
                for lin in self.srcs[s]:
                    assert lin.src is s
                    #print("%s.get_slice(%s)"%(src, lin.src))
                    cols = src.get_slice(lin.src)
                    rows = tgt.get_slice(lin.tgt)
                    A[rows, cols] = lin.A
            #print(shortstr(A))
            lin = Lin(tgt, src, A)
            #print(repr(lin))
            lins.append(lin)
            #print()
        self.lhs = lhs
        self.rhs = rhs
        Chain.__init__(self, lins)

    def get_swap(self):
        #print("MulChain.get_swap")
        ring = self.ring
        lookup = self.grades
        idxs = list(lookup.keys())
        idxs.sort(reverse=True)
        #print(idxs)
        items = [] # construct ChainMap
        for grade in idxs:
            #print(grade, lookup[grade])
            src = lookup[grade]
            N = len(src)
            src = AddSpace(ring, *src)
            perm = list(reversed(range(N)))
            swap = src.get_add_swap(perm)
            fs = [space.get_mul_swap() for space in swap.tgt.get_add_items()]
            #print([f.homstr() for f in fs])
            f = reduce(Lin.direct_sum, fs)
            #print("mul swap:", f.homstr())
            f = f*swap
            #print(f.homstr())
            items.append(f)
        tgt = self.rhs @ self.lhs
        return ChainMap(tgt, self, items)


class ChainMap(object):
    def __init__(self, tgt, src, lins=None, check=True):
        assert isinstance(tgt, Chain)
        assert isinstance(src, Chain)
        n = len(src)
        assert len(tgt) == n
        if lins is None:
            # zero map
            lins = [Lin(tgt.get(i), src.get(i)) for i in range(n+1)]
        assert len(lins) == n+1
        for i in range(n):
            assert src[i].src == lins[i].src
            assert tgt[i].src == lins[i].tgt
        if n:
            assert src[n-1].tgt == lins[n].src
            assert tgt[n-1].tgt == lins[n].tgt
        self.ring = tgt.ring
        self.tgt = tgt
        self.src = src
        self.hom = (tgt, src) # yes it's backwards, just like shape is.
        self.lins = list(lins)
        if check:
            self._check()

    def _check(self):
        tgt, src = self.hom
        lins = self.lins
        n = len(src)
        for i in range(n):
            assert tgt[i]*lins[i] == lins[i+1]*src[i], "not a chain map"

    @classmethod
    def rand(cls, tgt, src):
        # XXX only works for small Chain's !!
        assert isinstance(tgt, Chain)
        assert isinstance(src, Chain)
        assert len(tgt) == len(src)
        n = len(tgt)
        homs = []
        for f, g in zip(tgt, src):
            assert f.src.grade == g.src.grade
            assert f.tgt.grade == g.tgt.grade
            homs.append((f.src, g.src))
        homs.append((f.tgt, g.tgt))
        while 1:
            lins = [Lin.rand(*hom) for hom in homs]
            for i in range(n):
                if tgt[i]*lins[i] != lins[i+1]*src[i]:
                    break
            else:
                break
        return ChainMap(tgt, src, lins)

    def __str__(self):
        return "ChainMap(%s<---%s)"%(self.tgt, self.src)
    __repr__ = __str__

    def __len__(self):
        return len(self.lins)

    def __getitem__(self, idx):
        return self.lins[idx]

    def __eq__(self, other):
        assert self.hom == other.hom
        n = len(self)
        for i in range(n):
            if self[i] != other[i]:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        assert self.hom == other.hom
        n = len(self)
        lins = [self[i]+other[i] for i in range(n)]
        chain = ChainMap(self.tgt, self.src, lins)
        return chain

    def __sub__(self, other):
        assert self.hom == other.hom
        n = len(self)
        lins = [self[i]-other[i] for i in range(n)]
        chain = ChainMap(self.tgt, self.src, lins)
        return chain

    def __rmul__(self, r):
        r = self.ring.promote(r)
        lins = [r*lin for lin in self.lins]
        return ChainMap(self.tgt, self.src, lins)

    def __mul__(lhs, rhs):
        assert lhs.src == rhs.tgt
        lins = [f*g for (f,g) in zip(lhs, rhs)]
        return ChainMap(lhs.tgt, rhs.src, lins) 

    def __matmul__(lhs, rhs):
        assert isinstance(rhs, ChainMap)
        tgt = lhs.tgt @ rhs.tgt
        src = lhs.src @ rhs.src
        #print("__matmul__")
        #print(tgt)
        #print(src)
        idxs = list(src.grades.keys())
        idxs.sort(reverse=True)
        lookup = dict((idx, []) for idx in idxs)
        for r in rhs:
          for l in lhs:
            lr = l@r
            idx = lr.src.grade
            assert idx == lr.tgt.grade
            lookup[idx].append(lr)
        #for k,v in lookup.items():
        #    print(k, len(v))
        lins = []
        for idx in idxs:
            lin = reduce(Lin.direct_sum, lookup[idx])
            lins.append(lin)

        #srcs = src.get_grades()
        #tgts = tgt.get_grades()
        #for (s,t) in zip(srcs, tgts):
        #    print("%s <--- %s" % (t,s))

        return ChainMap(tgt, src, lins)


    def coequalizer(self, other):
        assert self.hom == other.hom
        n = len(self)
        tgt = [] # construct new chain
        lins = [] # construct chain map
        idx = n-1 # working backwards down to idx=0
        while idx >= 0:
            if not lins:
                lin = self[idx].coequalizer(other[idx])
                lins.insert(0, lin)
            else:
                ref = lins[0] * self.tgt[idx]
                lin, f = self[idx].coequalizer(other[idx], ref=ref)
                lins.insert(0, lin)
                tgt.insert(0, f)
            idx -= 1
        tgt = Chain(tgt)
        cmap = ChainMap(tgt, self.tgt, lins)
        return cmap

    def cokernel(self):
        zero = ChainMap(self.tgt, self.src)
        return self.coequalizer(zero)



def test_gf():
    from bruhat.frobenius import GF
    ring = GF(4)
    x = ring.x
    space = Space(ring)

    m = argv.get("m", 3)
    n = argv.get("n", 4)

    one = ring.one
    A = elim.zeros(ring, m, n)
    for i in range(m):
      for j in range(n):
        A[i, j] = choice(ring.elements)

    H = Lin(Space(ring, m, 0), Space(ring, n, 1), A)
    #print(H)

    c = Chain([H])
    #cc = c@c # XXX fix fix fix
    #for f in cc.lins:
    #    print(f)


def randchain(ring, n, m, stem="C"):
    V = Space(ring, n, 2, "%s_2"%stem)
    U = Space(ring, m, 1, "%s_1"%stem)
    B = Lin.rand(U, V)

    A = B.kernel()
    A = Lin(A.tgt, A.src.asgrade(3, "%s_3"%stem), A.A) # yikes
    assert (B*A).is_zero()

    C = B.cokernel()
    C = Lin(C.tgt.asgrade(0, "%s_0"%stem), C.src, C.A)
    assert (C*B).is_zero()

    c = Chain([A, B, C])
    return c


def test_chain():
    ring = element.Q

    m = argv.get("m", 3)
    n = argv.get("n", 3)

    c = randchain(ring, n, m)
    d = randchain(ring, n, m)

    print("c =")
    print(c.longstr())
    print("d =")
    print(d.longstr())

    cd = c+d
    print("c+d =")
    print(cd)
    print(cd.longstr())


def test_chain_tensor():
    #p = argv.get("p", 3)
    #ring = element.FiniteField(p)
    ring = element.Q

    space = Space(ring)

    m = argv.get("m", 3)
    n = argv.get("n", 4)

    one = ring.one
    if argv.toric:
        V = Space(ring, n, 1, "V")
        U = Space(ring, m, 0, "U")

        A = Lin(U, V)
        for i in range(m):
            A[i, i] = one
            A[i, (i+1)%m] = -one
    elif argv.surface:
        V = Space(ring, m, 1, "V")
        U = Space(ring, m-1, 0, "U")

        A = Lin(U, V)
        for i in range(m-1):
            A[i, i] = one
            A[i, (i+1)%m] = -one
    else:
        V = Space(ring, n, 1, "V")
        U = Space(ring, m, 0, "U")
        A = Lin.rand(U, V)

    c = Chain([A])
    #print(c)

    #c.dump()

    #cc = c @ c
    #print(cc)
    #ccc = c@cc
    #print(ccc)


def test_chainmap():

    p = argv.get("p", 2)
    ring = element.FiniteField(p)
    #ring = element.Q

    one = ring.one

    space = Space(ring)

    def mk_ising(m, n):
        U = Space(ring, m, 0, "U")
        V = Space(ring, n, 1, "V")
    
        A = Lin(U, V) # U <--- V
        for i in range(m):
            A[i, i] = one
            if i+1<n:
                A[i, i+1] = -one
        return A

    A = mk_ising(3, 4)
    c = Chain([A])

    f = c.identity()
    zero = ChainMap(c, c)
    assert f != zero
    assert f == f
    assert f != 2*f
    assert f+f == 2*f
    assert f-f == zero

    B = mk_ising(2, 2)
    d = Chain([B])

    fn = Lin(A.src, B.src)
    for i in range(len(B.src)):
        fn[i, i] = one

    fm = Lin(A.tgt, B.tgt)
    for i in range(len(B.tgt)):
        fm[i, i] = one

    f = ChainMap(c, d, [fn, fm])

    # -----------------------------------

    m, n = 8, 10
    V1 = Space(ring, n, 1, "V1")
    V0 = Space(ring, m, 0, "V0")
    A = Lin.rand(V0, V1) # V0 <--- V1

    c = Chain([A])

    U0 = Space(ring, m, 0, "U0")
    f0 = Lin.rand(V0, U0)

    f1, B = A.pullback(f0)
    d = Chain([B])
    f = ChainMap(c, d, [f1, f0])
    
    # -----------------------------------
    # construct a chain map (and a chain) from a 
    # pullback of a grade zero map.

    m, n, p = 5, 6, 1
    V1 = Space(ring, n, 1, "V1")
    V0 = Space(ring, m, 0, "V0")
    #A = Lin.rand(V0, V1) # V0 <--- V1
    A = Lin(V0, V1)
    for i in range(m):
        A[i, i] = one
        A[i, (i+1)] = one

    a = Chain([A])
    #print("A:")
    #print(A)

    U0 = Space(ring, p, 0, "U0")
    f0 = Lin(V0, U0)
    for i in range(p):
        f0[i,i] = one

    f1, B = A.pullback(f0)
    b = Chain([B])
    #print("B:")
    #print(B)
    f = ChainMap(a, b, [f1, f0])
    #print(f0)
    #print(f1)

    g = ChainMap(a, b)
    h = f.coequalizer(g)
    #print(h)
    #print(h[0])
    #print(h[1])
    c = h.tgt
    C = c[0]
    #print(C.shape, C.rank())
    #print(C)
    
    # -----------------------------------
    # construct a 'puncture' of a repitition code 
    # as a cokernel 

    m, n, p = 5, 6, 1
    V1 = Space(ring, n, 1, "V1")
    V0 = Space(ring, m, 0, "V0")
    #A = Lin.rand(V0, V1) # V0 <--- V1
    A = Lin(V0, V1)
    for i in range(m):
        A[i, i] = one
        A[i, (i+1)] = one

    a = Chain([A])

    U1 = Space(ring, 1, 1, "U1")
    U0 = Space(ring, 2, 0, "U0")
    B = Lin(U0, U1)
    B[0, 0] = one
    B[1, 0] = one
    b = Chain([B])

    offset = 2 # where to puncture

    f0 = Lin(V0, U0)
    f0[0+offset,0] = one
    f0[1+offset,1] = one

    f1 = Lin(V1, U1)
    f1[offset+1,0] = one

    f = ChainMap(a, b, [f1, f0])
    h = f.cokernel()
    c = h.tgt

    # -----------------------------------
    # Here we puncture a bit from a parity check matrix
    # by using a cokernel. This deletes that column,
    # as well as all the rows with support therein.

    m, n, p = 5, 6, 1
    V1 = Space(ring, n, 1, "V1")
    V0 = Space(ring, m, 0, "V0")
    A = Lin.rand(V0, V1) # V0 <--- V1

    a = Chain([A])
    #print(A)
    col = 0
    rows = []
    for i in range(V0.n):
        if A[i, col]:
            rows.append(i)

    U1 = Space(ring, 1, 1, "U1") 
    U0 = Space(ring, len(rows), 0, "U0")
    B = Lin(U0, U1)
    for i in range(U0.n):
        B[i, 0] = one
    b = Chain([B])
    #print(B)

    f0 = Lin(V0, U0)
    for i, row in enumerate(rows):
        f0[row, i] = one
    #print(f0)

    f1 = Lin(V1, U1)
    f1[col,0] = one

    f = ChainMap(a, b, [f1, f0])
    h = f.cokernel()
    c = h.tgt

    #print(c[0])


#def reassociate(ring, _C, _D, _E):
#    "C@(D@E) <---- (C@D)@E"

def chain_tensor(_C, _D):
    ring = _C.ring

    rig = Rig(ring)

    bdyC = _C.encode(rig)
    XC, C = bdyC.hom
    n = C.src
    assert C.tgt == n
    assert XC.hom == (n, n)

    N, K = rig.zero, rig.one
    X = [[(K if i-1==j else N) for i in n] for j in n]
    X = Cell1(n, n, X)
    assert XC == (X<<C).normalized
    assert XC == (C<<X).normalized # commutes !

    bdyD = _D.encode(rig)
    XD, D = bdyD.hom
    assert D.hom == (n, n)
    assert XD.hom == (n, n)

    CD = C<<D
    CD = CD.normalized

    front = bdyC << Cell2.identity(D)
    front = front.normalized
    assert front.src == CD

    back = Cell2.identity(C) << bdyD
    back = back.normalized
    assert back.src == CD

    XCD = (XC << D).normalized
    assert front.tgt == XCD
    assert back.tgt == XCD

    mid = front + back
    left = Cell1.codiagonal(n)
    right = Cell1.diagonal(n)
    mid = left << mid << right
    mid = mid.normalized

    bot = Cell2.diagonal(CD)

    cell = mid * bot

    top = Cell2.codiagonal(XCD)
    assert top.tgt == XCD
    
    cell = top * cell
    return cell


def test_tensor():

    p = 2
    ring = element.FiniteField(p)

    U = Space(ring, 4, 1, "U")
    V = Space(ring, 4, 0, "V")

    A = elim.parse(ring, "11.. .11. ..11 1..1")
    f = Lin(V, U, A)
    c = Chain([f])
    #print(c.longstr())

    U0 = Space(ring, 2, 1, "U_0")
    U1 = Space(ring, 2, 1, "U_1")
    U2 = Space(ring, 2, 1, "U_2")

    V0 = Space(ring, 2, 0, "V_0")
    V1 = Space(ring, 2, 0, "V_1")
    V2 = Space(ring, 2, 0, "V_2")

    W0 = Space(ring, 2, 0, "W_0")
    W1 = Space(ring, 2, 0, "W_1")
    W2 = Space(ring, 2, 0, "W_2")

    C = randchain(ring, 2, 2, "C")
    D = randchain(ring, 2, 2, "D")
    E = randchain(ring, 2, 2, "E")
    F = randchain(ring, 2, 2, "F")

    assert C.identity() * C.identity() == C.identity()

    f = ChainMap.rand(E, C)
    assert f * C.identity() == f
    assert E.identity() * f == f

    CD = C@D
    DC = D@C

    s = CD.get_swap()
    assert s.tgt == DC
    assert s.src == CD
    si = DC.get_swap()
    assert CD != DC
    assert s * si == DC.identity()
    assert si * s == CD.identity()

    g = ChainMap.rand(F, D)
    fg = f@g
    EF = E@F
    assert fg.tgt == EF

    lhs = EF.get_swap() * fg * DC.get_swap()
    rhs = g@f
    assert lhs == rhs

    cell = chain_tensor(C, D)
    #print(cell.homstr())
    #CC = Chain([cell[0,2], cell[0,1], cell[0,0]])
    _, n = cell.shape
    CC = Chain([cell[0,n-i-1] for i in range(n)])

    C1 = Space(ring, 3, 1, "C_1")
    C0 = Space(ring, 3, 0, "C_0")
    C = Chain([Lin(C0, C1, [[1,1,0],[0,1,1],[1,0,1]])])

    cell = chain_tensor(C, C)

    CC = Chain([cell[0,2], cell[0,1], cell[0,0]])

    return

    lhs = (c0@c1)@c2
    rhs = c0@(c1@c2)

    #print(lhs)
    #print(rhs)

    f = lhs[0]
    #print(f.homstr())
    tgt = f.tgt
    u = tgt.get_normal()
    #print(u.homstr())


def test_all():
    test_gf()
    #test_chain()
    test_chain_tensor()
    test_chainmap()
    test_tensor()



if __name__ == "__main__":


    fn = argv.next() or "test_all"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("OK\n")


