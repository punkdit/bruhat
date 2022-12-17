#!/usr/bin/env python


"""


"""


import sys, os
from time import time
start_time = time()
import random
from functools import reduce, lru_cache
from operator import add, mul
from string import ascii_lowercase

import numpy

from bruhat.poset import Poset, PreOrder
from bruhat.action import Perm, Group, Coset, mulclose, mulclose_hom, Action
from bruhat.algebraic import Matrix, Algebraic, Figure
from bruhat.solve import shortstr
from bruhat.poset import Poset
from bruhat.smap import SMap
from bruhat.argv import argv

cache = lru_cache(maxsize=None)





class Relation(object):
    "a matrix relation between left and right sets"
    def __init__(self, A, left, right):
        m = len(left)
        n = len(right)
        assert A.shape == (m, n), "%s != %s"%(A.shape, (m,n))
        assert 0 <= A.min()
        assert A.max() <= 1

        #assert A.dtype == numpy.int64
        #A = A.astype(numpy.int16) # faster but take care with overflow !!!
        self.A = A
        self.left = left
        self.right = right
        self.hom = (left, right)
        self.shape = A.shape

    @cache
    def get_llookup(self):
        llookup = {fig:i for (i, fig) in enumerate(self.left)}
        return llookup

    @cache
    def get_rlookup(self):
        rlookup = {fig:i for (i, fig) in enumerate(self.right)}
        return rlookup

    @classmethod
    def top(cls, left, right):
        m = len(left)
        n = len(right)
        A = numpy.zeros((m, n), dtype=int)
        A[:] = 1
        return cls(A, left, right)

    @classmethod
    def bot(cls, left, right):
        m = len(left)
        n = len(right)
        A = numpy.zeros((m, n), dtype=int)
        return cls(A, left, right)

    @classmethod
    def identity(cls, items):
        m = len(items)
        A = numpy.identity(m, dtype=int)
        return cls(A, items, items)

    def __str__(self):
        return shortstr(self.A)+"\n"

    def __eq__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        return numpy.alltrue(self.A==other.A)

    def __hash__(self):
        key = (
            id(self.left),
            id(self.right),
            hash(self.A.tobytes()))
        return hash(key)

    def __lt__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        return numpy.alltrue(self.A<=other.A) and not numpy.alltrue(self.A==other.A)

    def __le__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        return numpy.alltrue(self.A<=other.A)

    def __add__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        A = self.A + other.A
        A = numpy.clip(A, 0, 1)
        return Relation(A, self.left, self.right)
    __or__ = __add__

    def __sub__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        A = self.A - other.A
        A = numpy.clip(A, 0, 1)
        return Relation(A, self.left, self.right)

    def __mul__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        A = self.A * other.A
        return Relation(A, self.left, self.right)
    __and__ = __mul__

    def __matmul__(self, other):
        if isinstance(other, Relation):
            assert self.right is other.left
            A = numpy.dot(self.A, other.A)
            A = numpy.clip(A, 0, 1)
            return Relation(A, self.left, other.right)
        else:
            return self.get_left(other)

    def __rmatmul__(self, left):
        return self.get_right(left)

    @cache
    def transpose(self):
        A = self.A.transpose()
        A = A.copy() # clean up tobytes
        op = Relation(A, self.right, self.left)
        return op

    @property
    def op(self):
        return self.transpose()

    # because we use object identity for __eq__, __hash__
    the_star = ["*"]

    def get_right(self, left):
        llookup = self.get_llookup()
        lidx = llookup[left]
        row = self.A[lidx, :]
        #print(row)
        ridxs = numpy.where(row)[0]
        #print(ridxs)
        right = self.right
        A = numpy.zeros((1, len(right)), dtype=int)
        A[0, ridxs] = 1
        return Relation(A, Relation.the_star, self.right)
        #figs = [right[ridx] for ridx in ridxs]
        #return figs

    def get_left(self, fig):
        op = self.transpose()
        return op.get_right(fig).transpose()

    def nnz(self):
        return self.A.sum()

    def items(self):
        A = self.A
        for idx in zip(*numpy.where(A)):
            i, j = idx
            yield self.left[i], self.right[j]


def hecke_gen(G, left, right, hom=None, verbose=argv.get("verbose", False)):
    "build all the Hecke operators between the left and right figures"

    print("hecke_gen: %sx%s" % (len(left), len(right)))

    if hom is None:
        hom = {g:g for g in G.gen}

    llookup = dict((i, fig) for (fig, i) in enumerate(left))
    rlookup = dict((i, fig) for (fig, i) in enumerate(right))
    m = len(left)
    n = len(right)

    #H = numpy.zeros((m, n))
    remain = set(numpy.ndindex((m, n)))

    ops = []
    while remain:
        if verbose:
            print("[%s:%s]"%(len(ops),len(remain)), end="", flush=True)
        i, j = iter(remain).__next__()
        #assert H[i, j] == 0

        #fig = left[i] + right[j]
        J = numpy.zeros((m, n), dtype=int)

        bdy = set([(i, j)])
        while bdy:
            _bdy = set()
            #print(bdy)
            for (i, j) in bdy:
              l, r = left[i], right[j]
              for g in G.gen:
                i = llookup[g*l]
                j = rlookup[hom[g]*r]
                if J[i, j]:
                    continue
                _bdy.add((i, j))
                J[i, j] = 1
                remain.remove((i, j))
            bdy = _bdy
        #print(shortstr(J))
        #print()
        op = Relation(J, left, right)
        ops.append(op)

    if verbose:
        print()

    ops.sort(key = lambda op : op.A.sum())
    return ops


class Quantale(object):
    def __init__(self, basis):
        assert basis
        a = basis[0]
        assert a.left is a.right
        figs = a.left
        T = Relation.top(figs, figs)
        B = Relation.bot(figs, figs)
        I = Relation.identity(figs)
        for a in basis:
         for b in basis:
            assert a is b or a*b==B
        assert reduce(add, basis) == T

        self.figs = figs
        self.basis = basis
        self.n = len(basis)
        self.T = T
        self.B = B
        self.I = I

    def ldiv(self, q, s):
        basis = self.basis
        assert q in basis
        assert s in basis
        op = self.B
        for r in basis:
            if q@r <= s:
                op = op + r
        return op

    def rdiv(self, q, s):
        basis = self.basis
        assert q in basis
        assert s in basis
        op = self.B
        for r in basis:
            if r@q <= s:
                op = op + r
        return op

    def trinary(self, r, q, s):
        " r@q <= s ? "
        idxs = zip(*numpy.where(s.A==0))
        #result = (r@q <= s) # __matmul__ gets expensive.
        for (i, j) in idxs:
            assert s.A[i,j] == 0
            rq = numpy.dot(r.A[i, :], q.A[:, j])
            if rq>0:
                #assert not result
                return False
        #assert result
        return True

    def compare(self, q, s):
        basis = self.basis
        assert q in basis
        assert s in basis
        #return self.ldiv(q).nnz() or self.rdiv(q).nnz()
        trinary = self.trinary
        for r in basis:
            if trinary(r, q, s):
                return True
            if trinary(q, r, s):
                return True
            #if r@q <= s:
            #    return True
            #if q@r <= s:
            #    return True
        return False


def generate_quantale(gen, verbose=False, maxsize=None):
    if verbose:
        print("generate_quantale:", end=" ", flush=True)
    els = set(gen)
    bdy = list(els)
    changed = True
    while bdy:
        if verbose:
            print(len(els), end=" ", flush=True)
        _bdy = []
        for A in gen:
            for B in bdy:
                for C in [A@B, A+B]:
                  if C not in els:
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        return els
        bdy = _bdy
    if verbose:
        print()
    return els


def check_quantale(basis, check=True):
    # endo Hecke operators form a quantale

    assert basis
    assert len(basis) <= 14, "too big...?"
    op = basis[0]
    assert op.left is op.right
    flag = op.left
    T = Relation.top(flag, flag)
    B = Relation.bot(flag, flag)
    I = Relation.identity(flag)
    for op in basis:
        assert B < op < T
        assert B@op == B
        assert T@op == T
        assert I@op == op # I is monoidal unit 
        assert T & op == op

    assert reduce(add, basis) == T

    Q = generate_quantale(basis, verbose=True)
    Q.add(B)
    assert len(Q) == 2**len(basis)
    assert T in Q

    # A quantale is a partially ordered set:
    #   a<=b is the order
    #   a+b is join, or sup, or colimit
    #   a*b is meet, or inf, or limit
    #   a@b is the monoidal product, or relational composition

    # See: 
    # ON SOME BASIC CONSTRUCTIONS IN CATEGORIES OF QUANTALE-VALUED SUP-LATTICES
    # RADEK SLESINGER
    for a in Q:
     for b in Q:
      #assert a@b == b@a # nope !
      for c in Q:
        assert a @ (b + c) == a@b + a@c # preserves colimit

    quantale = Quantale(basis)
    for q in basis:
     for s in basis:
        lhs = quantale.compare(q, s)
        rhs = quantale.ldiv(q, s).nnz() or quantale.rdiv(q, s).nnz()
        assert lhs == bool(rhs)

    # internal hom, ldiv[q,_] is right adjoint to q@_
    ldiv = {}
    for q in Q:
     for s in Q:
        op = B
        for r in basis:
            if q@r <= s:
                op = op+r
        ldiv[q,s] = op

    # check hom-tensor adjunction
    for r in Q:
     for q in Q:
      for s in Q:
        lhs = q@r <= s
        rhs = r <= ldiv[q,s]
        assert lhs == rhs

    for r in Q:
     for s in Q:
      for t in Q:
        lhs = ldiv[r+s, t]
        rhs = ldiv[r, t] & ldiv[s, t]
        assert lhs == rhs

    # internal hom, rdiv[q,_] is right adjoint to _@q
    rdiv = {}
    for q in Q:
     for s in Q:
        op = B
        for r in basis:
            if r@q <= s:
                op = op+r
        rdiv[q,s] = op

    # check hom-tensor adjunction
    for r in Q:
     for q in Q:
      for s in Q:
        lhs = r@q <= s
        rhs = r <= rdiv[q,s]
        assert lhs == rhs


def test_hecke_GL3():

    print("test_hecke_GL3")

    n = 3
    G = Algebraic.GL(n)

    tps = [
        (G.all_figs([1])),
        (G.all_figs([2])),
        (G.all_figs([2,1])),
    ]
    points, lines, flags = tps
    assert len(points) == 7
    assert len(lines) == 7
    assert len(flags) == 21

    FF = hecke_gen(G, flags, flags)
    FL = hecke_gen(G, flags, lines)
    LL = hecke_gen(G, lines, lines)

    assert len(FL) == 3
    L, P, R = FL

    assert len(LL) == 2
    assert len(FF) == 6
    #I, L, P, LP, PL, LPL = FF

    flag = flags[0]
    print(flag)

    for op in FL:
        rel = flag @ op
        print(rel)

    bruhats = set()
    for op in FF:
        figs = op @ flag
        #print(figs)
        figs = figs.op @ L
        bruhats.add(figs)
    assert len(bruhats) == 3

    I, L, P, LP, PL, LPL = FF
    names = "I L P LP PL LPL".split()
    lookup = dict(kv for kv in zip(FF, names))

    LL = L@L
    assert LL == I+L
    assert LL >= I
    assert LL > I
    assert LL >= L
    assert LL > L
    assert not LL <= I
    assert not LL < I
    assert not LL <= L
    assert not LL < L

    LPLP = LP@LP
    assert LPLP >= PL
    assert LP @ LP == PL + LPL

    # much faster
    Q = Quantale(FF)
    pairs = []
    for a in FF:
     for b in FF:
        if Q.ldiv(a,b).nnz() or Q.rdiv(a,b).nnz():
            pairs.append((a, b))

    pairs = set(pairs)
    print("pairs:", len(pairs))

    for a, b in pairs:
        print(lookup[a], ">=", lookup[b])
    assert (L,I) not in pairs
    assert (I,L) in pairs
    assert (LP, PL) not in pairs
    assert (L, LP) in pairs
    assert (L, PL) in pairs

    order = PreOrder.generate(pairs, FF, check=True)
    #order.show() # _looks good!

    #s = order.get_dot(labels=False)
    #print(s)


def test_hecke_GL4():

    print("test_hecke_GL4")

    n = 4
    G = Algebraic.GL(n)

    tps = [
        (G.all_figs([1])),
        (G.all_figs([2])),
        (G.all_figs([3])),
        (G.all_figs([3,2,1])),
    ]
    points, lines, planes, flags = tps

    print(len(points), len(lines), len(planes), ":", len(flags))

    FF = hecke_gen(G, flags, flags)
    print("FF:", len(FF))

    I, P, L, A = FF[:4]

    FL = hecke_gen(G, flags, lines)
    print("FL:", len(FL))

    flag = flags[0]
    print(flag)

    Q = Quantale(FF)

    bruhats = set()
    lookup = {}
    for a in FF:
        figs = a @ flag
        figs = figs.op @ FL[0]
        bruhats.add(figs)
        lookup[a] = figs
    print( len(bruhats) )
    for figs in bruhats:
        print("bruhat")
        for _,fig in figs.items():
            print('\t', fig)

    # reverse
    lookup = dict((k,v) for (v,k) in lookup.items())

    pairs = []
    for a in bruhats:
     for b in bruhats:
        left = lookup[a]
        right = lookup[b]
        if Q.compare(left, right):
            pairs.append((a, b))
            print("*", end="", flush=True)
    print()

    order = PreOrder.generate(pairs, check=True)
    order.show(labels=False) # _looks good!
    return

    pairs = []
    for a in FF:
     for b in FF:
        if Q.ldiv(a,b).nnz() or Q.rdiv(a,b).nnz():
            pairs.append((a, b))
     print(len(pairs))

    pairs = set(pairs)
    print("pairs:", len(pairs))

    order = PreOrder.generate(pairs, FF, check=True)
    order.show(labels=False) # _looks good!

    s = order.get_dot(labels=False)
    print(s)


def bruhat_str(figs):
    items = []
    for fig in figs:
        assert len(fig.items)==1, "TODO"
        items.append(Matrix(fig.items[0]))
    G = Algebraic(items, G=items) # bit of a hack ?!?
    return G.show()

def fig_str(figs):
    items = []
    for _,fig in figs.items():
        items.append(fig)
    return bruhat_str(items)

def fig_key(a, b):
    return fig_str(a), fig_str(b)


def test_hecke_Sp2():

    print("test_hecke_Sp2")

    n = 2
    nn = 2*n
    G = Algebraic.Sp(nn)

    tps = [
        (G.all_figs([1])),
        (G.all_figs([2])),
        (G.all_figs([2,1])),
    ]
    points, lines, flags = tps

    print(len(points), len(lines), ":", len(flags))

    FF = hecke_gen(G, flags, flags, verbose=True)
    print("FF:", len(FF))

    for op in FF:
        print(op.nnz(), end=' ')
    print()

    I, P, L = FF[:3]

    FL = hecke_gen(G, flags, lines, verbose=True)
    print("FL:", len(FL))

    flag = flags[0]
    print(flag)

    Q = Quantale(FF)

    bruhats = set()
    lookup = {}
    for a in FF:
        figs = a @ flag
        figs = figs.op @ FL[0]
        bruhats.add(figs)
        lookup[a] = figs
    print( len(bruhats) )

    for figs in bruhats:
        print("bruhat")
        print(fig_str(figs))

    # reverse
    lookup = dict((k,v) for (v,k) in lookup.items())

    pairs = [fig_key(a,a) for a in bruhats]
    count = 0
    for a in bruhats:
     for b in bruhats:
        left = lookup[a]
        right = lookup[b]
        order = PreOrder.generate(pairs, check=True)
        assert fig_key(a, a) in order.pairs
        key = fig_key(a, b)
        if key in order.pairs:
            continue
        count += 1
        if Q.compare(left, right):
            pairs.append(key)
            print("*", end="", flush=True)
    print()
    print("count:", count)

    order = PreOrder.generate(pairs, check=True)
    order.show(labels=True) # _looks good!


def test_hecke_Sp3():

    print("test_hecke_Sp3")

    n = 3
    nn = 2*n
    G = Algebraic.Sp(nn)

    tps = [
        (G.all_figs([1])),
        (G.all_figs([2])),
        (G.all_figs([3])),
        (G.all_figs([3,2,1])),
    ]
    points, lines, planes, flags = tps

    print(len(points), len(lines), len(planes), ":", len(flags))

    FF = hecke_gen(G, flags, flags, verbose=True)
    print("FF:", len(FF))

    for op in FF:
        print(op.nnz(), end=' ')
    print()

    I, P, L, A = FF[:4]

    FL = hecke_gen(G, flags, lines, verbose=True)
    print("FL:", len(FL))

    flag = flags[0]
    print(flag)

    Q = Quantale(FF)

    bruhats = set()
    lookup = {}
    for a in FF:
        figs = a @ flag
        figs = figs.op @ FL[0]
        bruhats.add(figs)
        lookup[a] = figs
    print( len(bruhats) )
    for figs in bruhats:
        print("bruhat")
        print(fig_str(figs))
    print()

    # reverse
    lookup = dict((k,v) for (v,k) in lookup.items())

    pairs = [fig_key(a,a) for a in bruhats]
    for a in bruhats:
     for b in bruhats:
        left = lookup[a]
        right = lookup[b]
        order = PreOrder.generate(pairs, check=True)
        assert fig_key(a, a) in order.pairs
        key = fig_key(a, b)
        if key in order.pairs:
            pairs.append(key)
            print("/\\", end="", flush=True)
            continue
        print("_", end="", flush=True)
        if Q.compare(left, right):
            pairs.append(key)
            print("*", end="", flush=True)
        else:
            print(".", end="", flush=True)
     print()
    print()

    order = PreOrder.generate(pairs, check=True)
    order.show(labels=True)



def test_hecke_monoid():

    print("test_hecke_monoid")

    n = 3
    G = Algebraic.GL(n)

    flag = G.all_figs([2,1])
    assert len(flag) == 21

    gen = hecke_gen(G, flag, flag)
    assert len(gen) == 6
    check_quantale(gen)

    T = Relation.top(flag, flag)
    B = Relation.bot(flag, flag)

    I, L, P, LP, PL, LPL = gen

    assert I@L == L
    assert L@L == I+L
    assert P@P == I+P
    assert P < P@P
    assert LP == L@P
    assert LP@LP == PL + LPL
    assert LP@LP > LPL
    assert PL == P@L
    assert LPL == L@P@L
    assert LPL == P@L@P

    assert LPL @ LPL == T



def test_hecke_induce():

    n = 4
    G0 = Algebraic.GL(n-1)
    hom = {}
    gen = []
    for g0 in G0.gen:
        A = numpy.identity(n, dtype=scalar)
        A[:n-1, :n-1] = g0.A
        g = Matrix(A)
        gen.append(g)
        hom[g0] = g
    G01 = Algebraic(gen)
    G1 = Algebraic.GL(n)

    print(len(G01), ">-->", len(G1))
    for g in G01.gen:
        assert g in G1

    print("(3 1):")
    for fig in G0.all_figs([1]):
        print('\t', fig)
    print("(3 2):")
    for fig in G0.all_figs([2]):
        print('\t', fig)

    print(len((G1.all_figs([1]))))
    print(len((G1.all_figs([2]))))
    print(len((G1.all_figs([3]))))

    left = (G0.all_figs([2,1]))
    right = (G1.all_figs([3,2,1]))
    ops = hecke_gen(G0, left, right, hom)
    print("ops:", len(ops))
    for op in ops[:7]:
        A = op.A
        #print(shortstr(A))
        print(A.sum(0))
        print(A.sum(1))

    #B = reduce(add, [op.A for op in ops[:7]])
    #print(shortstr(B))


def test_symplectic():

    n = argv.get("n", 2)
    nn = n*2
    m = argv.get("m", 1)
    assert m<=n

    p = argv.get("p", 2)

    G = Algebraic.Sp(nn, p)

    print("|G| =", len(G))

    if n==2:
        left = (G.all_figs([2, 1]))
        right = (G.all_figs([m]))
    elif n==3:
        left = (G.all_figs([3, 2, 1]))
        #right = (G.all_figs([3, 2, 1]))
        right = (G.all_figs([m]))

    print("left:", len(left))
    print("right:", len(right))

    W = G.get_weyl()
    B = G.get_borel()
    fig = right[0]
    P = set()
    for g in G:
        if g*fig == fig:
            P.add(g)
    print(len(G)//len(P))

    PW = [w for w in W if w in P]
    #gen = PW + list(B)
    #P1 = Algebraic(gen)
    #P = Algebraic(list(P))
    #assert P1==P
    #for g in P1:
    #    assert g in P # yes

    ops = hecke_gen(G, left, right)
    print("left*right:", len(ops))

    for op in ops:
        print("bruhat class")
        for flag in left:
            print(flag.shortstr(), "->", end=" ")
            for fig in op.get_right(flag):
                print(fig.shortstr(), end=" ")
            print()
            break


def test_hecke_Sp():

    print("test_hecke_Sp")

    n = 3
    nn = 2*n
    G = Algebraic.Sp(nn)

    tps = [
        G.all_figs([1]),
        G.all_figs([2]),
        G.all_figs([3]),
        G.all_figs([3,2,1]),
    ]
    point, line, plane, flag = tps
    desc = "point line plane flag".split()

    print(len(line), len(flag))

    ops = hecke_gen(G, flag, line)
    print(len(ops))
    fig = flag[0]
    for op in ops:
        #print(shortstr(op.A))
        #print(op.A.sum())
        for other in op.get_right(fig):
            print(other)
        print()


def test_hecke_eigs():

    n = argv.get("n", 3)
    G = Algebraic.SL(n)
    print("|G| =", len(G))

    left = argv.get("left", [n,1]) 
    right = argv.get("right", left)

    left = list(Figure.qchoose(left))
    right = list(Figure.qchoose(right))

    ops = hecke_gen(G, left, right)
    print("Hecke operators:", len(ops))

    if argv.eigvals:
      for op in ops:

        J = op.A
        print(J.shape, int(round(J.sum())))
        vals = numpy.linalg.eigvals(J)
        #print(vals)
        ss = []
        for x in vals:
            if abs(x.imag)>EPSILON and abs(x.real)>EPSILON:
                ss.append(str(x))
            elif abs(x.imag)>EPSILON:
                ss.append("%.4f"%(x.imag)+"j")
            elif abs(x.real - int(round(x.real))) < EPSILON:
                ss.append("%.0f"%x.real)
            else:
                ss.append("%.4f"%x.real)
        ss = set(ss)
        print("eigs:", ' '.join(ss), '\n')



if __name__ == "__main__":
    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))

    

    




