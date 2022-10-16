#!/usr/bin/env python

"""
rewrite of chain.py, reversing order of Chain lins

"""


from time import time
start_time = time()

from random import choice
from functools import reduce
from operator import matmul

import numpy

from bruhat import elim
from bruhat.lin import Lin, Space, AddSpace, MulSpace, element
from bruhat.smap import SMap
from bruhat.argv import argv
from bruhat.util import cross, allperms



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
            #if maxsize and len(ops) >= maxsize:
            #    break
        bdy = _bdy
        if verbose:
            print("mulclose:", len(ops))
        if maxsize and len(ops) >= maxsize:
            break
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

    def __mul__(self, other):
        assert isinstance(other, ChainMap)
        assert other.tgt == self.src
        lins = [l*r for (l,r) in zip(self, other)]
        return ChainMap(self.tgt, other.src, lins)

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



def test_kagome():

    p = 2
    ring = element.FiniteField(p)

    C1 = Space(ring, 3, 1, "C_1")
    C0 = Space(ring, 3, 0, "C_0")
    f = Lin(C0, C1, [[1,1,0],[0,1,1],[1,0,1]])
    f = Chain([f])

    send = ChainMap(f, f, [
        Lin(C0, C0, [[0,1,0],[0,0,1],[1,0,0]]),
        Lin(C1, C1, [[0,1,0],[0,0,1],[1,0,0]]),
    ])
    inv = ChainMap(f, f, [lin.transpose() for lin in send])

    f3 = Chain.product(ring, [f]*3)

    i = f.identity()
    assert i*i == i
    assert send != i
    assert send*send != i
    assert send*send*send == i

    assert send != inv
    assert inv != i
    assert inv*inv != i
    assert inv*inv*inv == i

    assert send*inv == i
    assert inv*send == i

    a = ChainMap.product(ring, [send, inv, i])
    b = ChainMap.product(ring, [send, i, inv])
    c = ChainMap.product(ring, [i, send, inv])

    G = mulclose([a, b, c], verbose=True)
    print("|G| =", len(G))

    for g in G:
      for h in G:
        assert g*h == h*g



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



