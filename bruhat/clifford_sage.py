#!/usr/bin/env python

"""
build clifford/pauli groups in sage
"""

from time import time
start_time = time()
from random import choice

from bruhat.argv import argv

from sage.all_cmdline import FiniteField, CyclotomicField
from sage import all_cmdline 

#Matrix = lambda *args : all_cmdline.Matrix(*args, immutable=True)

#def Matrix(*args):
#    m = all_cmdline.Matrix(*args)
#    m.set_immutable()
#    return m


class Matrix(object):
    def __init__(self, ring, rows):
        m = all_cmdline.Matrix(ring, rows)
        m.set_immutable()
        self.m = m
        self.ring = ring

    def __eq__(self, other):
        assert isinstance(other, Matrix)
        assert self.ring == other.ring
        return self.m == other.m

    def __hash__(self):
        return hash(self.m)

    def __mul__(self, other):
        assert isinstance(other, Matrix)
        assert self.ring == other.ring
        m = self.m * other.m
        return Matrix(self.ring, m)

    def __rmul__(self, r):
        m = r*self.m
        return Matrix(self.ring, m)

    def __matmul__(self, other):
        assert isinstance(other, Matrix)
        assert self.ring == other.ring
        m = self.m.tensor_product(other.m)
        return Matrix(self.ring, m)
    tensor_product = __matmul__

    def inverse(self):
        m = self.m.inverse()
        return Matrix(self.ring, m)

    def transpose(self):
        m = self.m.transpose()
        return Matrix(self.ring, m)



def mulclose(gen, verbose=False, maxsize=None):
#    for g in gen:
#        g.set_immutable()
    els = set(gen)
    bdy = list(els)
    changed = True
    while bdy:
        if verbose:
            print(len(els), end=" ", flush=True)
        _bdy = []
        for A in gen:
            for B in bdy:
                C = A*B
#                C.set_immutable()
                if C not in els:
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        if verbose:
                            print()
                        return els
        bdy = _bdy
    if verbose:
        print()
    return els


class Coset(object):
    def __init__(self, group, items):
        self.group = group
        self.items = list(items)
        self.g = items[0] # pick one

    def __mul__(self, other):
        assert isinstance(other, Coset)
        gh = self.g * other.g
#        gh.set_immutable()
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
#                gh.set_immutable()
                remain.remove(gh)
                coset.append(gh)
            coset = Coset(self, coset) # promote
            cosets.append(coset)
            if len(remain)%1000 == 0:
                print(len(remain), end=" ", flush=True)
        print()
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


def main():
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
    
    Cliff2 = mulclose([SI, IS, HI, IH, w2II, CZ], verbose=True)
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
#        h.set_immutable()
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
#        gh.set_immutable()
#        g.set_immutable()
#        h.set_immutable()
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




if __name__ == "__main__":

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()


