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


class Matrix(object):
    def __init__(self, ring, rows):
        M = all_cmdline.Matrix(ring, rows)
        M.set_immutable()
        self.M = M
        self.ring = ring

    def __eq__(self, other):
        assert isinstance(other, Matrix)
        assert self.ring == other.ring
        return self.M == other.M

    def __hash__(self):
        return hash(self.M)

    def __str__(self):
        return str(self.M)

    def __mul__(self, other):
        assert isinstance(other, Matrix)
        assert self.ring == other.ring
        M = self.M * other.M
        return Matrix(self.ring, M)

    def __rmul__(self, r):
        M = r*self.M
        return Matrix(self.ring, M)

    def __matmul__(self, other):
        assert isinstance(other, Matrix)
        assert self.ring == other.ring
        M = self.M.tensor_product(other.M)
        return Matrix(self.ring, M)
    tensor_product = __matmul__

    def inverse(self):
        M = self.M.inverse()
        return Matrix(self.ring, M)

    def transpose(self):
        M = self.M.transpose()
        return Matrix(self.ring, M)

    def is_diagonal(self):
        M = self.M
        return M.is_diagonal()
        print(' '.join(dir(M)))
        n = len(M) # XX
        for i in range(n):
          for j in range(n):
            if i==j:
                continue
            if M[i,j] != 0:
                return False
        return True



def mulclose(gen, verbose=False, maxsize=None):
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
    
    CNOT = Matrix(K,
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
    # XXX this _appears to be the projective Clifford group only :-(

    n = 2
    W = mulclose([HI, IH, SWAP]) # Weyl group
    assert len(W) == 8

    #B = mulclose([XI, IX, CNOT, CZ, SI, IS])
    B = mulclose([XI, IX, CNOT, CZ, SI, IS]+gen, verbose=True)
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


