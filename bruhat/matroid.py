#!/usr/bin/env python

"""
Use z3 sat solver to find all matroids (and symplectic matroids)
on an n-element set.

"""

import os
from random import shuffle
from math import sin, cos, pi
from functools import cache, reduce
from operator import add

import numpy

from z3 import Bool, And, Or, Xor, Not, Implies, Sum, If, Solver

from bruhat.action import get_orbits
from bruhat.gset import Group, Perm
from bruhat.algebraic import qchoose, Matrix, row_reduce_p, Algebraic
from bruhat.util import all_subsets
from bruhat.argv import argv



def le(n, a, b):
    for i in range(n):
        if a[i]>b[i]:
            return False
    return True


#@cache
#def all_sp_masks(nn): # uturn order
#    n = nn//2
#    masks = []
#    for mask in numpy.ndindex((2,)*nn):
#        for i in range(n):
#            if mask[i] and mask[nn-i-1]:
#                #print("\tskip", mask)
#                break
#        else:
#            masks.append(mask)
#    return masks


@cache
def all_sp_masks(nn): # zip order
    n = nn//2
    masks = []
    for mask in numpy.ndindex((2,)*nn):
        for i in range(n):
            if mask[2*i] and mask[2*i+1]:
                #print("\tskip", mask)
                break
        else:
            masks.append(mask)
    return masks



class Matroid:
    def __init__(self, n, masks):
        "n:number of elements, masks: the independent sets as bit-vectors"
        self.n = n
        masks = list(masks)
        masks.sort()
        masks = tuple(masks)
        self.masks = masks
        self.items = set(masks)
        self.key = (n, masks)
        self.rank = max(sum(a) for a in masks)

    @classmethod
    def from_basis(cls, n, basis):
        masks = []
        for mask in numpy.ndindex((2,)*n):
            for b in basis:
                if le(n, mask, b):
                    masks.append(mask)
                    break
        M = cls(n, masks)
        return M

    @classmethod
    def from_lin(cls, H, q):
        m, n = H.shape
    
        masks = []
        #for idxs in all_subsets(n):
        for mask in numpy.ndindex((2,)*n):
            idxs = [i for i in range(n) if mask[i]]
            #H1 = H[:, idxs].transpose()
            H1 = H[:, idxs].t
            #H1 = row_reduce_p(H1, q, True)
            H1 = H1.row_reduce()
            if len(H1) == len(idxs):
            #if Matrix(H1,q).t.rank() == len(idxs):
                masks.append(mask)
        M = cls(n, masks)
        #M.check()
        return M
    
    
    #def mul(self, a, b):
    #    return tuple(ai

    def get_basis(self):
        rank = self.rank
        items = []
        for a in self.masks:
            if sum(a)==rank:
                items.append(a)
        return items

    def get_dual(self):
        basis = []
        n = self.n
        for b in self.get_basis():
            basis.append(tuple(1-i for i in b))
        M = self.__class__.from_basis(n, basis)
        return M

    def le(self, mask, nask):
        return le(self.n, mask, nask)

    def restrict(self, mask):
        assert len(mask) == self.n
        masks = {m for m in self.masks if self.le(m, mask)}
        return self.__class__(self.n, masks)

    def all_masks(self):
        return list(numpy.ndindex((2,)*self.n))

    def rankfunc(self):
        #masks = [tuple(mask) for mask in self.all_masks()]
        func = {mask:self.restrict(mask).rank for mask in self.all_masks()}
        return func

    def __str__(self):
        masks = [{i for (i,ii) in enumerate(mask) if ii} for mask in self.get_basis()]
        return "%s(%d, %s)"%(self.__class__.__name__, self.n, masks)
    __repr__ = __str__

    def __rmul__(self, g):
        #print("__rmul__")
        masks = [tuple(m[g[i]] for i in range(self.n)) for m in self.masks]
        return self.__class__(self.n, masks)

    def check(self):
        # check the Matroid axioms
        n = self.n
        items = self.items

        # (1) the empty set is independent
        assert (0,)*n in items

        # (2) independence is down-closed
        for a in self.all_masks():
            for b in items:
                if self.le(a, b):
                    assert a in items

        # (3) independent set exchange 
        for a in items:
          for b in items:
            if sum(a) <= sum(b):
                continue
            # A has more elements than B
            for i in range(n):
                if a[i] and not b[i]:
                    break
            else:
                assert 0


    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return self.key == other.key

    def __le__(self, other):
        assert self.n == other.n
        return self.items.issubset(other.items)


class SpMatroid(Matroid):
    def all_masks(self):
        # all admissable sets as bitvectors on 2*n elements
        nn = self.n
        assert nn%2 == 0
        n = nn//2
        masks = []
        for vec in numpy.ndindex((3,)*n):
            mask = [0]*nn
            for i in range(n):
                mask[2*i:2*i+2] = [[0,0], [0,1], [1,0]][vec[i]]
            mask = tuple(mask)
            masks.append(mask)
        return masks

    def get_dual(self):
        assert 0, "SpMatroid"
    
    def check(self):
        # check the symplectic matroid axioms
        Matroid.check(self)
        nn = self.n
        items = self.items

        masks = set(all_sp_masks(nn))
        for mask in items:
            assert mask in masks, mask

        # TODO: add symplectic exchange condition:
        r"""
        A subset-closed family I of admissible subsets of {0,...,nn-1} is
        the collection of independent sets of a symplectic matroid
        if and only if it has the following property:
            If X and Y are members of I such that |X| < |Y|, then either
            1. there exists y in Y-X such that XU{y} in I, or
            2. X U Y is inadmissible, and there exists x not in XUY such that both
                X U {x} in I and (X-Y*) U {x*} in I.
        """

    @classmethod
    def from_sp(cls, H, q): # zip-order
        m, nn = H.shape
        assert nn%2 == 0
        n = nn//2
    
        masks = []
        for mask in all_sp_masks(nn):
            #print("\t", mask)
            idxs = [i for i in range(nn) if mask[i]]
            #H1 = H[:, idxs].transpose()
            H1 = H[:, idxs].t
            #H1 = row_reduce_p(H1, q, True)
            H1 = H1.row_reduce()
            if len(H1) == len(idxs):
            #if cls(H1,q).t.rank() == len(idxs):
                masks.append(mask)
        #print(masks)
        M = cls(nn, masks)
        #M.check()
        return M




def all_matroids(n):

    bits = list(numpy.ndindex((2,)*n))

    # is this subset an independent set?
    vs = []
    lookup = {}
    for a in bits:
        name = "v"+''.join(str(i) for i in a)
        v = Bool(name)
        vs.append(v)
        lookup[a] = v

    solver = Solver()
    Add = solver.add

    def le(a, b):
        for i in range(n):
            if a[i]>b[i]:
                return False
        return True

    pairs = [(a,b) for a in bits for b in bits]
    for a in bits:
        if sum(a)==0:
            Add(lookup[a])

    for (a,b) in pairs:
        if le(a,b):
            Add( If(lookup[b], lookup[a], True) )

    for (a,b) in pairs:
        if sum(a) <= sum(b):
            continue
        terms = []
        for i in range(n):
            if a[i] and not b[i]:
                c = list(b)
                c[i] = 1
                c = tuple(c)
                terms.append(lookup[c])
        assert terms
        Add( If(And(lookup[a], lookup[b]), Or(*terms), True) )
    
    count = 0
    matroids = set()
    while 1:
        result = solver.check()
        if str(result) != "sat":
            break

        model = solver.model()
    
        sol = {a:model.get_interp(lookup[a]) for a in bits}
        found = []
        for (k,v) in sol.items():
            if v:
                found.append(k)
        M = Matroid(n, found)
        M.check()
        assert M not in matroids
        matroids.add(M)
        yield M

        term = []
        for a in bits:
            if sol[a]:
                term.append(Not(lookup[a]))
            else:
                term.append((lookup[a]))
        term = Or(*term)
        Add(term)
        count += 1
    #print("all_matroids(%d) = %d" % (n, count))


def coxeter_bc(n):
    # We use ziporder: (2*i,2*i+1) are the symplectic pairs
    nn = 2*n
    gen = []
    for i in range(n):
        perm = list(range(nn))
        perm[2*i:2*i+2] = 2*i+1, 2*i # coin flip
        gen.append(Perm(perm))
    if n > 1:
        perm = [(i+2)%nn for i in range(nn)] # rotate
        gen.append(Perm(perm))
    if n > 2:
        perm = [2,3,0,1] + list(range(4, nn))
        gen.append(Perm(perm))

    return gen


def all_sp_matroids(n, rank=None):
    # We use ziporder: (2*i,2*i+1) are the symplectic pairs
    # See: Theorem 3.8.1 in [Borovik,Gelfand,White]

    nn = 2*n

    # all admissable sets as bitvectors on 2*n elements
    bits = []
    for vec in numpy.ndindex((3,)*n):
        mask = [0]*nn
        for i in range(n):
            mask[2*i:2*i+2] = [[0,0], [0,1], [1,0]][vec[i]]
        mask = tuple(mask)
        bits.append(mask)
        
    # is this subset an independent set?
    lookup = {} # variables
    for a in bits:
        name = "v"+''.join(str(i) for i in a)
        v = Bool(name)
        lookup[a] = v


    solver = Solver()
    Add = solver.add
    #def Add(term):
    #    print("Add: %s"%(term))
    #    solver.add(term)

    def le(a, b):
        for i in range(nn):
            if a[i]>b[i]:
                return False
        return True

    pairs = [(a,b) for a in bits for b in bits]
    for a in bits:
        if sum(a)==0:
            Add(lookup[a]) # empty set

    # down-closed
    for (a,b) in pairs:
        if a!=b and le(a,b):
            Add( If(lookup[b], lookup[a], True) )

    if rank is not None:
        #print("\nrank =", rank)
        terms = []
        for a in bits:
            if sum(a) > rank:
                Add(Not(lookup[a]))
                #print(Not(lookup[a]), end=' ')
            if sum(a) == rank:
                terms.append(lookup[a])
        #assert len(terms) == 2**n
        #print()
        #print(terms)
        if len(terms):
            Add(Or(*[terms]))
    
    star = lambda i : [i+1, i-1][i%2]
    for i in range(nn):
        assert 0<=star(i)<nn

    #print("symplectic exchange:")
    for (X,Y) in pairs:
        #print("pair:", X, Y)
        if sum(X) >= sum(Y):
            continue
        # |X| < |Y|
        src = And(lookup[X], lookup[Y]) # "if X and Y are independent sets"
        terms = []
        for i in range(nn):
            if Y[i] and not X[i]:
                Z = list(X)
                Z[i] = 1
                Z = tuple(Z)
                if Z in lookup: # Z admissable
                    terms.append(lookup[Z])
        assert terms
        Z = tuple(max(i,j) for (i,j) in zip(X,Y)) # Z = XuY
        if Z not in lookup:
            #print("XuY inadmissable", X, Y, Z)
            #items = []
            for i in range(nn):
                if Z[i]:
                    continue
                Xi = list(X)
                Xi[i] = 1
                Xi = tuple(Xi)
                if Xi not in lookup:
                    continue
                Yi = [0]*nn
                for j in range(n):
                    if X[j] and not Y[star(j)]:
                        Yi[j] = 1
                assert Yi[i] == 0
                Yi[i] = 1
                Yi = tuple(Yi)
                assert Yi in lookup
                t = And(lookup[Xi], lookup[Yi])
                #print("\t", t)
                terms.append(t)
            #terms.append
        Add( If(src, Or(*terms), True) )

    #return
    
    count = 0
    matroids = set()
    while 1:
        result = solver.check()
        if str(result) != "sat":
            break

        model = solver.model()
    
        sol = {a:model.get_interp(lookup[a]) for a in bits}
        found = []
        for (k,v) in sol.items():
            if v:
                found.append(k)
        M = SpMatroid(nn, found)
        M.check()
        assert M not in matroids
        matroids.add(M)
        yield M

        term = []
        for a in bits:
            if sol[a]:
                term.append(Not(lookup[a]))
            else:
                term.append((lookup[a]))
        term = Or(*term)
        Add(term)
        count += 1
    #print("all_matroids(%d) = %d" % (n, count))


def all_orbits(n):


    G = Group.symmetric(n)

    remain = set(all_matroids(n))

    orbits = []
    while remain:
        m = remain.pop()
        orbit = {g*m for g in G}
        remain.difference_update(orbit)
        orbit = list(orbit)
        orbits.append(orbit)
    return orbits


def show_sp():
    for rank in [2]:
        print("rank:", rank)
        for M in all_sp_matroids(rank, rank):
            if (1,0,1,0) in M.masks and (0,1,1,0) in M.masks and (0,1,0,1) in M.masks:
                print(M)


def test_lag():
    n = 2
    items = set(all_sp_matroids(n, n))
    for M in items:
        N = M.get_dual() 
        assert N in items
        print(M==N, end=' ')
    print()


def test_sp():

    for n in range(4):
      gen = coxeter_bc(n)
      print("n=%d:"%n, end=' ', flush=True)
      for rank in range(0, n+1):
        count = 0
        items = list(all_sp_matroids(n, rank))
        orbits = get_orbits(gen, items)
        print(len(orbits), end=' ', flush=True)
      print()

    
def test_sp_labelled():

    for n in range(4):
      print("n=%d:"%n, end=' ', flush=True)
      for rank in range(0, n+1):
        count = 0
        for M in all_sp_matroids(n, rank):
            M.check()
            assert M.rank == rank, M
            count += 1
        print(count, end=' ', flush=True)
      print()

    
def test_matroids():

    n = 4
    for M in all_matroids(n):
        N = M.from_basis(n, M.get_basis())
        N.check()
        #print(M, M.get_basis())
        #print(N, N.get_basis())
        assert M == N
        K = M.get_dual()
        K.check()
    #return
    
    # https://oeis.org/A058669
    for n in range(6):
        found = {i:[] for i in range(n+1)}
        for M in all_matroids(n):
            #print(M, M.rank)
            found[M.rank].append(M)
        print(n, [len(found[i]) for i in range(n+1)])


    # https://oeis.org/A058673
    assert len(list(all_matroids(0))) == 1
    assert len(list(all_matroids(1))) == 2
    assert len(list(all_matroids(2))) == 5
    assert len(list(all_matroids(3))) == 16
    assert len(list(all_matroids(4))) == 68
    assert len(list(all_matroids(5))) == 406
    #assert len(list(all_matroids(6))) == 3807
    #assert len(list(all_matroids(7))) == 75164

    

def test_repr():
    "representable (linear) matroids"

    q = argv.get("q", 2)

    n = argv.get("n", 3)
    m = argv.get("m", 2)

    #G = Algebraic.GL(n, p=q)

    count = 0
    found = set()
    for H in qchoose(n, m, q):
        H = Matrix(H, q) # <--- XXX Expensive constructor
        M = Matroid.from_lin(H, q)
        #print(H, M)
        found.add(M)
        count += 1
    print(count, len(found))

    
def test_repr_sp(nn=6, m=3, q=3):
    "representable (linear) symplectic matroids"

    q = argv.get("q", q)

    nn = argv.get("nn", nn)
    assert nn%2 == 0
    n = nn//2
    m = argv.get("m", m)

    remain = set( all_sp_matroids(n, m) )
    print("remain:", len(remain))

    G = Algebraic.Sp(nn, p=q)

    perm = []
    for i in range(n):
        perm.append(i)
        perm.append(nn-i-1)
    perm = Perm(perm)
    #print(perm)

    count = 0
    found = set()
    #for H in qchoose(n, m, q):
    for H in G.qchoose(m): # qchoose is bottleneck here...
        A = H.A[:, perm.perm]
        H = Matrix(A, q)
        M = SpMatroid.from_sp(H, q)
        M.check()
        #print(H, M)
        #M = perm*M # convert uturn to zip-order
        found.add(M)
        if M in remain:
            remain.remove(M)
        count += 1
    print(count, len(found))

    print("remain:")
    for M in remain:
        print(M)
    print(len(remain))
    #print("found:")
    #for M in found:
    #    print(M)
    #print(len(found))

    return remain
    
    




if __name__ == "__main__":
    from time import time
    start_time = time()
    fn = argv.next() or "test_matroids"

    print("%s()"%fn)

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))


