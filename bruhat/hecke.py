#!/usr/bin/env python

"""
these are G-sets where G is the clifford group
acting on various codes & tuples of codes, etc.
orbit counting, etc....

"""

from operator import mul, matmul, add
from functools import reduce

from sage.all import ZZ, QQ

import numpy

from bruhat import matrix_sage
from bruhat.argv import argv
from bruhat.gset import mulclose, Perm, Group, GL
from bruhat.action import mulclose_names
from bruhat.solve import shortstr
from bruhat.util import all_subsets
from bruhat.tom import Tom
from bruhat.todd_coxeter import Schreier
from bruhat.repr_sage import GL32

from bruhat.matrix_sage import CyclotomicField, Matrix
from bruhat.clifford_sage import Clifford, K, w4, get_pauli_gens, get_clifford_gens




def test_double_cosets():

    n = argv.get("n", 2)
    k = argv.get("k", 0)

    gen = get_pauli_gens(n)

    Pauli = mulclose(gen)
    Pauli = list(Pauli)
    lookup = dict((g,i) for (i,g) in enumerate(Pauli))
    print("Pauli:", len(Pauli))
    assert len(Pauli) == 4**(n+1)

    gen = get_clifford_gens(n)

    perms = []
    for g in gen:
        ig = ~g
        idxs = [lookup[ig*h*g] for h in Pauli]
        perm = Perm(idxs)
        perms.append(perm)
    gen = perms

    stab = []
    I, X = Clifford(1).I, Clifford(1).X()
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


def test_double_cosets_faster(): 
    n = argv.get("n", 2)
    k = argv.get("k", 0)

    gen = get_pauli_gens(n)
    Pauli = mulclose(gen)
    Pauli = list(Pauli)
    lookup = dict((g,i) for (i,g) in enumerate(Pauli))
    print("Pauli:", len(Pauli))
    assert len(Pauli) == 4**(n+1)

    gen = get_clifford_gens(n)

    perms = []
    for g in gen:
        ig = ~g
        idxs = [lookup[ig*h*g] for h in Pauli]
        perm = Perm(idxs)
        perms.append(perm)
    gen = perms

    I = Clifford(1).I
    X = Clifford(1).X()
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
    print()
    codes = list(codes)
    codes.sort()
    lookup = {c:idx for (idx,c) in enumerate(codes)}

    perms = []
    for g in gen:
        perm = []
        for (i,c) in enumerate(codes):
            j = lookup[mkpoint([g[k] for k in c])]
            perm.append(j)
        perm = Perm(perm)
        perms.append(perm)
    gen = perms

    find_double_cosets(gen)


def find_double_cosets(gen):
    N = gen[0].rank
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
        bdy = [(row,col)]
        count = 1
        while bdy:
            if argv.verbose:
                print("%d:%d"%(count,len(bdy)), end=" ", flush=True)
            _bdy = []
            while bdy:
                i0,j0 = bdy.pop()
                for g in gen:
                    i,j = g[i0], g[j0]
                    if mask[i,j]:
                        _bdy.append((i,j))
                        mask[i,j] = 0
                        count += 1
            bdy = _bdy
        if argv.verbose:
            print()
        print("[%s]" % (count//N), end="", flush=True)
        counts.append(count//N)
    print()

    counts.sort()
    print("orbits:", counts, len(counts))
    assert sum(counts) == N


def get_flags(Pauli, gen, n, ks):

    assert Pauli[0].shape == (2**n, 2**n)
    assert gen[0].shape == (2**n, 2**n)

    print("get_flags", n, ks, end=" ... ", flush=True)

    lookup = dict((g,i) for (i,g) in enumerate(Pauli))
    perms = []
    for g in gen:
        ig = ~g
        idxs = [lookup[ig*h*g] for h in Pauli]
        perm = Perm(idxs)
        perms.append(perm)
    gen = perms

    def mkpoint(stab):
        stab = list(stab)
        stab.sort()
        stab = tuple(stab)
        return stab

    def get_stab(n, k):
        I = Clifford(1).I
        X = Clifford(1).X()
        stab = []
        for i in range(n-k):
            op = [I]*n
            op[i] = X
            stab.append(reduce(matmul, op))
        stab = mulclose(stab)
        stab = [lookup[g] for g in stab]
        stab = mkpoint(stab)
        return stab


    flag = tuple(get_stab(n, k) for k in ks)
    del lookup

    flags = {flag}
    bdy = list(flags)
    while bdy:
        _bdy = []
        for g in gen:
            for flag in bdy:
                glag = tuple(mkpoint([g[i] for i in code]) for code in flag)
                if glag in flags:
                    continue
                flags.add(glag)
                _bdy.append(glag)
        bdy = _bdy
    flags = list(flags)
    flags.sort()
    lookup = {c:idx for (idx,c) in enumerate(flags)}

    perms = []
    for g in gen:
        perm = []
        for (i,flag) in enumerate(flags):
            glag = tuple(mkpoint([g[i] for i in code]) for code in flag)
            j = lookup[glag]
            perm.append(j)
        perm = Perm(perm)
        perms.append(perm)
    gen = perms

    print(gen[0].rank)

    return gen



def big_double_cosets(gen):
    # find_double_cosets without the big mask array
    assert isinstance(gen, list)
    N = gen[0].rank

    remain = set(range(N)) # rows
    counts = []
    while remain:
        row = remain.pop()

        # now find the orbit of (row,0)
        bdy = [(row,0)]
        found = set(bdy)
        count = 1
        while bdy:
            if argv.verbose:
                print("%d:%d"%(len(found),len(bdy)), end=" ", flush=True)
            _bdy = []
            while bdy:
                i,j = bdy.pop()
                for g in gen:
                    tgt = g[i], g[j]
                    if tgt in found:
                        continue
                    found.add(tgt)
                    _bdy.append(tgt)
                    if tgt[1] == 0:
                        remain.remove(tgt[0])
                        count += 1
            bdy = _bdy
        if argv.verbose:
            print()
        print("[%s]" % count, end="", flush=True)
        counts.append(count)
    print()

    counts.sort()
    print("orbits:", counts, len(counts))
    assert sum(counts) == N


def count_hecke(lgen, rgen):
    assert len(lgen) == len(rgen)
    assert isinstance(lgen, list)
    assert isinstance(rgen, list)
    M = lgen[0].rank # rows
    N = rgen[0].rank # cols

    pairs = list(zip(lgen, rgen))

    remain = set(range(M)) # rows
    counts = []
    while remain:
        row = remain.pop()

        # now find the orbit of (row,0)
        bdy = [(row,0)]
        found = set(bdy)
        count = 1
        while bdy:
            if argv.verbose:
                print("%d:%d"%(len(found),len(bdy)), end=" ", flush=True)
            _bdy = []
            while bdy:
                i,j = bdy.pop()
                for l,r in pairs:
                    tgt = l[i], r[j]
                    if tgt in found:
                        continue
                    found.add(tgt)
                    _bdy.append(tgt)
                    if tgt[1] == 0:
                        remain.remove(tgt[0])
                        count += 1
            bdy = _bdy
        if argv.verbose:
            print()
        #print("[%s]" % count, end="", flush=True)
        counts.append(count)
    #print()

    counts.sort()
    #print("orbits:", counts, len(counts))
    assert sum(counts) == M

    return len(counts)


def get_operators(lgen, rgen):
    "generate all hecke operators"
    assert len(lgen) == len(rgen)
    assert isinstance(lgen, list)
    assert isinstance(rgen, list)
    M = lgen[0].rank # rows
    N = rgen[0].rank # cols

    pairs = list(zip(lgen, rgen))

    remain = set(range(M)) # rows
    counts = []
    while remain:
        row = remain.pop()

        # now find the orbit of (row,0)
        bdy = [(row,0)]
        found = set(bdy)
        op = numpy.zeros((M, N), dtype=int)
        op[row, 0] = 1
        count = 1
        while bdy:
            if argv.verbose:
                print("%d:%d"%(len(found),len(bdy)), end=" ", flush=True)
            _bdy = []
            while bdy:
                i,j = bdy.pop()
                for l,r in pairs:
                    tgt = l[i], r[j]
                    if tgt in found:
                        continue
                    op[tgt] = 1
                    found.add(tgt)
                    _bdy.append(tgt)
                    if tgt[1] == 0:
                        remain.remove(tgt[0])
                        count += 1
            bdy = _bdy
        if argv.verbose:
            print()
        #print("[%s]" % count, end="", flush=True)
        yield op


def get_product(lgen, rgen):
    assert len(lgen) == len(rgen)
    assert isinstance(lgen, list)
    assert isinstance(rgen, list)
    M = lgen[0].rank # rows
    N = rgen[0].rank # cols

    gen = list(zip(lgen, rgen))

    remain = set(range(M)) # rows
    counts = []
    while remain:
        row = remain.pop()

        # now find the orbit of (row,0)
        bdy = [(row,0)]
        found = set(bdy)
        op = numpy.zeros((M, N), dtype=int)
        op[row, 0] = 1
        count = 1
        while bdy:
            if argv.verbose:
                print("%d:%d"%(len(found),len(bdy)), end=" ", flush=True)
            _bdy = []
            while bdy:
                i,j = bdy.pop()
                for l,r in gen:
                    tgt = l[i], r[j]
                    if tgt in found:
                        continue
                    op[tgt] = 1
                    found.add(tgt)
                    _bdy.append(tgt)
                    if tgt[1] == 0:
                        remain.remove(tgt[0])
                        count += 1
            bdy = _bdy
        if argv.verbose:
            print()
        #print("[%s]" % count, end="", flush=True)
        #yield op
        found = list(found)
        found.sort()
        lookup = {v:i for i,v in enumerate(found)}
        perms = [Perm([lookup[l[i],r[j]] for (i,j) in found]) for l,r in gen]
        G = Group(None, perms)
        yield G


def count_hecke_injections(lgen, rgen, max_count=None):
    "this counts (i believe/hope) a single entry in the table of marks"
    assert len(lgen) == len(rgen)
    assert isinstance(lgen, list)
    assert isinstance(rgen, list)
    M = lgen[0].rank # rows
    N = rgen[0].rank # cols

    pairs = list(zip(lgen, rgen))

    remain = set(range(M)) # rows
    counts = []
    while remain:
        row = remain.pop()

        # now find the orbit of (row,0)
        bdy = [(row,0)]
        found = set(bdy)
        count = 1
        while bdy and count == 1:
            if argv.verbose:
                print("%d:%d"%(len(found),len(bdy)), end=" ", flush=True)
            _bdy = []
            while bdy and count == 1:
                i,j = bdy.pop()
                for l,r in pairs:
                    tgt = l[i], r[j]
                    if tgt in found:
                        continue
                    found.add(tgt)
                    _bdy.append(tgt)
                    if tgt[1] == 0:
                        if tgt[0] in remain:
                            remain.remove(tgt[0])
                        count += 1
                        break
            bdy = _bdy
        if argv.verbose:
            print()
        #print("[%s]" % count, end="", flush=True)
        if count == 1:
            counts.append(count)
            if max_count:
                return 1
    return len(counts)


def test_hecke_GL32():

    #G = Group.alternating(5)
    #G = GL(3,2)
    G = GL32()

    gens = G.gens

    #count_hecke(G.gens, G.gens)

    #Hs = G.conjugacy_subgroups()
    Hs = list(reversed(G.parabolics))
    Hs.append(Group([G.identity]))

#    for H in Hs:
#        print(H)
#        if len(H) != 8:
#            continue
#        for h in H:
#            print(h)
    print(Hs)

    Xs = [G.action_subgroup(H) for H in Hs]

    dump_tom(G, Xs)



def test_hecke():

    #n = 4
    #G = Group.coxeter_bc(n)

    if 0:
        #s = Schreier.make_D(n)
        #s = Schreier.make_F(); n=4
    
        #H_3 = Schreier.make_reflection(3, {(0,1):5, (1,2):3}); n=3

        H_4 = Schreier.make_reflection(4, {(0,1):5, (1,2):3, (2,3):3}); n=4
        #s = H_4 # order 14400

        s = Schreier.make_reflection(2, {(0,1):16})
    
        print(len(s))
        #G = s.get_group()
        #print(G)
        #perms = G.gens
        perms = s.get_gens()
        print("gens:", len(perms))
        N = perms[0].n
        gens = [Perm([g[i] for i in range(N)]) for g in perms]
        print("gens:", len(gens))
        G = Group.generate(gens, verbose=True)
        #G = object()
        #G.gens = gens

    #G = Group.symmetric(n)

    #G = Group.alternating(6)
    G = GL32()

    gens = G.gens
    n = len(gens)
    print(G, "gens:", n)
    assert n < 10

    if 0:
        Hs = G.conjugacy_subgroups()
    elif 1:
        Hs = list(G.parabolics)
        Hs.append(Group([G.identity]))
    else:
    
        Hs = []
        for idxs in all_subsets(len(G.gens)):
            if len(idxs)==0:
                #Hs.append(Group([G.identity]))
                H = Group([G.identity])
            else:
                H = Group.generate([gens[i] for i in idxs])
    
            for J in Hs:
                if G.is_conjugate_subgroup(H, J):
                    break
            else:
                Hs.append(H)

    Hs.sort(key = len, reverse=True)
    print("Hs:", Hs)
    Xs = []
    for H in Hs:
        X = G.act_subgroup(H)
        print(X)
        Xs.append(X)

    dump_tom(G, Xs)


def dump_tom(G, Xs):
    tgts = []
    for X in Xs:

        if isinstance(X, Group):
            Xgens = X.gens
        else:
            lgens = []
            lookup = X.src.lookup
            #print(X.send_perms)
            Xgens = []
            for g in G.gens:
                i = lookup[g]
                j = X.send_perms[i]
                g = X.tgt[j]
                Xgens.append(g)
        tgts.append(Xgens)

    N = len(Xs)
    for i in range(N):
      for j in range(N):
        lgens = tgts[i]
        rgens = tgts[j]
        c = count_hecke(lgens, rgens)
        print("%2s"%c, end=" ")
      print()

    rows = []

    print()
    print("table of marks:")
    for i in range(N):
      row = []
      for j in range(N):
        lgens = tgts[i]
        rgens = tgts[j]
        c = count_hecke_injections(lgens, rgens)
        row.append(c)
        print("%2s"%(c or '.'), end=" ")
      print()
      rows.append(row)

    #return

    tom = Tom(rows)
    print(tom)
    ops = tom.names
    for a in ops:
        row = tom[a]
        for b in ops:
            print("%s*%s=%s"%(a,b,tom.get_desc(tom[a]*tom[b])), end=" ")
        print()




class Builder:
    "Incrementally _build parts of the table of marks (TOM) for a group"
#    def __init__(self, G):
#        self.G = G
#        self.gens = G.gens
#        e = G.identity
#        G0 = Group([e], [e for g in self.gens])
#        Xs = [G.act_subgroup(H) for H in [G, G0]]
#        self.Xs = Xs

    def __init__(self, ngens):
        self.ngens = ngens # number of gens
        self.Xs = []

    def __getitem__(self, idx):
        return self.Xs[idx]

    def __len__(self):
        return len(self.Xs)

    def __contains__(self, X):
        assert isinstance(X, Group)
        assert len(X.gens) == self.ngens
        Xs = self.Xs
        #print("Builder.add", [Y.rank for Y in Xs])
        for Y in Xs:
            if Y.rank != X.rank:
                continue
            #print("\t", Y.rank, "?")
            a = X.count_hecke_injections(Y, max_count=1)
            if not a:
                continue
            #print("\t", a)
            b = Y.count_hecke_injections(X, max_count=1)
            #print("\t", b)
            if b:
                return True
        return False

    def add(self, X):
        assert isinstance(X, Group)
        assert len(X.gens) == self.ngens
        Xs = self.Xs
        if X in self:
            return

        idx = 0
        while idx < len(Xs):
            if X.rank <= Xs[idx].rank:
                Xs.insert(idx, X)
                return
            idx += 1
        Xs.append(X)

    def get_product(self, i, j):
        lgens = self[i].gens
        rgens = self[j].gens
        for X in get_product(lgens, rgens):
            #self.insert(X)
            yield X

    def get_tom(self, names=None, hecke=False, augment=True):

        N = len(self)

        print()
        print("get_tom():")
        if hecke:
            print("hecke:")
            for i in range(N):
              for j in range(N):
                lgens = self[i].gens
                rgens = self[j].gens
                c = count_hecke(lgens, rgens)
                print("%2s"%c, end=" ")
              print()
    
        rows = []

        if not N:
            return Tom(rows, names)
    
        #print([self[i].rank for i in range(N)])
        print("hecke injections:")
        for i in range(N):
          row = []
          for j in range(N):
            lgens = self[i].gens
            rgens = self[j].gens
            c = count_hecke_injections(lgens, rgens)
            row.append(c)
          if augment:
            row.append(self[i].rank)
          print(' '.join("%3s"%(c or '.') for c in row))
          rows.append(row)
        #if self[j].rank > c:
        if augment:
            row = [0]*(N+1)
            row[-1] = self[j].rank
            rows.append(row)
            print(' '.join("%3s"%(c or '.') for c in row))
        #print(rows)
    
        tom = Tom(rows, names)
        return tom

    def find_missing(self):
        tom = self.get_tom()
        ops = tom.names
        for i,a in enumerate(ops):
            for j,b in enumerate(ops):
                vec = tom[a]*tom[b]
                desc = tom.get_desc(vec)
                if desc is None:
                    return (i,j)

    def check(self):
        tom = self.get_tom()
        #print(tom)
        ops = tom.names
        for a in ops:
            row = tom[a]
            for b in ops:
                desc = tom.get_desc(tom[a]*tom[b])
                assert desc is not None
                #print("%s*%s=%s"%(a,b,desc), end=" ")
            #print()

    def dump(self, names=None):
        tom = self.get_tom(names)
        print()
        print(tom)
        print()
        ops = tom.names
        N = len(ops)
        for i in range(N):
            for j in range(i,N):
                a, b = ops[i], ops[j]
                vec = tom[a]*tom[b]
                desc = tom.get_desc(vec)
                print("%s*%s=%s"%(a,b,desc), end=" ")
            print()



def test_builder():

    #G = GL32()
    G = Group.symmetric(4)
    n = len(G.gens)

    builder = Builder(n)
    builder.check()

    Hs = G.conjugacy_subgroups()
    for H in Hs:
        X = G.act_subgroup(H)

        builder.add(X)
        #builder.check()

        builder.get_tom()

        idx = builder.find_missing()
        print("missing:", idx)
        if idx is None:
            continue

        Xs = list(builder.get_product(*idx))
        for X in Xs:
            #print("add", X.rankstr())
            builder.add(X)
            builder.get_tom()
            print()


        #break

    builder.check()

    #names = "N P L F PPP PP LL LPP x".split()
    names = None
    #builder.dump(names)

    #for X in builder:
    #    print(len(G) // X.rank, X.rank)
        

def test_flags(): 

    n = argv.get("n", 2)
    if n==1:
        flags = [[], [0]]
    elif n==2:
        flags = [[], [1], [0], [0,1]]
    elif n==3:
        flags = [
            [],
            [2], [1], [0],
            [1,2], [0,2], [0,1],
            [0,1,2]]
    elif n==4:
        flags = [
            [],
            [3], [2], [1], [0],
            [2,3], [1,3], [0,3], [1,2], [0,2], [0,1],
            [1,2,3], [0,2,3], [0,1,3], [0,1,2],
            [0,1,2,3]]
    else:
        assert 0

    gen = get_pauli_gens(n)
    Pauli = mulclose(gen)
    Pauli = list(Pauli)
    assert len(Pauli) == 4**(n+1)
    print("Pauli:", len(Pauli))

    cliff_gens = get_clifford_gens(n)

    genss = [get_flags(Pauli, cliff_gens, n, flag) for flag in flags]
    Xs = [Group(gens=gens, build=False) for gens in genss]

    method = argv.get("meth", "count_hecke")

    for X in Xs:
      for Y in Xs:
        f = getattr(X, method)
        c = f(Y)
        print("%3s"%(c or '.'), end=' ', flush=True)
      print()

    print()

#    builder = Builder()
#    for X in Xs:
#        builder.add(X)

    if n>2:
        return

    gens = genss[-1]
    total = None
    for op in get_operators(gens, gens):
        #print(shortstr(op))
        #print(op.sum(0)[0], op.sum(1)[0])
        total = op if total is None else (total + op)
    #print(total)
    rhs = numpy.zeros(total.shape, dtype=int)
    rhs[:] = 1
    assert numpy.all(total == rhs)

    if n>2:
        return

    #G = mulclose(cliff_gens, verbose=True)
    #G = list(G)
    #G.sort(key = str)
    #print(len(G))

    #lookup = {g:i for (i,g) in enumerate(G)}
    #gens = [Perm([lookup[gen*g] for g in G]) for gen in cliff_gens]
    #G = Group(gens=gens, verbose=True) # XXX

    builder = Builder(len(gens))
    builder.check()

    for X in Xs:

        builder.add(X)
        #builder.check()

        builder.get_tom()

        idx = builder.find_missing()
        print("missing:", idx)
        if idx is None:
            continue

        Xs = list(builder.get_product(*idx))
        for X in Xs:
            #print("add", X.rankstr())
            builder.add(X)
            builder.get_tom()
            print()

    for X in builder.get_product(3,3):
        print("add", X.rank)
        builder.add(X)

    builder.get_tom(hecke=True)

        #break
    #names = "N P L F PPP PP LL LPP x".split()
    builder.check()
    names = None
    builder.dump(names)
        


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

