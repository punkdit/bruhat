#!/usr/bin/env python
"""
"""

import string, os
from random import randint, choice
from time import sleep, time
from functools import reduce
from operator import matmul

import numpy
from numpy import alltrue, zeros, dot

from bruhat.util import cross
from qupy.ldpc.solve import (
    array2, zeros2, shortstr, dot2, linear_independent, row_reduce, find_kernel,
    span, intersect, rank, enum2)
from bruhat.action import Perm, Group, Coset, mulclose, close_hom, is_hom
from bruhat.todd_coxeter import Schreier
from bruhat.argv import argv
from bruhat.smap import SMap


def make_dot(graph, cmd="neato"):
    n = len(graph)
    gens = graph.get_gens()

    f = open("schreier.dot", "w")
    labels = string.ascii_lowercase
    if len(labels) < n:
        labels = []
        i = 0
        while len(labels) < n:
            for c in string.ascii_lowercase:
                labels.append("%s%d"%(c,i))
            i += 1
    print("graph\n{\n", file=f)
    #cls = "green blue red".split()
    cls = "red blue green".split()
    for k in range(3):
        gen = gens[k]
        for i in gen.items:
            if i <= gen(i):
                print("    %s -- %s [color=%s];" % (
                    labels[i], labels[gen(i)], cls[k]), file=f)
    print("}\n", file=f)
    f.close()

    os.system("%s -Tpdf schreier.dot > schreier.pdf" % cmd)
    #data = os.popen("dot -Tpdf schreier.dot").read()
    #print("data:", data)
    #open("schreier.pdf", "w").write(data)



def build_group(idx=0, halve=False):
    # start with hyperbolic Coxeter reflection group: a--5--b--5--c
    # Then we add another generator "d" that halves these triangles.

    ngens = 3
    a, b, c, = range(ngens)
    rels_552 = [
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*5, (b,c)*5, (a,c)*2,
    ]

    if halve:
        d = ngens
        ngens += 1
        rels_552 += [ (d,)*2, (c,d,a,d), (a,d)*4, (b,d)*2, (d,c)*4, ]
        # Now we get:
        # d--4--a--5--b or,
        # b--5--c--4--d
        # where a equals c.

    rels = [
        (0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1),
        (0, 1, 0, 2, 1, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 0, 1,
            0, 2, 0, 1),
        (1, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 1, 2, 0, 1,
            0, 1, 0, 1, 0, 2),
        (0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 2, 1, 2, 1, 2, 0, 1,
            0, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 1, 0, 2, 0, 1,
            0, 2, 1, 2, 0, 2, 1, 2),
        (0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1, 0, 2, 0, 1,
            0, 2, 1, 2, 1, 2, 1, 2),
        (1, 2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 1, 0, 1, 0, 2, 1, 2,
            0, 1, 0, 1, 0, 2, 1, 2, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1),
        (1, 2, 1, 2, 0, 1, 0, 1, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2,
            1, 2, 0, 2, 1, 2, 0, 1, 0, 1, 0, 2, 1, 2, 0, 1, 0, 1,
            0, 1, 0, 2, 1, 2, 0, 2, 1, 2, 0, 1, 0, 1, 0, 1),
    ]

    for rel in rels:
        assert len(rel)%2 == 0

    rel = rels[idx]
    graph = Schreier(ngens, rels_552 + [rel])
    graph.build()
    G = graph.get_group()
    return G

gap = """
LoadPackage("recog");
LoadPackage("LINS");
F := FreeGroup("a", "b", "c");;
AssignGeneratorVariables(F);;
G := F/[a^2,b^2,c^2,(a*b)^5,(a*c)^2,(b*c)^4];;
gr:=LowIndexNormalSubgroupsSearchForAll(G, 600);
L := List(gr);
S := Grp(L[9]);
GeneratorsOfGroup(S);
"""

def parse(s, a, b, c):
    s = s.replace("^-1", "")
    s = s.replace("*", "+")
    s = s.replace("^", "*")
    s = s.replace(" ", "")
    return eval(s)


def build_group_524():
    ngens = 3
    a, b, c, = range(ngens)
    rels_524 = [
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*5, (b,c)*4, (a,c)*2,
    ]
    a, b, c = (a,), (b,), (c,)

    relss = []

    # order 160, non-checkerboard
    rels = [ 
        a+(b+a+c)*2+b+(c+a+b)*2+c, (a+c+b)*2+a+c+(b+c+a)*2+b, 
        (b+c+b+a)*2+(b+c+b+a)*2, (c+b+a+b)*2+(c+b+a+b)*2, 
        b+a+(b+a+c)*2+b+(c+a+b)*2+c+b, b+c+(b+a)*2+c+b+a+c+b+c+(a+b)*2+c,
    ]
    assert rels == parse("""
    [ a*(b*a*c)^2*b*(c^-1*a^-1*b^-1)^2*c^-1, (a*c*b)^2*a*c*(b^-1*c^-1*a^-1)^2*b^-1,
    (b*c*b*a)^2*(b^-1*c^-1*b^-1*a^-1)^2, (c*b*a*b)^2*(c^-1*b^-1*a^-1*b^-1)^2,
    b*a*(b*a*c)^2*b*(c^-1*a^-1*b^-1)^2*c^-1*b^-1,
    b*c*(b*a)^2*c*b*a*c^-1*b^-1*c^-1*(a^-1*b^-1)^2*c^-1 ]
    """, a, b, c)
    relss.append(rels)

    # order 240, checkerboard, [[30,8,3]] Bring's curve
    relss.append(parse("""
    [ (b*a*c)^3*(b^-1*c^-1*a^-1)^3, (c*b*a)^3*c^-1*(b^-1*c^-1*a^-1)^2*b^-1*a^-1, 
    ((b*a)^2*c)^2*(b^-1*a^-1*b^-1*c^-1*a^-1)^2, 
    b*a*b*c*(b*a)^2*c*b*a^-1*b^-1*c^-1*(a^-1*b^-1)^2*c^-1*b^-1*a^-1, 
    (b*a*c*b*a)^2*(b^-1*c^-1*a^-1*b^-1*a^-1)^2, 
    b*c*(b*a)^2*c*b*a*b*c^-1*(a^-1*b^-1)^2*c^-1*(b^-1*a^-1)^2, 
    b*(c*b*a)^3*c^-1*(b^-1*c^-1*a^-1)^2*b^-1*a^-1*b^-1, 
    a*b*(c*b*a)^3*c^-1*(b^-1*c^-1*a^-1)^2*(b^-1*a^-1)^2 ]
    """, a, b, c))

    idx = argv.get("idx", 0)
    graph = Schreier(ngens, rels_524 + relss[idx])
    graph.build(maxsize=100000)
    n = len(graph)
    print("|G| =", n)
    G = graph.get_group()
    return G


def build_group_624():
    ngens = 3
    a, b, c, = range(ngens)
    rels_624 = [
        (a,)*2, (b,)*2, (c,)*2,
        (a,b)*6, (b,c)*4, (a,c)*2,
    ]
    a, b, c = (a,), (b,), (c,)

    relss = []

    rels = parse("""
        [b*c*b*a*b*c*b^-1*a^-1*b^-1*c^-1*b^-1*a^-1, 
        c*b*a*b*c*b*a^-1*b^-1*c^-1*b^-1*a^-1*b^-1,
        a*c*b*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2, 
        b*a*c*b*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1,
        a*b*a*c*b*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^3,
        (b*a)^2*c*b*a*b*c*b^-1*a^-1*b^-1*c^-1*(a^-1*b^-1)^2*a^-1,
        (c*b*a)^2*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1,
        (a*c*b)^2*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1*a^-1,
        b*(c*b*a)^2*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1*b^-1,
        (b*a*c)^2*b*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1*a^-1*b^-1,
        a*(b*a*c)^2*b*a*b*c*b*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1*a^-1*b^-1*a^-1,
        c*b*a*(b*a*c)^2*b*a*b*c^-1*b^-1*(a^-1*b^-1*c^-1)^2*(a^-1*b^-1)^2 ]
    """, a, b, c)
    relss.append(rels) # [[30,11,3]]

    rels = parse("""
[ ((b*a)^2*c)^2*(b^-1*a^-1*b^-1*c^-1*a^-1)^2, (b*a*c*b*a)^2*(b^-1*c^-1*a^-1*b^-1*a^-1)^2,
  (c*(b*a)^2)^2*c^-1*b^-1*a^-1*b^-1*c^-1*(a^-1*b^-1)^2*a^-1,
  (b*a)^2*b*c*(b*a)^2*c*b^-1*a^-1*b^-1*c^-1*(a^-1*b^-1)^2*c^-1*b^-1*a^-1,
  b*a*b*c*(b*a)^2*c*b*a*b^-1*c^-1*(a^-1*b^-1)^2*c^-1*(b^-1*a^-1)^2,
  (b*a*b*c)^2*b*a*c*b^-1*c^-1*(a^-1*b^-1*c^-1*b^-1)^2*a^-1,
  b*c*(b*a)^2*b*c*b*a*b*c^-1*(b^-1*c^-1*a^-1)^2*b^-1*c^-1*b^-1*a^-1,
  b*(c*(b*a)^2)^2*(c^-1*b^-1*a^-1*b^-1)^2*a^-1*b^-1*a^-1,
  (b*c*b*a)^2*b*a*b*c^-1*b^-1*c^-1*a^-1*b^-1*c^-1*b^-1*a^-1*b^-1*c^-1*a^-1,
  (b*c*b*a)^2*c*b*a*c^-1*(b^-1*c^-1*b^-1*a^-1)^2*b^-1*a^-1,
  b*(c*b*a)^2*b*c*b*a*(c^-1*b^-1)^2*(a^-1*b^-1)^2*c^-1*b^-1*a^-1,
  c*b*a*(b*a*b*c)^2*(b^-1*c^-1*a^-1)^2*b^-1*c^-1*b^-1*a^-1*b^-1,
  (c*b*a*b)^2*a*b*c*b^-1*c^-1*a^-1*b^-1*c^-1*b^-1*a^-1*b^-1*c^-1*a^-1*b^-1,
  (c*b*a*b)^2*c*b*a*(c^-1*b^-1*a^-1*b^-1)^2*c^-1*b^-1*a^-1,
  c*(b*a*c*b*a)^2*b^-1*c^-1*a^-1*b^-1*(a^-1*b^-1*c^-1)^2*a^-1,
  c*b*a*(c*b*a*b)^2*c^-1*b^-1*c^-1*(a^-1*b^-1)^2*c^-1*b^-1*a^-1*b^-1,
  a*c*b*a*(b*a*b*c)^2*(b^-1*c^-1*a^-1)^2*b^-1*c^-1*(b^-1*a^-1)^2,
  (b*a*c)^2*b*a*b*c*b*a*(b^-1*c^-1*a^-1)^2*b^-1*a^-1*b^-1*c^-1*b^-1*a^-1,
  b*c*b*a*(b*a*c)^2*b*a*b^-1*c^-1*b^-1*(a^-1*b^-1*c^-1)^2*a^-1*b^-1*a^-1,
  c*b*a*(b*a*c)^2*b*a*b*c^-1*b^-1*(a^-1*b^-1*c^-1)^2*(a^-1*b^-1)^2,
  (c*b*a)^3*b*a*c*b^-1*a^-1*b^-1*c^-1*(b^-1*a^-1)^2*b^-1*c^-1*a^-1*b^-1,
  a*c*b*a*(b*a*c)^2*b*a*b*c^-1*b^-1*(a^-1*b^-1*c^-1)^2*(a^-1*b^-1)^2*a^-1 ]
    """, a, b, c)
    relss.append(rels) # [[30,11,3]]

    idx = argv.get("idx", 0)
    graph = Schreier(ngens, rels_624 + relss[idx])
    graph.build(maxsize=100000)
    n = len(graph)
    print("|G| =", n)
    G = graph.get_group()
    return G


class QCode(object):
    def __init__(self, H):
        m, n, k = H.shape
        assert k == 2
        assert H.max() <= 1
        self.H = H
        self.m = m
        self.n = n
        self.shape = m, n
        self.check()

    @property
    def flatH(self):
        H = self.H.view()
        m, n = self.shape
        H.shape = m, 2*n
        return H

    def check(self):
        H = self.H
        m, n = self.shape
        J = H.copy()
        J[:, :, :] = J[:, :, [1,0]]
        H = H.view()
        H.shape = J.shape = (m, 2*n)
        HJt = dot2(H, J.transpose())
        if HJt.sum() == 0:
            return
        print("not symplectic:")
        print(self)
        print(HJt)
        assert 0

    def dual(self):
        D = array2([[0,1],[1,0]])
        H = dot(self.H, D) % 2
        return QCode(H)

    def apply(self, idx, gate):
        H = self.H.copy()
        H[:, idx] = dot(self.H[:, idx], gate) % 2
        return QCode(H)

    def apply_H(self, idx):
        # swap X<-->Z on bit idx
        H = array2([[0,1],[1,0]])
        return self.apply(idx, H)

    def apply_S(self, idx):
        # swap X<-->Y
        S = array2([[1,1],[0,1]])
        return self.apply(idx, S)

    def apply_SH(self, idx):
        # X-->Z-->Y-->X 
        SH = array2([[0,1],[1,1]])
        return self.apply(idx, SH)

    def row_reduce(self):
        H = self.H.copy()
        m, n = self.shape
        H.shape = m, 2*n
        H = row_reduce(H)
        m, nn = H.shape
        H.shape = m, nn//2, 2
        return QCode(H)

    def get_logops(self):
        self = self.row_reduce()
        m, n = self.shape
        nn = 2*n
        H = self.flatH
        #print(shortstr(H))
        pivots = []
        row = col = 0
        while row < m:
            while col < nn and H[row, col] == 0:
                #print(row, col, H[row, col])
                pivots.append(col)
                col += 1
            row += 1
            col += 1
        while col < nn:
            pivots.append(col)
            col += 1
        #print(len(pivots))
        #s = ['.']*nn
        #for i in pivots:
        #    s[i] = "*"
        #print(''.join(s))
        #print(pivots)
        W = zeros2(len(pivots), nn)
        for i, ii in enumerate(pivots):
            W[i, ii] = 1
        #print()
        #print(shortstr(W))
        Ht = self.dual().flatH
        K = find_kernel(Ht)
        L = intersect(W, K)
        #print("L")
        #print(shortstr(L))
        kk = len(L)
        assert kk%2 == 0
        k = kk//2
        assert dot2(Ht, L.transpose()).sum() == 0 # commutes w stabilizers
        HL = numpy.concatenate((H, L))
        assert rank(HL) == m + 2*k # linearly independant
        return L

    def get_params(self):
        L = self.get_logops()
        kk = len(L)
        print("get_params: k=", kk//2)
        if kk > 22:
            return self.n, kk//2, None
        #Lt = L.transpose()
        #print("L:", L.shape)
        d = self.n
        #for v in span(L):
        #for v in enum2(kk):
        for v in numpy.ndindex((2,)*kk):
            v = dot2(v, L)
            count = 0
            for i in range(self.n):
              if v[2*i] or v[2*i+1]:
                count += 1
            if count:
                d = min(count, d)
        return self.n, kk//2, d

    def __str__(self):
        smap = SMap()
        m, n = self.shape
        H = self.H
        for i in range(m):
          for j in range(n):
            x, z = H[i, j]
            c = '.'
            if x and z:
                c = 'Y'
            elif x:
                c = 'X'
            elif z:
                c = 'Z'
            smap[i,j] = c
        return str(smap)

    def shortstr(self):
        H = self.H.view()
        H.shape = (self.m, 2*self.n)
        return shortstr(H)

    
def find_code(H):
    "build a QCode from face--vert adjacency matrix"
    import z3
    from z3 import Bool, And, Or, Not, Implies, If, Solver

    m, n = H.shape
    vs = {}
    xvars = {}
    zvars = {}
    yvars = {}

    clauses = []
    solver = Solver()

    add = solver.add

    for i in range(m):
      for j in range(n):
        if H[i, j] == 0:
            continue
        X = Bool("X_%d_%d"%(i,j))
        Z = Bool("Z_%d_%d"%(i,j))
        vs[i, j] = (X, Z)
        xvars[i,j] = X
        zvars[i,j] = Z
        yvars[i,j] = And(X, Z)
        # X or Z (or both)
        add(Or(X, Z))

    # each qubit has X & Z check
    for j in range(n): # cols
        idxs = [i for i in range(m) if H[i, j]] # rows
        assert len(idxs) > 1, idxs
        add(Or(*[xvars[i,j] for i in idxs]))
        add(Or(*[zvars[i,j] for i in idxs]))
        add(Not(And(*[yvars[i,j] for i in idxs])))
        # yes?

#    # commuting stabilizers
#    for i0 in range(m):
#      for i1 in range(i0+1, m):
#        bits = [j for j in range(n) if H[i0,j]==H[i1,j]>0]
#        if not bits:
#            continue
#        print(bits)
#        for lhs in cross([(0,1)]*len(bits)):
#          lvars = [vs[i0, k][l] for (k,l) in zip(bits, lhs)]
#          #print(lhs, lvars)
#          for rhs in cross([(0,1)]*len(bits)):
#            rvars = [vs[i1, k][r] for (k,r) in zip(bits, rhs)]
#            syndrome = len([idx for idx in range(len(bits)) if lhs[idx]!=rhs[idx]]) % 2
#            #print('\t', rhs, rvars, syndrome)
#            if syndrome:
#                clause = Not(And(*lvars, *rvars))
#                #print('\t', clause)
#                solver.add(clause)

    # commuting stabilizers
    for i0 in range(m):
      for i1 in range(i0+1, m):
        bits = [j for j in range(n) if H[i0,j]>0 and H[i1,j]>0]
        if not bits:
            continue
        #print(bits)
        for bits0 in cross([((1,0),(0,1),(1,1))]*len(bits)):
         for bits1 in cross([((1,0),(0,1),(1,1))]*len(bits)):
            syndrome = sum(l[0]*r[1]+l[1]*r[0] for (l,r) in zip(bits0, bits1))
            #print(bits0, bits1, syndrome)
            if syndrome%2 == 0:
                continue # allowed
            clause = []
            for j, jj in enumerate(bits):
                X0, Z0 = vs[i0, jj]
                X1, Z1 = vs[i1, jj]
                clause.append(X0 if bits0[j][0] else Not(X0))
                clause.append(Z0 if bits0[j][1] else Not(Z0))
                clause.append(X1 if bits1[j][0] else Not(X1))
                clause.append(Z1 if bits1[j][1] else Not(Z1))
            clause = Not(And(*clause)) # not allowed
            #print(str(clause).replace(" ", "").replace("\n", ""))
            solver.add(clause)

    result = solver.check()
    if result != z3.sat:
        print(result)
        return None

    H = zeros2(m, n, 2)
    model = solver.model()
    for (i, j), (X, Z) in vs.items():
        X = model.evaluate(X)
        Z = model.evaluate(Z)
        assert X or Z
        if X:
            H[i, j, 0] = 1
        if Z:
            H[i, j, 1] = 1
    code = QCode(H)
    return code


def test_surface():

    sides = argv.get("sides", 5)
    if sides==5:
        G = build_group_524()
    elif sides==6:
        G = build_group_624()
    else:
        assert 0

    if argv.dot:
        graph = Schreier(ngens, rels_524 + relss[idx])
        #graph.build([b+c+b+c]) # halve a vert
        graph.build([a+c]) # halve an edge
        make_dot(graph, "neato")
        return

    a, b, c = G.gens
    red, blue, green = G.gens

#    count = 0
#    for g in G:
#        if g.order() == 2:
#            count += 1
#    print("order 2 elements:", count)

    if 0:
        assert g.order() == 2
        H = [g, g*g]
        pairs = G.left_cosets(H)


    def get_adj(left, right):
        A = zeros2((len(left), len(right)))
        for i, l in enumerate(left):
          for j, r in enumerate(right):
            lr = l.intersection(r)
            A[i, j] = len(lr)>0
        return A

    if argv.fold_vert:
        g = blue*green*blue*green # divide a vert in half
        assert g.order() == 2
    elif argv.fold_edge:
        g = red*green # divide an edge in half
        assert g.order() == 2
    else:
        g = G.identity

    H = Group.generate([g])
    #pairs = G.left_cosets(H)
    act = G.action_subgroup(H)
    faces = act.orbits(Group.generate([blue, red]))
    print("faces:", len(faces))
    print([len(f) for f in faces])
    edges = act.orbits(Group.generate([green, red]))
    print("edges:", len(edges))
    print([len(f) for f in edges])
    verts = act.orbits(Group.generate([green, blue]))
    print("verts:", len(verts))
    print([len(f) for f in verts])

#    def setkey(items):
#        items = list(items)
#        items.sort(key=str)
#        return str(items)
#
#    def canonical(items):
#        # XXX FAIL: items is a list of set's
#        items = list(items)
#        items.sort(key=setkey)
#        return items
#
#    faces = canonical(faces)
#    edges = canonical(edges)
#    verts = canonical(verts)

    A = get_adj(faces, edges)
    B = get_adj(verts, edges)

    AA = dot(A, A.transpose())
    BB = dot(B, B.transpose())

    print("face -- face")
    print(AA.shape)
    print(AA)

    AB = dot(A, B.transpose())
    print("face -- vert")
    print(AB.shape)
    print(AB)

    if argv.find_code:
        code = find_code(AB)
        if code is None:
            return
        print(code)
        n,k,d = code.get_params()
        print("[[%s,%s,%s]]"%(n,k,d))
        if d is None:
            print(shortstr(code.get_logops()))
        return

    labels = string.ascii_lowercase

    def show_incident(A):
        print("graph")
        print("{")
        for (i,j) in numpy.ndindex(A.shape):
            if A[i,j] and i<j:
                print("\t", labels[i], "--", labels[j], ";")
        print("}")
    #print("vert -- edge")

#    print("face -- face")
#    show_incident(AA)
#    print("vert -- vert")
#    show_incident(BB)
#    print("face -- edge")
#    show_incident(A)

    from bruhat.isomorph import Tanner, search
    U = Tanner.build2(A, B)
    V = Tanner.build2(A, B)
    perms = []
    for f in search(U, V):
        perms.append(f)
    print("autos:", len(perms))
    
    assert dot2(A, B.transpose()).sum() == 0
    Hx, Hz = A, B
    Hx = linear_independent(Hx)
    Hz = linear_independent(Hz)

    return

    # checkerboard subgroup
    H = mulclose([red, blue, blue*green*blue*green])
    is_checkerboard = 2*len(H)==len(G)
    print("is_checkerboard:", is_checkerboard)
    checkerboard = G.left_cosets(H)

    if 1:
        # orientation subgroup
        L = Group.generate([red*blue, blue*green])
        print("L:", len(L))
        orients = G.left_cosets(L)
        is_orientable = len(orients) == 2
        print("is_orientable:", is_orientable)
        orients = list(orients)

    # vert'ex subgroup
    H = mulclose([blue, green])
    verts = G.left_cosets(H)
    print("verts:", len(verts))
    n = len(verts) # qubits

    # face subgroup
    H = mulclose([red, blue])
    assert len(H) == 10
    faces = G.left_cosets(H)
    print("faces:", len(faces))

    # edge subgroup
    H = mulclose([green, red])
    assert len(H) == 4
    edges = G.left_cosets(H)
    print("edges:", len(edges))

    if not is_checkerboard:
        return

    mx = len(faces)//2 - 1
    mz = len(faces)//2 - 1
    print("k =", n-mx-mz)

#    for face in faces:
#        print(face)
#        for black in checkerboard:
#            print('\t', face.intersect(black))


    Hx = zeros2(mx, n)
    Hz = zeros2(mz, n)
    ix = iz = 0
    black, white = checkerboard
    black, white = [], []
    for face in faces:
        if face.intersect(checkerboard[0]):
            black.append(face)
        else:
            white.append(face)

    Hx = get_adj(black, verts)
    Hz = get_adj(white, verts)
    print(shortstr(Hx))
    print()
    print(shortstr(Hz))
    assert dot2(Hx, Hz.transpose()).sum() == 0


    


def make_hyperbolic():

    G = build_group(2)
    H = build_group(1)
    print("|G| =", len(G))
    print("|G.gens| =", len(G.gens))
    #a, b, c, d = G.gens

    count = 0
    for g in G:
        if g.order() != 2:
            continue
        for h in G:
            if g*h != h*g:
                break
        else:
            print(g)
            print()
            count += 1
    print(count)

    hom = dict(zip(G.gens, H.gens))
    hom = close_hom(hom)
    #assert is_hom(hom)

    kernel = []
    for g,h in hom.items():
        if h.is_identity():
            kernel.append(g)
    print(len(kernel))
    for g in kernel:
        if not g.is_identity():
            print(g)
        else:
            print("is_identity")


if __name__ == "__main__":

    start_time = time()


    profile = argv.profile
    name = argv.next()
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
        #test_real_pauli()
        #render_1()
        #render_2()
        #main()
        #test_dessins()
        #make_hyperbolic()



    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")


