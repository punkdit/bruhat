#!/usr/bin/env python


#  copied & modified from quantale.py 

discussion = """

I was calling this bicategory a "relational algebroid"
because Tarski's axioms easily generalize to the
many-object setting. There's also this Hoffnung
paper, https://arxiv.org/abs/1007.1931 where he
calls it the Hecke bicategory. 

The objects are the finite G-sets, and the morphisms
are the Hecke relations. (This could be generalized
to groupoids and spans, which is what we talked about
back in January of 2024.)

Anyway, back to our old friend the categorified Gram-Schmidt process
https://ncatlab.org/nlab/show/Gram-Schmidt+process#CategorifiedGramSchmidtProcess

the matrix of hom's there is:

1 1 1  1  1
1 2 2  3  4
1 2 3  4  6
1 3 4  7 12
1 4 6 12 24

these numbers count the primitive (atomic) Hecke relations between S_4-sets
corresponding to the parabolic subgroups of S_4. 

Let's name the objects A,D,E,I,K :

 A  D  E  I  K   
 1  1  1  1  1  A
 1  2  2  3  4  D
 1  2  3  4  6  E
 1  3  4  7 12  I
 1  4  6 12 24  K

For example, hom(D,I) has 3 atomic Hecke relations,
and so the entire hom(D,I) is then a boolean algebra with 2^3=8 elements.

Of particular interest are the Hecke injections.
These are the f in hom(D,I) such that f^t.f = 1_D.
It turns out that 2 of the 3 Hecke relations are injections.
Here is a table showing just the Hecke injections:

 A  D  E  I  K   
 1  1  1  1  1  A
 .  1  .  2  4  D
 .  .  2  2  6  E
 .  .  .  2 12  I
 .  .  .  . 24  K

When you do this for all the subgroups of S_4, not just the
parabolics, you get 11 objects (up to isomorphism), and the
corresponding table is called the table of marks.
Here it is:

 A  B  C  D  E  F  G  H  I  J  K   
 1  1  1  1  1  1  1  1  1  1  1  A *
 .  2  .  .  .  2  .  2  .  2  2  B
 .  .  1  .  1  3  1  .  1  3  3  C
 .  .  .  1  .  .  .  1  2  .  4  D *
 .  .  .  .  2  .  .  .  2  2  6  E *
 .  .  .  .  .  6  .  .  .  6  6  F
 .  .  .  .  .  .  2  .  .  2  6  G
 .  .  .  .  .  .  .  2  .  .  8  H
 .  .  .  .  .  .  .  .  2  . 12  I *
 .  .  .  .  .  .  .  .  .  4 12  J
 .  .  .  .  .  .  .  .  .  . 24  K *

Or rather, this is a conjecture that I don't know how to prove.
(The usual definition of this table is different, but it's
probably just a three line proof to show equivalence..)

I haven't yet explained why the table of marks is so important...

One question I might ask: do these Hecke injections form a sub-bicategory
of the Hecke bicategory somehow? 

from zelim.py:
    D_8:
      | A B  C  D  E    F    G  H
    --+----------------------------
    A | A B  C  D  E    F    G  H
    B | B 2B G  G  2E   H    2G 2H
    C | C G  2C G  H    2F   2G 2H
    D | D G  G  2D H    H    2G 2H
    E | E 2E H  H  2E+H 2H   2H 4H
    F | F H  2F H  2H   2F+H 2H 4H
    G | G 2G 2G 2G 2H   2H   4G 4H
    H | H 2H 2H 2H 4H   4H   4H 8H

    not the same letters as follows >>>
    swap G<->F

    
[Space(1, 'A'), Space(2, 'B'), Space(2, 'C'), Space(2, 'D'),
Space(4, 'E'), Space(4, 'F'), Space(4, 'G'), Space(8, 'H')] 8

 A  B  C  D  E  F  G  H
 1  1  1  1  1  1  1  1  A
 1  2  1  1  2  2  1  2  B
 1  1  2  1  1  2  2  2  C
 1  1  1  2  1  2  1  2  D
 1  2  1  1  3  2  2  4  E
 1  2  2  2  2  4  2  4  F
 1  1  2  1  2  2  3  4  G
 1  2  2  2  4  4  4  8  H

the table of marks (tom):
 A  B  C  D  E  F  G  H
 1  1  1  1  1  1  1  1  A
 .  2  .  .  2  2  .  2  B
 .  .  2  .  .  2  2  2  C
 .  .  .  2  .  2  .  2  D
 .  .  .  .  2  .  .  4  E
 .  .  .  .  .  4  .  4  F
 .  .  .  .  .  .  2  4  G
 .  .  .  .  .  .  .  8  H

Multiplication in the burnside ring *is* 
just pointwise multiplying the rows in the tom

So:
E*F = 
 .  .  .  .  2  .  .  4  E
 .  .  .  .  .  4  .  4  F
=
 .  .  .  .  .  .  . 16  
= 2*H

C*D =
 .  .  2  .  .  2  2  2  C
 .  .  .  2  .  2  .  2  D
=
 .  .  .  .  .  4  .  4  F

E*E =
 .  .  .  .  2  .  .  4  E
 .  .  .  .  2  .  .  4  E
=
 .  .  .  .  4  .  . 16
= 2*E + H


Back to the parabolic subgroups of S_4:

 A  D  E  I  K   
 1  1  1  1  1  A
 .  1  .  2  4  D
 .  .  2  2  6  E
 .  .  .  2 12  I
 .  .  .  . 24  K

we can read off products:
D*D = 
 .  1  .  4 16 = D+I
D*E = 2*I
D*I = 2*I+K
etc.

somehow we get a sub-burnside ring from just ADEIK:
 *        *  *           *     * 
 A  B  C  D  E  F  G  H  I  J  K   
 1  1  1  1  1  1  1  1  1  1  1  A *
 .  .  .  1  .  .  .  1  2  .  4  D *
 .  .  .  .  2  .  .  .  2  2  6  E *
 .  .  .  .  .  .  .  .  2  . 12  I *
 .  .  .  .  .  .  .  .  .  . 24  K *
the missing columns (BCFGHJ) become redundant...

"""


import sys, os
import random
from functools import reduce, lru_cache
cache = lru_cache(maxsize=None)

from operator import add, mul
from string import ascii_uppercase, ascii_lowercase
ascii_letters = ascii_uppercase + ascii_lowercase

import numpy

from bruhat.gset import Group, Perm, GL
from bruhat.solve import shortstr
from bruhat.smap import SMap
from bruhat.argv import argv


class Space:
    "just a finite set of points"
    def __init__(self, n, name="?", **kw):
        assert n >= 0
        self.n = n
        self.name = name
        self.__dict__.update(kw)

    def __str__(self):
        return self.name
    #__repr__ = __str__

    def __repr__(self):
        return "Space(%s, %r)"%(self.n, self.name)

    def __len__(self):
        return self.n

    def __getitem__(self, i):
        assert 0<=i
        if i < self.n:
            return i
        raise IndexError

    def get_identity(self):
        A = numpy.identity(self.n, dtype=int)
        return Relation(A, self, self)

assert list(Space(3)) == [0,1,2]


class Relation:
    "a (matrix) boolean relation between left and right spaces"
    def __init__(self, A, left, right):
        assert isinstance(left, Space)
        assert isinstance(right, Space)
        m = len(left)
        n = len(right)
        assert A.shape == (m, n), "%s != %s"%(A.shape, (m,n))
        assert 0 <= A.min()
        assert A.max() <= 1

        self.A = A
        self.left = left
        self.right = right
        self.hom = (left, right)
        self.shape = A.shape
        self.homstr = "%s<--%s"%(left, right)

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
        return shortstr(self.A, keepshape=True)+(" (%s <-- %s)"%(self.left, self.right))

    def __eq__(self, other):
        assert self.left is other.left, (self.left, other.left)
        assert self.right is other.right, (self.right, other.right)
        return numpy.all(self.A==other.A)

    def __hash__(self):
        key = (
            id(self.left),
            id(self.right),
            hash(self.A.tobytes()))
        return hash(key)

    def __lt__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        return numpy.all(self.A<=other.A) and not numpy.all(self.A==other.A)

    def __le__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        return numpy.all(self.A<=other.A)

    def __add__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        A = self.A + other.A
        A = numpy.clip(A, 0, 1)
        return Relation(A, self.left, self.right)
    __or__ = __add__

    def __xor__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        A = (self.A + other.A) % 2
        return Relation(A, self.left, self.right)

    def __sub__(self, other):
        assert self.left is other.left
        assert self.right is other.right
        A = self.A - other.A
        A = numpy.clip(A, 0, 1)
        return Relation(A, self.left, self.right)

    def __neg__(self):
        A = 1 - self.A
        A = numpy.clip(A, 0, 1)
        return Relation(A, self.left, self.right)

    def __and__(self, other): # ??
        assert self.left is other.left
        assert self.right is other.right
        A = self.A * other.A
        return Relation(A, self.left, self.right)

    def __mul__(self, other):
        if isinstance(other, Relation):
            assert self.right is other.left
            A = numpy.dot(self.A, other.A)
            A = numpy.clip(A, 0, 1)
            return Relation(A, self.left, other.right)
        else:
            return self.get_left(other)

    def __rmul__(self, left):
        return self.get_right(left)

    @cache
    def transpose(self):
        A = self.A.transpose()
        A = A.copy() # clean up tobytes
        op = Relation(A, self.right, self.left)
        return op
    __invert__ = transpose

    @property
    def op(self):
        return self.transpose()

    # because we use object identity for __eq__, __hash__
    #the_star = ["*"]
    the_star = Space(1)

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




class Geometry:
    """
    The Hecke algebroid for a finite group G.
    This is a category with objects the transitive G-sets G/H,
    and morphisms the Hecke operators G/K <-- G/H.
    """
    def __init__(self, G, cgy=True, verbose=False):
        if cgy:
            Hs = G.conjugacy_subgroups(verbose=verbose, sort=True)
        else:
            Hs = list(G.subgroups(verbose=verbose))
            Hs.sort(key = lambda H:-len(H))
        if verbose:
            print("Geometry.__init__:", len(Hs))
    
        Xs = []
        for i,H in enumerate(Hs):
            #print(H, end=" ")
            X = G.action_subgroup(H)
            X.name = ascii_letters[i]
            X.space = Space(X.rank, X.name, X=X, H=H)
            Xs.append(X)
        self.G = G
        self.Hs = Hs
        self.Xs = Xs
        self.spaces = [X.space for X in Xs]
        #self.build_homs()

    def get_X(self, Y):
        for X in self.Xs:
            if X.isomorphic(Y):
                return X
        assert 0

    def get_name(self, XY):
        Xs = self.Xs
        names = []
        for nat in XY.get_atoms():
            Z = nat.src
            name = None
            for X in Xs:
                if nat.src.isomorphic(X):
                    assert name is None
                    name = X.name
            assert name is not None
            names.append(name)
        unique = list(set(names))
        unique.sort()
        name = '+'.join("%d*%s"%(names.count(name),name) for name in unique)
        #return "+".join(names)
        name = name.replace("+1*", "+")
        if name.startswith("1*"):
            name = name[2:]
        return name

    def get_hecke(self, cone):
        left, right = cone.legs
        XY = cone.apex
        m = left.tgt.rank
        n = right.tgt.rank
        lspace = left.tgt.space
        rspace = right.tgt.space
        #if lspace is None:
        #    lspace = self.get_X(left.tgt).space
        #if rspace is None:
        #    rspace = self.get_X(right.tgt).space
        rels = []
        for nat in XY.get_atoms():
            Z = nat.src
            M = numpy.zeros((m,n), dtype=int)
            for u in range(Z.rank):
                v = nat.send_items[u]
                i = left.send_items[v]
                j = right.send_items[v]
                M[i,j] = 1
            rel = Relation(M, lspace, rspace)
            rels.append(rel)
        return rels

    @property
    @cache
    def homs(self):
        homs = {}
        Xs = self.Xs
        N = len(Xs)
        for i in range(N):
          for j in range(N):
            X, Y = Xs[i], Xs[j]
            cone = X.mul(Y)
            rels = self.get_hecke(cone)
            homs[X.space, Y.space] = rels
        #self.homs = homs
        return homs




def check_geometry(geometry):
    homs = geometry.homs
    pairs = lambda hom : [(r,s) for r in hom for s in hom]
    for a in geometry.spaces:
        ai = a.get_identity()
        a_a = homs[a,a]
        for r in a_a:
            assert r*ai == r
            assert ~ ~r == r
        for r,s in pairs(a_a):
            assert ~(r*s) == (~s)*(~r)
            assert (~r)*(-(r*s)) + -s == -s
            for t in a_a:
                assert r+(s+t) == (r+s)+t
                assert r*(s*t) == (r*s)*t
                assert (r+s)*t == r*t + s*t

    for a in geometry.spaces:
      for b in geometry.spaces:
        a_b = homs[a,b]
        for r,s in pairs(a_b):
            assert r+s == s+r
            assert -(-r + s) + -(-r+ -s) == r
            assert ~(r+s) == (~r)+(~s)
            assert r & s == -(-r + -s)
            assert r - s == r & (-s)
            assert r ^ s == (r-s) + (s-r)

    for a in geometry.spaces:
      for b in geometry.spaces:
        for c in geometry.spaces:
            b_a = homs[b,a]
            c_b = homs[c,b]
            for r in c_b:
              for s in b_a:
                assert ~(r*s) == ~s * ~r
                assert (~r)*(-(r*s)) + -s == -s
                for t in b_a:
                    assert r*(s+t) == r*s + r*t
            for d in geometry.spaces:
                d_c = homs[d,c]
                for t in d_c:
                    assert t*(r*s) == (t*r)*s

    #print("check_geometry:", geometry.G, "OK")

            
def test_hecke():
    G = Group.alternating(4) # 5 spaces
    geometry = Geometry(G)
    check_geometry(geometry)

    #return

    #G = Group.alternating(5) # 9 spaces
    #G = Group.symmetric(4)   # 11 spaces
    #G = Group.symmetric(3)   # 4 spaces
    G = GL(3,2)
    #G = Group.dihedral(4)

    geometry = Geometry(G, True)
    homs = geometry.homs
    spaces = geometry.spaces
    #spaces = [s for s in spaces if s.name in "ADEIK"]
    #Hs = geometry.Hs
    Hs = [space.H for space in spaces]
    Xs = [space.X for space in spaces]
    lookup = {(X.name, Y.name):homs[X,Y] for X in spaces for Y in spaces}
    print(spaces, len(spaces))
    print()

    M = []

    tom = SMap()
    hecke = SMap()
    nu = SMap()
    for j,space in enumerate(spaces):
        tom[0, 3*j] = "%2s"%space.name
        tom[1+j, 3*len(spaces)+1] = space.name
    hecke = tom.copy()
    nu = tom.copy()
    table = tom.copy()
    for i,u in enumerate(spaces):
      for j,v in enumerate(spaces):
        hecke[1+i, 3*j] = "%2s"%(len(homs[u,v]))
        count = len([f for f in homs[u,v] if f.op*f==v.get_identity()])
        table[1+j, 3*i] = "%2s"%(count if count else ".")
    #for j,H in enumerate(Hs):
    #    conjs = {H.conjugate(g) for g in G}
    #    for i,K in enumerate(Hs):
    #        count = len([h for h in conjs if K.is_subgroup(h)])
    #        nu[1+i, 3*j] = "%2s"%(count if count else ".")
    for i,X in enumerate(Xs):
        sig = X.signature(Hs)
        #print("".join([" %2s"%(i if i else ".") for i in sig]))
        for j,u in enumerate(sig):
            tom[1+i, 3*j] = "%2s"%(u if u else ".")
        M.append(list(sig))
    
    print(tom)
    print()
    print(hecke)
    print()
    #print(nu)
    #print()
    print(table)
    print()

    assert tom == table

    from sage.all_cmdline import matrix
    M = matrix(M)
    #U = M.inverse()
    #print(len(G)*U)
    #print(M*M)

    return

    #A, B, C, D, E, F, G, H, I = spaces
    A,D,E,I,K = spaces
    for f in homs[I,D]:
        print(f)
        print(f.op * f == D.get_identity())
        print()






if __name__ == "__main__":
    from time import time
    start_time = time()

    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))

    

    
