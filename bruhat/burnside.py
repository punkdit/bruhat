#!/usr/bin/env python3

"""
Find the rank of the burnside ring for cyclic, dihedral and dicyclic groups.
"""

import string

from bruhat.action import mulclose, Perm, Group, conjugacy_subgroups, CFunc
from bruhat.action import burnside as _burnside
from bruhat.smap import tabulate
from bruhat.util import cross, write
from bruhat.argv import argv
from bruhat.element import cyclotomic, PolynomialRing, Linear, cayley, Z


def burnside(G):

    def print(*args, **kw):
        pass
    write = print

    Hs = conjugacy_subgroups(G)
    letters = list(string.ascii_uppercase + string.ascii_lowercase)
    letters.remove("O")
    letters.remove("o")
    letters = letters + [l+"'" for l in letters] + [l+"''" for l in letters]
    assert len(letters) >= len(Hs)
    letters = letters[:len(Hs)]

    homs = []
    rows = []
    for i, H in enumerate(Hs):
        cosets = G.left_cosets(H)
        assert len(G) == len(cosets)*len(H)

        hom = G.left_action(cosets)
        assert hom.src is G
        hom.name = letters[i]
        H.name = hom.name
        homs.append(hom)
        assert len(hom.components())==1 # transitive
        #print hom.tgt.perms

        row = [hom.name, len(H), len(cosets), len(H.conjugates)]
        print("%s subgroup order = %d, number of cosets = %d, conjugates = %d" %
            tuple(row), end="")
        if H.is_cyclic():
            print(", cyclic")
            row.append(r"\checkmark")
        else:
            print()
            row.append(r"")
        rows.append(row)

    T = []
    for hom in homs:
        chi = CFunc.from_action(hom)
        #print(hom.name, chi)
        T.append(chi)

    table = {}
    width = 0

    for i in range(len(homs)):
      for j in range(i, len(homs)):
        A = homs[i]
        B = homs[j]
        C = A.pushout(B)
        assert C.src is G
        write("%s*%s ="%(A.name, B.name))
        names = []
        for hom in C.components():
            assert hom.src is G
            name = '?'

            # We know it must be one of these possibilities:
            possible = [hom1 for hom1 in homs if not hom.refute_isomorphism(hom1, Hs)]
            assert possible
            if len(possible)==1:
                name = possible[0].name
                #write('.')

            else:
                assert 0
                return
                #write('/')
                for hom1 in possible:
                    if hom.isomorphic(hom1):
                        name = hom1.name
                        break

            names.append(name)
            if name == '?':
                assert 0

        #print len(C.components())
        uniq = set(names)
        counts = dict((name, 0) for name in uniq)
        for name in names:
            counts[name] += 1
        #print "+".join(names)
        names = list(uniq)
        names.sort()
        ss = []
        for name in names:
            c = counts[name]
            ss.append(name if c==1 else "%s*%s"%(c, name))
        s = '+'.join(ss)
        write(s+' ')
        width = max(width, len(s)+2)
        table[A.name, B.name] = s
        table[B.name, A.name] = s # commutative
      print()

    space = 1
    table = dict((k, v.replace("*", "")) for (k, v) in table.items())

    rows = cols = [hom.name for hom in homs]

    if argv.latex:
        print()
        print("$$")
        print(latex_table(table, rows, cols, upper=argv.get("upper")))
        print("$$")

    s = str(tabulate(table, rows, cols, space))
    print(s)

    import zelim
    import numpy
    A = zelim.parse(s)
    print(repr(A))

    if argv.latex:
        print()
        print(r"\noindent Table of multiplicities:")
        print()
        s = latex_dump(cols, A, sider=cols, sep=False)
        print(s)

    L, B = zelim.zelim(A)
    if argv.latex:
        print()
        print(r"\noindent Upper triangular form:")
        print()
        s = latex_dump(cols, B, sep=False)
        s = s.replace("0 ", ". ")
        print(s)
    else:
        print("H~:")
        print(B)
#        print("UMU:")
#        UMU = numpy.dot(L, numpy.dot(A, L.transpose()))
#        print(UMU)

    print("rank:", len(B))

    # diagonal matrix
    LT = numpy.dot(L, T)
    for chi1 in LT:
      for chi2 in LT:
        print(chi1.dot(chi2), end=" ")
      print()

    #return len(B)
    return LT


def make_dicyclic(n, debug=False):
    # Binary dihedral group
    # https://people.maths.bris.ac.uk/~matyd/GroupNames/dicyclic.html

    p = cyclotomic(PolynomialRing(Z), 4*n)
    ring = PolynomialRing(Z) / p
    rzeta = ring.x
    zeta = rzeta**2
    izeta = zeta**(2*n-1) # inverse

    for k in range(1, 6*n):
        if zeta**k == 1:
            break
    assert k==2*n

    i = rzeta**n
    assert i != 1
    assert i**2 == -1
    assert i**3 == -i
    assert i**4 == 1

    assert zeta**(2*n) == 1
    assert izeta != 1
    assert zeta*izeta == 1

    # we use the natural 2-dim representation
    GL = Linear(2, ring)
    e = GL.get([[1, 0], [0, 1]])
    r = GL.get([[zeta, 0], [0, izeta]])
    s = GL.get([[0, -1], [1, 0]])

    Di = mulclose([r, s])
    Di = list(Di)
    assert len(Di) == 4*n

    eidx = Di.index(e)
    ridx = Di.index(r)
    sidx = Di.index(s)

    G = cayley(Di)
    e = G[eidx]
    r = G[ridx]
    s = G[sidx]

    assert r**(2*n) == e
    assert s**2 == r**n
    assert r*s*r == s

    # complex characters ----------------------------------------------------

    # build the transpose of the character table.
    cols = []
    elems = [] # group elements with order the cols
    cgy = [] # element of each cgy class
    kidxs = range(2, n)

    desc = []

    # class: e, size 1
    elems.append(e)
    cgy.append(e)
    desc.append("e")
    cols.append([1, 1, 1, 1, 2] + [2 for k in kidxs])

    # class: s^2, size 1
    cols.append([1, 1, (-1)**n, (-1)**n, -2] + [((-1)**k)*2 for k in kidxs])
    elems.append(s*s)
    cgy.append(s*s)
    desc.append("s^2")

    for idx in range(1, n):
        # class: r**idx, size 2
        col = [1, 1, (-1)**idx, (-1)**idx, zeta**idx+izeta**idx]+[zeta**(idx*k) + izeta**(idx*k) for k in kidxs]
        cols.append(col)
        cols.append(col)
        elems.append(r**idx)
        cgy.append(r**idx)
        desc.append("r^%d"%idx)
        elems.append(r**(2*n-idx))

    for idx in range(n):
        # class: s, size n
        cols.append([1, -1, i**n, -i**n, 0] + [0 for k in kidxs])
        elems.append(s*(r**(2*idx)))
    cgy.append(s)
    desc.append("s")

    for idx in range(n):
        # class: s*r, size n
        cols.append([1, -1, -i**n, i**n, 0] + [0 for k in kidxs])
        elems.append(s*(r**(2*idx+1)))
    cgy.append(s*r)
    desc.append("sr")

    #print(i, -i, i**n, -i**n)

    assert len(set(elems)) == len(elems) == len(G)
    assert set(elems) == set(G)

    for col in cols:
        assert len(col) == n+3, len(col)

    if debug:
        print()
        for s in desc:
            s = "%10s"%s
            print(s, end=" ")
        print()
        c_chars = []
        for k in range(n+3):
            chi = dict((elems[i], cols[i][k]) for i in range(4*n))
            f = CFunc(G, chi)
            for g in cgy:
                s = "%10s"%(f[g])
                print(s, end=" ")
            print()
            c_chars.append(f)
        print()
    
        for f in c_chars:
          for g in c_chars:
            val = f.dot(g, normalize=False)
            #assert (f == g) == (val == 1)
            print(val, end=" ")
          print()

    else:
        c_chars = []
        for k in range(n+3):
            chi = dict((elems[i], cols[i][k]) for i in range(4*n))
            f = CFunc(G, chi)
            c_chars.append(f)
    
        for f in c_chars:
          for g in c_chars:
            val = f.dot(g)
            assert (f == g) == (val == 1)

    G.c_chars = c_chars

    # real characters ----------------------------------------------------

    # XXX NOT FINISHED XXX

    # build the transpose of the character table.
    cols = []
    elems = [] # group elements with order the cols
    cgy = [] # element of each cgy class
    kidxs = range(2, n)

    desc = []

    # class: e, size 1
    elems.append(e)
    cgy.append(e)
    desc.append("e")
    cols.append([1, 1, 1, 1, 2] + [2 for k in kidxs])

    # class: s^2, size 1
    cols.append([1, 1, (-1)**n, (-1)**n, -2] + [((-1)**k)*2 for k in kidxs])
    elems.append(s*s)
    cgy.append(s*s)
    desc.append("s^2")

    for idx in range(1, n):
        # class: r**idx, size 2
        col = [1, 1, (-1)**idx, (-1)**idx, zeta**idx+izeta**idx]+[zeta**(idx*k) + izeta**(idx*k) for k in kidxs]
        cols.append(col)
        cols.append(col)
        elems.append(r**idx)
        cgy.append(r**idx)
        desc.append("r^%d"%idx)
        elems.append(r**(2*n-idx))

    for idx in range(n):
        # class: s, size n
        cols.append([1, -1, i**n, -i**n, 0] + [0 for k in kidxs])
        elems.append(s*(r**(2*idx)))
    cgy.append(s)
    desc.append("s")

    for idx in range(n):
        # class: s*r, size n
        cols.append([1, -1, -i**n, i**n, 0] + [0 for k in kidxs])
        elems.append(s*(r**(2*idx+1)))
    cgy.append(s*r)
    desc.append("sr")

    #print(i, -i, i**n, -i**n)

    assert len(set(elems)) == len(elems) == len(G)
    assert set(elems) == set(G)

    for col in cols:
        assert len(col) == n+3, len(col)

    # XXX NOT FINISHED XXX

    if debug:
        print()
        for s in desc:
            s = "%10s"%s
            print(s, end=" ")
        print()
        r_chars = []
        for k in range(n+3):
            chi = dict((elems[i], cols[i][k]) for i in range(4*n))
            f = CFunc(G, chi)
            for g in cgy:
                s = "%10s"%(f[g])
                print(s, end=" ")
            print()
            r_chars.append(f)
        print()
    
        for f in r_chars:
          for g in r_chars:
            r = f.dot(g, normalize=False)
            print(r, end=" ")
          print()

    else:
        r_chars = []
        for k in range(n+3):
            chi = dict((elems[i], cols[i][k]) for i in range(4*n))
            f = CFunc(G, chi)
            r_chars.append(f)
    
        for f in r_chars:
          for g in r_chars:
            r = f.dot(g)
            assert (f == g) == (r != 0)

    G.r_chars = r_chars

    return G


def make_cyclic(n, check=False):

    R = PolynomialRing(Z)

    p = cyclotomic(R, n)
    #print(p)
    #print(p(R.x*2))

    ring = R / p
    zeta = ring.x
    #print("zeta: %r"%zeta)
    izeta = zeta**(n-1) # inverse

    #print(zeta, izeta)
    #print(zeta.conj())

    #G = mulclose([zeta])
    #for g in G:
    #    print(g, end=" ")
    #print()
    #G = cayley(G)

    items = list(range(n))
    perms = []
    e = Perm(dict((i, i) for i in range(n)), items)
    a = Perm(dict((i, (i+1)%n) for i in range(n)), items)
    perms = [e, a]
    for i in range(2, n):
        b = perms[-1] * a
        perms.append(b)

    G = Group(perms, items, check=True)

    # complex characters  --------------------------------------------------
    c_chars = []
    for k in range(n):
        chi = dict((perms[i], zeta**(i*k)) for i in range(n))
        chi = CFunc(G, chi)
        c_chars.append(chi)

    G.c_chars = c_chars

    for f in c_chars:
      for g in c_chars:
        r = f.dot(g)
        assert (f == g) == (r == 1)
        #print(f.dot(g), end=" ")
      #print()

    # real characters  --------------------------------------------------
    r_chars = [c_chars[0]] # triv
    if n%2:
        for k in range(1, (n+1)//2):
            a = c_chars[k] + c_chars[n-k]
            r_chars.append(a)
    else:
        chi = dict((perms[i], (-1)**i) for i in range(n))
        chi = CFunc(G, chi)
        r_chars.append(chi)
        for k in range(1, n//2):
            #print("a = c_chars[k] + c_chars[n-k]")
            a = c_chars[k] + c_chars[n-k]
            r_chars.append(a)

    for f in r_chars:
      for g in r_chars:
        r = f.dot(g)
        assert (f == g) == (r != 0)
        #print(f.dot(g), end=" ")
      #print()
    G.r_chars = r_chars

    return G


def main():

    n = argv.get("n", 5)

    #G = make_cyclic(n)
    G = make_dicyclic(n)
    print(len(G))

    q_chars = burnside(G)
    for f in q_chars:
        #print(f)
        for g in G.c_chars:
            val = f.dot(g)
            print(val, end=" ")
        print()



def do_rank():

#    for n in range(20, 40):
#        G = Group.cyclic(range(n))
#        rank = burnside(G)

    for n in range(25, 35):

        #print(n, end=" ")
        G = Group.cyclic(range(n))
        rank_C = burnside(G)
        #print(rank, end=" ")

        if n != 2:
            G = Group.dihedral(range(n))
            rank_D = burnside(G)
            #print(rank, end=" ")
        #else:
            #print("?", end=" ")

        G = make_dicyclic(n)
        rank_Dic = burnside(G)
        #print(rank, end=" ")
        #print()

        print(r"    %d  & %d   & %d  & %d  \\"%(n, rank_C, rank_D, rank_Dic))


if __name__ == "__main__":

    main()



