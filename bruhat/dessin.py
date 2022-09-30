#!/usr/bin/env python
"""
build dessins d'enfants & their homology (CSS quantum codes),
and draw some pictures..
"""

import string, os
from random import randint, choice
from time import sleep, time
from functools import reduce
from operator import matmul

import numpy
from numpy import alltrue, zeros, dot

from bruhat.solve import array2, shortstr, dot2, linear_independent
from bruhat.action import Perm, Group, Coset, mulclose, close_hom, is_hom
from bruhat.todd_coxeter import Schreier
from bruhat.argv import argv
from bruhat.solve import shortstr, zeros2, dot2, array2, solve

from huygens.namespace import *
from huygens.pov import Mat

from qupy.dev._algebra import Algebra, build_algebra, Tensor


EPSILON=1e-8

pauli = build_algebra("IXZY", "X*X=I Z*Z=I Y*Y=-I X*Z=Y Z*X=-Y X*Y=Z Y*X=-Z Z*Y=-X Y*Z=X")
I = pauli.I
Z = pauli.Z
X = pauli.X

def make_op(h, op):
    n = len(h)
    ops = [I]*n
    for i, hi in enumerate(h):
        if hi:
            ops[i] = op
    op = reduce(matmul, ops)
    return op


def test_real_pauli():
    I = pauli.I
    X = pauli.X
    Y = pauli.Y
    Z = pauli.Z

    op = Tensor(pauli)
    assert op.get_keys() == []
    op[(0,)] = 1.
    assert op[(0,)] == 1.
    assert op[(1,)] == 0.
    assert op.get_keys() == [(0,)]

    assert str(op) == "I", repr(str(op))
    op = I*X
    assert str(op) == "X"

    zero = Tensor(pauli)
    assert zero.norm() == 0.
    assert (X-X) == zero
    assert str(X-X) == "0", str(X-X)
    assert (X-X).norm() < EPSILON, (X-X).norm()
    assert X==X

    assert I*X == X
    assert X*X==I
    assert Z*Z==I
    assert Y*Y==-I
    assert X*Z==Y
    assert Z*X==-Y
    assert X*Y==Z
    assert Y*X==-Z
    assert Z*Y==-X
    assert Y*Z==X

    II = I@I
    XI = X@I
    IX = I@X
    XX = X@X
    assert II+II == 2*II



def cycleperm(orbits):

    items = list(orbits.replace(" ", ""))
    items.sort()
    orbits = [list(orbit) for orbit in orbits.split()]

    perm = {}
    for orbit in orbits:
        n = len(orbit)
        for i in range(n):
            perm[orbit[i]] = orbit[(i+1)%n]
    perm = Perm(perm, items)
    return perm


def test():

    g = cycleperm("aceg bf d h")
    r = cycleperm("abcdefgh")
    b = cycleperm("adeh bg cf")

    assert r*g == b

    G = mulclose([g, r])
    assert len(G) == 64


#class Point(object):
#    def __init__(self, x, y):
#        self.v = 

def render_1(radius=0.07):

    cvs = Canvas()
    cvs.append(Scale(5., 5.))

    point = lambda x, y : Mat([x, y])

    dx, dy = 2, 1
    ps = {}
    for row in range(5):
      for col in range(3):
        ps[row, col] = point((col-0.25*row)*dx, row*dy)

    st = [grey]+st_Thick
    circle = lambda p : cvs.fill(path.circle(*p, 0.04))

    paths = []
    def stroke(p0, p1, rev=False): 
        if rev:
            p0, p1 = p1, p0
        pth = path.line(*p0, *p1)
        #cvs.stroke(pth, st) #+[deco.marrow])
        p01 = 0.5*(p0+p1)
        #cvs.text(*p01, "%d"%len(paths), st_center+[text.size.tiny])
        paths.append(pth)
        return len(paths)-1

    def curve(p0, p1, dx=0., dy=0., rev=False):
        if rev:
            p0, p1 = p1, p0
        p2 = 0.5*(p0 + p1) + Mat([dx, dy])
        pth = path.curve(*p0, *p2, *p2, *p1)
        #cvs.stroke(pth, st) #+[deco.marrow]) 
        #cvs.text(*p2, "%d"%len(paths), st_center+[text.size.tiny])
        paths.append(pth)
        return len(paths)-1

    def fill(i0, i1, i2):
        #pths = [paths[idx] for idx in idxs]
        p = paths[i0] + paths[i1] + paths[i2]
        i = 1
        while i < len(p.items):
            if isinstance(p.items[i], MoveTo):
                p.items.pop(i)
            else:
                i += 1
        cvs.fill(p, st)

    for i in range(5):
        stroke(ps[0, 1], ps[i, 2], i%2)
        stroke(ps[4, 1], ps[i, 0], i%2)
    for i in range(4):
        stroke(ps[i, 0], ps[i+1, 0], i%2)
        stroke(ps[i, 1], ps[i+1, 1])
        stroke(ps[i, 2], ps[i+1, 2], i%2)
    stroke(ps[0,0], ps[0,1])
    stroke(ps[4,1], ps[4,2])
    stroke(ps[0,0], ps[2,1], True)
    stroke(ps[2,1], ps[4,2], True)

    curve(ps[0,1], ps[2,1], -0.4, -0.4)
    curve(ps[0,1], ps[2,1], +0.3, +0.4, True)
    curve(ps[2,1], ps[4,1], -0.3, -0.4, True)
    curve(ps[2,1], ps[4,1], +0.4, +0.4)

    fill(7, 5, 16)
    fill(3, 1, 10)
    fill(24, 22, 26)
    fill(27, 11, 14)
    fill(8, 21, 6)
    fill(28, 17, 20)
    fill(25, 29, 23)
    fill(4, 15, 2)

    for pth in paths:
        cvs.stroke(pth, st_thin)

    for row in range(5):
      for col in range(3):
        if (row, col) in [(0, 1), (4, 1)]:
            cl = red
        elif (row, col) in [(1,0),(1,2),(2,1),(3,0),(3,2)]:
            cl = blue
        else:
            cl = green
        p = ps[row, col]
        pth = path.circle(*p, radius)
        cvs.fill(pth, [cl])

    cvs.writePDFfile("dessin_1.pdf")
    cvs.writeSVGfile("dessin_1.svg")


def render_2(radius=0.07):

    cvs = Canvas()
    cvs.append(Scale(5., 5.))

    point = lambda x, y : Mat([x, y])

    dx, dy = 2, 2
    p00 = 0*dx, 0*dy
    p01 = 1*dx, 0*dy
    p02 = 2*dx, 0*dy

    p10 = 0*dx, 1*dy
    p11 = (1/3)*dx, 1*dy
    p12 = (2/3)*dx, 1*dy
    p13 = 1*dx, 1*dy
    p14 = (4/3)*dx, 1*dy
    p15 = (5/3)*dx, 1*dy
    p16 = 2*dx, 1*dy

    p20 = 0*dx, 2*dy
    p21 = 1*dx, 2*dy
    p22 = 2*dx, 2*dy

    st = [grey]+st_Thick
    circle = lambda p : cvs.fill(path.circle(*p, 0.04))

    debug = False
    paths = []
    def stroke(p0, p1, rev=False): 
        p0, p1 = point(*p0), point(*p1)
        if rev:
            p0, p1 = p1, p0
        pth = path.line(*p0, *p1)
        if debug:
            cvs.stroke(pth, st+[deco.marrow])
        p01 = 0.5*(p0+p1)
        if debug:
            cvs.text(*p01, "%d"%len(paths), st_center+[text.size.tiny])
        paths.append(pth)
        return len(paths)-1

    def curve(p0, p1, dx=0., dy=0., rev=False):
        p0, p1 = point(*p0), point(*p1)
        if rev:
            p0, p1 = p1, p0
        p2 = 0.5*(p0 + p1) + Mat([dx, dy])
        pth = path.curve(*p0, *p2, *p2, *p1)
        if debug:
            cvs.stroke(pth, st+[deco.marrow]) 
            cvs.text(*p2, "%d"%len(paths), st_center+[text.size.tiny])
        paths.append(pth)
        return len(paths)-1

    def fill(i0=None, i1=None, i2=None):
        if i0 is None:
            n = len(paths)
            i0, i1, i2 = n-3, n-2, n-1
        p = paths[i0] + paths[i1] + paths[i2]
        i = 1
        while i < len(p.items):
            if isinstance(p.items[i], MoveTo):
                p.items.pop(i)
            else:
                i += 1
        cvs.fill(p, st)

    stroke(p00, p01)
    stroke(p01, p13)
    curve(p13, p00, 0., -0.2)
    fill()

    stroke(p02, p16)
    curve(p16, p13, 0., -0.8)
    curve(p13, p02, 0., -0.2)
    fill()

    curve(p13, p14)
    curve(p14, p15)
    curve(p15, p13, 0., -0.4)
    fill()

    curve(p13, p10, 0., -0.8)
    curve(p10, p11)
    curve(p11, p13, 0., -0.4)
    fill()

    curve(p11, p12)
    curve(p12, p13)
    curve(p13, p11, 0., +0.4)
    fill()

    curve(p13, p15, 0., +0.4)
    curve(p15, p16)
    curve(p16, p13, 0., +0.8)
    fill()

    curve(p13, p22, 0., 0.2)
    curve(p22, p21)
    curve(p21, p13)
    fill()

    curve(p13, p20, 0., +0.2)
    curve(p20, p10)
    curve(p10, p13, 0., +0.8)
    fill()

    stroke(p01, p02)
    stroke(p20, p21)
    stroke(p00, p10)
    stroke(p16, p22)

    for pth in paths:
        cvs.stroke(pth, st_thin)

    for p in [p00, p02, p11, p15, p20, p22]:
        pth = path.circle(*p, radius)
        cvs.fill(pth, [blue])
    for p in [p01, p10, p12, p14, p16, p21]:
        pth = path.circle(*p, radius)
        cvs.fill(pth, [green])
    pth = path.circle(*p13, radius)
    cvs.fill(pth, [red])

    cvs.writePDFfile("dessin_2.pdf")
    cvs.writeSVGfile("dessin_2.svg")

    
def build_code_1():
    faces = [
        (2, 6, 3),   #  0
        (0, 3, 4),   #  1
        (1, 4, 5),   #  2
        (5, 8, 2),   #  3
        (6, 10, 7),   #  4
        (7, 11, 12),   #  5
        (8, 9, 15),   #  6
        (9, 13, 14),   #  7
        (16, 19, 17),   #  8
        (17, 18, 10),   #  9
        (11, 18, 12),   # 10
        (19, 0, 20),   # 11
        (20, 1, 21),   # 12
        (22, 21, 16),   # 13
        (23, 22, 15),   # 14
        (13, 23, 14),   # 15
    ]
    verts = [
        (0, 19, 16, 21, 1, 5, 2, 3),
        (0, 4, 1, 20),
        (2, 8, 15, 22, 16, 17, 10, 6),
        (10, 18, 11, 7),
        (11, 12),
        (4, 3, 6, 7, 12, 18, 17, 19, 20, 21, 22, 23, 13, 9, 8, 5),
        (13, 14),
        (14, 23, 15, 9),
    ]

    return build_code(faces, verts, 24)
    

def build_code_2():
    faces = [
        (0, 3, 2),  # 0
        (3, 4, 12), # 1
        (4, 5, 15), # 2
        (5, 6, 20), # 3 
        (6, 7, 16),# 4 
        (7, 8, 13),# 5 
        (8, 9, 13),# 6 
        (1, 10, 9),# 7 
        (10, 11, 19),# 8 
        (11, 2, 14),  #  9
        (14, 12, 21),  #  10
        (16, 17, 0), # 11
        (17, 18, 23), # 12
        (18, 19, 23), # 13
        (21, 15, 22),  #  14
        (22, 20, 1), # 15
    ]
    verts = [
        (0, 2, 11, 10, 1, 20, 6, 16), # 0
        (0, 17, 23, 19, 11, 14, 21, 22, 1, 9, 8, 7, 6, 5, 4, 3), #  1
        (2, 3, 12, 14), #  2
        (8, 13), #  3
        (12, 4, 15, 21), #  4
        (7, 13, 9, 10, 19, 18, 17, 16), #  5
        (15, 5, 20, 22), #  6
        (18, 23), #  7
    ]
    return build_code(faces, verts, 24)


def build_code(faces, verts, n=24):
    Hz = []
    for face in faces:
        h = [0]*n
        for i in face:
            assert 0<=i<n, i
            h[i] = 1
        Hz.append(h)

    Hx = []
    for vert in verts:
        h = [0]*n
        for i in vert:
            assert 0<=i<n, i
            h[i] = 1
        Hx.append(h)
        
    Hz = array2(Hz)
    Hx = array2(Hx)
    #print(shortstr(Hz))
    #print()
    #print(shortstr(Hx))

    #print(shortstr(dot2(Hz, Hx.transpose())))
    assert dot2(Hz, Hx.transpose()).sum() == 0

    Hz = linear_independent(Hz)
    Hx = linear_independent(Hx)

    return Hz, Hx


def build_proj(Hz, Hx):
    n = Hz.shape[1]
    zgens = []
    for hz in Hz:
        op = make_op(hz, Z)
        zgens.append(op)

    xgens = []
    for hx in Hx:
        op = make_op(hx, X)
        xgens.append(op)

    # build projector
    In = reduce(matmul, [I]*n) 

    def mkproj(gens):
        P = In
        #print(P)
        for g in gens:
            #print(g)
            #P = (0.5)*P*(In + g)
            P = P*(In+g)
        return P
    #print()

    Px = mkproj(xgens)
    Pz = mkproj(zgens)
    P = Px*Pz
    return P



def get_wenum(P):
    counts = {}
    for k in P.i_keys():
        value = P[k].real
        nz = k.count(2) + k.count(3)
        nx = k.count(1) + k.count(3)
        key = nz, nx
        counts[key] = counts.get(key, 0) + int(value)
    ks = list(counts.keys())
    ks.sort()
    #print(len(counts))    
    for k in ks:
        print("%s:%s "%(k, counts[k]), end="")
    print()
    return counts


def main():

    Hz, Hx = build_code_1()
    P = build_proj(Hz, Hx)
    c1 = get_wenum(P)
    del P

    Hz, Hx = build_code_2()
    P = build_proj(Hz, Hx)
    c2 = get_wenum(P)
    del P

    for k in c2.keys():
        assert k in c1
    for k in c1.keys():
        assert k in c2

    ks = list(c1.keys())
    ks.sort()
    for k in ks:
        print(c1[k] - c2[k], end=" ")
    print()


def test_dessins():

    # red, green, blue reflections
    ngens = 3
    r, g, b = (0, 1, 2)
    rels = [(r, r), (g, g), (b, b)] # self-inverse

    def make(hgens, maxsize=10000, rev=False):
        hgens = [tuple('rgb'.index(item) for item in gen) for gen in hgens]
        if rev:
            hgens = [tuple(reversed(gen)) for gen in hgens]
        graph = Schreier(ngens, rels)
        graph.build(hgens, maxsize)
        return graph

    # The dessin for the Gaussian elliptic curve
    hgens = ["rgrb", "rb"*4, "rbgbrb", "grbr", "bgbg"]
    graph = make(hgens, rev=True)
    assert len(graph) == 8

    rev = lambda gen : tuple(reversed(gen))

    # this is dessin/code #1 above:
    hgens = [
        "bg"*8,  # around red
        "rb"*4, # around green
        "grgr", # around blue
        "grbrbg", # green 4
        "bgrgrgrgrb", # blue 5
        rev("bgbrgb"), # green 13
        rev("rbrgbrgrbr"), # green 5
        "gbgbrb",  # homology cycle
        "rbrgbgbgbg", # homology cycle
        #"gbgrgbgbgbgb", # another homology cycle.. not needed
    ]
    graph = make(hgens, rev=True)
    assert len(graph) == 16

    # make some random codes
    count = 0
    while count < 0:
        hgens = []
        for n in range(10):
            rel = []
            #for m in range(randint(8, 12)):
            for m in range(16):
                rel.append(choice('rgb'))
            hgens.append(tuple(rel))
        graph = make(hgens)
        N = len(graph)
        if N<5:
            continue
        elif N<1000:
            print(len(graph), end=" ", flush=True)
        else:
            print(".", end="", flush=True)
        count += 1
    print()

    # Klein Quartic: it's a homogeneous space
    red, green, blue = r, g, b
    a = (green, blue)
    b = (blue, red)
    bi = (red, blue)
    klein = rels + [a*3, b*7, (a+b)*2, (a+bi+bi)*4]
    graph = Schreier(ngens, klein)
    graph.build()
    assert len(graph) == 168*2


    # Bring's curve 
    ngens = 3
    a, b, c = range(ngens)
    rels = [
        (a,a), (b,b), (c,c),
        (a,b)*5, (b,c)*5, (a,c)*2,
        (3*(b,a)+(a,c))*3,
    ]
    graph = Schreier(ngens, rels)
    graph.build()
    assert len(graph) == 120

    graph = Schreier(ngens, rels)
    #graph.build([(a,b,a,b)]) # index 5 subgroup
    #graph.build([(a,b,c,b)]) # index 3
    graph.build([(a,b,c,a)]) # index 3
    n = len(graph)
    print("|G/H| =", n)


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
    print(n)

    if argv.dot:
        graph = Schreier(ngens, rels_524 + relss[idx])
        #graph.build([b+c+b+c]) # halve a vert
        graph.build([a+c]) # halve an edge
        make_dot(graph, "neato")
        return

    G = graph.get_group()
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

    g = blue*green*blue*green # divide a vert in half
    #g = red*green # divide an edge in half
    assert g.order() == 2
    H = Group.generate([g, g*g])
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

    def canonical(items):
        items = list(items)
        items.sort(key=str)
        return items

    faces = canonical(faces)
    edges = canonical(edges)
    verts = canonical(verts)

    A = get_adj(faces, edges)
    B = get_adj(verts, edges)

    AA = dot(A, A.transpose())
    BB = dot(B, B.transpose())
    print(BB.shape)
    print(BB)

    labels = string.ascii_lowercase

    def show_incident(A):
        print("graph")
        print("{")
        for (i,j) in numpy.ndindex(A.shape):
            if A[i,j] and i<j:
                print("\t", labels[i], "--", labels[j], ";")
        print("}")
    #print("face -- edge")
    #print("vert -- edge")

    #show_incident(AA)
    show_incident(BB)

    assert dot2(A, B.transpose()).sum() == 0
    Hx, Hz = A, B
    Hx = linear_independent(Hx)
    Hz = linear_independent(Hz)
    print(Hx.shape)
    print(Hz.shape)

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

    G = build_group_524()
    return

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

    #test()
    #test_real_pauli()
    #render_1()
    #render_2()
    #main()
    #test_dessins()
    make_hyperbolic()

    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")


