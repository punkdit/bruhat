#!/usr/bin/env python

from time import sleep, time
from functools import reduce
from operator import matmul

from bruhat.solve import array2, shortstr, dot2, linear_independent

from bruhat.action import Perm, Group, Coset, mulclose, close_hom

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

def render_1():

    cvs = Canvas()
    cvs.append(Scale(5., 5.))

    point = lambda x, y : Mat([x, y])

    dx, dy = 2, 1
    ps = {}
    for row in range(5):
      for col in range(3):
        ps[row, col] = point((col-0.25*row)*dx, row*dy)

    st = [grey]+st_Thick
    circle = lambda p : cvs.fill(path.circle(*p, 0.03))

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
        pth = path.circle(*p, 0.05)
        cvs.fill(pth, [cl])

    cvs.writePDFfile("dessin_1.pdf")


def render_2():

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
    circle = lambda p : cvs.fill(path.circle(*p, 0.03))

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
        pth = path.circle(*p, 0.05)
        cvs.fill(pth, [blue])
    for p in [p01, p10, p12, p14, p16, p21]:
        pth = path.circle(*p, 0.05)
        cvs.fill(pth, [green])
    pth = path.circle(*p13, 0.05)
    cvs.fill(pth, [red])

    cvs.writePDFfile("dessin_2.pdf")

    
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




if __name__ == "__main__":

    start_time = time()

    test()
    test_real_pauli()
    #render_1()
    #render_2()

    #main()

    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")

