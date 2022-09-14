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

def render():

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
    stroke = lambda p0, p1: cvs.stroke(path.line(*p0, *p1), st)
    def curve(p0, p1, dx=0., dy=0.):
        p2 = 0.5*(p0 + p1) + Mat([dx, dy])
        cvs.stroke(path.curve(*p0, *p2, *p2, *p1), st) 

    for i in range(5):
        stroke(ps[0, 1], ps[i, 2])
        stroke(ps[4, 1], ps[i, 0])
    for i in range(4):
        stroke(ps[i, 0], ps[i+1, 0])
        stroke(ps[i, 1], ps[i+1, 1])
        stroke(ps[i, 2], ps[i+1, 2])
    stroke(ps[0,0], ps[0,1])
    stroke(ps[4,1], ps[4,2])
    stroke(ps[0,0], ps[2,1])
    stroke(ps[2,1], ps[4,2])

    curve(ps[0,1], ps[2,1], -0.4)
    curve(ps[0,1], ps[2,1], +0.4)
    curve(ps[2,1], ps[4,1], -0.4)
    curve(ps[2,1], ps[4,1], +0.4)

    for row in range(5):
      for col in range(3):
        p = ps[row, col]
        pth = path.circle(*p, 0.03)
        if row and col<2:
            cvs.fill(pth)
        else:
            cvs.stroke(pth)

    #cvs.writePDFfile("dessin.pdf")

    
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
    #render()

    main()

    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")


