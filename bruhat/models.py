#!/usr/bin/env python

import sys, os
from random import seed

import numpy
from numpy import concatenate
import scipy.sparse.linalg
from scipy import sparse

from solve import shortstr, shortstrx, parse, eq2, dot2, zeros2, array2, identity2
from solve import row_reduce, RowReduction, span, get_reductor
from solve import u_inverse, find_logops, solve, find_kernel, linear_independent
from solve import rand2, find_stabilizers, find_errors, enum2

from argv import Argv
argv = Argv()


def write(s):
    sys.stdout.write(str(s)+' ')
    sys.stdout.flush()


def genidx(shape):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx(shape[1:]):
                yield (idx,)+_idx


def check_conjugate(A, B):
    if A is None or B is None:
        return
    assert A.shape == B.shape
    I = numpy.identity(A.shape[0], dtype=numpy.int32)
    assert eq2(dot2(A, B.transpose()), I)


def check_commute(A, B):
    if A is None or B is None:
        return
    C = dot2(A, B.transpose())
    assert C.sum() == 0, "\n%s"%shortstr(C)



def build_gcolor(size):

    from qupy.ldpc import gcolor

    lattice = gcolor.Lattice(size)

    n = len(lattice.qubits)
    print lattice

    code = lattice.build_code(check=False)
    #Ex = lattice.Ex
    Gx, Gz = code.Gx, code.Gz
    Hx, Hz = code.Hx, code.Hz

    return Gx, Gz, Hx


def build_compass(li, lj=None):

    if lj is None:
        lj = li

    n = li*lj

    keys = [(i, j) for i in range(li) for j in range(lj)]
    coords = {}  
    for i, j in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            coords[i+di, j+dj] = keys.index(((i+di)%li, (j+dj)%lj))

    m = n 
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)

    idx = 0 
    for i in range(li):
      for j in range(lj):
        Gx[idx, coords[i, j]] = 1 
        Gx[idx, coords[i, j+1]] = 1 

        Gz[idx, coords[i, j]] = 1 
        Gz[idx, coords[i+1, j]] = 1 
        idx += 1

    assert idx == m

    mx = lj-1
    Hx = zeros2(mx, n)
    for idx in range(mx):
      for i in range(li):
        Hx[idx, coords[i, idx]] = 1
        Hx[idx, coords[i, idx+1]] = 1

    mz = li-1
    Hz = zeros2(mz, n)
    for idx in range(mz):
      for j in range(lj):
        Hz[idx, coords[idx, j]] = 1
        Hz[idx, coords[idx+1, j]] = 1

    assert dot2(Hx, Hz.transpose()).sum() == 0

    return Gx, Gz, Hx, Hz


def build_compass3(li, lj=None, lk=None):

    if lj is None:
        lj = li

    if lk is None:
        lk = li

    n = li*lj*lk

    keys = [(i, j, k) for i in range(li) for j in range(lj) for k in range(lk)]
    coords = {}  
    for i, j, k in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            for dk in range(-lk, lk+1):
              coords[i+di, j+dj, k+dk] = keys.index(((i+di)%li, (j+dj)%lj, (k+dk)%lk))

    m = 2*n 
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)

    idx = 0 
    for i in range(li):
      for j in range(lj):
       for k in range(lk):
        Gx[idx, coords[i, j, k]] = 1 
        Gx[idx, coords[i+1, j, k]] = 1 

        Gz[idx, coords[i, j, k]] = 1 
        Gz[idx, coords[i, j+1, k]] = 1 
        idx += 1

        Gx[idx, coords[i, j, k]] = 1 
        Gx[idx, coords[i, j+1, k]] = 1 

        Gz[idx, coords[i, j, k]] = 1 
        Gz[idx, coords[i, j, k+1]] = 1 
        idx += 1

    assert idx == m

#    mx = lj-1
#    Hx = zeros2(mx, n)
#    for idx in range(mx):
#      for i in range(li):
#        Hx[idx, coords[i, idx]] = 1
#        Hx[idx, coords[i, idx+1]] = 1
#
#    mz = li-1
#    Hz = zeros2(mz, n)
#    for idx in range(mz):
#      for j in range(lj):
#        Hz[idx, coords[idx, j]] = 1
#        Hz[idx, coords[idx+1, j]] = 1
#
#    assert dot2(Hx, Hz.transpose()).sum() == 0

    Hx = Hz = None

    return Gx, Gz, Hx, Hz


def build_random(n):

    weight = argv.get("weight", 3)
    coweight = argv.get("coweight")

    p = argv.get("p", 0.3)
    m = argv.get("m", n)
    mx = argv.get("mx", m)
    mz = argv.get("mz", m)

    if coweight is not None:
        Gx = rand2(n, mx, weight=coweight).transpose()
        Gz = rand2(n, mz, weight=coweight).transpose()

    else:
        Gx = rand2(mx, n, p=p, weight=weight)
        Gz = rand2(mz, n, p=p, weight=weight)

    Hx = Hz = None

    Gx = Gx[[i for i in range(m) if Gx[i].sum()], :]
    Gz = Gz[[i for i in range(m) if Gz[i].sum()], :]

    li = argv.get("li", True)

    if li:
        Gx = linear_independent(Gx)
        Gz = linear_independent(Gz)

    return Gx, Gz, Hx, Hz


def build_random_selfdual(n):

    weight = argv.get("weight", 3)
    m = argv.get("m", n)
    h = argv.get("h", 0)

    while 1:
        Gx = rand2(m, n, weight=weight)
        Gz = Gx.copy()
    
        Hx = Hz = None
    
        Gx = Gx[[i for i in range(m) if Gx[i].sum()], :]
        Gx = linear_independent(Gx)
    
        if len(Gx)<m:
            write("m")
            continue

        Gz = Gx.copy()

        Hx = find_stabilizers(Gz, Gx)
        Hz = find_stabilizers(Gx, Gz)

        if len(Hx)==h and len(Hz)==h:
            break

        write("H(%d,%d)"%(len(Hx), len(Hz)))

    print
    return Gx, Gz, Hx, Hz


def build_random_nostabs(n):

    m = argv.get("m", n)
    mx = argv.get("mx", m)
    mz = argv.get("mz", m)
    h = argv.get("h", 0)
    hx = argv.get("hx", h)
    hz = argv.get("hz", h)
    while 1:

        Gx, Gz, Hx, Hz = build_random(n)

        if len(Gx)<mx or len(Gz)<mz:
            write("m")
            continue

        Hx = find_stabilizers(Gz, Gx)
        Hz = find_stabilizers(Gx, Gz)

        if len(Hx)==hx and len(Hz)==hz:
            break

        write("H(%d,%d)"%(len(Hx), len(Hz)))

    print

    return Gx, Gz, Hx, Hz


def build_pauli(n):

    m = 2**n
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)
    for i, idx in enumerate(genidx((2,)*n)):
        for j in idx:
            Gx[i,j] = 1
            Gz[i,j] = 1

    Hx = zeros2(0, n)
    Hz = zeros2(0, n)

    return Gx, Gz, Hx, Hz


def build_ising(n):

    assert n>=2

    mx = mz = n
    if n==2:
        mz = 1
    Gx = zeros2(mx, n)
    Gz = zeros2(mz, n)
    for i in range(mx):
        Gx[i, i] = 1 # transverse field

    for i in range(mz):
        Gz[i, i] = 1
        Gz[i, (i+1)%n] = 1

    Hx = zeros2(1, n)
    Hz = zeros2(0, n)

    Hx[:] = 1

    return Gx, Gz, Hx, Hz


def build_isingdual(n):

    Gx, Gz, Hx, Hz = build_ising(n)

    return Gz, Gx, Hz, Hx



#def build_ising2():
#
#    l = argv.get("l", 6)
#    assert l%2 == 0
#
#    li = lj = l
#    n = l**2
#
#    keys = [(i, j) for i in range(li) for j in range(lj)]
#    coords = {}  
#    for i, j in keys:
#        for di in range(-li, li+1):
#          for dj in range(-lj, lj+1):
#            coords[i+di, j+dj] = keys.index(((i+di)%li, (j+dj)%lj))
#
#    assert n%4==0
#    m = n/4
#
#    Gx = zeros2(m, n)
#    Gz = zeros2(m, n)
#
#    idx = 0
#    for i1 in range(l//2):
#      for j1 in range(l//2):
#        i = 2*i1
#        j = 2*j1
#        Gx[idx, coords[i, j]] = 1
#        Gx[idx, coords[i+1, j]] = 1
#        Gx[idx, coords[i, j+1]] = 1
#        Gx[idx, coords[i+1, j+1]] = 1
#
#        Gz[idx, coords[i+1, j+1]] = 1
#        Gz[idx, coords[i+2, j+1]] = 1
#        Gz[idx, coords[i+1, j+2]] = 1
#        Gz[idx, coords[i+2, j+2]] = 1
#
#        idx += 1
#
#    return Gx, Gz, None, None


def build_hex(li, lj=None):
    if lj is None:
        lj = li
    n = li*lj

    keys = [(i, j) for i in range(li) for j in range(lj)]
    coords = {}  
    for i, j in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            coords[i+di, j+dj] = keys.index(((i+di)%li, (j+dj)%lj))

    Gx = []
    if argv.open:
        idxs = range(li-1)
        jdxs = range(lj-1)
    else:
        idxs = range(li)
        jdxs = range(lj)

    for i in idxs:
      for j in jdxs:
        g = zeros2(n)
        g[coords[i,   j]] = 1 
        g[coords[i,   j+1]] = 1 
        g[coords[i+1, j+1]] = 1 
        Gx.append(g)

        g = zeros2(n)
        g[coords[i,   j]] = 1 
        g[coords[i+1, j]] = 1 
        g[coords[i+1, j+1]] = 1 
        Gx.append(g)
    Gx = array2(Gx)

    Gz = Gx.copy()

    return Gx, Gz, None, None


def build_hex2(li, lj=None):
    if lj is None:
        lj = li
    n = li*lj

    keys = [(i, j) for i in range(li) for j in range(lj)]
    coords = {}  
    for i, j in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            coords[i+di, j+dj] = keys.index(((i+di)%li, (j+dj)%lj))

    Gx = []
    Gz = []
    if argv.open:
        idxs = range(li-1)
        jdxs = range(lj-1)
    else:
        idxs = range(li)
        jdxs = range(lj)

    for i in idxs:
      for j in jdxs:
        g = zeros2(n)
        g[coords[i,   j]] = 1 
        g[coords[i,   j+1]] = 1 
        g[coords[i+1, j+1]] = 1 
        Gx.append(g)

        g = zeros2(n)
        g[coords[i,   j]] = 1 
        g[coords[i+1, j]] = 1 
        g[coords[i+1, j+1]] = 1 
        Gz.append(g)
    Gx = array2(Gx)
    Gz = array2(Gz)

    return Gx, Gz, None, None


def build_xy(n):

    m = n
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)
    for i in range(m):
        Gx[i, i] = 1
        Gx[i, (i+1)%n] = 1
        
        Gz[i, i] = 1
        Gz[i, (i+1)%n] = 1

    if n%2 == 0:
        Hx = zeros2(1, n)
        Hz = zeros2(1, n)
    
        Hx[:] = 1
        Hz[:] = 1

    else:
        Hx = Hz = None


    return Gx, Gz, Hx, Hz


def build_xy2(li, lj=None):
    if lj is None:
        lj = li
    n = li*lj

    keys = [(i, j) for i in range(li) for j in range(lj)]
    coords = {}  
    for i, j in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            coords[i+di, j+dj] = keys.index(((i+di)%li, (j+dj)%lj))

    Gx = []
    if argv.open:
        idxs = range(li-1)
        jdxs = range(lj-1)
    else:
        idxs = range(li)
        jdxs = range(lj)

    for i in idxs:
      for j in jdxs:
        g = zeros2(n)
        g[coords[i,   j]] = 1 
        g[coords[i,   j+1]] = 1 
        g[coords[i+1, j]] = 1 
        g[coords[i+1, j+1]] = 1 
        Gx.append(g)

    Gx = array2(Gx)

    Gz = Gx.copy()

    return Gx, Gz, None, None


def build_xy21(li, lj=None):
    if lj is None:
        lj = li
    n = li*lj

    keys = [(i, j) for i in range(li) for j in range(lj)]
    coords = {}  
    for i, j in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            coords[i+di, j+dj] = keys.index(((i+di)%li, (j+dj)%lj))

    Gx = []
    if argv.open:
        idxs = range(li-1)
        jdxs = range(lj-1)
    else:
        idxs = range(li)
        jdxs = range(lj)

    for i in idxs:
      for j in jdxs:
        g = zeros2(n)
        g[coords[i,   j]] = 1 
        g[coords[i,   j+1]] = 1 
        Gx.append(g)

        g = zeros2(n)
        g[coords[i, j]] = 1 
        g[coords[i+1, j]] = 1 
        Gx.append(g)

    Gx = array2(Gx)

    Gz = Gx.copy()

    return Gx, Gz, None, None


def build_xy3(li, lj=None, lk=None):
    if lj is None:
        lj = li
    if lk is None:
        lk = li
    n = li*lj*lk

    keys = [(i, j, k) for i in range(li) for j in range(lj) for k in range(lk)]
    coords = {}  
    for i, j, k in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            for dk in range(-lk, lk+1):
              coords[i+di, j+dj, k+dk] = keys.index(((i+di)%li, (j+dj)%lj, (k+dk)%lk))

    Gx = []
    if argv.open:
        idxs = range(li-1)
        jdxs = range(lj-1)
        kdxs = range(lk-1)
    else:
        idxs = range(li)
        jdxs = range(lj)
        kdxs = range(lk)

    for i in idxs:
     for j in jdxs:
      for k in kdxs:
        g = zeros2(n)
        g[coords[i,   j,   k]] = 1 
        g[coords[i,   j+1, k]] = 1 
        g[coords[i+1, j,   k]] = 1 
        g[coords[i+1, j+1, k]] = 1 
        g[coords[i,   j,   k+1]] = 1 
        g[coords[i,   j+1, k+1]] = 1 
        g[coords[i+1, j,   k+1]] = 1 
        g[coords[i+1, j+1, k+1]] = 1 
        Gx.append(g)

    Gx = array2(Gx)

    Gz = Gx.copy()

    return Gx, Gz, None, None


def build_xy32(li, lj=None, lk=None):

    # TOO BIG ...

    if lj is None:
        lj = li
    if lk is None:
        lk = li
    n = li*lj*lk

    keys = [(i, j, k) for i in range(li) for j in range(lj) for k in range(lk)]
    coords = {}  
    for i, j, k in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            for dk in range(-lk, lk+1):
              coords[i+di, j+dj, k+dk] = keys.index(((i+di)%li, (j+dj)%lj, (k+dk)%lk))

    Gx = []
    if argv.open:
        idxs = range(li-1)
        jdxs = range(lj-1)
        kdxs = range(lk-1)
    else:
        idxs = range(li)
        jdxs = range(lj)
        kdxs = range(lk)

    for i in idxs:
     for j in jdxs:
      for k in kdxs:
        g = zeros2(n)
        g[coords[i,   j,   k]] = 1 
        g[coords[i,   j+1, k]] = 1 
        g[coords[i+1, j,   k]] = 1 
        g[coords[i+1, j+1, k]] = 1 
        Gx.append(g)

        g = zeros2(n)
        g[coords[i,   j,   k+1]] = 1 
        g[coords[i,   j+1, k+1]] = 1 
        g[coords[i+1, j,   k+1]] = 1 
        g[coords[i+1, j+1, k+1]] = 1 
        Gx.append(g)

        g = zeros2(n)
        g[coords[i,   j,   k]] = 1 
        g[coords[i,   j,   k+1]] = 1 
        g[coords[i+1, j,   k]] = 1 
        g[coords[i+1, j,   k+1]] = 1 
        Gx.append(g)

        g = zeros2(n)
        g[coords[i,   j+1, k]] = 1 
        g[coords[i,   j+1, k+1]] = 1 
        g[coords[i+1, j+1, k]] = 1 
        g[coords[i+1, j+1, k+1]] = 1 
        Gx.append(g)

#        g = zeros2(n)
#        g[coords[i,   j,   k]] = 1 
#        g[coords[i,   j,   k+1]] = 1 
#        g[coords[i,   j+1, k]] = 1 
#        g[coords[i,   j+1, k+1]] = 1 
#        Gx.append(g)
#
#        g = zeros2(n)
#        g[coords[i+1, j,   k]] = 1 
#        g[coords[i+1, j,   k+1]] = 1 
#        g[coords[i+1, j+1, k]] = 1 
#        g[coords[i+1, j+1, k+1]] = 1 
#        Gx.append(g)

    Gx = array2(Gx)

    Gz = Gx.copy()

    return Gx, Gz, None, None



def mkop(n, ops):
    A = zeros2(len(ops), n)
    for i, op in enumerate(ops):
        for j in op:
            A[i, j] = 1
    return A


def build_gcolor2():

    n = 39
    m = 10

    delta = 19

    top = n-1 # top qubit

    # bottom faces: points must be adjacent in each face
    bfaces = [
        [0, 1, 2, 3],
        [1, 4, 5, 6, 7, 2],
        [3, 2, 7, 8],
        [4, 9, 10, 5],
        [8, 7, 6, 11, 12, 13],
        [9, 14, 15, 10],
        [5, 10, 15, 16, 11, 6],
        [11, 16, 17, 12],
        [13, 12, 17, 18]]

    faces = list(bfaces) + [[i+delta for i in face] for face in bfaces]

    def check_faces():
        items = [list(face) for face in faces]
        for face in items:
            assert len(face)%2 == 0, face
            face.sort()
        assert len(set([tuple(face) for face in items]))==len(items)
    check_faces()

    # bottom edges
    bedges = []
    for face in bfaces:
        f = len(face)
        for i in range(f):
            bedges.append([face[i], face[(i+1)%f]])

    # edges are not yet unique..
    for edge in bedges:
        edge.sort()
    bedges = list(set([tuple(e) for e in bedges]))

    # extrude bottom edges to make a face
    for edge in bedges:
        edge = list(edge)
        a, b = edge
        face = edge + [a+delta, b+delta]
        faces.append(face)
    check_faces()

    stabs = []
    for face in bfaces:
        stabs.append(face + [i+delta for i in face])

    # top faces
    for face in [
        [0, 1, 4, 9, 14],
        [0, 3, 8, 13, 18],
        [14, 15, 16, 17, 18]]:
        face = [i+delta for i in face] + [top]
        faces.append(face)
    check_faces()

    stabs.append([i+delta for i in range(19)] + [top])

    g = len(faces)
    #print "faces:", g

    for stab in stabs:
        assert len(stab)%2 == 0, stab

    #faces.sort()
    #for face in faces:
    #    print face

    Gx = mkop(n, faces)
    Gz = Gx.copy()

    rows = [shortstr(g) for g in Gx]
    #rows.sort()
    #for i, row in enumerate(rows):
    #    print row, faces[i]
    assert len(set(rows))==len(rows)

    Hz = mkop(n, stabs)
    Hx = Hz.copy()

    # bottom face
    Lx = mkop(n, [range(19)])
    Lz = Lx.copy()

    check_commute(Hx, Hz)
    check_commute(Hx, Gz)
    check_commute(Hz, Gx)
    check_commute(Gx, Lz)
    check_commute(Gz, Lx)
    check_commute(Hx, Lz)
    check_commute(Hz, Lx)

    #code = CSSCode(Hx=Hx, Gx=Gx, Hz=Hz, Gz=Gz, build=False)

    Lx = find_logops(Gz, Hx)

    #print "Lx:", shortstr(Lx)

    return Gx, Gz, Hx


def build_projective(n, dim=2):

    from bruhat import incidence as geometry
    g = geometry.projective(n, dim)

    P = g.types[0]
    L = g.types[1]
    if dim==3:
        L = g.types[2]

    points = g.tplookup[P]
    lines = g.tplookup[L]

    #lines = lines[:-4] # throw one out
    #points = points[:-1] # throw one out

    n = len(points)
    m = len(lines)
    Gx = zeros2(m, n)
    for i, line in enumerate(lines):
        for j, point in enumerate(points):
            if (line, point) in g.incidence:
                Gx[i, j] = 1

    #print shortstr(Gx)

    Gz = Gx.copy()

    Hx = None
    Hz = None

    return Gx, Gz, Hx, Hz




def build(name=""):

    if name:
        setattr(argv, name, True) # hack this

    _seed = argv.get("seed")
    if _seed is not None:
        numpy.random.seed(_seed)
        seed(_seed)

    size = argv.get("size", 1)
    l = argv.get('l', 4)
    li = argv.get('li', l)
    lj = argv.get('lj', l)
    lk = argv.get('lk', l)

    if argv.gcolor2 or (argv.gcolor and size==1.5):
        Gx, Gz, Hx = build_gcolor2()
        Hz = Hx.copy()

    elif argv.gcolor:
        Gx, Gz, Hx = build_gcolor(size)
        Hz = Hx.copy()

    elif argv.compass:
        Gx, Gz, Hx, Hz = build_compass(li, lj)

    elif argv.compass3:
        Gx, Gz, Hx, Hz = build_compass3(li, lj, lk)

    elif argv.hex:
        Gx, Gz, Hx, Hz = build_hex(li, lj)

    elif argv.hex2:
        Gx, Gz, Hx, Hz = build_hex2(li, lj)

    elif argv.xy:
        Gx, Gz, Hx, Hz = build_xy(l)

    elif argv.xy2:
        Gx, Gz, Hx, Hz = build_xy2(li, lj)

    elif argv.xy21:
        Gx, Gz, Hx, Hz = build_xy21(li, lj)

    elif argv.xy3:
        Gx, Gz, Hx, Hz = build_xy3(li, lj, lk)

    elif argv.xy32:
        Gx, Gz, Hx, Hz = build_xy32(li, lj, lk)

    elif argv.ising:
        Gx, Gz, Hx, Hz = build_ising(l)

    elif argv.isingdual:
        Gx, Gz, Hx, Hz = build_isingdual(l)

    elif argv.random:
        Gx, Gz, Hx, Hz = build_random(l)

    elif argv.random_nostabs:
        Gx, Gz, Hx, Hz = build_random_nostabs(l)

    elif argv.random_selfdual:
        Gx, Gz, Hx, Hz = build_random_selfdual(l)

    elif argv.pauli:
        Gx, Gz, Hx, Hz = build_pauli(l)

    elif argv.projective:
        n = argv.get('n', 3)
        dim = argv.get('dim', 2)
        Gx, Gz, Hx, Hz = build_projective(n, dim)

    elif argv.test:
        Gx, Gz, Hx, Hz = build_test()

    else:

        name = argv.next()
        try:
            fn = eval("build_%s"%name)
        except NameError:
            print "no model found"
            raise

        Gx, Gz, Hx, Hz = fn()

    if Hx is None:
        Hx = find_stabilizers(Gz, Gx)
    if Hz is None:
        Hz = find_stabilizers(Gx, Gz)

    if argv.flip:
        Gz, Gx = Gx, Gz
        Hz, Hx = Hx, Hz

    if argv.show:
        print "Gx Gz:"
        print shortstrx(Gx, Gz)
        if len(Hx):
            print "Hx Hz:"
            print shortstrx(Hx, Hz)

    return Gx, Gz, Hx, Hz


def build_reduced():

    Gx, Gz, Hx, Hz = build()

    Px = get_reductor(Hx) # projector onto complement of rowspan of Hx
    Pz = get_reductor(Hz)
    Rz = [dot2(Pz, g) for g in Gz]
    Rz = array2(Rz)
    Rz = row_reduce(Rz, truncate=True)

    Rx = [dot2(Px, g) for g in Gx]
    Rx = array2(Rx)
    Rx = row_reduce(Rx, truncate=True)

    return Rx, Rz



class Model(object):
    def __init__(self, attrs):
        self.__dict__.update(attrs)
        self.Qx = self.Rz.transpose() # backwards compat
        self.cache = {}

    def __str__(self):
        return "Model(n=%d, Lx/z: %d, Gx: %d, Gz: %d, Hx: %d, Hz: %d, Rx/z: %d)" % (
            self.n,
            len(self.Lx), len(self.Gx), len(self.Gz), 
            len(self.Hx), len(self.Hz), len(self.Rx))

    attrs = "Gz Gx Rz Rx Hz Hx Tz Tx Pz Px Lz Lx".split()
    def get_dual(self):
        Gz, Gx = self.Gx, self.Gz        
        Rz, Rx = self.Rx, self.Rz        
        Hz, Hx = self.Hx, self.Hz        
        Tz, Tx = self.Tx, self.Tz        
        Pz, Px = self.Px, self.Pz        
        Lz, Lx = self.Lx, self.Lz        
        #Qz, Qx = self.Qx, self.Qz        
        n = self.n
        return Model(locals())

    def build_ham(self, excite=None, weights=None, Jx=1., Jz=1.):
        Gx, Gz = self.Gx, self.Gz        
        Rx, Rz = self.Rx, self.Rz        
        Hx, Hz = self.Hx, self.Hz        
        Tx, Tz = self.Tx, self.Tz        
        gz = len(Gz)
        r = len(Rx)
        n = self.n

        if type(excite) is int:
            _excite = [0]*len(Tx)
            _excite[excite] = 1
            excite = tuple(_excite)

        if excite is not None:
            assert len(excite)==len(Tx)
    
            t = zeros2(n)
            for i, ex in enumerate(excite):
                if ex:
                    t = (t + Tx[i])%2
            #print "t:", shortstr(t)
            Gzt = dot2(Gz, t)

        else:
            Gzt = 0

        if weights is None:
            weights = [1.]*len(Gx)
        assert len(weights) == len(Gx), len(weights)

        H = numpy.zeros((2**r, 2**r))
        for i, v in enumerate(genidx((2,)*r)):
            v = array2(v)
            syndrome = (dot2(Gz, Rx.transpose(), v) + Gzt)%2
            value = gz - 2*syndrome.sum()
            #print shortstr(dot2(Rx.transpose(), v)), value
            H[i, i] = Jz*value
            #U.append(value)

        Pxt = self.Px.transpose()
        Qx = Rz.transpose()
        #print dot2(Rx, Qx)
        PxtQx = dot2(Pxt, Qx)
        for i, v in enumerate(genidx((2,)*r)):
            v = array2(v)
            #print shortstr(v),
            #for g in Gx:
            for j, g in enumerate(Gx):
                u = (v + dot2(g, PxtQx))%2
                k = eval('0b'+shortstr(u, zero='0'))
                H[i, k] += Jx*weights[j]
                #A[i, k] = A.get((i, k), 0) + 1

        return H

    def sparse_ham_eigs(self, excite=None, weights=None, Jx=1., Jz=1.):

        key = str((excite, weights, Jx, Jz))
        if key in self.cache:
            return self.cache[key]

        Gx, Gz = self.Gx, self.Gz        
        Rx, Rz = self.Rx, self.Rz        
        Hx, Hz = self.Hx, self.Hz        
        Tx, Tz = self.Tx, self.Tz        
        Px, Pz = self.Px, self.Pz

        gz = len(Gz)
        r = len(Rx)
        n = self.n

        if type(excite) is int:
            _excite = [0]*len(Tx)
            _excite[excite] = 1
            excite = tuple(_excite)

        if excite is not None:
            assert len(excite)==len(Tx)
    
            t = zeros2(n)
            for i, ex in enumerate(excite):
                if ex:
                    t = (t + Tx[i])%2
            #print "t:", shortstr(t)
            Gzt = dot2(Gz, t)

        else:
            Gzt = 0

        verts = []
        lookup = {}
        for i, v in enumerate(span(Rx)): # XXX does not scale well
            #if v0 is not None:
            #    v = (v+v0)%2
            #    v = dot2(Px, v)
            lookup[v.tostring()] = i
            verts.append(v)
        print "span:", len(verts)
        assert len(lookup) == len(verts)
    
        mz = len(Gz)
        n = len(verts)

        print "building H",
        H = {} # adjacency
        U = [] # potential

        #if offset is None:
        offset = mz + 1 # make H positive definite

        for i, v in enumerate(verts):
            if i%1000==0:
                write('.')
            #count = dot2(Gz, v).sum()
            syndrome = (dot2(Gz, v) + Gzt) % 2
            count = syndrome.sum()
            #syndrome = (dot2(Gz, Rx.transpose(), v) + Gzt)%2
            #H[i, i] = mz - 2*count
            U.append(offset + mz - 2*count)
            for g in Gx:
                v1 = (g+v)%2
                v1 = dot2(Px, v1)
                j = lookup[v1.tostring()]
                H[i, j] = H.get((i, j), 0) + 1

        print "\nnnz:", len(H)
        for i in range(len(U)):
            H[i, i] = H.get((i, i), 0) + U[i] 
        N = len(U)
        del U
        #H1 = sparse.lil_matrix(N, N)
        keys = H.keys()
        keys.sort()
        data = [] 
        rows = [] 
        cols = [] 
        for idx in keys:
            #H1[idx] = H[idx]
            data.append(H[idx])
            rows.append(idx[0])
            cols.append(idx[1])
        del H
        H1 = sparse.coo_matrix((data, (rows, cols)), (N, N))
        H1 = sparse.csr_matrix(H1, dtype=numpy.float64)
    
        #print "do_lanczos: eigsh"
        vals, vecs = sparse.linalg.eigsh(H1, k=min(N-5, 40), which="LM")
    
        vals -= offset
        self.cache[key] = vals
        return vals

    def do_slepc(self, excite=None, weights=None, Jx=1., Jz=1.):
        key = str((excite, weights, Jx, Jz))
        if key in self.cache:
            #print "CACHE HIT"
            return self.cache[key]

        from slepc import slepc
        vals = slepc(excite=excite, **self.__dict__)
        self.cache[key] = vals
        return vals

    def do_lp(self):
        # so far a failed experiment to apply LP...
        Gx, Gz = self.Gx, self.Gz        
        Rx, Rz = self.Rx, self.Rz        
        Hx, Hz = self.Hx, self.Hz        
        Tx, Tz = self.Tx, self.Tz        
        Px, Pz = self.Px, self.Pz

        gz = len(Gz)
        r = len(Rx)
        n = self.n

        assert len(Hz)

        import pulp
        
        prob = pulp.LpProblem("test1", pulp.LpMinimize) 
        
        # Variables 
        #x = pulp.LpVariable("x", 0, 4) 
        #y = pulp.LpVariable("y", -1, 1) 
        #z = pulp.LpVariable("z", 0) 

        points = list(enum2(r))

        lookup = {}
        #ps = []
        for u in points:
            #name = "u"+shortstr(u)
            #var = pulp.LpVariable(name, 0., 1.)
            #lookup[u.tostring()] = var
            #ps.append(var)
            for v in points:
                name1 = "u%sv%s" % (shortstr(u), shortstr(v))
                lookup[u.tostring(), v.tostring()] = pulp.LpVariable(name1, 0., 1.)
        
        # Objective 
        #prob += x + 4*y + 9*z 
        
        ps = [lookup[u.tostring(), u.tostring()] for u in points]
        prob += sum(ps)==1.

        if 0:
            for t in enum2(len(Tx)):
                txt = dot2(Tx.transpose(), t)
                items = []
                for u in points:
                    Rxtu = dot2(Rx.transpose(), u)
                    coeff = dot2(Gz, Rxtu + txt).sum() - dot2(Gz, Rxtu).sum()
                    items.append(coeff * lookup[shortstr(u)])
                prob += sum(items) > 0

        ham = []
        #pairs = []
        for u in points:
            w = dot2(Gz, Rx.transpose(), u).sum()
            ham.append((len(Gz) - 2*w) * lookup[u.tostring(), u.tostring()])
            for gx in Gx:
                v = (u + dot2(gx, Rz.transpose()))%2
                key = u.tostring(), v.tostring()
                ham.append(lookup[key])
                #pairs.append(key)
                
            #print w, shortstr(v)
        print "ham", len(ham)

        #pairs = set(pairs) # uniq
        #for u, v in pairs:
        spoints = [u.tostring() for u in points]
        for u in spoints:
          for v in spoints:
            # 1/2(x**2 + y**2) >= xy
            prob += 0.5*(lookup[u, u] + lookup[v, v]) >= lookup[u, v]
            for w in spoints:
                prob += 0.5*(lookup[u,u]+lookup[v,v]+lookup[w,w])>=\
                    lookup[u,v]+lookup[u,w]-lookup[v,w]

            prob += (lookup[u, v]==lookup[v, u])

        # Objective
        prob += -sum(ham)
        
        print "solving..."
        pulp.GLPK().solve(prob) 
        
        # Solution 
        for v in prob.variables(): 
            if v.varValue > 0.:
                print v.name, "=", v.varValue 
        
        print "objective=", pulp.value(prob.objective)  
    

    def show_stabx(self, sx):
        gxs = []
        for gx in self.Gx:
            if eq2(gx*sx, gx):
                gxs.append(gx)
        Gx = array2(gxs)
        #print "Gx:", Gx.shape
        #print shortstr(Gx)
        print "sx.sum() =", sx.sum(),
        Gxt = Gx.transpose()
        K = find_kernel(Gxt)
        #print "kernel:", K
        K = array2(K)
        #print "kernel:", len(K)
        #print shortstr(K)
        #print
        best = None
        ubest = None
        u = solve(Gxt, sx)
        for w in enum2(len(K)):
            u2 = (u+dot2(w, K))%2
            if best is None or u2.sum() < best:
                best = u2.sum()
                ubest = u2
            print "u.sum() =", u2.sum(),
        print

#        print shortstr(sx)
#        print "u.sum() =", ubest.sum()
        #print shortstr(sx)
        #print #"-"*len(sx)
#        for i in range(len(ubest)):
#            if ubest[i]:
#                print shortstr(Gx[i])
#        print



def check_sy(Lx, Hx, Tx, Rx, Lz, Hz, Tz, Rz, **kw):

    check_conjugate(Lx, Lz)
    check_commute  (Lx, Hz)
    check_commute  (Lx, Tz)
    check_commute  (Lx, Rz)

    check_commute  (Hx, Lz)
    check_conjugate(Hx, Tz)
    check_commute  (Hx, Hz)
    check_commute  (Hx, Rz)

    check_commute  (Tx, Lz)
    check_commute  (Tx, Tz)
    check_conjugate(Tx, Hz)
    check_commute  (Tx, Rz)

    check_commute  (Rx, Lz)
    check_commute  (Rx, Hz)
    check_commute  (Rx, Tz)
    check_conjugate(Rx, Rz)



def build_model(Gx=None, Gz=None, Hx=None, Hz=None):

    if Gx is None:
        Gx, Gz, Hx, Hz = build()

    n = Gx.shape[1]

    if Hx is None:
        Hx = find_stabilizers(Gz, Gx)
    if Hz is None:
        Hz = find_stabilizers(Gx, Gz)

    check_commute(Hx, Hz)
    check_commute(Gx, Hz)
    check_commute(Hx, Gz)

    #Px = get_reductor(concatenate((Lx, Hx)))
    #Pz = get_reductor(concatenate((Lz, Hz)))
    Px = get_reductor(Hx)
    Pz = get_reductor(Hz)

    # Lz = find_logops( Hx            , Hz            )
    #      find_logops( ............. , ............. )
    #                 ( commutes with , orthogonal to )
    #                 ( ............. , ............. )

    Lz = find_logops(Gx, Hz)
    assert Lz.shape[1] == n

    if 0:
        PGz = get_reductor(Gz)
        Lz = dot2(Lz, PGz.transpose())
        Lz = row_reduce(Lz)
    
        print shortstrx(Lz, Gz, Hz)

    if len(Lz):
        #print Lz.shape, Hz.shape
        assert len(row_reduce(concatenate((Lz, Hz))))==len(Lz)+len(Hz)
        assert len(row_reduce(concatenate((Lz, Gz))))==len(Lz)+len(row_reduce(Gz))

    # Tz = find_errors( Hx            , Lx            )
    #      find_errors( ............. , ............. )
    #                 ( conjugate to  , commutes with )
    #                 ( ............. , ............. )

    Lx = find_errors(Lz, Gz) # invert Lz, commuting with Gz

    check_commute  (Lx, Gz)
    check_commute  (Lx, Hz)
    check_conjugate(Lx, Lz)
    check_commute  (Lz, Gx)
    check_commute  (Lz, Hx)


    # Lx | Lz
    # Hx | ?
    # ?  | Hz
    # ?  | ?
    #Rz = find_logops(concatenate((Lx, Hx)), Hz)
    Rz = dot2(Gz, Pz.transpose())
    Rz = row_reduce(Rz)

    check_commute  (Rz, Lx)
    check_commute  (Rz, Hx)

    Rx = dot2(Gx, Px.transpose())
    Rx = row_reduce(Rx)

    check_commute  (Rx, Lz)
    check_commute  (Rx, Hz)

    # Lx | Lz
    # Hx | ?
    # ?  | Hz
    # Rx'| Rz'

    Tz = find_errors(Hx, concatenate((Lx, Rx)))
    Tx = find_errors(Hz, concatenate((Lz, Rz, Tz)))

    assert len((concatenate((Lx, Hx, Tx, Rx)))) == n
    assert len((concatenate((Lz, Hz, Tz, Rz)))) == n
    assert len(row_reduce(concatenate((Lx, Hx, Tx, Rx)))) == n
    assert len(row_reduce(concatenate((Lz, Hz, Tz, Rz)))) == n

    check_commute  (Rz, Tx)

    Rx = find_errors(Rz, concatenate((Lz, Hz, Tz)))

    check_conjugate(Rx, Rz)
    check_commute  (Rx, Hz)
    check_commute  (Rx, Tz)
    check_commute  (Rx, Lz)

    Rxt = Rx.transpose()
    Rzt = Rz.transpose()

    Pxt = Px.transpose()
    Pzt = Pz.transpose()

    check_sy(Lx, Hx, Tx, Rx, Lz, Hz, Tz, Rz)

    assert eq2(dot2(Gz, Rxt), dot2(Gz, Pzt, Rxt))
    assert eq2(dot2(Gx, Rzt), dot2(Gx, Pxt, Rzt))

#    print shortstrx(dot2(Rx, Pz), Rx)

    assert eq2(dot2(Rx, Pz), Rx) 
    assert eq2(dot2(Rz, Px), Rz) 

    assert len(find_kernel(dot2(Gz, Rx.transpose())))==0

    model = Model(locals())

    if argv.dual:
        print("get_dual")
        model = model.get_dual()
        argv.dual = False # HACK !!

    return model



if __name__ == "__main__":

    Gx, Gz, Hx, Hz = build()

    model = build_model(Gx, Gz, Hx, Hz)


    if argv.extend:
        k = len(model.Lx)
        n = model.n + k
        mx = len(model.Hx)
        mz = len(model.Hz)
            
        Hx = zeros2(mx+k, n+k)
        Hz = zeros2(mz+k, n+k)

        Hx[:mx, :n] = model.Hx
        Hz[:mz, :n] = model.Hz
        Hx[mx:, :n] = model.Lx
        Hz[mz:, :n] = model.Lz

        for i in range(k):
            Hx[mx+i, n+i] = 1
            Hz[mz+i, n+i] = 1

        model = build_model() # um....
            
    print model

    if argv.show:
        print "Hx/Hz:"
        print shortstrx(model.Hx, model.Hz)
        print
        print "Gx/Gz:"
        print shortstrx(Gx, Gz)
        print
        print "Lx/Lz:"
        print shortstrx(model.Lx, model.Lz)

    if len(model.Lx) and argv.distance:
        w = min([v.sum() for v in span(model.Lx) if v.sum()])
        print "distance:", w
        
    if argv.do_lp:
        model.do_lp()

    if argv.do_slepc:
        model.do_slepc()

    if argv.solve:
        vals = model.sparse_ham_eigs()
        print vals

#    for g in Gx:
#        print g.sum(),
#    print

    Rx = model.Rx
    HR = numpy.concatenate((Hx, Rx))
#    print shortstr(Hx)
#    print
#    print shortstr(SR)

#    U = solve(Gx.transpose(), HR.transpose())
#    print U.shape
#    h = len(Hx)
#    #GH = dot2(Gx.transpose(), U)
#    U = U[:, :h]

    #K = find_kernel(Gx.transpose())
    #print "kernel:", len(K)

    #for g in Hx:
    #    model.show_stabx(g)

#        best = None
#        vbest = None
#        for u in enum2(h):
#            v = dot2(U, u)
#            vsum = v.sum()
#            if vsum == 0:
#                continue
#            if best is None or v.sum() < best:
#                best = v.sum()
#        print "best:", best

        
            
        



