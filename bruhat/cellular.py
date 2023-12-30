#!/usr/bin/env python

"""
Build some 
2-dimensional cellular complexes, mostly by mutation.
"""


from time import time
from random import choice, randint, seed, shuffle

import numpy

from bruhat.argv import argv


class Cell(object):
    def __init__(self, dim=0, children=[]):
        self.dim = dim
        if type(children) in [list, tuple, set]:
            children = {child:1 for child in children}
        self.children = dict(children)

    def __str__(self):
        name = "Vert Edge Face".split()[self.dim]
        return "%s(%d)"%(name, len(self.children))
    __repr__ = __str__

    def __delitem__(self, child):
        del self.children[child]

    def __setitem__(self, child, r):
        assert child not in self.children
        self.children[child] = r

    def __iter__(self):
        return iter(self.children)

    def get(self, child):
        return self.children.get(child)

    def __getitem__(self, child):
        return self.children[child]

    def items(self):
        return self.children.items()

    def __contains__(self, cell):
        return cell in self.children

    def __len__(self):
        return len(self.children)

    def intersect(self, other):
        return [child for child in self if child in other]

    @property
    def src(self):
        assert self.dim == 1, "not an edge"
        assert len(self.children) == 2
        for (cell, r) in self.children.items():
            if r == -1:
                return cell
        assert 0, self.children

    @property
    def tgt(self):
        assert self.dim == 1, "not an edge"
        assert len(self.children) == 2
        for (cell, r) in self.children.items():
            if r == +1:
                return cell
        assert 0, self.children


class Complex(object):
    def __init__(self, grades=None):
        if grades is None:
            grades = [[], [], []]
        self.grades = [list(g) for g in grades]
        self.lookup = {}

    def check(self):
        grades, lookup = self.grades, self.lookup
        for (dim, grade) in enumerate(grades):
            for i,cell in enumerate(grade):
                assert lookup[cell] == (dim, i)
        for (cell, k) in lookup.items():
            dim, i = k
            assert grades[dim][i] is cell
        for e in self.edges:
            assert len(e) == 2
            vals = []
            for v,r in e.items():
                vals.append(r)
            vals.sort()
            assert vals == [-1, 1], e.children

    def is_built(self):
        H0 = self.bdy(0)
        H1 = self.bdy(1)
        #print(H0)
        #print(H1)
        for e in self.edges:
            vals = []
            for f in self.faces:
                if f.get(e):
                    vals.append(f[e])
            vals.sort()
            assert vals == [-1, 1], e
        return True

    @property
    def sig(self):
        return [len(g) for g in self.grades]

    def __str__(self):
        return "Complex(%s)"%(self.sig,)

    @property
    def verts(self):
        return self.grades[0]

    @property
    def edges(self):
        return self.grades[1]

    @property
    def faces(self):
        return self.grades[2]

    @property
    def euler(self):
        r = 0
        sgn = 1
        for grade in self.grades:
            r += sgn * len(grade)
            sgn = -sgn
        return r

    def bdy(self, dim=None):
        if dim is None:
            return [self.bdy(d) for d in range(len(self.grades)-1)]
        tgt, src = self.grades[dim:dim+2]
        lookup = self.lookup
        m, n = len(tgt), len(src)
        H = numpy.zeros((m, n), dtype=int)
        #print("bdy", dim, src)
        for j,cell in enumerate(src):
            #print(j, cell)
            for dell,r in cell.children.items():
                _, i = lookup[dell]
                #i = tgt.index(dell)
                #print('\t', (i, j), "<--", r)
                H[i, j] = r
        return H

    @classmethod
    def frombdy(cls, H0, H1):
        #print("frombdy")
        #print(H0)
        #print(H1)
        HH = numpy.dot(H0, H1)
        assert numpy.alltrue(HH==0)
        cx = cls()
        m, n = H0.shape
        verts = [cx.vertex() for i in range(m)]
        edges = []
        for j in range(n):
            src, tgt = None, None
            for i in range(m):
                r = H0[i,j]
                if r == -1:
                    assert src is None
                    src = verts[i]
                if r == +1:
                    assert tgt is None
                    tgt = verts[i]
            assert src is not None
            assert tgt is not None
            edges.append(cx.edge(src, tgt))
        m, n = H1.shape
        faces = []
        for j in range(n):
            rs = {}
            for i in range(m):
                r = H1[i,j]
                if r != 0:
                    rs[edges[i]] = r
            faces.append(cx.cell(2, rs))
        return cx

    def dual(self):
        H0 = self.bdy(0)
        H1 = self.bdy(1)
        #print("dual")
        #print(H0)
        #print(H1)
        cx = Complex.frombdy(H1.transpose(), H0.transpose())
        return cx

    def is_homology(self):
        assert self.is_built()
        H0, H1 = self.bdy(0), self.bdy(1)
        H = numpy.dot(H0, H1)
        for idx in numpy.ndindex(H.shape):
            if H[idx] != 0:
                return False
        return True

    def cell(self, dim=0, children=[]):
        cell = Cell(dim, children)
        grade = self.grades[dim]
        self.lookup[cell] = (dim, len(grade))
        grade.append(cell)
        self.check()
        return cell

    def remove(self, cell):
        grades, lookup = self.grades, self.lookup
        dim, idx = lookup[cell]
        grade = grades[dim]
        assert grade.pop(idx) == cell
        del lookup[cell]
        for jdx in range(idx, len(grade)):
            cell = grade[jdx]
            lookup[cell] = (dim, jdx)
        self.check()

    def face(self, children):
        cell = self.cell(2, children)
        return cell

    def edge(self, v0, v1):
        cell = self.cell(1, {v0:-1, v1:+1})
        return cell

    def vertex(self):
        cell = self.cell(0)
        return cell

    def cut_edge(self, e):
        assert e.dim == 1
        self.remove(e)
        v0 = e.src
        v2 = e.tgt
        v1 = self.vertex()
        e0 = self.edge(v0, v1)
        e1 = self.edge(v1, v2)
        for face in self.grades[2]:
            r = face.get(e)
            if r is None:
                continue
            del face[e]
            face[e0] = r
            face[e1] = r
        self.check()
        return e0, e1

    def cut_face(self, face):
        assert face.dim == 2
        self.remove(face)
        cv = self.vertex() # center
        es = list(face) # edges
        vs = [e.tgt for e in es] # _vertices
        vs += [e.src for e in es] # _vertices
        vs = set(vs)
        assert len(es) == len(vs)
        edges = {v:self.edge(v, cv) for v in vs}
        for e,r in face.items():
            e0 = edges[e.src]
            e1 = edges[e.tgt]
            self.face({e:r, e1:r, e0:-r})
        self.check()

    def barycenter(self, face):
        assert face.dim == 2
        for e in list(face):
            e0, e1 = self.cut_edge(e)
        self.cut_face(face)

    def bevel(self, v0):
        assert v0.dim == 0
        faces = []
        for face in self.faces:
          for e in face:
            if v0 in e:
                faces.append(face)
                break
        #print(len(faces))
        edges = []
        for e in list(self.edges):
            if v0 not in e:
                continue
            e0, e1 = self.cut_edge(e)
            if v0 in e0:
                edges.append(e0)
            elif v0 in e1:
                edges.append(e1)
            else:
                assert 0
        assert len(edges) == len(faces)
        face = {} # new face
        for f in faces:
            #print(len(f), end=" ")
            src = tgt = None
            for e in edges:
                if e not in f:
                    continue
                for v in e:
                    if v != v0:
                        break
                r = e[v]
                if f[e] == r:
                    assert tgt is None
                    tgt = v
                elif f[e] == -r:
                    assert src is None
                    src = v
                del f[e]
            assert tgt is not None
            assert src is not None
            #print(len(f), end=" ")
            e = self.edge(src, tgt)
            f[e] = 1
            #print(len(f))
            face[e] = -1
        for e in edges:
            self.remove(e)
        self.remove(v0)
        f = self.face(face)
        self.check()
        return f


def make_ball():
    cx = Complex()
    vertex, edge, face = cx.vertex, cx.edge, cx.face
    v0 = vertex()
    v1 = vertex()
    e0 = edge(v0, v1)
    e1 = edge(v1, v0)
    f0 = face({e0, e1})
    f1 = face({e0:-1, e1:-1})
    return cx


def make_torus(rows=2, cols=2):
    cx = Complex()
    vertex, edge, face = cx.vertex, cx.edge, cx.face
    verts = {}
    for i in range(rows):
      for j in range(cols):
        verts[i, j] = vertex()
    hedge = {}
    vedge = {}
    for i in range(rows):
      for j in range(cols):
        hedge[i, j] = edge(verts[i,j], verts[i,(j+1)%cols])
        vedge[i, j] = edge(verts[i,j], verts[(i+1)%rows,j])
    for i in range(rows):
      for j in range(cols):
        face({hedge[i,j]:1, vedge[i,(j+1)%cols]:1, 
            hedge[(i+1)%rows,j]:-1, vedge[i,j]:-1})
    return cx


def make_tetrahedron():
    cx = make_ball()
    cx.cut_edge(cx.edges[0])
    cx.cut_face(cx.faces[0])
    return cx


def make_octahedron():
    cx = make_ball()
    e0, e1 = cx.edges
    cx.cut_edge(e0)
    cx.cut_edge(e1)
    f0, f1 = cx.faces
    cx.cut_face(f0)
    cx.cut_face(f1)
    return cx


def make_cube():
    cx = make_octahedron()
    cx = cx.dual()
    return cx


def test():

    cx = make_ball()
    H0 = cx.bdy(0)
    H1 = cx.bdy(1)
    HH = numpy.dot(H0, H1)
    assert numpy.alltrue(HH==0)
    assert cx.is_homology()
    assert cx.euler == 2

    e0 = cx.edges[0]
    cx.cut_edge(e0)
    assert cx.is_homology()
    assert cx.euler == 2

    f = cx.faces[0]
    cx.cut_face(f)
    assert cx.is_homology()
    assert cx.euler == 2

    f = cx.faces[0]
    cx.barycenter(f)
    assert cx.is_homology()
    assert cx.euler == 2

    cx = make_tetrahedron()
    assert cx.is_homology()
    assert cx.euler == 2
    cx.cut_face(cx.faces[0])
    assert cx.is_homology()
    assert cx.euler == 2

    cx = make_octahedron()
    assert cx.is_homology()
    assert cx.euler == 2

    cx = make_cube()
    assert cx.is_homology()
    assert cx.euler == 2
    cx.cut_face(cx.faces[0])
    assert cx.is_homology()
    assert cx.euler == 2
    
    cx = make_torus()
    assert cx.is_homology()
    assert cx.euler == 0

    d = cx.dual()
    assert d.euler == 0

    cx.cut_face(cx.faces[0])
    assert cx.is_homology()
    assert cx.euler == 0
    
    d = cx.dual()
    assert d.euler == 0

    cx = make_tetrahedron()
    for v in list(cx.verts):
        f = cx.bevel(v)
    assert cx.is_homology()
    assert cx.euler == 2, cx.euler
    assert cx.sig == [12, 18, 8]

    cx = make_cube()
    for v in list(cx.verts):
        f = cx.bevel(v)
    assert cx.is_homology()
    assert cx.euler == 2, cx.euler
    assert cx.sig == [24, 36, 14]


def get_code(cx):
    H0, H1 = cx.bdy()
    H0 %= 2
    H1 %= 2
    A = numpy.dot(H0, H1) # vert -- face
    #print(A)
    vdeg = (cx.bdy(0)%2).sum(1)
    assert vdeg.min() >= 3
    assert vdeg.max() <= 4

    n = len(cx.verts)
    m = len(cx.faces)
    H = numpy.empty((m, n), dtype=object)
    H[:] = "."
    def swap(j0,j1):
        assert H[j0,i] != H[j1,i]
        H[j0,i], H[j1,i] = H[j1,i], H[j0,i]
    for i in range(n):
        #for j in range(m):
        #    if A[i, j]:
        #        H[j, i] = 'x'
        if vdeg[i] == 3:
            labels = list('XYZ')
            shuffle(labels)
            for j in range(m):
                if A[i,j]:
                    H[j,i] = labels.pop()
        elif vdeg[i] == 4:
            labels = list('XZXZ')
            shuffle(labels)
            jdxs = [j for j in range(m) if A[i,j]]
            for j in jdxs:
                H[j,i] = labels.pop()
            j0, j1, j2, j3 = jdxs
            f0 = cx.faces[j0]
            f1 = cx.faces[j1]
            f2 = cx.faces[j2]
            f3 = cx.faces[j3]
            #print(j0, j1, j2, j3, ''.join(H[:,i].transpose()))

            # 0 is next to 1,3 or 2,3 or 1,2:
            if f0.intersect(f1) and f0.intersect(f3):
                if H[j0,i] == H[j1,i]: swap(j0,j3)
                elif H[j0,i] == H[j3,i]: swap(j0,j1)
            elif f0.intersect(f2) and f0.intersect(f3):
                if H[j0,i] == H[j2,i]: swap(j0,j3)
                elif H[j0,i] == H[j3,i]: swap(j0,j2)
            elif f0.intersect(f1) and f0.intersect(f2):
                if H[j0,i] == H[j1,i]: swap(j0,j2)
                elif H[j0,i] == H[j2,i]: swap(j0,j1)
            else:
                assert 0
            #print('       ', ''.join(H[:,i].transpose()))
        else:
            assert 0
    shortstr = lambda H : ('\n'.join(''.join(row) for row in H))
    print(shortstr(H))
    print()
    #H = H[:-1, :]

    from qumba.qcode import QCode, fromstr, linear_independent
    H = fromstr(H)
    H = linear_independent(H)
    code = QCode(H)
    return code


def mutate(cx):
    ch = cx.euler

    for _ in range(3):
        i = randint(0, 2)
        if i==0:
            v = choice(cx.verts)
            cx.bevel(v)
        #elif i==1:
        #    e = choice(cx.edges)
        #    #cx.cut_edge(e)
        elif i==1:
            f = choice(cx.faces)
            cx.barycenter(f)
        elif i==2:
            f = choice(cx.faces)
            cx.cut_face(f)
        assert cx.euler == ch
        print(cx)
    H0, H1 = cx.bdy()
    H0 %= 2
    H1 %= 2
    vdeg = H0.sum(1)
    fdeg = H1.sum(0)
    print(vdeg)
    print(fdeg)
    verts = list(cx.verts)
    for i,v in enumerate(verts):
        if vdeg[i] > 4:
            cx.bevel(v)
    vdeg = (cx.bdy(0)%2).sum(1)
    assert vdeg.min() == 3
    assert vdeg.max() == 4
    v3 = (vdeg==3).sum()
    v4 = (vdeg==4).sum()
    assert len(cx.faces) - v3//2 - v4 == ch
    return cx


def main():
    cx = make_octahedron()
    code = get_code(cx)
    assert code.k == 0

    cx = make_tetrahedron()
    code = get_code(cx)
    assert code.get_params() == (4, 1, 2)

    cx = make_cube()
    code = get_code(cx)
    assert code.get_params() == (8, 3, 2)

    cx = mutate(cx)
    code = get_code(cx)
    print(code.get_params())

    from qumba.distance import distance_z3
    d = distance_z3(code)
    print("d =", d)


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
        main()

    t = time() - start_time
    print("finished in %.3f seconds"%t)
    print("OK!\n")

