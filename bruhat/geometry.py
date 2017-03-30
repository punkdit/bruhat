#!/usr/bin/env python

import os, sys

import numpy
import networkx as nx

from solve import zeros2, enum2, row_reduce, span, shortstr, shortstrx, solve, rank, find_kernel, find_logops, identity2
import isomorph
from isomorph import Bag, Point, write
from bruhat import BruhatMonoid
from action import Perm, Group
from util import all_subsets, factorial, choose

from argv import Argv
argv = Argv()



#
# http://mathworld.wolfram.com/q-BinomialCoefficient.html
#
def qbinomial(n, m, q=2):
    "n choose_q m"
    assert n>=m>=0, (n, m)

    if n==m or m==0:
        return 1

    if m==1:
        top = 1-q**n
        bot = 1-q
        assert top%bot == 0
        return top//bot

    return q**m * qbinomial(n-1, m) + qbinomial(n-1, m-1)

assert qbinomial(4, 2) == 35


class Geometry(object):

    def __init__(self, incidence, tpmap):
        """
            incidence: list of item pairs (i, j)
            tpmap: dict mapping each item to it's type
        """
        items = set() # all items
        nbd = {} # map item -> list of incident items
        tplookup = {} # map type -> list of items of that type
        for i, j in incidence:
            items.add(i)
            items.add(j)
        for i in items:
            nbd[i] = []
        for i, j in incidence:
            if j not in nbd[i]:
                nbd[i].append(j)
            if i not in nbd[j]:
                nbd[j].append(i)
        for jtems in nbd.values():
            assert len(set(jtems))==len(jtems), "nbd:%s"%nbd # uniq

        for i in items:
          for j in items:
            try:
                i==j
            except ValueError:
                print "cant compare %r and %r" % (i, j)
                raise
        for i in items:
            if i not in nbd[i]:
                nbd[i].append(i)
            for j in nbd[i]:
                if i==j:
                    continue

                assert tpmap[i] != tpmap[j]
                if i not in nbd[j]:
                    nbd[j].append(i)

        for t in tpmap.values():
            tplookup[t] = []
        for item, t in tpmap.items():
            tplookup[t].append(item)
        for jtems in nbd.values():
            assert len(set(jtems))==len(jtems), "nbd:%s"%nbd # uniq

        incidence = [] # rebuild this
        for item, jtems in nbd.items():
            for jtem in jtems:
                incidence.append((item, jtem))
        incidence.sort()

        self.incidence = set(incidence) # incidence relation: list of pairs (i, j)
        self.tpmap = dict(tpmap) # map item -> type
        self.tplookup = tplookup # map type -> list of items
        self.types = list(set(tpmap.values()))
        self.rank = len(self.types)
        self.items = items
        self.nbd = nbd
        self.check_geometry()

    def __str__(self):
        incidence = list(self.incidence)
        incidence.sort()
        ss = []
        for i, j in incidence:
            if i<j:
                ss.append("%s--%s" % (i, j))
        return "Geometry([%s], %s)"%(', '.join(ss), self.tpmap)
    __repr__ = __str__

    def __contains__(self, pair):
        return pair in self.incidence

    @classmethod
    def polygon(cls, n):
        verts = ['p%d'%i for i in range(n)]
        edges = ['l%d'%i for i in range(n)]
        incidence = []
        tpmap = {}
        for i in range(n):
            incidence.append((verts[i], edges[i]))
            incidence.append((verts[i], edges[(i+1)%n]))
            tpmap[verts[i]] = 'p'
            tpmap[edges[i]] = 'l'
        #print incidence
        #print tpmap
        return cls(incidence, tpmap)

    @classmethod
    def simplex(cls, dim):
        tps = range(dim+1)
        items = [tuple(idxs) for idxs in all_subsets(dim+1) if idxs] # skip empty set
        incidence = []
        tpmap = {}
        for i in items:
            tpmap[i] = len(i)-1
            for j in items:
                if set(i).intersection(set(j)) == set(i):
                    incidence.append((i, j))
        #print "simplex: %d"%dim
        #print "incidence:", incidence
        #print "tpmap:", tpmap
        return cls(incidence, tpmap)

    @classmethod
    def cube(cls):
        "the verts, edges and faces of a cube"
        verts = []
        edges = []
        faces = []
        incidence = []
        for i in [0, 1]:
         for j in [0, 1]:
          for k in [0, 1]:
            verts.append((i, j, k))
        for edge in choose(verts, 2):
            i, j = edge
            if sum([abs(i[ii]-j[ii]) for ii in range(3)])==1:
                e = (i, j)
                edges.append(e)
                incidence.append((i, e))
                incidence.append((j, e))
        assert len(edges)==12
        assert len(incidence)==24
        for idx in range(3):
          for face in choose(verts, 4):
            r = face[0][idx]
            for v in face:
                if v[idx] != r:
                    break
            else:
                faces.append(face)
                for v in face:
                    incidence.append((v, face))
                for edge in choose(face, 2):
                    if edge in edges:
                        incidence.append((edge, face))
        assert len(faces)==6
        #for (i,j) in incidence[24:]:
        #    print i, '---', j
        assert len(incidence)==24 + 6*4 + 6*4
        tpmap = {}
        for v in verts:
            tpmap[v] = 'v'
        for e in edges:
            tpmap[e] = 'e'
        for f in faces:
            tpmap[f] = 'f'
        g = Geometry(incidence, tpmap)
        return g

    def is_flag(self, flag):
        if not len(set(flag))==len(flag):
            print "X"
            return False
        for i in flag:
            if not i in self.items:
                print "I"
                return False
            for j in flag:
                if not (i, j) in self.incidence:
                    print "C"
                    return False
        #print "."
        return True

    def ordered_flags(self, flag=[]):
        #yield flag
        nbd = set(self.items)
        for i in flag:
            nbd = nbd.intersection(self.nbd[i])
        #print "all_flags:", flag, nbd
        for i in nbd:
            if i in flag:
                continue
            yield flag+[i]
            for _flag in self.ordered_flags(flag+[i]):
                yield _flag

    def all_flags(self):
        flags = [flag for flag in self.ordered_flags()]
        for flag in flags:
            flag.sort()
        flags = [tuple(flag) for flag in flags]
        flags = list(set(flags))
        flags.sort()
        return flags

    def maximal_flags(self):
        "all maximal flags"
        flags = self.all_flags()
        maximal = []
        for f in flags:
            #print "is maximal %s ?"%(f,)
            for g in flags:
                if f is g:
                    continue
                if set(g).intersection(f) == set(f):
                    #print "%s contains %s" % (g, f)
                    break
            else:
                #print "maximal:", f
                maximal.append(f)
        return maximal

    def flags_type(self, tps):
        "all flags of type tps"
        tps = list(tps)
        if len(tps)==0:
            return [()]
        if len(tps)==1:
            return [(i,) for i in self.tplookup[tps[0]]]
        flags = []
        incidence = self.incidence
        for flag in self.flags_type(tps[:-1]):
            #print "flags_type", flag
            tp = tps[-1]
            for i in self.tplookup[tp]:
                for j in flag:
                    #print "?", (i, j)
                    if (i, j) not in incidence:
                        break
                else:
                    flags.append(flag + (i,))
        return flags

    def check_geometry(self):
        # every maximal flag is a chamber
        for flag in self.maximal_flags():
            if len(flag)!=self.rank:
                print "flag %s is not a chamber (rank=%d)"%(flag, self.rank)
                assert 0

    def is_digon(self):
        if self.rank != 2:
            return False
        tpi = self.types[0]
        tpj = self.types[1]
        for i in self.tplookup[tpi]:
          for j in self.tplookup[tpj]:
            if (i, j) not in self.incidence:
                return False
        return True

    def residue(self, flag):
        items = []
        for i in flag:
            assert i in self.items
        for i in self.items:
            if i in flag:
                continue
            for j in flag:
                if (i, j) not in self.incidence:
                    break
            else:
                items.append(i)
        items = list(set(items))
        incidence = [(i, j) for i in items for j in items if (i, j) in self.incidence]
        assert len(set(incidence))==len(incidence), str(incidence)
        #print "residue:", incidence
        tpmap = dict((i, self.tpmap[i]) for i in items)
        g = Geometry(incidence, tpmap)
        return g

    def get_diagram(self):
        diagram = []
        for (i, j) in choose(self.types, 2):
            tps = [k for k in self.types if k!=i and k!=j]
            assert len(tps)<20
            #print "flag_type:", tps
            for flag in self.flags_type(tps):
                h = self.residue(flag)
                assert h.types == [i, j] or h.types == [j, i]
                if h.is_digon():
                    break
            else:
                diagram.append((i, j))
        diagram.sort()
        return diagram

    def is_projective_plane(self):
        if len(self.types) != 2:
            return False
        ptp, ltp = self.types
        points = self.tplookup[ptp]
        lines = self.tplookup[ltp]
        nbd = self.nbd
        for line in lines:
            for other in lines:
                if other is line:
                    continue
                ps = set(nbd[line]).intersection(set(nbd[other]))
                if len(ps)!=1:
                    #print line
                    #print other
                    #print ps
                    return False
                p = list(ps)[0]
                assert self.tpmap[p]==ptp
        for i in points:
          for j in points:
            if i==j:
                continue
            ilines = set(nbd[i]).intersection(set(nbd[j]))
            if len(ilines)!=1:
                print "L"
                return False
            line = list(ilines)[0]
            assert self.tpmap[line]==ltp
        return True

    def get_graph(self):
        graph = nx.Graph()
        for item in self.items:
            graph.add_node(item)
        for i, j in self.incidence:
            if i!=j and i<j:
                graph.add_edge(i, j)
        return graph

    def rel(self, a, b):
        incidence = self.incidence
        A = zeros2(len(a), len(b))
        for i, aa in enumerate(a):
          for j, bb in enumerate(b):
            A[i, j] = (aa, bb) in incidence
        return A

    def get_bag(self):
        items = list(self.items)
        items.sort()
        tpmap = self.tpmap
        points = []
        lookup = {}
        for idx, item in enumerate(items):
            tp = tpmap[item]
            p = Point(str(tp), idx)
            #p = Point("XX", idx)
            lookup[item] = p
            points.append(p)
        #print points
        for a, b in self.incidence:
            lookup[a].nbd.append(lookup[b])
        bag = Bag(points)
        #for p in points:
        #    print p, p.nbd
        return bag
        
    def get_symmetry(self):
        bag0 = self.get_bag()
        bag1 = self.get_bag()
        return isomorph.search(bag0, bag1)


class System(object):
    """
        Chamber (maximal flag) system
    """

    def __init__(self, geometry):
        flags = list(geometry.maximal_flags())
        tpmap = geometry.tpmap

        graphs = {} # map tp -> graph
        for tp in geometry.types:
            graph = nx.Graph()
            for flag in flags:
                graph.add_node(flag)
                #print "node:", tp, flag
            graphs[tp] = graph

        for flag in flags:
            for _flag in flags:
                assert len(flag)==len(_flag)
                rank = len(flag)
                idxs = [idx for idx in range(rank) if flag[idx]!=_flag[idx]]
                if len(idxs)!=1:
                    continue
                idx = idxs[0]
                assert tpmap[flag[idx]]==tpmap[_flag[idx]]
                tp = tpmap[flag[idx]]
                graphs[tp].add_edge(flag, _flag)

        # union graph
        graph = nx.Graph()
        for i in flags:
            graph.add_node(i)

        for tp, equ in graphs.items():
            for flag, _flag in nx.edges(equ):
                graph.add_edge(flag, _flag, tp=tp)

                assert len(flag)==len(_flag)
                rank = len(flag)
                idxs = [idx for idx in range(rank) if flag[idx]!=_flag[idx]]
                assert len(idxs)==1
                idx = idxs[0]
                assert tpmap[flag[idx]]==tpmap[_flag[idx]]
                assert tp == tpmap[flag[idx]]

        self.flags = flags
        self.tpmap = tpmap
        self.types = geometry.types
        self.graphs = graphs
        self.graph = graph
        self.geometry = geometry
        self.monoid = None
        self._geodesics = {} # cache

    def get_geodesic(self, c, d):
        """
            Find a minimal gallery from c to d, return as a list of types.
        """
        key = (c, d)
        if key in self._geodesics:
            return list(self._geodesics[c, d])
        graph = self.graph
        tpmap = self.tpmap
        path = nx.shortest_path(graph, c, d)
        assert path[0] == c
        assert path[-1] == d
        word = []
        if c==d:
            self._geodesics[key] = list(word)
            return word # <---- return
        assert len(c)==len(d)
        rank = len(c)
        i = 0
        while i+1 < len(path):
            idx = None
            for j in range(rank):
                if path[i][j] != path[i+1][j]:
                    assert idx is None, (idx, j)
                    idx = j
            tp = tpmap[path[i][idx]]
            assert tp == tpmap[path[i+1][idx]]
            word.append(tp)
            if i:
                assert word[i-1] != word[i]
            i += 1
        self._geodesics[key] = list(word)
        return word

    def build_weyl(self, monoid):
        mul = monoid.mul
        words = monoid.words
        print "words:", len(words)
        print words
        self.monoid = monoid
        self.mul = monoid.mul
        self.words = monoid.words
        self.I = I = ""
        weyl = {} # map each pair of flags to the Weyl distance
        flags = self.flags
        gen = monoid.gen
        for i in flags:
          for j in flags:
            a = I
            path = self.get_geodesic(i, j)
            for b in path:
                a = mul[a, gen[b]]
            weyl[i, j] = a
        self.weyl = weyl

    def build(self, letters):
        geometry = self.geometry
        gen = letters[:geometry.rank]
        diagram = geometry.get_diagram()
        desc = {}
        for i, j in diagram:
            desc[gen[i], gen[j]] = 3
        print "BruhatMonoid(%s, %s)" % (gen, desc)
        monoid = BruhatMonoid(gen, desc)
        words = monoid.build()
        self.build_weyl(monoid)

    def __str__(self):
        return "System(%d)"%len(self)
        #return "System(%s, %s)"%(self.flags, self.equs)

    def __len__(self):
        return len(self.flags)

    def __getitem__(self, idx):
        return self.flags[idx]

    def XXget_bag(self):
        flags = self.flags
        equs = self.equs
        lookup = {}
        points = []
        for flag in flags:
            p = Point("f", len(lookup))
            lookup[flag] = p
            points.append(p)
        for tp, flagss in equs.items():
            p = Point("t", len(lookup))
            lookup[tp] = p
            points.append(p)
            for i, flags in enumerate(flagss):
                e = Point("e", len(lookup))
                lookup[tp, i] = e
                points.append(e)
                p.nbd.append(e)
                e.nbd.append(p)
                for flag in flags:
                    e.nbd.append(lookup[flag])
                    lookup[flag].nbd.append(e)
        bag = Bag(points)
        return bag

    def get_symmetry(self):
        bag0 = self.get_bag()
        bag1 = self.get_bag()
        count = 0
        for f in isomorph.search(bag0, bag1):
            #print f
            count += 1
        return count

    def length(self, chain):
        assert len(chain)
        if len(chain)==1:
            return self.I  # <--- return
        monoid = self.monoid
        mul = self.mul
        gen = monoid.gen
        A = self.get_geodesic(chain[0], chain[1])
        for i in range(1, len(chain)-1):
            B = self.get_geodesic(chain[i], chain[i+1])
            A = A+B
        if not len(A):
            return self.I  # <--- return
        if len(A)==1:
            A = gen[A[0]]
        elif len(A)==2:
            A = mul[gen[A[0]], gen[A[1]]]
        elif len(A)==3:
            A = mul[mul[gen[A[0]], gen[A[1]]], gen[A[2]]]
        else:
            a = gen[A[0]]
            for b in A[1:]:
                a = mul[a, gen[b]]
            A = a
        assert A in self.words
        return A

    def get_nchains_0(self, accept):
        flags = self.flags
        if self.I in accept:
            return [(flag,) for flag in flags]
        return []

    def get_nchains_1(self, accept):
        flags = self.flags
        weyl = self.weyl
        chains = []
        for key, value in weyl.items():
            if value in accept:
                chains.append(key)
        return chains

    def get_nchains_2(self, accept):
        flags = self.flags
        weyl = self.weyl
        chains = []
        mul = self.mul
        for k, v in weyl.items():
            a, b = k
            for c in flags:
                if mul[v, weyl[b, c]] in accept:
                    chains.append((a, b, c))
        return chains

    def get_nchains_3(self, accept):
        flags = self.flags
        weyl = self.weyl
        chains = []
        mul = self.mul
        for k, v in weyl.items():
            a, b = k
            for c in flags:
              w = mul[v, weyl[b, c]]
              for d in flags:
                if mul[w, weyl[c, d]] in accept:
                    chains.append((a, b, c, d))
        return chains

    def get_nchains(self, n, accept):
        if n==0:
            return self.get_nchains_0(accept) # <--- return
        if n==1:
            return self.get_nchains_1(accept) # <--- return
        if n==2:
            return self.get_nchains_2(accept) # <--- return
        if n==3:
            return self.get_nchains_3(accept) # <--- return
        flags = self.flags
        nchains = []
        for chain in alltuples((flags,)*(n+1)):
            if self.length(chain) in accept:
                nchains.append(chain)

        assert len(set(nchains))==len(nchains)
        return nchains


class MonoidSystem(System):
    def __init__(self, monoid):
        flags = list(monoid.words)
        
        graphs = {} # map tp -> graph
        for tp in monoid.gen:
            graph = nx.Graph()
            for flag in flags:
                graph.add_node(flag)
                #print "node:", tp, flag
            graphs[tp] = graph

        # this is a thin geometry
        for flag in flags:
            for tp in monoid.gen:
                _flag = monoid.mul[tp, flag]
                graphs[tp].add_edge(flag, _flag)
                #graphs[tp].add_edge(_flag, flag)

        # union graph
        graph = nx.Graph()
        for i in flags:
            graph.add_node(i)

        for tp, equ in graphs.items():
            for flag, _flag in nx.edges(equ):
                graph.add_edge(flag, _flag, tp=tp)

        self.flags = flags
        #self.tpmap = tpmap
        self.types = monoid.gen
        self.graphs = graphs
        self.graph = graph
        self.geometry = None
        self.monoid = monoid
        self._geodesics = {} # cache

        self.build_weyl(monoid)

    def get_geodesic(self, c, d):
        """
            Find a minimal gallery from c to d, return as a list of types.
        """
        key = (c, d)
        if key in self._geodesics:
            return list(self._geodesics[c, d])
        graph = self.graph
        monoid = self.monoid
        word = []
        if c==d:
            self._geodesics[key] = list(word)
            return word # <---- return
        path = nx.shortest_path(graph, c, d)
        assert path[0] == c
        assert path[-1] == d
        print "get_geodesic %20s" % ((c, d),),
        print "path:", path
        i = 0
        while i+1 < len(path):
            a, b = path[i:i+2]
            data = graph.get_edge_data(a, b)
            tp = data['tp']
            #for tp in monoid.gen:
            #    if monoid.mul[tp, a] == b or monoid.mul[tp, b] == a:
            #        break
            #else:
            #    assert 0, (a, b)
            word.append(monoid.gen.index(tp))
            if i:
                assert word[i-1] != word[i], "word=%s, path=%s" % (word, path)
            i += 1
        self._geodesics[key] = list(word)
        return word

    def build_weyl(self, monoid):
        mul = monoid.mul
        words = monoid.words
        print "words:", len(words)
        print words
        self.monoid = monoid
        self.mul = monoid.mul
        self.words = monoid.words
        self.I = I = ""
        weyl = {} # map each pair of flags to the Weyl distance
        flags = self.flags
        gen = monoid.gen
        for i in flags:
          for j in flags:
            a = I
            path = self.get_geodesic(i, j)
            for b in path:
                a = mul[a, gen[b]]
            weyl[i, j] = a
        self.weyl = weyl



def genidx(*shape):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx(*shape[1:]):
                yield (idx,)+_idx


def tesellate_3d(n=3):
    "cubic tesselation of 3-torus by n*n*n cubes and verts"

    assert n<10
    cubes = {}
    verts = {}
    incidence = []
    tpmap = {}
    for i,j,k in genidx(n,n,n):
        cubes[i,j,k] = 'c%d%d%d'%(i,j,k)
        verts[i,j,k] = 'v%d%d%d'%(i,j,k)
    for i,j,k in genidx(n,n,n):
        c = cubes[i,j,k]
        for di,dj,dk in genidx(2,2,2):
            v = verts[(i+di)%n, (j+dj)%n, (k+dk)%n]
            incidence.append((c, v))
    for c in cubes.values():
        tpmap[c] = 'c'
    for v in verts.values():
        tpmap[v] = 'v'
    g = Geometry(incidence, tpmap)
    return g



def fano():
    # https://www.maa.org/sites/default/files/pdf/upload_library/22/Ford/Lam305-318.pdf
    points = range(1, 8)
    lines = [(1,2,4), (2,3,5), (3,4,6), (4,5,7), (5,6,1), (6,7,2), (7,1,3)]
    assert len(lines)==7
    incidence = []
    tpmap = {}
    for point in points:
        tpmap[point] = 'p'
    for line in lines:
        tpmap[line] = 'l'
        for point in line:
            incidence.append((point, line))
    g = Geometry(incidence, tpmap)
    assert g.is_projective_plane()
    return g


def projective(n, dim=2):
    # Take n-dim F_2-vector space
    # points are subspaces of dimension 1
    # lines are subspaces of dimension 2
    # etc.

    def get_key(L):
        vs = [str(v) for v in span(L) if v.sum()]
        vs.sort()
        key = ''.join(vs)
        return key

    assert n>1

    points = []
    for P in enum2(n):
        if P.sum()==0:
            continue
        points.append(P)
    #print "points:", len(points)

    lines = []
    lookup = {}
    for L in enum2(2*n):
        L.shape = (2, n)
        L = row_reduce(L)
        if len(L)!=2:
            continue
        key = get_key(L)
        if key in lookup:
            continue
        lines.append(L)
        lookup[key] = L
    #print "lines:", len(lines)

    spaces = []
    if n>3 and dim>2:
        m = 3
        lookup = {}
        for A in enum2(m*n):
            A.shape = (m, n)
            A = row_reduce(A)
            if len(A)!=m:
                continue
            key = get_key(A)
            if key in lookup:
                continue
            spaces.append(A)
            lookup[key] = A
        #print "spaces:", len(spaces)

    incidence = []
    tpmap = {} 
    for point in points:
        point = str(point)
        tpmap[point] = 0
        #print point

    for L in lines:
        line = str(tuple(tuple(row) for row in L))
        tpmap[line] = 1
        for P in span(L):
            if P.sum():
                incidence.append((str(P), line))

    for A in spaces:
        space = freeze(A)
        tpmap[space] = 2
        for P in span(A):
            if P.sum()==0:
                continue
            incidence.append((str(P), space))
        for L in lines:
            B = solve(A.transpose(), L.transpose())
            if B is not None:
                line = str(tuple(tuple(row) for row in L))
                incidence.append((space, line))

    g = Geometry(incidence, tpmap)
    if dim==2:
        assert g.get_diagram() == [(0, 1)]
    elif dim==3:
        assert n>3
        assert g.get_diagram() == [(0, 1), (1, 2)]
    return g


def test():

    # triangle:
    g = Geometry(
        "aB aC bA bC cB cA".split(), 
        {'a':'p', 'b':'p', 'c':'p', 'A':'l', 'B':'l', 'C':'l'})

    assert len(g.items) == 6
    #print g.flags_type('l') 
    #print g.flags_type('p') 

    #g = g.residue(["a"])

    for n in range(2, 7):
        g = Geometry.polygon(n)
        assert len(g.items)==2*n
        assert len(g.maximal_flags())==2*n

        assert len(g.flags_type(['p'])) == n
        assert len(g.flags_type(['l'])) == n
        assert len(g.flags_type(['p', 'l'])) == 2*n
        if n>3:
            assert not g.is_projective_plane()

    g = Geometry.simplex(2)
    h = g.residue([(0,), (0,1)])
    assert len(h.items)==1
    assert len(g.flags_type([0, 1])) == 6

    for dim in range(4):
        g = Geometry.simplex(dim)
        assert len(g.items)==2**(dim+1)-1
        assert len(g.maximal_flags())==factorial(dim+1)

        h = g.residue([(0,)])

        if dim>0:
            g.residue([(0,), (0,1)])

    g = Geometry.cube()
    d = g.get_diagram()
    assert d == [('e', 'v'), ('f', 'e')] # woo !

    assert Geometry.simplex(4).get_diagram() == [(0, 1), (1, 2), (2, 3)] # A_4 Dynkin diagram !

    #g = tesellate_3d()

    g = fano()
    assert g.get_diagram() == [('p', 'l')]

    n = argv.get("n", 3)
    assert n>=3
    g = projective(n)

    assert len(g.tplookup[0]) == qbinomial(n, 1)
    assert len(g.tplookup[1]) == qbinomial(n, 2)

    if argv.cube:
        g = Geometry.cube()
    elif argv.simplex:
        g = Geometry.simplex(2)
    elif argv.triangle:
        g = Geometry.polygon(3)
    elif argv.projective:
        dim = argv.get("dim", 2)
        g = projective(n, dim)

        if dim>2:
            assert len(g.tplookup[2]) == qbinomial(n, 3)
    elif argv.bruhat:
        if argv.A_2:
            monoid = BruhatMonoid("LP", {("L", "P") : 3})
        elif argv.A_3:
            monoid = BruhatMonoid("LPS", {("L", "P") : 3, ("L", "S"):3})
        elif argv.A_4:
            monoid = BruhatMonoid("LPSH", {("L", "P") : 3, ("L", "S"):3, ("S", "H"):3})
        elif argv.B_2:
            monoid = BruhatMonoid("LP", {("L", "P"):4})
        elif argv.B_3:
            monoid = BruhatMonoid("LPS", {("L", "P"):3, ("P", "S"):4})
        elif argv.D_4:
            monoid = BruhatMonoid("LPSH", {("L", "P"):3, ("L", "S"):3, ("L", "H"):3})
        else:
            return
        monoid.build()
        system = MonoidSystem(monoid)
        magnitude_homology(system=system)

    else:
        return

    if argv.system:
        system = System(g)
        return

    if argv.magnitude:
        magnitude_homology(g)

    if argv.hecke:
        hecke(g)
        #print g.get_symmetry()
        #s = g.get_system()
        #s.get_symmetry()
        return

    print "OK"


def hecke(self):

    bag0 = self.get_bag()
    bag1 = self.get_bag()

    points = [p for p in bag0 if p.desc=='0']
    lines = [p for p in bag0 if p.desc=='1']
    print points
    print lines
    flags = []
    for p in points:
        ls = [l for l in p.nbd if l!=p]
        print p, ls
        for l in ls:
            flags.append((p, l))

    print "flags:", len(flags)

#    for point in bag0:
#        print point, repr(point.desc)

    #items = [(point.idx, line.idx) for point in points for line in lines]
    #items = [(a.idx, b.idx) for a in points for b in points]
    items = [(a, b) for a in flags for b in flags]

    # fixing two points gives the Klein four group, Z_2 x Z_2:
    #for f in isomorph.search(bag0, bag1):
    #    if f[0]==0 and f[1]==1:
    #        print f

    perms = []
    for f in isomorph.search(bag0, bag1):
        #print f
        g = dict((bag0[idx], bag0[f[idx]]) for idx in f.keys())
        perm = {}
        for a, b in items:
            a1 = (g[a[0]], g[a[1]])
            b1 = (g[b[0]], g[b[1]])
            assert (a1, b1) in items
            perm[a, b] = a1, b1
        perm = Perm(perm, items)
        perms.append(perm)
        write('.')
    print

    g = Group(perms, items)
    print "group:", len(g)

    orbits = g.orbits()
    #print orbits
    print "orbits:", len(orbits)

    eq = numpy.allclose
    dot = numpy.dot

    n = len(flags)
    basis = []
    gen = []
    I = identity2(n)
    for orbit in orbits:
        A = zeros2(n, n)
        #print orbits
        for a, b in orbit:
            A[flags.index(a), flags.index(b)] = 1
        print shortstr(A)
        print A.sum()
        #print

        basis.append(A)
        if A.sum()==42:
            gen.append(A)

    assert len(gen)==2
    L = gen[0]
    P = gen[1]
    q = 2

    assert eq(dot(L, L), L+q*I)
    assert eq(dot(P, P), P+q*I)
    assert eq(dot(L, dot(P, L)), dot(P, dot(L, P)))

    #print shortstr(dot(L, dot(P, L)))
    #print
    #print shortstr(dot(P, dot(L, P)))




def uniqtuples(itemss):
    if len(itemss)==0:
        yield ()
    elif len(itemss)==1:
        items = itemss[0]
        for head in items:
            yield (head,)
        return
    for head in itemss[0]:
        for tail in uniqtuples(itemss[1:]):
            if head != tail[0]:
                yield (head,)+tail


def alltuples(itemss):
    if len(itemss)==0:
        yield ()
    elif len(itemss)==1:
        items = itemss[0]
        for head in items:
            yield (head,)
        return
    for head in itemss[0]:
        for tail in alltuples(itemss[1:]):
            yield (head,)+tail

"""
References:
https://golem.ph.utexas.edu/category/2011/06/the_magnitude_of_an_enriched_c.html
https://golem.ph.utexas.edu/category/2015/05/categorifying_the_magnitude_of.html
https://golem.ph.utexas.edu/category/2016/09/magnitude_homology.html
https://golem.ph.utexas.edu/category/2016/08/a_survey_of_magnitude.html#c050913
https://golem.ph.utexas.edu/category/2016/08/monoidal_categories_with_proje.html#c051440
"""

def magnitude_homology(geometry=None, system=None):

    letters = "PLSXYZ"
    I, P, L, S = "", "P", "L", "S"

    if system is None:
        print "find_building"
        assert geometry is not None
        assert geometry.rank <= 3, "really?"
        system = System(geometry)

        system.build(letters)

    flags = system.flags
    print "flags:", len(flags)

    chains = {}

    # accept is a tuple of lengths
    accept = argv.get("accept")
    if accept is not None:
        accept = accept.replace("I", '')
        accept = accept.split(',')
    else:
        accept = (I,L,P)
    if accept in system.words:
        accept = (accept,) # tuple-ize
    for l in accept:
        assert l in system.words
    print "accept:", repr(accept)

    N = argv.get("N", 2)
    for n in range(N):
        chains[n] = system.get_nchains(n, accept)
        print "|%d-chains|=%d" % (n, len(chains[n]))

    bdys = {} # map 
    for n in range(N-1):
        # bdy maps n+1 chains to n chains
        print "bdy: chains[%d] -> chains[%d]" % (n+1, n),
        bdy = {}
        source = chains[n+1]
        target = chains[n]
        for col, chain in enumerate(source):
            assert len(chain)==n+2
            #for i in range(1, n+1): # skip first & last, like in the graph homology
            for i in range(n+2):
                chain1 = chain[:i] + chain[i+1:]
                assert len(chain1)==n+1
                if system.length(chain1) in accept:
                    #print chain1
                    #assert len(set(chain1))==n+1 # uniq NOT !
                    row = target.index(chain1)
                    bdy[row, col] = bdy.get((row, col), 0) + (-1)**i
        print "nnz:", len(bdy), "range:", set(bdy.values())
        bdys[n+1] = bdy

    # bdys[n]: map n-chains to (n-1)-chains

    if not argv.homology:
        return

    for i in range(1, N-1):
        b1, b2 = bdys[i], bdys[i+1]
        b12 = compose(b1, b2)
        #print "len(bdy*bdy):", len(b12.values())
        for value in b12.values():
            assert value == 0, value

        if argv.Z2:
            find_homology_2(b1, b2, len(chains[i-1]), len(chains[i]), len(chains[i+1]))
        else:
            find_homology(b1, b2, len(chains[i-1]), len(chains[i]), len(chains[i+1]))

        #if i==2 and names[l]=="L":
        #    return



def find_homology(g, f, *dims):

    print "find_homology"
    print dims[2], "--f-->", dims[1], "--g-->", dims[0]

    F = numpy.zeros((dims[1], dims[2]))
    for row, col in f.keys():
        v = f[row, col]
        F[row, col] = v
    #print shortstr(F)

    G = numpy.zeros((dims[0], dims[1]))
    for row, col in g.keys():
        v = g[row, col]
        G[row, col] = v
    #print shortstr(G)

    GF = numpy.dot(G, F)
    #print shortstr(GF)
    assert numpy.abs(GF).sum() == 0

    import homology
    if f and g:
        d = homology.bettiNumber(G, F)
        print "MH:", d

# See also:
# http://stackoverflow.com/questions/15638650/is-there-a-standard-solution-for-gauss-elimination-in-python


def find_homology_2(g, f, *dims):

    print "find_homology"
    print dims[2], "--f-->", dims[1], "--g-->", dims[0]

    F = numpy.zeros((dims[1], dims[2]))
    for row, col in f.keys():
        v = f[row, col]
        F[row, col] = v % 2
    #print shortstr(F.transpose())
    #print

    G = numpy.zeros((dims[0], dims[1]))
    for row, col in g.keys():
        v = g[row, col]
        G[row, col] = v % 2
    #print shortstr(G)

    if argv.ldpc:
        from qupy.ldpc.css import CSSCode
        code = CSSCode(Hx=F.transpose(), Hz=G)
        code.build()
        code = code.prune_deadbits()
        print code
        name = argv.get("name")
        if name.endswith(".ldpc"):
            print "saving", name
            code.save(name)

    if argv.logops:
        L = find_logops(G, F.transpose())
        w = min([v.sum() for v in L])
        print "logops weight:", w
        L = find_logops(F, G.transpose())
        w = min([v.sum() for v in L])
        print "logops weight:", w

    GF = numpy.dot(G, F) % 2
    #print shortstr(GF)
    assert numpy.abs(GF).sum() == 0

    F = F.astype(numpy.int32)
    G = G.astype(numpy.int32)

    print "rank(f)=%d, rank(g)=%d" % (rank(F), rank(G))
    if g:
        kerg = find_kernel(G)
        print "ker(g):", len(kerg)
    #if f:
    #    kerf = find_kernel(F)
    #    print "ker(f):", len(kerf)




def compose(g, f):
    # first f then g
    h = {}
    #print g.keys()
    for row, c in g.keys():
        for r, col in f.keys():
            if c==r:
                #print row, col
                h[row, col] = h.get((row, col), 0) + g[row, c] * f[r, col]

    return h




def freeze(A):
    items = [A.shape]
    for idx in genidx(*A.shape):
        items.append((idx, A[idx]))
    items = tuple(items)
    return items

def thaw(items):
    #idx, value = items[0]
    shape = items[0]
    A = numpy.zeros(shape, dtype=numpy.int32)
    for (idx, value) in items[1:]:
        A[idx] = value
    return A


def all_chains(n, above, chain=()):
    elements = above.keys()

    if n==1:
        for a in elements:
            yield (a,)
    elif n==2:
        for a in elements:
          for b in above[a]:
            yield (a, b)
        return

    # build a nested for-loop
    lines = [
        "def get_all(elements, above):",
        "  for a0 in elements:",
    ]
    indent = '   '
    for i in range(n-1):
        lines.append(indent + "for a%d in above[a%d]:" % (i+1, i))
        indent += ' '
    tpl = "(" + ','.join('a%d'%i for i in range(n)) + ")"
    lines.append(indent + "yield %s"%tpl)
    #print '\n'.join(lines)
    exec ''.join(l+'\n' for l in lines)

    for chain in get_all(elements, above):
        yield chain


def poset_homology(elements, order):

    above = {}
    for a in elements:
        above[a] = [b for b in elements if order[a, b] and b!=a] # a <= b

    r = len(elements) # 0-th homology
    homology = [r]
    n = 1
    while 1:

        chains = set(all_chains(n+1, above))
        for chain in chains:
            assert len(set(chain))==len(chain)==n+1

        size = len(chains)
        if size==0:
            break
        homology.append(size)
        r += (-1)**n * size

        n += 1

    print "poset_homology:", homology, "=", r

    return r


"""
TWF links:

ADE classifications and Simple Lie Algebras:
http://math.ucr.edu/home/baez/week64.html

http://math.ucr.edu/home/baez/week162.html # jordan algebras

Lie groups & geometry, Borel groups, parabolic subgroups:
http://math.ucr.edu/home/baez/week178.html
http://math.ucr.edu/home/baez/week180.html

the finite case & q-polynomials:
http://math.ucr.edu/home/baez/week186.html

Projective planes, relation to finite fields:
http://math.ucr.edu/home/baez/week145.html

"""

"""
buildings as enriched categories:
https://golem.ph.utexas.edu/category/2010/02/3000_and_one_thing_to_think_ab.html#c032006
https://ncatlab.org/nlab/show/building#buildings_as_metric_spaces

magnitude of enriched category:
https://golem.ph.utexas.edu/category/2011/06/the_magnitude_of_an_enriched_c.html

https://golem.ph.utexas.edu/category/2016/08/monoidal_categories_with_proje.html#c051238

magnitude homology of a graph:
https://golem.ph.utexas.edu/category/2015/05/categorifying_the_magnitude_of.html

magnitude homology in general:
https://golem.ph.utexas.edu/category/2016/09/magnitude_homology.html#more
"""



if __name__ == "__main__":

    if argv.profile:
       import cProfile as profile
       profile.run("test()")

    else:
       test()



