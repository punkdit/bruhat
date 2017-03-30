#!/usr/bin/env python

import sys
#from heapq import heappush, heappop, heapify
from random import randint, choice, seed

try:
    import numpy
    #import scipy.sparse.linalg as la
except ImportError:
    print "numpy not found"


from util import write


class Point(object):
    """
        A Point is a vertex in a Graph (an undirected graph).
        Each Point has a "desc", this is any distinguishing 
        characteristic (colour/type, etc.)
        as respected by isomorphisms of Graph's.
        The "desc" can be any string.
    """
    def __init__(self, desc, idx, nbd=None, colour="", **kw):
        self.desc = desc
        self._colour = colour
        self._desc = {} # cache get_desc
        self.idx = idx
        if nbd is None:
            nbd = []
        self.nbd = nbd
        self.__dict__.update(kw)

    def __str__(self):
        return "Point(desc=%s, idx=%s, nbd=%s)"%(
            self.desc, self.idx, [p.idx for p in self.nbd])
    __repr__ = __str__

    def get_colour(self):
        return self._colour

    def set_colour(self, colour):
        self._desc = {} # clear cache
        self._colour = colour
    colour = property(get_colour, set_colour)

    def get_desc(self, depth=1, source=None):
        assert self.nbd is not None
        assert depth>=0
        assert depth<=1
        #_desc = self._desc.get(depth)
        #if _desc:
        #    return _desc
        desc = self.desc+str(self._colour)
        if depth==0:
            #self._desc = desc
            return desc
        if source is None:
            source = []
        else:
            assert self not in source
        descs = [a.get_desc(depth-1, source+[self]) for a in self.nbd if a not in source]
        descs.sort()
        desc = "%s[%s]"%(desc, ' '.join(descs))
        #self._desc = desc
        return desc

    #def __str__(self):
    #    return "Point(%s: %s)"%(self.desc, descs)


class Graph(object):
    """
        Undirected graph.
    """
    def __init__(self, points=[], **attrs):
        self.__dict__.update(attrs)
        self.descs = {}  # cache, map point -> desc
        self.deps = None # map point -> list of points
        self.attrs = dict(attrs)
        self.points = list(points)
        for i, point in enumerate(points):
            assert point.idx == i

    def add(self, desc='', **kw):
        "add a Point"
        assert not self.descs
        assert self.deps is None
        i = len(self.points)
        point = Point(desc, i, **kw)
        self.points.append(point)
        return point

    def add_directed(self, pi, pj, desc='directed'):
        "encode a directed edge using a path with two extra (coloured) Point's"
        pa = self.add("%s_a"%desc)
        pb = self.add("%s_b"%desc)
        self.join(pi, pa)
        self.join(pa, pb)
        self.join(pb, pj)

    def __str__(self):
        return "Graph(%s)"%(self.points,)

    def __len__(self):
        return len(self.points)

    def __getitem__(self, idx):
        return self.points[idx]

    def join(self, pi, pj):
        points = self.points
        if type(pi) in [int, long]:
            pi = points[pi]
        if type(pj) in [int, long]:
            pj = points[pj]
        if pi not in pj.nbd:
            pj.nbd.append(pi)
        if pj not in pi.nbd:
            pi.nbd.append(pj)

    @classmethod
    def build(cls, Gx):
    
        m, n = Gx.shape
        points = []
        for i in range(m):
            g = Gx[i]
            assert g.sum()==4
            weights = []
            for j in numpy.where(g)[0]:
                weights.append(Gx[:, j].sum())
            weights.sort()
            desc = ''.join(str(w) for w in weights)
            a = Point(desc, i)
            points.append(a)
        #print [a.desc for a in points]
        
        for i in range(m):
            g = Gx[i]
            a = points[i]
            for j in numpy.where(g)[0]:
                for i1 in numpy.where(Gx[:, j])[0]:
                    if i1 != i:
                        a.nbd.append(points[i1])

        return cls(points, m=m, n=n)

    def map(self, fn):
        points = [None]*len(self)
        for p in self.points:
            p = Point(p.desc, fn[p.idx])
            points[p.idx] = p
        for p in self.points:
            for p1 in p.nbd:
                points[fn[p.idx]].nbd.append(points[fn[p1.idx]]) # whoops.. tricky
        
        return self.__class__(points, **self.attrs)

    def get_desc(self, depth=1):
        return [v.get_desc(depth) for v in self.points]

    def get_stats(self, depth=1):
        stats = {}
        for point in self:
            desc = point.get_desc(depth)
            stats[desc] = stats.get(desc, 0) + 1
        return stats

    # ---------- HOTSPOT -----------------------------> 
    def get_orbits(self, depth=1):
        orbits = {}
        assert depth==1
        if self.deps is None:
            deps = {}
            for p in self.points:
                deps[p] = [p]+p.nbd # 1-neighbours
            self.deps = deps
        descs = self.descs
        for p in self.points:
            desc = descs.get(p)
            if desc is None:
                desc = p.get_desc(depth)
                descs[p] = desc
            orbit = orbits.setdefault(desc, [])
            orbit.append(p)
        return orbits # map desc -> list of points

    def set_colour(self, p, colour=''):
        if colour:
            assert p.colour==''
        else:
            assert p.colour
        p.colour = colour

        for p in self.deps[p]:
            self.descs[p] = None # clear cache

Bag = Graph # backwards compat


class Tanner(Graph):
    # This is the Tanner graph

    @classmethod
    def build(cls, Gx, Gz=None):
        if Gz is not None:
            return cls.build2(Gx, Gz)

        m, n = Gx.shape
        checks = [Point('c', i) for i in range(m)]
        bits = [Point('b', i+m) for i in range(n)]
        for i in range(m):
            for j in range(n):
                if Gx[i, j]==0:
                    continue
                checks[i].nbd.append(bits[j])
                bits[j].nbd.append(checks[i])
        return cls(checks+bits, m=m, n=n)

    @classmethod
    def build2(cls, Gx, Gz):
        # This is the Tanner graph
        mx, n = Gx.shape
        mz, n = Gz.shape
        xchecks = [Point('x', i, row=i) for i in range(mx)]
        zchecks = [Point('z', i+mx, row=i) for i in range(mz)]
        bits = [Point('b', i+mx+mz, row=i) for i in range(n)]
        for i in range(mx):
            for j in range(n):
                if Gx[i, j]==0:
                    continue
                xchecks[i].nbd.append(bits[j])
                bits[j].nbd.append(xchecks[i])
        for i in range(mz):
            for j in range(n):
                if Gz[i, j]==0:
                    continue
                zchecks[i].nbd.append(bits[j])
                bits[j].nbd.append(zchecks[i])
        return cls(xchecks+zchecks+bits, mx=mx, mz=mz, n=n)

    def shortstr(self):
        m, n = self.m, self.n
        rows = []
        for i in range(m): # checks
            row = ['.']*n
            p = self.points[i]
            for p1 in p.nbd:
                row[p1.idx-m] = '1'
            row = ''.join(row)
            rows.append(row)
        return '\n'.join(rows)


def from_sparse_ham(n, H):
    points = []
    for i in range(n):
        p = Point('(%s)'%H[i, i], i)
        points.append(p)
    for i, j in H.keys():
        if i!=j:
            points[i].nbd.append(points[j])
    bag = Graph(points)
    return bag


def from_ham(H, syndromes=None):
    if syndromes is not None:
        return from_ham_syndromes(H, syndromes) # <------ return
    n = len(H)
    points = []
    for i in range(n):
        p = Point('(%s)'%H[i, i], i)
        points.append(p)
    for i in range(n):
      for j in range(n):
        if i==j:
            continue
        if H[i, j]:
            points[i].nbd.append(points[j])
    bag = Graph(points)
    return bag


def from_ham_syndromes(H, syndromes):
    n = len(H) # dimension of state space
    assert len(syndromes)==n # one syndrome for each basis vector
    m = len(syndromes[0]) # each syndrome has m check values
    points = []
    for i in range(n):
        p = Point('(%s)'%H[i, i], i)
        points.append(p)
    checks = []
    for i in range(m):
        c = Point('c', n+i)
        checks.append(c)
    for i in range(n):
      for j in range(n):
        if i==j:
            continue
        if H[i, j]:
            points[i].nbd.append(points[j])
      for j in range(m):
        if syndromes[i][j]:
            points[i].nbd.append(checks[j])
            checks[j].nbd.append(points[i])
    bag = Graph(points+checks)
    return bag




def get_perm(m, n, fn):

    U = numpy.zeros((m, m), dtype=int)
    for i in range(m):
        j = fn[i]
        U[i, j] = 1

    V = numpy.zeros((n, n), dtype=int)
    for i in range(n):
        j = fn[i+m]-m
        V[j, i] = 1

    return U, V


def search_recursive(bag0, bag1, fn=None, depth=1):

    assert depth>0

    if fn is None:
        fn = {}
        if len(bag0)!=len(bag1):
            return

        assert bag0 is not bag1

    orbits0 = bag0.get_orbits(depth)
    orbits1 = bag1.get_orbits(depth)

    if len(orbits0) != len(orbits1):
        return

    keys0 = orbits0.keys()
    keys1 = orbits1.keys()
    keys0.sort()
    keys1.sort()
    if keys0 != keys1:
        return

    idx = len(fn)

    # choose any uncoloured bag0 point
    p = bag0.points[idx]
    assert p.colour == ''

    key = p.get_desc(depth)
    orbit = orbits1[key]

    #p.colour = str(idx)
    bag0.set_colour(p, str(idx))

    # go through each candidate in bag1
    for p1 in orbit:
        assert p1.colour == ''
    
        #p1.colour = str(idx)
        bag1.set_colour(p1, str(idx))
    
        assert fn.get(idx) is None
        fn[idx] = p1.idx
    
        if len(fn) == len(bag0):
            yield dict(fn)

        else:

            for _fn in search_recursive(bag0, bag1, fn, depth):
                yield _fn

        del fn[idx]
        assert len(fn) == idx

        #p1.colour = ''
        bag1.set_colour(p1)

    #p.colour = ''
    bag0.set_colour(p, '')


class Backtrack(Exception):
    pass

class State(object):
    def __init__(self, bag0, bag1, depth=1):
        orbits0 = bag0.get_orbits(depth) # map: desc -> list of points
        orbits1 = bag1.get_orbits(depth) # map: desc -> list of points
    
        if len(orbits0) != len(orbits1):
            raise Backtrack() # <-------------- raise
    
        keys0 = orbits0.keys()
        keys1 = orbits1.keys()
        keys0.sort()
        keys1.sort()
        if keys0 != keys1:
            raise Backtrack() # <-------------- raise
        self.bags = bag0, bag1
        self.orbitss = orbits0, orbits1
        self.keyss = keys0, keys1
        self.idx0 = None
        self.depth = depth

    def choose(self, idx0):
        assert self.idx0 is None
        assert idx0 is not None

        bag0, bag1 = self.bags
        p0 = bag0.points[idx0]
        assert p0.colour == ''
    
        key0 = p0.get_desc(self.depth)
        self.orbit1 = self.orbitss[1][key0]
        assert self.orbit1 # otherwise: wtf?
        self.idx0 = idx0 # source index: this is constant
        self.idx1 = 0 # search target index
        self.p0 = p0
        self.p1 = None

    def choose_best(self):
        XXX
        orbits0 = self.orbitss[0]
        items = orbits0.items()
        items.sort(key = lambda item : len(item[1]))
        p = items[0][1][0] # first guy in smallest orbit
        self.choose(p.idx)
        return p.idx

    def do(self, fn):
        bag0, bag1 = self.bags
        # make assignment: idx0 -> idx1
        p0 = self.p0
        #assert p0.colour == ''
        #p0.colour = str(self.idx0)
        bag0.set_colour(p0, str(self.idx0))

        p1 = self.orbit1[self.idx1]
        #assert p1.colour == ''
        #p1.colour = str(self.idx0)
        bag1.set_colour(p1, str(self.idx0))
    
        assert fn.get(self.idx0) is None
        fn[self.idx0] = p1.idx
        assert self.p1 is None
        self.p1 = p1

    def undo(self, fn):
        bag0, bag1 = self.bags
        # undo assignment
        del fn[self.idx0]
        assert self.p1 is not None
        p0 = self.p0
        p1 = self.p1
        assert p1.colour==str(self.idx0)
        assert p0.colour==str(self.idx0)
        #p0.colour = ''
        #p1.colour = ''
        bag0.set_colour(p0)
        bag1.set_colour(p1)
        self.p1 = None

    def next(self):
        assert self.p1 is None
        self.idx1 += 1
        if self.idx1 >= len(self.orbit1):
            raise Backtrack() # <-------------- raise


def search(bag0, bag1, depth=1, fn=None, verbose=False):

    assert bag0 is not bag1
    if len(bag0) != len(bag1):
        return

    # doesn't help any:
    #if bag0.get_stats() != bag1.get_stats():
    #    return

    if fn is None:
        fn = {}

    remain = range(len(bag0))

    orbits = bag0.get_orbits(depth)
    bag1.get_orbits()

    keys = orbits.keys()
    keys.sort(key = lambda key : len(orbits[key]))
    remain = []
    for key in keys:
        for p in orbits[key]:
            if p.idx not in fn:
                remain.append(p.idx)

    #for idx in fn.keys():
    #    remain.remove(idx)
    remain.sort()

    for idx in fn:
        bag0.set_colour(bag0[idx], str(idx))
        bag1.set_colour(bag1[fn[idx]], str(idx))

    try:
        state = State(bag0, bag1, depth)
    except Backtrack:
        return

    idx = remain.pop(0)
    state.choose(idx)
    #idx = remain.pop(randint(0, len(remain)-1))
    #state.choose(idx)
    #idx = state.choose_best()
    #remain.remove(idx)

    stack = [state]

    while stack:

        if verbose:
            print "SEARCH", len(stack)

        for idx in remain:
            assert fn.get(idx) is None

        assert len(remain)+len(fn)+1==len(bag0)

        state = stack[-1]
        state.do(fn)

        assert len(remain)+len(fn)==len(bag0)

        if verbose:
            print fn

        if len(fn) == len(bag0):
            if verbose:
                print "FOUND"
            yield dict(fn)

        else:
            # try to add another state
            try:
                _state = State(bag0, bag1, depth)
                #idx = remain.pop(randint(0, len(remain)-1))
                idx = remain.pop(0)
                _state.choose(idx)
                #idx = _state.choose_best()
                #remain.remove(idx)
                stack.append(_state)
                if verbose: 
                    print "PUSH"
                continue
    
            except Backtrack:
                if verbose: 
                    print "BACK"
                # the above do() doesn't work
                pass

        # next
        while stack:
            state = stack[-1]
            if verbose:
                print "UNDO"
            assert len(remain)+len(fn)==len(bag0)
            state.undo(fn)
            assert len(remain)+len(fn)+1==len(bag0)
            try:
                if verbose:
                    print "NEXT"
                state.next()
                break # ok, finished backtracking
            except Backtrack:
                if verbose:
                    print "POP"
                state = stack.pop() # discard this guy
                #remain.append(state.idx0)
                remain.insert(0, state.idx0)



def all_autos(Gx):
    #Gx = parse(gcolor_gauge)
    m, n = Gx.shape

    bag0 = Tanner.build(Gx)
    bag1 = Tanner.build(Gx)

    for fn in search(bag0, bag1):
        U, V = get_perm(m, n, fn)
        yield U, V


def peterson_graph():

    inside = [Point('', i) for i in range(5)]
    outside = [Point('', i+5) for i in range(5)]

    bag = Graph(inside+outside)

    for i in range(5):

        bag.join(i, (i+2)%5)
        bag.join(i, (i+3)%5)
        bag.join(i, i+5)

        if i<4:
            bag.join(i+5, i+6)
        else:
            bag.join(i+5, i+1)

    return bag


def cyclic_graph():
    n = 5

    points = [Point('', i) for i in range(n)]
    bag = Graph(points)

#    for i in range(n):
#        points[i].nbd.append(points[(i+1)%n])
#        points[(i+1)%n].nbd.append(points[i])

    for i in range(n):
        bag.add_directed(points[i], points[(i+1)%n])

    return bag



gcolor_gauge = """
1111...........
11..11.........
1.1.1.1........
..11..11.......
.1.1.1.1.......
....1111.......
11......11.....
1.1.....1.1....
........1111...
..11......11...
.1.1.....1.1...
1...1...1...1..
........11..11.
.1...1...1...1.
....11......11.
........1.1.1.1
..1...1...1...1
....1.1.....1.1
"""

gcolor_stab = """
11111111.......
1111....1111...
11..11..11..11.
1.1.1.1.1.1.1.1
"""

cube_ham = """
6111....
14..11..
1.4.1.1.
1..4.11.
.11.2..1
.1.1.2.1
..11..21
....1110
"""


def parse(s):
    s = s.replace('.', '0') 
    lines = s.split()
    lines = [l.strip() for l in lines if l.strip()]
    rows = [list(int(c) for c in l) for l in lines]
    if rows:
        n = len(rows[0])
        for row in rows:
            assert len(row)==n, "rows have varying lengths"
    a = numpy.array(rows, dtype=numpy.int32)
    return a


def test():

    # Find rotation symmetry of the code. It's S_4 with order 24.

    Gx = parse(gcolor_gauge)
    m, n = Gx.shape

    bag0 = Tanner.build(Gx)
    bag1 = Tanner.build(Gx)

    #global search
    #search = search_recursive

    count = 0
    for fn in search(bag0, bag1):
        #print "iso", fn
        bag = bag0.map(fn)
        #print bag.shortstr()
        U, V = get_perm(m, n, fn)
        Gx1 = numpy.dot(U, numpy.dot(Gx, V))
        assert numpy.abs(Gx-Gx1).sum()==0
        count += 1
    #print "count:", count
    assert count == 24

    # S_3 symmetry of cubical hamiltonian
    depth = 1
    H = parse(cube_ham)
    bag0 = from_ham(H)
    bag1 = from_ham(H)
    count = 0 
    for fn in search(bag0, bag1, depth=depth):
        count += 1
    assert count == 6

    bag0 = peterson_graph()
    bag1 = peterson_graph()
    assert len(list(search(bag0, bag1, depth=1))) == 120

    # directed graph
    bag0 = cyclic_graph()
    bag1 = cyclic_graph()
    assert len(list(search(bag0, bag1))) == 5


from argv import Argv 
argv = Argv()
seed(0)


if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("test()")

    else:
        test()

    print "OK"



