#!/usr/bin/env python

"""
refs:
https://arxiv.org/pdf/math/0505518
https://arxiv.org/pdf/math/0102166

Toric Topology. Chapter 1: Geometry and combinatorics of polytopes
Victor M. Buchstaber
Taras E. Panov

"""

from bruhat.argv import argv

class Space(object):
    def __init__(self, n):
        self.n = n
        self.pairs = [(i,j) for i in range(n) for j in range(n) if i!=j]

    def le(self, a, b):
        "a <= b ?"
        n = self.n
        i0, j0 = a
        i1, j1 = b
        assert 0<=i0<n
        assert 0<=j0<n
        assert 0<=i1<n
        assert 0<=j1<n
        assert i0!=j0
        assert i1!=j1
        if i0>j0:
            j0 = j0+n
        if i1>j1:
            j1 = j1+n
        assert i0<j0
        assert i1<j1
        assert i0<n
        assert i1<n
        return (i0>=i1 and j0<=j1) or (i0+n>=i1 and j0+n<=j1) or (i0>=i1+n and j0<=j1+n)

    def disjoint(self, a0, a1):
        n = self.n
        i0, j0 = a0
        i1, j1 = a1
        if i0>j0:
            j0 = j0+n
        if i1>j1:
            j1 = j1+n
        assert i0<j0
        assert i1<j1
        assert i0<n
        assert i1<n
        for b0 in [(i0,j0), (i0+n,j0+n)]:
          for b1 in [(i1,j1), (i1+n,j1+n)]:
            ii0, jj0 = b0
            ii1, jj1 = b1
            if jj1>=ii0 and jj0>=ii1:
                return False
        #return j1 < i0 or j0 < i1
        return True

    def accept(self, a, b):
        return a!=b and (self.le(a, b) or self.le(b, a) or self.disjoint(a,b))

    def find(self, pairs=[]):
        if len(pairs) == self.n-1:
            yield list(pairs)
            return
        assert len(pairs) < self.n-1
        for a in self.pairs:
            for b in pairs:
                if not self.accept(a, b):
                    break
            else:
                for stack in self.find(pairs + [a]):
                    yield stack

    def get_verts(self):
        verts = set()
        for stack in self.find():
            #print(stack)
            stack.sort()
            stack = tuple(stack)
            verts.add(stack)
        verts = list(verts)
        verts.sort()
        return verts


def main():
    n = 3
    space = Space(n)

    pairs = space.pairs
    le = space.le
    assert le((1,2), (1,0))
    assert le((2,0), (1,0))
    assert not le((0,1), (1,0))

    assert le((0,1), (2,1))
    assert space.accept( (0,1), (2,1) )

    verts = space.get_verts()
    #for v in verts:
    #    print(v)
    assert len(verts) == 6, len(verts)
    #print(verts)

    assert Space(4).accept((0,1),(2,3))

    for n in [3,4,5,6,7]:
        space = Space(n)
        le = space.le
        accept = space.accept
        verts = space.get_verts()
    
        counts = []
        for pair in space.pairs:
            c = len([v for v in verts if pair in v])
            counts.append(c)
        counts.sort()
        print("n=%d"%n, counts, len(counts))
        print()

    #for n in range(2,10):
    #    print(len(Space(n).get_verts()))


if __name__ == "__main__":

    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()


