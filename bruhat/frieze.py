#!/usr/bin/env python3

"""
Coxeter's frieze patterns: 
https://arxiv.org/abs/1503.05049

"""

from random import choice

from bruhat.smap import SMap
from bruhat.argv import argv


class Frieze(object):


    def __init__(self, m):
        self.values = {}
        self.m = m
        self.jmax = 0

    def __getitem__(self, ij):
        i, j = ij
        assert (i+j)%2 == 0
        assert 0 <= i <= self.m+1
        if i==0 or i==self.m+1:
            return 1
        value = self.values.get((i, j))
        return value

    def __setitem__(self, ij, value):
        i, j = ij
        assert (i+j)%2 == 0
        assert 0 < i < self.m+1
        if value is None:
            del self.values[ij]
        else:
            self.values[ij] = value
            self.jmax = max(self.jmax, j)
        return value

    #def keys(self):

    def get_bdy(self):
        bdy = set()
        m = self.m
        for i in range(1, m+1):
          for j in range(self.jmax+2):
            if (i+j)%2==0:
                continue
            idxs = [(i+1, j), (i-1, j), (i, j+1), (i, j-1)]
            nbd = [self[i, j] for i, j in idxs]
            if nbd.count(None)==1:
                bdy.add(idxs[nbd.index(None)])
        return bdy

    def check(self):
        m = self.m
        for i in range(1, m+1):
          for j in range(self.jmax+2):
            if (i+j)%2==0:
                continue
            idxs = [(i+1, j), (i-1, j), (i, j+1), (i, j-1)]
            N, S, E, W = [self[i, j] for i, j in idxs]
            if None not in [N, S, E, W]:
                assert E*W - N*S == 1

    def fill1(self, i, j):
        N, S, W = self[i+1, j-1], self[i-1, j-1], self[i, j-2]
        if None not in [N, S, W]:
            assert (1+N*S)%W == 0
            self[i, j] = (1+N*S)//W
            return
        N, S, E = self[i+1, j+1], self[i-1, j+1], self[i, j+2]
        if None not in [N, S, E]:
            assert (1+N*S)%E == 0
            self[i, j] = (1+N*S)//E
            return
        if i>1:
            N, W, E = self[i+2, j], self[i+1, j-1], self[i+1, j+1]
            if None not in [N, W, E]:
                assert (E*W-1)%N == 0
                self[i, j] = (E*W-1)//N
                return
        if i<self.m:
            S, W, E = self[i-2, j], self[i-1, j-1], self[i-1, j+1]
            if None not in [S, W, E]:
                assert (E*W-1)%S == 0
                self[i, j] = (E*W-1)//S
                return
        assert 0

    def fill(self, count=50):
        for _ in range(count):
            bdy = self.get_bdy()
            if not bdy:
                break
            for (i, j) in bdy:
                self.fill1(i, j)

    def __str__(self):
        i0 = i1 = j0 = j1 = 0
        values = self.values
        s = SMap()
        m = self.m
        jmax = max(4*m, self.jmax)
        dj = 1
        val = max(list(values.values()))
        if val > 9:
            dj = 2
        if val > 99:
            dj = 3
        for i in range(m+2):
          for j in range(jmax):
            if (i+j)%2==0:
                value = self[i, j]
                c = str(value) if value is not None else '.'
                s[i, dj*j] = c
        for (key, value) in values.items():
            i, j = key
            s[i, dj*j] = str(value)
        for (i, j) in self.get_bdy():
            s[i, dj*j] = "*"
        return str(s)


def main():

    m = argv.get("m", 3)

    f = Frieze(m)
    j = m + (m%2) + 1
    for i in range(1, m+1):
        f[i, j] = 1
        j += choice([-1, +1])
    print(f)
    print()

    f.fill()
    s = str(f).replace(" 1 ", " . ")
    print(s)

    f.check()

if __name__ == "__main__":

    main()




