#!/usr/bin/env python3

from random import randint, shuffle

import numpy

from bruhat.argv import argv
from bruhat.smap import SMap
from bruhat.util import cross
from bruhat.action import mulclose


def array(m, n):
    return numpy.zeros((m, n), dtype=int)


class Network(object):
    def __init__(self, T):
        N = len(T)
        assert T.shape == (N, N)
        self.N = N
        self.T = T
        self.ceil = T.sum(0)
        for idx in range(1, N):
            assert self.ceil[idx] > 0
        gen = []
        for i in range(N-1):
            values = [0]*(N-1)
            values[i] = 1
            gen.append(self.make_state(values))
        self.gen = gen

    @classmethod
    def build_lattice(cls, m, n):
        N = 1 + m*n
        T = array(N, N)
        lookup = {}
        for i in range(m):
          for j in range(n):
            lookup[i, j] = i*n + j + 1
        get = lambda i, j : lookup.get((i, j), 0)
    
        for i in range(m):
          for j in range(n):
            T[get(i+1, j), get(i, j)] += 1
            T[get(i, j+1), get(i, j)] += 1
            T[get(i-1, j), get(i, j)] += 1
            T[get(i, j-1), get(i, j)] += 1
    
        self = cls(T)
        self.shape = (m, n)
        self.lookup = lookup
        return self

    @classmethod
    def build_line(cls, m):
        N = m+1
        T = array(N, N)
        lookup = {}
        for i in range(m):
            lookup[i] = i+1
        get = lambda i : lookup.get(i, 0)
    
        for i in range(m):
            T[get(i+1), get(i)] += 1
            T[get(i-1), get(i)] += 1
    
        self = cls(T)
        self.shape = (m, 1)
        self.lookup = lookup
        return self
    
    def make_state(self, values=None):
        v = numpy.zeros(self.N, dtype=int)
        if values is not None:
            assert len(values) == self.N-1
            v[1:] = values
        return State(self, v)

    def get_critical(self):
        values = list(self.ceil - 1)[1:]
        return self.make_state(values)

    def all_states(self):
        N = self.N
        ceil = self.ceil
        itemss = [list(range(ceil[idx])) for idx in range(1, N)]
        for values in cross(itemss):
            state = self.make_state(values)
            yield state

    def rand_state(self):
        N = self.N
        ceil = self.ceil
        values = [randint(0, ceil[idx]-1) for idx in range(1, N)]
        return self.make_state(values)

    def _topple(self, v):
        v = v.copy()
        N = self.N
        T = self.T
        ceil = self.ceil
        done = False
        while not done:
            done = True
            for idx in range(1, N):
                if v[idx] < ceil[idx]:
                    continue
                done = False
                # topple location idx
                #assert ceil[idx]
                v[idx] -= ceil[idx]
                v += T[:, idx]
        return v
    

class State(object):
    def __init__(self, network, v):
        v = network._topple(v)
        self.network = network
        self.v = v
        self.key = v[1:].tostring()

    def __str__(self):
        network = self.network
        shape = network.shape
        lookup = network.lookup
        smap = SMap()
        m, n = self.shape
        v = self.v
        for i in range(m):
          for j in range(n):
            smap[i, j+1] = str(v[i*n+j+1])
        smap[0,0] = "["
        smap[m-1,n+1] = "]"
        return str(smap)

    def __str__(self):
        return str(self.v[1:])

    #def __setitem__(self, idx, value): # um... mutable?
    #    self.v[self.lookup[idx]] = value

    #def __getitem__(self, idx):
    #    return self.v[self.lookup[idx]]

    def __add__(self, other):
        assert self.network is other.network
        v = self.v + other.v
        state = State(self.network, v)
        return state
    __mul__ = __add__

    def __eq__(self, other):
        #return str(self.v[1:]) == str(other.v[1:])
        return self.key == other.key

    def __ne__(self, other):
        #return str(self.v[1:]) != str(other.v[1:])
        return self.key != other.key

    def __hash__(self):
        #return hash(str(self.v[1:]))
        return hash(self.key)


def main():

    if 0:
        m = argv.get("m", 3)
        n = argv.get("n", 3)
        L = Network.build_lattice(m, n)
    
        crit = L.get_critical()
    
        zero = L.make_state([2,1,2,1,0,1,2,1,2])
        a = L.make_state([2]*9)
        nega = L.make_state([2,3,2,3,2,3,2,3,2])
    
        assert zero == zero
        assert zero != crit
        assert a+zero == a
        assert zero+zero == zero
        assert a+nega == zero

    else:
        L = Network.build_line(3)
        crit = L.get_critical()

    #states = mulclose([L.get_critical()] + L.gen, verbose=True)
    #print(len(states))

    counts = set()
    states = list(L.all_states())
    print(len(states))

    for b in states:
      for a in L.gen:
        print(a, ":", b, a+b)

    if 0:
        shuffle(states)
        for b in states[:100]:
            c = crit+b
            count = 1
            while c!=crit:
                c = c+b
                count += 1
                assert count < 1000
            if count not in counts:
                print(count, end=" ",flush=True)
                counts.add(count)

    
    if 0:
        group = set()
        count = 0
        for state in L.all_states():
            a = crit + state
            group.add(a)
            count += 1
        print(count)
        print(len(group))



if __name__ == "__main__":

    main()


