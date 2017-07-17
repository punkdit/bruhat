#!/usr/bin/env python3

import sys, os
from math import sqrt, floor


def isprime_slow(n):
    assert n>0
    assert int(n)==n
    if n==1:
        return False

    i = floor(sqrt(n))
    #print("isprime", i, n)

    for j in range(2, i+1):
        if n%j == 0:
            return False

    return True


def sieve(n, ps=None):
    items = [0]*n
    p = 2

    while p**2 < n:
        i = 2
        while p*i < n:
            items[p*i] = 1
            i += 1

        p += 1
        while p < n and items[p]:
            p += 1

    ps = [i for i in range(2, n) if items[i]==0]
    return ps


class Sieve(object):
    def __init__(self):
        self.ps = []
        self.ms = []
        self.n = 2

    def pull(self):
        n = self.n
        ps = self.ps
        ms = self.ms

        

        self.n = n
        


#class Form(object):
#    def __init__(self, 


def main():
    for i in [1, 4, 6, 8, 9]:
        assert not isprime_slow(i), i
    
    for i in [2, 3, 5, 7, 11]:
        assert isprime_slow(i), i


    N = 10000
    ps = sieve(N)
    for i in range(2, N):
        assert isprime_slow(i) == (i in ps)

    N = 10000000
    ps = sieve(N)


if __name__ == "__main__":

    main()

