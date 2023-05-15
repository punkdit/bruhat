#!/usr/bin/env python

import numpy

from random import randint
from time import time


N = 34
A = numpy.zeros((2**N,), dtype=int)

K = 2**20
print("warmup..")
count = 0
while count < 2**N:
    A[count : count + K] = randint(0, 2**63)
    count += K

def pump(total):
    count = 0
    while count < total:
        i = randint(0, len(A)-K)
        A[i:i+K] = randint(0, 2**63)
        count += 1

print("go")
M = 2**16
start = time()
pump(M)
t = time() - start
speed = M * K * 8 / t
print("speed:", speed/(1024*1024), "MB/s")


