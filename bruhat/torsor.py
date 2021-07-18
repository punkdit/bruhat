#!/usr/bin/env python3

"""
Find monoid, group, torsor using abstract tables.
"""

from bruhat.util import cross
from bruhat.argv import argv

import numpy


def table(*shape):
    return numpy.zeros(shape, dtype=int)


def all_mul(n):
    itemss = [tuple(range(n))]*(n-1)*(n-1)
    for idxs in cross(itemss):
        mul = table(n, n)
        for i in range(n):
            mul[0, i] = i # identity is always zero
            mul[i, 0] = i # identity is always zero
        for i in range(n-1):
          for j in range(n-1):
            mul[1+i, 1+j] = idxs[i + (n-1)*j]
        yield mul


def is_assoc(mul):
    n = mul.shape[0]
    for a in range(n):
      for b in range(n):
        for c in range(n):
            if mul[mul[a, b], c] != mul[a, mul[b, c]]:
                return False
    return True


def all_act(n, m):
    itemss = [tuple(range(m))]*(n-1)*m
    for idxs in cross(itemss):
        act = table(n, m)
        for i in range(m):
            act[0, i] = i # identity is always zero
        for i in range(n-1):
          for j in range(m):
            act[1+i, j] = idxs[i + (n-1)*j]
        yield act

def is_action(mul, act):
    n = mul.shape[0]
    m = act.shape[1]
    for g in range(n):
      for h in range(n):
        for x in range(m):
            if act[g, act[h, x]] != act[mul[g, h], x]:
                return False # not an action
    return True 


def is_group(mul):
    n = mul.shape[0]
    for a in range(n):
        for b in range(n):
            if mul[a, b] == mul[b, a] == 0:
                break # found inverse
        else:
            return False # did not find inverse
    return True


def is_torsor(mul, act):
    n = mul.shape[0]
    m = act.shape[1]
    found = set()
    for g in range(n):
      for x in range(m):
        found.add((act[g, x], x))
    return len(found) == m*m
    

def is_self_torsor(mul):
    n = mul.shape[0]
    send = set()
    for a in range(n):
      for b in range(n):
        send.add( (mul[a, b], b) )
    return len(send) == n*n


def main():

    n = argv.get("n", 3)
    m = argv.get("m", n)

    monoids, groups = 0, 0
    for mul in all_mul(n):
        if not is_assoc(mul):
            continue
        print(".", flush=True, end="")
        lhs = is_group(mul)
        rhs = is_self_torsor(mul)
        assert lhs == rhs
        if lhs:
            print("*", end="")
        #print()
        monoids += 1
        groups += int(lhs)
        count = 0
        for act in all_act(n, m):
            if not is_action(mul, act):
                continue
            if is_torsor(mul, act):
                print("T", end="")
            count += 1
        print(count)
    print("monoids:", monoids)
    print("groups:", groups)



if __name__=="__main__":

    main()

