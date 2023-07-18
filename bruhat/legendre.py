#!/usr/bin/env python

from bruhat.action import Group, Perm
from bruhat.util import all_primes
from bruhat.smap import SMap
from bruhat.argv import argv


def dump(table):
    smap = SMap()
    for (i,j) in table.keys():
        smap[i, 3*j] = str(table[i,j]+1).rjust(3)
    #return str(smap)
    print(smap)
    print()

def qresidue(p, q):
    assert 0<p<q
    for i in range(1, q):
        if (i**2)%q == p:
            return True
    return False


# https://www.youtube.com/watch?v=X63MWZIN3gM&t=1656s
def legendre_mathologer(p, q):
    #print("legendre", p, q)

    items = [(i, j) for i in range(p) for j in range(q)]
    #print(items)
    table = dict((ij, idx) for (idx,ij) in enumerate(items))
    #dump(table)

    sable = {}
    count = 0
    for idx in range(p*q):
        sable[idx%p, idx%q] = table[idx%p, idx//p]
        count += 1
    #dump(sable)
    perm = [sable[ij] for ij in items]
    perm = dict(enumerate(perm))
    return Perm(perm, list(range(len(perm)))).sign()


def rowcol(p, q):
    #print("rowcol", p, q)
    items = list(range(p*q))
    #print(items)
    jtems = [i*q+j for j in range(q) for i in range(p)]
    #print(jtems)
    perm = dict(zip(jtems, items))
    perm = Perm(perm, list(range(p*q)))
    return perm.sign()


# https://en.wikipedia.org/wiki/Zolotarev%27s_lemma
def legendre(q, p):
    assert 0<q<p
    items = list(range(1, p))
    jtems = [(i*q)%p for i in items]
    perm = Perm(dict(zip(items,jtems)), items)
    return perm.sign()


def main():

    #sgn = legendre(3,5)
    #print("sgn =", sgn)
    ps = [p for p in all_primes(100) if p>2]

    #for (p,q) in [(3,5), (3,7)]:
    for p in ps:
      for q in ps:
        i = rowcol(p, q)
        ii = (-1)**( (p-1)*(q-1) // 4 )
        assert i == ii
        

    for p in ps:
      for q in ps:
        if q>=p:
            break
        lhs = legendre(q,p) == 1
        rhs = qresidue(q,p)
        assert lhs == rhs

    


if __name__ == "__main__":

    main()

    print("OK\n")


