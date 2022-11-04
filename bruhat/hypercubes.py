#!/usr/bin/env python

from bruhat.smap import SMap
from bruhat.argv import argv
from bruhat.util import cross, allperms, distinct


def nbd(v):
    for i, vv in enumerate(v):
        if vv==0:
            u = list(v)
            u[i] = 1
            u = tuple(u)
            yield u

    
def search(n):

    #print("search(%d)"%n)
    verts = list(cross([(0,1)]*n))

    start = (0,)*n
    found = set([start])
    bdy = set([start])
    while bdy:

        edges = []
        _bdy = set()
        for v in bdy:
          for u in nbd(v):
            if u not in found:
                _bdy.add(u)
                edges.append((v, u))

        if edges:
            print(len(edges), end=" ")
        
        bdy = _bdy
        found.update(bdy)

    print()



if __name__ == "__main__":

    for n in range(8):
        search(n)


    print("OK\n")



