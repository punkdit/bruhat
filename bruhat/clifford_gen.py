#!/usr/bin/env python
"""
Build the clifford group on two qubits
from generators and relations ( https://arxiv.org/pdf/1310.6813.pdf ).
It has order 92160.

See also: word.py

"""

from bruhat.argv import argv
from bruhat.todd_coxeter import Schreier


def parse(gens, rels):
    gens = gens.split()
    n = len(gens)
    ngens = 2*n
    get = {gen:i for (i, gen) in enumerate(gens)}
    geti = {gen:(i+n) for (i, gen) in enumerate(gens)}
    rels = rels.split('\n')
    _rels = []
    for gen in gens:
        _rels.append( (get[gen], geti[gen]) )
    for rel in rels:
        rel = rel.strip()
        if not rel:
            continue
        lhs, rhs = rel.split("=")
        lhs = lhs.strip().split()
        rhs = rhs.strip().split()
        if '1' in rhs:
            rhs.remove('1')
        print(lhs, rhs)
        lhs = [get[g] for g in lhs]
        rhs = [geti[g] for g in reversed(rhs)]
        _rels.append( tuple(rhs + lhs) )
    print(_rels)

    graph = Schreier(ngens, _rels)
    #graph.build(maxsize=4000000)
    graph.build()
    print(len(graph))
    return graph


gens = "w HI IH SI IS CZ"
rels = """
w w w w w w w w = 1

w HI = HI w
w IH = IH w
w SI = SI w
w IS = IS w
w CZ = CZ w

HI IH = IH HI
HI IS = IS HI
SI IS = IS SI
SI IH = IH SI

IH IH = 1
IS IS IS IS = 1
IS IH IS IH IS IH  = w

HI HI = 1
SI SI SI SI = 1
SI HI SI HI SI HI  = w

CZ CZ = 1
SI CZ = CZ SI
IS CZ = CZ IS
HI SI SI HI CZ = CZ HI IS SI IS SI HI
IH IS IS IH CZ = CZ SI IH SI IS IS IH
CZ HI CZ = SI HI CZ SI IS HI SI w w w w w w w
CZ IH CZ = IS IH CZ IS SI IH IS w w w w w w w
"""

graph = parse(gens, rels)


print("OK\n\n")

