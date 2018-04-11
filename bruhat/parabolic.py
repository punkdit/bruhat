#!/usr/bin/env python2


from __future__ import print_function

import sys, os
import string
from fractions import Fraction
from operator import mul


from bruhat.util import all_subsets
from bruhat.weyl import Weyl
from bruhat.action import Perm, Group, conjugacy_subgroups, burnside

from argv import argv


def parabolic(roots, simple):
    "use _generators from reflections of simple roots"
    for root in simple:
        assert root in roots
    n = len(roots[0])
    idxs = range(n)
    lookup = dict((root, i) for (i, root) in enumerate(roots))
    gen = [] 
    for alpha in simple:
        #print "alpha:", alpha
        r0 = sum(alpha[i]*alpha[i] for i in idxs)
        perm = [] 
        for root in roots:
            #print "    root:", root
            r = sum(alpha[i]*root[i] for i in idxs)
            _root = tuple(root[i] - 2*Fraction(r, r0)*alpha[i] for i in idxs)
            #print "    _root:", _root
            perm.append(lookup[_root])
        perm = Perm(perm, roots)
        gen.append(perm)
    return Group.generate(gen, items=roots)


def main():

    if argv.A_3:
        W = Weyl.build_A(3)
    elif argv.A_4:
        W = Weyl.build_A(4)
    elif argv.B_2:
        W = Weyl.build_B(2)
    elif argv.B_3:
        W = Weyl.build_B(3)
    elif argv.D_4:
        W = Weyl.build_D(4)
    else:
        return

    G = Group.generate(W.gen)
#    for g in W.gen:
#        print(g)
#
#    for root in W.roots:
#        print(root, end=" ")
#        if root in W.simple:
#            print("*")
#        else:
#            print("")

    groups = []
    for subset in all_subsets(len(W.simple)):
        #print("subset:", subset)
        simple = tuple(W.simple[i] for i in subset)
        H = parabolic(W.roots, simple)
        #print("H:", len(H))
        assert G.is_subgroup(H)
        groups.append(H)

    assert len(groups)==2**(len(W.simple))
    print("parabolic subgroups:", len(groups))
    print("orders:", [len(H) for H in groups])

    Hs = conjugacy_subgroups(G, groups)
    print("conjugacy_subgroups:", len(Hs))

    burnside(G, Hs)


if __name__ == "__main__":

    main()



