#!/usr/bin/env python3

import numpy


def mulclose(gen, verbose=False, maxsize=None):
    els = set(gen)
    bdy = list(els)
    changed = True
    while bdy:
        #if verbose:
        #    print "mulclose:", len(els)
        _bdy = []
        for A in gen:
            for B in bdy:
                C = A*B
                if C not in els:
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
        bdy = _bdy
    return els


class Perm(object):

    def __init__(self, perm):
        perm = numpy.array(perm)
        self.perm = perm.copy()

    def __str__(self):
        return "Perm(%s)"%(self.perm,)
    __repr__ = __str__

    def __hash__(self):
        return hash(self.perm.tostring())

    def __eq__(self, other):
        return numpy.alltrue(self.perm == other.perm)

    def __ne__(self, other):
        return not numpy.alltrue(self.perm == other.perm)

    def __lt__(self, other):
        return self.perm.tostring() < other.perm.tostring()

    def __le__(self, other):
        return self.perm.tostring() <= other.perm.tostring()

    def __mul__(self, other):
        perm = self.perm[other.perm]
        return Perm(perm)


class Group(object):
    def __init__(self, perms=None, gen=None):
        if perms is None:
            assert gen is not None
            perms = list(mulclose(gen))
        else:
            perms = list(perms)
        perms.sort()
        #gen = list(gen)
        #self.gen = gen
        self.perms = perms
        self.n = len(self.perms)

    def __str__(self):
        return "Group(order=%s)"%(self.n,)
    __repr__ = __str__

    def __len__(self):
        return self.n

    def __getitem__(self, idx):
        return self.perms[idx]

    def __eq__(self, other):
        return self.perms == other.perms

    def regular_rep(self):
        "the left _regular gset"
        lookup = dict((perm, idx) for (idx, perm) in enumerate(self.perms))
        perms = []
        for idx, p in enumerate(self):
            perm = []
            for jdx, q in enumerate(self):
                r = p*q # left action on self
                kdx = lookup[r]
                perm.append(kdx)
            perm = Perm(perm)
            perms.append(perm)
        tgt = Group(perms)
        send_perms = list(range(self.n))
        hom = Hom(self, tgt, send_perms)
        return hom




class Hom(object):
    def __init__(self, src, tgt, send_perms):
        assert isinstance(src, Group)
        assert isinstance(tgt, Group)
        self.src = src
        self.tgt = tgt
        self.send_perms = send_perms


def test():

    a = Perm([1, 0, 2])
    b = Perm([0, 2, 1])

    G = Group(gen=[a, b])
    assert len(G) == 6
    
    hom = G.regular_rep()
    R = hom.tgt
    assert R != G
    assert R == R.regular_rep().tgt


if __name__ == "__main__":

    test()


