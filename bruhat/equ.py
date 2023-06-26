#!/usr/bin/env python3

class Equ(object):
    "Equivalence class for an _equivalance relation"
    def __init__(self, item):
        self.items = [item]
        self.parent = None 

    def __str__(self):
        return "Equ(%s, top=%s)" % (
            self.items, self.top if self.parent else None)

    @property
    def top(self):
        top = self
        while top.parent:
            top = top.parent
        if top is not self:
            self.parent = top
        return top

#    def add(self, item):
#        assert item not in self.items
#        self.items.append(item)

    def merge(self, other):
        top = self
        while top.parent:
            top = top.parent
        while other.parent:
            other = other.parent
        if top is other:
            return top
        other.parent = top
        top.items += other.items
        return top

    def eq(e1, e2):
        while e1.parent:
            e1 = e1.parent
        while e2.parent:
            e2 = e2.parent
        return e1 is e2


def quotient(items, relation=None, verbose=False):
    """ return a dict:item->items that sends each item 
        to the list of equivalant items 
        under the equivalance relation.
    """
    if relation is None:
        relation = lambda a,b : (a==b)
    equs = [Equ(item) for item in items]
    n = len(items)
    for i in range(n):
        ei = equs[i]
        if verbose:
            print(".", end="", flush=True)
        for j in range(i+1, n):
            ej = equs[j]
            if ei.eq(ej):
                continue
            if relation(items[i], items[j]):
                ei.merge(ej)
        if verbose:
            found = set(equ.top for equ in equs)
            print("[%d]"%(len(found)))
    try:
        hom = {}
        for i in range(n):
            hom[items[i]] = list(equs[i].top.items)
        if verbose:
            print()
        return hom
    except TypeError:
        return equs


def quotient_rep(items, relation=None):
    """ return a dict:item->item that sends each item 
        to a representative equivalant item
        under the equivalance relation.
    """
    hom = quotient(items, relation)
    hom = dict((item, items[0]) for (item, items) in hom.items()) # _confused ?
    return hom


assert quotient(range(7), (lambda i,j : (i-j)%3==0)) \
    == {0: [0, 3, 6], 1: [1, 4], 2: [2, 5], 3: [0, 3, 6], 4: [1, 4], 5: [2, 5], 6: [0, 3, 6]}

assert quotient_rep(range(7), (lambda i,j : (i-j)%3==0)) \
    == {0:0, 1:1, 2:2, 3:0, 4:1, 5:2, 6:0}


