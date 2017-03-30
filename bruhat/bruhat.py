#!/usr/bin/env python

import sys, os

from argv import argv

"""
Coxeter groups:
http://math.ucr.edu/home/baez/week62.html
"""


def pairs(gen):
    for s in gen:
        for t in gen:
            yield (s, t)

def distinct_pairs(gen):
    for s in gen:
        for t in gen:
            if s != t:
                yield (s, t)

def rewrite(word, src, tgt): # 40% hotspot
    if src not in word:
        return
    assert src != tgt
    stem = ''
    while src in word:
        idx = word.index(src)
        n = len(src)
        yield stem + word[:idx] + tgt + word[idx+n:]
        stem, word = word[:idx+1], word[idx+1:]

assert list(rewrite("zxxxw", 'xx', 'y')) == ['zyxw', 'zxyw']


def cancel_pairs(word):
    for i in range(len(word)-1):
        if word[i]==word[i+1]:
            yield word[:i] + word[i+2:]

def idempotent_pairs(word):
    for i in range(len(word)-1):
        if word[i]==word[i+1]:
            yield word[:i] + word[i+1:]

assert list(cancel_pairs('xyzzz')) == ['xyz', 'xyz']

class Coxeter(object):
    """ Coxeter group
    """

    def __init__(self, gen, rel, identity='', bruhat=False, build=False):
        """
        gen : list of generators
        rel : map pairs (i,j) of generators to m_{ij} (this is the Coxeter matrix)
        """
        for i in gen:
          for j in gen:
            if rel.get((i, j)) and not rel.get((j, i)):
                rel[j, i] = rel[i, j]
        for i in gen:
          for j in gen:
            m = rel.get((i, j))
            if m is None:
                rel[i, j] = 2 # default 
        for i in gen:
          for j in gen:
            assert rel[i, j] == rel[j, i]
            #assert rel[i, j] in (2, 3, 4, 6)
        self.gen = gen
        self.rel = rel

        reduced = {'':('',)} # map word -> sorted tuple of equivalent reduced words
        lookup = {('',):''} # map sorted tuple -> word
        for g in gen:
            reduced[g] = (g,)
            lookup[(g,)] = g
            if bruhat:
                reduced[g+g] = reduced[g] # Bruhat monoid
            else:
                reduced[g+g] = reduced[''] # Coxeter group
            for h in gen:
                if g==h:
                    continue
                m = rel[g, h]
                ghm = ''.join(([g, h]*m)[:m])
                hgm = ''.join(([h, g]*m)[:m])
                r = [ghm, hgm]
                r.sort()
                r = tuple(r)
                reduced[ghm] = reduced[hgm] = r
                lookup[r] = ghm
        self.reduced = reduced
        self.lookup = lookup
        self.bruhat = bruhat
        if build:
            self.build()

    def get_canonical(self, orig): # 30% hotspot
        reduced = self.reduced
        if orig in reduced:
            return reduced[orig]
        gen = self.gen
        rel = self.rel
        items = set([orig])
        done = False
        #print "E"
        pair_iter = idempotent_pairs if self.bruhat else cancel_pairs
        while not done:
            done = True
            for word in list(items):
                for word1 in pair_iter(word):
                    items = set([word1])
                    assert len(word1)<len(word)
                    #print "X"
                    done = False
                    break
                else:
                    continue
                break

            for s, t in distinct_pairs(gen):
                m = rel[s, t]
                stm = ''.join(([s, t]*m)[:m])
                tsm = ''.join(([t, s]*m)[:m])
                for word in list(items):
                    if stm in word:
                        for word1 in rewrite(word, stm, tsm):
                            if word1 not in items:
                                items.add(word1)
                                #print "Y"
                                done = False
        #print "F"
        items = list(items)
        items.sort()
        items = tuple(items)
        reduced[orig] = items
        lookup = self.lookup
        if items not in lookup:
            lookup[items] = orig
        elif len(lookup[items]) > orig:
            lookup[items] = orig
        return items

    def is_equal(self, lhs, rhs):
        if lhs==rhs:
            return True
        left = self.get_canonical(lhs)
        right = self.get_canonical(rhs)
        return left == right

    def build(self, max_size=None):
        lookup = self.lookup # map canonical -> word
        words = set(lookup.keys()) # canonical forms
        #print "words:", words
        pairs = [(i, j) for i in words for j in words]
        mul = {} # map (i,j) -> i*j
        while 1:
            newpairs = []
            #print "pairs:", pairs
            for i, j in pairs:
                g = lookup[i]
                h = lookup[j]
                gh = g+h # multiply words
                k = self.get_canonical(gh) # updates lookup aswell
                gh = lookup[k]
                mul[g, h] = gh
                if k not in words:
                    newpairs += [(k, g) for g in words]
                    newpairs += [(g, k) for g in words]
                    newpairs += [(k, k) for g in words]
                    words.add(k)
                    if max_size is not None and len(words)>=max_size:
                        break
            if not newpairs:
                break
            pairs = newpairs
        self.mul = mul
        self.words = lookup.values()
        self.check()
        return self.words

    def check(self):
        words, mul = self.words, self.mul
        for i in words:
          for j in words:
            k = mul.get((i, j))
            assert k is not None
            assert k in words


class BruhatMonoid(Coxeter):
    def __init__(self, *args, **kw):
        kw["bruhat"] = True
        Coxeter.__init__(self, *args, **kw)


def main():

    bruhat = argv.bruhat

    # |A_n| = (n+1)!
    A_2 = Coxeter("LP", {("L", "P") : 3}, bruhat=bruhat)
    words = A_2.build()
    assert len(words)==6

    # |I_n| = 2*n
    I_5 = Coxeter("LP", {("L", "P") : 5})
    assert len(I_5.build())==10
    # B_2 : 8

    # |B_n| = 2^n n!
    B_3 = Coxeter("ABC", {("A", "B"):3, ("B", "C"):4})
    assert len(B_3.build())==48

    G_2 = Coxeter("LP", {("L", "P") : 6})
    assert len(G_2.build())==12

    #H_3 = Coxeter("ABC", {("A", "B") : 5, ("B", "C") : 3})
    #assert len(H_3.build())==120

    I_4 = Coxeter("XZ", {("X", "Z") : 4})
    assert len(I_4.build()) == 8

    I_44 = Coxeter("UVXZ", {
        ("U", "X"):2, ("U", "Z"):2, ("V", "X"):2, ("V", "Z"):2, 
        ("U", "V"):4, ("X", "Z") : 4})
    #print len(I_44.build()) # 64 ???

    if argv.F_4:
        F_4 = Coxeter("ABCD", {("A", "B"):3, ("A", "C"):4, ("A", "D"):3})
        g = F_4.build()
        print len(g)
        assert len(g)==1152

    if argv.D_4:
        D_4 = Coxeter("ABCD", {("A", "B"):3, ("A", "C"):3, ("A", "D"):3})
        g = D_4.build()
        assert len(g)==192

    if argv.A_3:
        A_3 = Coxeter("LPS", {("L", "P") : 3, ("L", "S"):3}, bruhat=bruhat)
        assert len(A_3.build())==24
        #print A_3.lookup
        #print A_3.reduced

    if argv.A_4:
        A_4 = Coxeter("LPSH", {("L", "P") : 3, ("L", "S"):3, ("S", "H"):3})
        assert len(A_4.build())==120

    print "OK"



if __name__ == "__main__":

     if argv.profile:
        import cProfile as profile
        profile.run("main()")

     else:
        main()


