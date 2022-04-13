#!/usr/bin/env python

"""
Object oriented version of todd_coxeter.py
Somewhat slower...

See also: https://jdbm.me/posts/2021-07-11-todd-coxeter/

"""

import heapq

from bruhat.argv import argv

#class Gen(object):
#    def __init__(self, group, idx):
#        self.group = group
#        self.idx = idx
#        self.lookup = [None]*group.ngens
#        self.send = self
#
#    def __mul__(self, other):
        

class Word(object):
    def __init__(self, group, gens):
        self.group = group
        self.gens = gens
        self._send = self
        self.mul = [None]*group.ngens
        group.append(self)

    @property
    def send(self):
        item = self._send
        while item._send != item:
            item = item._send
        self._send = item
        return item

    def dosend(self, other):
        #print("dosend: %s --> %s" % (self, other))
        assert self._send == self
        assert other._send == other
        self._send = other
        #print("dosend", repr(self))

    def __getitem__(self, idx):
        group = self.group
        return group.names[self.gens[idx]]

    def __hash__(self):
        return id(self)
        
    def __str__(self):
        names = self.group.names
        return '*'.join(names[i] for i in self.gens) or "1"

    def __repr__(self):
        return "Word(%s --> %s, [%s])" % (
            self, self.send, ','.join(str(g) for g in self.mul))

    def __lt__(self, other):
        return self.group.lt_word(self, other)

    def __le__(self, other):
        return self==other or self < other

    def __mul__(self, other):
        assert isinstance(other, Word)
        #assert len(other.gens) == 1
        for gen in other.gens:    
            word = self.mul[gen]
            if word is None:
                gens = self.gens + (gen,)
                word = Word(self.group, gens)
                self.mul[gen] = word
            self = word.send
        return self

    def __pow__(self, n):
        if n==0:
            return self.group.I
        assert n>0
        op = self
        while n > 1:
            op = self*op
            n -= 1
        return op


class Group(object):
    def __init__(self, names):
        self.names = names
        self.ngens = len(names)
        self.words = []
        self.I = Word(self, ())
        self.gens = [Word(self, (idx,)) for idx in range(self.ngens)]
        for idx, g in enumerate(self.gens):
            self.I.mul[idx] = g

    def append(self, word):
        #for w in self.words:
        #    assert str(w) != str(word)
        self.words.append(word)

    def lt_word(self, lhs, rhs):
        lhs, rhs = lhs.gens, rhs.gens
        return len(lhs)<len(rhs) or len(lhs)==len(rhs) and lhs<rhs

    def follow_path(self, g, rel):
        for jdx in reversed(rel.gens):
            gen = self.gens[jdx]
            g = g*gen
        return g

    def unify(self, g, h):
        items = [(g, h)]
        while items:
            g, h = items.pop()
            #print("unify", g, h)
            g, h = g.send, h.send
            if g==h:
                continue
            if h < g:
                g, h = h, g
            h.dosend(g)
            for idx in range(self.ngens):
                gn = g.mul[idx]
                hn = h.mul[idx]
                if gn is None:
                    g.mul[idx] = hn
                elif hn is not None:
                    items.append((gn, hn))

    def build(self, rels, maxsize=None):
        words = self.words
        #print(" ".join(str(g) for g in words))
        idx = 0
        while idx < len(words):
            #print("idx =", idx)
            #print(self.words)
            g = words[idx]
            if g.send == g:
                for rel in rels:
                    h = self.follow_path(g, rel)
                    self.unify(h, g)
                    #assert g.send == h.send
            idx += 1
            if maxsize is not None and idx>=maxsize:
                return False
            else:
                assert idx < 10000000
        return True

    def get(self):
        G = [g for g in self.words if g==g._send]
        return G



class PriorityGroup(Group):
    def __init__(self, names, lt_word=None):
        Group.__init__(self, names)
        if lt_word is not None:
            self.lt_word = lt_word

    def lt_word(self, lhs, rhs):
        lhs, rhs = lhs.gens, rhs.gens
        return len(lhs)<len(rhs) or len(lhs)==len(rhs) and lhs<rhs

    def append(self, word):
        words = self.words
        heapq.heappush(words, word)

    def build(self, rels, maxsize=None, maxcount=None):
        I = self.I
        for rel in rels:
            self.unify(I, rel)
        words = self.words
        remain = []
        count = 0
        while words:
            g = heapq.heappop(words)
            count += 1
            if g._send == g:
                for rel in rels:
                    h = self.follow_path(g, rel)
                    self.unify(h, g)
                if g._send == g:
                    remain.append(g)
                if maxsize and len(remain)>=maxsize:
                    break
            #if count % 10000 == 0:
            #    print("[%d] -->" % len(remain), end="")
            #    remain = [g for g in remain if g._send==g]
            #    print("[%d]"%len(remain))
            #assert len(words) < 100000
            #if maxsize and len(words)>=maxsize:
            #    break
            if maxcount and maxcount <= count:
                break
        #else:
        #    print("no words")
        heapq.heapify(remain) # ?
        self.words = remain
        print("PriorityGroup.build: count =", count)

    def filter(self, accept):
        words = [g for g in self.words if accept(g)]
        heapq.heapify(words)
        self.words = words


def test(builder):

    G = builder(list("a"))
    a, = G.gens
    G.build([a*a])
    assert len(G.get()) == 2, len(G.get())

    G = builder(list("ab"))
    a, b = G.gens
    G.build([a*a, b*b, a*b*a*b])
    assert len(G.get()) == 4, len(G.get())

    G = builder("a ai b bi".split())
    a, ai, b, bi = G.gens
    G.build([a*ai, b*bi, a*a, b*b, a*b*a*b*a*b])
    assert len(G.get()) == 6, len(G.get())

    G = builder("a ai b bi c ci".split())
    a, ai, b, bi, c, ci = G.gens
    G.build([
        a*ai, b*bi, c*ci, 
        a*a, b*b, c*c, 
        (a*b)**3, (a*c)**2, (b*c)**3])
    assert len(G.get()) == 24, len(G.get())

    G = builder("a ai b bi c ci d di".split())
    a, ai, b, bi, c, ci, d, di = G.gens
    rels = [ 
        (ai* a), (bi* b), (c*ci), (d*di), 
        (a*a), (b*b), (c*c), (d*d), (a*c)**2, (a*d)**2, 
        (b*d)**2, (a*b)**3, (b*c)**4, (c*d)**3 ]
    G.build(rels)
    assert len(G.get()) == 1152, len(G.get()) # F_4

    if argv.slow:
        G = builder("a ai b bi c ci d di".split())
        a, ai, b, bi, c, ci, d, di = G.gens
        rels = [ 
            (ai* a), (bi* b), (c*ci), (d*di),
            (a*a), (b*b), (c*c), (d*d), 
            (a*c)**2, (a*d)**2, (b*d)**2, (a*b)**5, (b*c)**3, (c*d)**3 ]
        G.build(rels)
        assert len(G.get()) == 14400 # H_4

    G = builder("s si t ti".split())
    s, si, t, ti = G.gens
    rels = [(si* s), (s* s), (ti* t), (s*t)**3]
    G.build(rels, maxsize=10000)
    print(len(G.get()))


def main():
    #test(Group)
    test(PriorityGroup)




if __name__ == "__main__":

    profile = argv.profile
    name = argv.next() or "main"

    if profile:
        import cProfile as profile
        profile.run("%s()"%name)
    else:
        fn = eval(name)
        fn()

    print("OK")


