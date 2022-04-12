#!/usr/bin/env python

"""
Object oriented version of todd_coxeter.py
Somewhat slower...
"""

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
        
    def __str__(self):
        names = self.group.names
        return '*'.join(names[i] for i in self.gens) or "1"

    def __repr__(self):
        return "Word(%s --> %s, [%s])" % (
            self, self.send, ','.join(str(g) for g in self.mul))

    def __lt__(self, other):
        lhs, rhs = self.gens, other.gens
        return len(lhs)<len(rhs) or len(lhs)==len(rhs) and lhs<rhs

    def __mul__(self, other):
        assert isinstance(other, Word)
        #assert len(other.gens) == 1
        for gen in other.gens:    
            word = self.mul[gen]
            if word is None:
                gens = self.gens + (gen,)
                word = Word(self.group, gens)
                self.mul[gen] = word
            self = word
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

    def follow_path(self, g, rel):
        for jdx in reversed(rel.gens):
            gen = self.gens[jdx]
            g = g*gen
            g = g.send
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
            assert idx < 1000000
            if maxsize is not None and idx>=maxsize:
                return False
        return True

    def get(self):
        G = [g for g in self.words if g==g._send]
        return G


def main():

    G = Group(list("a"))
    a, = G.gens
    G.build([a*a])
    assert len(G.get()) == 2

    G = Group(list("ab"))
    a, b = G.gens
    G.build([a*a, b*b, a*b*a*b])
    assert len(G.get()) == 4, len(G.get())

    G = Group("a ai b bi".split())
    a, ai, b, bi = G.gens
    G.build([a*ai, b*bi, a*a, b*b, a*b*a*b*a*b])
    assert len(G.get()) == 6, len(G.get())
    print(len(G.words))

    G = Group("a ai b bi c ci".split())
    a, ai, b, bi, c, ci = G.gens
    G.build([
        a*ai, b*bi, c*ci, 
        a*a, b*b, c*c, 
        (a*b)**3, (a*c)**2, (b*c)**3])
    assert len(G.get()) == 24, len(G.get())
    print(len(G.words))

    G = Group("a ai b bi c ci d di".split())
    a, ai, b, bi, c, ci, d, di = G.gens
    rels = [ (ai* a), (bi* b), (c*ci), (d*di)]
    rels += [ (a*a), (b*b), (c*c), (d*d), (a*c)**2, (a*d)**2, 
        (b*d)**2, (a*b)**3, (b*c)**4, (c*d)**3 ]
    G.build(rels)
    assert len(G.get()) == 1152, len(G.get()) # F_4
    print(len(G.words))




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


