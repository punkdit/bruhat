#!/usr/bin/env python3

from bruhat.species import all_subsets
from bruhat.action import all_functions
from bruhat import equ


def closure(pairs, els=None):
    homs = set((a,b) for (a,b) in pairs) # set of (a,b) where a<=b
    ident = set()
    changed = True
    while changed:
        changed = False
        pairs = list(homs)
        for (a,b) in pairs:
          for (c,d) in pairs:
            if b!=c:
                continue
            if (a,d) not in homs:
                homs.add((a,d))
                changed = True
            if a!=b and a==d:
                ident.add((a,b))
            #if a!=b:
            #    assert a!=d, "loop found: %s = %s" % (a, b)

    # At this point we have a PreOrder, now we
    # identify elements to form a Poset.
    relation = lambda a,b : (a,b) in ident
    els = list(set([a for (a,b) in homs]+[b for (a,b) in homs]))
    send = equ.quotient_rep(els, relation)
    get = send.get
    homs = set((get(a, a), get(b, b)) for (a, b) in homs)

    els = set([a for (a,b) in homs]+[b for (a,b) in homs])
    for a in els:
        homs.add((a,a))
    return homs, els


class Poset(object):
    def __init__(self, pairs, els=None, check=True):
        pairs, els = closure(pairs, els)
        items = list(pairs)
        items.sort() # canonical
        els = list(els)
        els.sort() # canonical

        ups = {}
        dns = {}
        for (a, b) in pairs:
            up = ups.get(a, [])
            up.append(b)
            ups[a] = up
            dn = dns.get(b, [])
            dn.append(a)
            dns[b] = dn

        self.items = items
        self.pairs = pairs # set(items)
        self.els = els
        self.ups = ups
        self.dns = dns
        self._sup = {} # cache
        if check:
            self.check()

    def check(self):
        pairs = self.pairs
        ups = self.ups
        dns = self.dns
        for (a, b) in pairs:
            assert a in ups
            assert b in ups
            assert b in ups[a]
            assert a in dns[b]
            for (c, d) in pairs:
                if b!=c:
                    continue
                assert (a, d) in pairs
                assert a==b or a!=d, self # irreflexive

    def add_pairs(self, pairs):
        "add relations to form a new Poset"
        ups = self.ups
        for (a, b) in pairs:
            assert a in ups
            assert b in ups
        pairs = self.items + pairs
        return Poset(pairs)

    def __str__(self):
        return "%s(%s, %s)"%(self.__class__.__name__, self.els, self.pairs)
    __repr__ = __str__

    def __eq__(self, other):
        return self.items == other.items

    def __ne__(self, other):
        return self.items != other.items

    def __hash__(self):
        return hash(tuple(self.items))

    def __getitem__(self, i):
        return self.els[i]

    def __contains__(self, el):
        return el in self.ups

    def __len__(self):
        return len(self.els)

    def min(self, items):
        if not items:
            return [] # <----- return
        #print("min", self)
        pairs = self.pairs
        up = set(items)
        while len(up) > 1:
            items = list(up)
            #print("\t", items)
            n = len(up)
            for a1 in items:
              for b1 in items:
                if a1==b1:
                    continue
                if (a1, b1) in pairs and b1 in up:
                    up.remove(b1)
            #print("\t", up)
            if len(up) == n:
                break
        assert up
        return up

    def max(self, items):
        if not items:
            return [] # <----- return
        pairs = self.pairs
        dn = set(items)
        while len(dn) > 1:
            items = list(dn)
            n = len(dn)
            for a1 in items:
              for b1 in items:
                if a1==b1:
                    continue
                if (b1, a1) in pairs and b1 in dn:
                    dn.remove(b1)
            if len(dn) == n:
                break
        assert dn
        return dn

    def uniq_min(self, items):
        up = self.min(items)
        if len(up)==1:
            return list(up)[0]

    def uniq_max(self, items):
        dn = self.max(items)
        if len(dn)==1:
            return list(dn)[0]

    _top = None
    @property
    def top(self):
        if self._top is None:
            self._top = self.uniq_max(self.els)
        return self._top

    _bot = None
    @property
    def bot(self):
        if self._bot is None:
            self._bot = self.uniq_min(self.els)
        return self._bot

    def sup(self, a, b): 
        "return colimit of a & b, if it exists, otherwise None"
        if (a, b) in self._sup:
            return self._sup[a, b]
        upa = set(self.ups[a])
        upb = set(self.ups[b])
        up = upa.intersection(upb)
        value = self.uniq_min(up)
        self._sup[a, b] = value
        return value

    def get_id(self):
        "build _identity Hom"
        send = dict((a, a) for a in self.els)
        return Hom(self, self, send)

    def hom_iter(self, other):
        left = self.els
        right = other.els
        for f in all_functions(left, right):
            for (a,b) in self.pairs:
                fa, fb = f[a], f[b]
                if (fa, fb) not in other.pairs:
                    break
            else:
                yield Hom(self, other, f)

    def hom(self, other):
        "poset of hom's"
        hom = list(self.hom_iter(other))
        pairs = [(f,g) for f in hom for g in hom if f<=g]
        return Poset(pairs)

    def tensor(L, R):
        gen = [(l, r) for l in L for r in R if l!=L.bot and r!=R.bot]
        pairs = []
        for l in L:
            for (r1, r2) in R.pairs:
                pairs.append(((l,r1), (l,r2)))
        for r in R:
            for (l1, l2) in L.pairs:
                pairs.append(((l1, r), (l2, r)))

        F = Poset(pairs, els)
        assert F.bot == (L.bot, R.bot)
        for (l1, l2) in L.pairs:
          for r in R:
            lsup = L.sup(l1, l2)
            assert F.sup((l1, r), (l2, r)) == (lsup, r)
        for l in L:
          for (r1, r2) in R.pairs:
            rsup = R.sup(r1, r2)
            assert F.sup((l, r1), (l, r2)) == (l, rsup)
        return F


    __matmul__ = tensor


class SupPoset(Poset):

    def __init__(self, pairs, els=None, check=True):
        Poset.__init__(self, pairs, els, check=check)
        if check:
            self.check()

    def check(self):
        Poset.check(self)
        assert self.bot is not None
        for a in self.els:
          for b in self.els:
            assert self.sup(a, b) is not None

    @classmethod
    def promote(cls, item, check=True):
        if isinstance(item, SupPoset):
            return item
        if not isinstance(item, Poset):
            raise TypeError
        return cls(item.pairs, check=check)

    @classmethod
    def free(cls, gens, check=True):
        items = list(all_subsets(gens))
        pairs = [(a, b) for a in items for b in items if set(a).issubset(set(b))]
        try:
            # try to simplify this...
            assert "+" not in gens
            pairs = [('+'.join(a), '+'.join(b)) for (a,b) in pairs]
        except:
            pass
        return cls(pairs, check=check)

    def add_pairs(self, pairs, check=True):
        "add relations to form a new SupPoset"
        ups = self.ups
        for (a, b) in pairs:
            assert a in ups
            assert b in ups
        pairs = set(self.items + pairs)

        while 1:
            self = Poset(pairs, check=check)
            assert self.bot is not None
            ident = set()
            for a in self:
              for b in self:
                upa = set(self.ups[a])
                upb = set(self.ups[b])
                up = upa.intersection(upb)
                assert len(up)
                up = self.min(up)
                assert len(up)
                if len(up) == 1:
                    continue
                for a1 in up:
                  for b1 in up:
                    ident.add((a1, b1))
            print("ident:", ident)
            if not ident:
                return SupPoset.promote(self)
            pairs = self.items + list(ident)

    def hom(self, other):
        "poset of hom's that preserve sup's, including empty sup == bot"
        els = self.els
        pairs = [(a, b) for a in els for b in els]
        assert self.bot is not None, "self is not sup-complete"
        assert other.bot is not None, "other is not sup-complete"
        hom = []
        for f in self.hom_iter(other):
            #print(f, "?")
            if f[self.bot] != other.bot:
                #print("\t fail 1")
                continue
            for (a, b) in pairs:
                c = self.sup(a, b)
                assert c is not None, "sup(%s, %s) == None" % (a, b)
                fa, fb = f[a], f[b]
                d = other.sup(fa, fb)
                if d!=f[c]:
                    #print("\t fail 2", a, b, c, fa, fb, d)
                    break
            else:
                #print("\tappend")
                hom.append(f)
        pairs = [(f,g) for f in hom for g in hom if f<=g]
        return Poset(pairs)


class Hom(object):
    "Morphism between Poset's"
    def __init__(self, src, tgt, send, check=True):
        assert isinstance(src, Poset)
        assert isinstance(tgt, Poset)
        send = dict(send)
        els = set(send.keys())
        items = list(send.items())
        items.sort()
        self.items = tuple(items)
        assert len(els) == len(src), (els, len(src))
        for k,v in send.items():
            assert k in src
            assert v in tgt
        self.src = src
        self.tgt = tgt
        self.send = send
        if check:
            self.check()

    def check(self):
        src = self.src
        tgt = self.tgt
        send = self.send
        for a,b in src.pairs:
            p = send[a], send[b]
            assert p in tgt.pairs

    def __str__(self):
        return "Hom(%s)"%(self.items,)
    __repr__ = __str__

    def __mul__(self, other):
        # other then self
        assert other.tgt == self.src
        pairs = {}
        for (k,v) in other.send.items():
            v = self.send[v]
            pairs[k] = v
        return Hom(other.src, self.tgt, pairs)

    def __getitem__(self, k):
        return self.send[k]

    def __eq__(self, other):
        assert self.src == other.src
        assert self.tgt == other.tgt
        return self.send == other.send

    def __ne__(self, other):
        assert self.src == other.src
        assert self.tgt == other.tgt
        return self.send != other.send
    
    _hash = None
    def __hash__(self):
        if self._hash is None:
            self._hash = hash(self.items)
        return self._hash

    def __le__(self, other):
        assert self.src == other.src
        assert self.tgt == other.tgt
        for a in self.src.els:
            lhs = self.send[a]
            rhs = other.send[a]
            if lhs > rhs:
                return False
        return True
    
    def __lt__(self, other):
        assert self.src == other.src
        assert self.tgt == other.tgt
        return not self==other and self<=other
    
        

def main():

    I = Poset(['01']) # 0 <= 1
    #print(I)

    P = Poset(['ab', 'ba'])
    assert len(P) == 1

    P = Poset(['ab', 'bc', 'ca'])
    assert len(P) == 1, P

#    ok = False
#    try:
#        P = Poset(['ab', 'bc', 'ca'])
#    except AssertionError:
#        ok = True
#    assert ok
    

    P = Poset('aa bb'.split())
    assert P.sup('a', 'b') == None

    P = Poset('a1 b1 a0 b0'.split())
    assert P.sup('a', 'b') == None

    P = Poset('a1 b1 a0 b0 12 02'.split())
    assert P.sup('a', 'b') == None

    P = Poset('a1 b1 a0 b0 a2 b2 20 21'.split())
    assert P.sup('a', 'b') == '2'

    P = Poset('0a 0b 0c a1 b1 c1'.split())
    assert len(P) == 5

    assert P.sup('a', 'b') == '1'
    assert P.sup('a', '0') == 'a'
    assert P.sup('a', '1') == '1'

    i = P.get_id()
    assert i*i == i

    F = SupPoset.free([])
    assert len(F)==1
    assert F.bot is not None

    F = SupPoset.free('ABC')
    assert len(F) == 8
    assert F.sup("A", "A") == "A"
    assert F.sup("", "A") == "A"
    assert F.sup("", "A+B") == "A+B"
    assert ("A", "A+B") in F.pairs
    assert F.sup("A", "A+B") == "A+B"

    if 1:
        sup = F.sup
        G = F.add_pairs([
            ("A", sup("B", "C")),
            ("B", sup("A", "C")),
            ("C", sup("B", "A")),
        ])
        print(G)
    
        return

    assert len(list(f for f in I.hom_iter(I))) == 3

    P = Poset('0a 0b 0c a1 b1 c1'.split())
    for f in P.hom_iter(I):
        f.check()

    Q = P.hom(I)
    for f in Q:
        f.check()
    assert len(Q) == 10
    assert len(list(P.hom_iter(P))) == 178

    P = SupPoset.promote(P)
    I = SupPoset.promote(I)

    Q = P.hom(I)
    for f in Q:
        f.check()
    assert len(Q) == 5

    #assert len(P.hom(P))==50 # slow...

#    print("II =", len(I@I))
#
#    PP = P@P
#    print(PP.els)
#    print(len(P@P))

    F = SupPoset.free("ABC")
    assert len(F.hom(I)) == len(F)

    print("OK")


if __name__ == "__main__":

    main()

