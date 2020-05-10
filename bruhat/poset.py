#!/usr/bin/env python3


from bruhat.species import all_subsets
from bruhat.action import all_functions
from bruhat import equ
from bruhat.argv import argv

#CHECK = True # much slower...
CHECK = False


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
    els = set([a for (a,b) in homs]+[b for (a,b) in homs])
    for a in els:
        homs.add((a,a))
    return homs, els



class PreOrder(object):
    def __init__(self, els, pairs, check=CHECK):
        if els is None:
            els = set([a for (a,b) in pairs]+[b for (a,b) in pairs])
        els = list(els)
        els.sort() # canonical
        items = list(pairs)
        items.sort() # canonical

        ups = dict((a, set()) for a in els)
        dns = dict((a, set()) for a in els)
        for (a, b) in pairs:
            ups[a].add(b)
            dns[b].add(a)

        self.items = items
        self.pairs = pairs # set(items)
        self.els = els
        self.ups = ups
        self.dns = dns
        self._sup = {} # cache
        self._inf = {} # cache
        if check:
            self.check()

    @classmethod
    def generate(cls, pairs, els=None, check=CHECK):
        pairs, els = closure(pairs, els)
        P = cls(els, pairs, check)
        return P

    def add_pairs(self, pairs, check=CHECK):
        pairs = list(pairs) + self.items
        P = PreOrder.generate(pairs, check=check)
        send = dict((a,a) for a in self.els)
        hom = Hom(self, P, send)
        return hom

    def check(self):
        pairs = self.pairs
        els = self.els
        for a in els:
            assert (a, a) in pairs, self
        ups = self.ups
        dns = self.dns
        for (a, b) in pairs:
            assert a in ups, self
            assert b in ups, self
            assert b in ups[a], self
            assert a in dns[b], self
            for (c, d) in pairs:
                if b!=c:
                    continue
                assert (a, d) in pairs, self

    def __str__(self):
        return "%s(%s, %s)"%(self.__class__.__name__, self.els, self.pairs)
    __repr__ = __str__

    def __eq__(self, other):
        return self.items == other.items

    def __ne__(self, other):
        return self.items != other.items

    def __hash__(self):
        return hash(tuple(self.items)) # XX cache me

    def __getitem__(self, i):
        return self.els[i]

    def __contains__(self, el):
        return el in self.ups

    def __len__(self):
        return len(self.els)


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
                #yield Hom(self, other, f)
                yield f

    def hom(self, other):
        "poset of hom's"
        hom = [Hom(self, other, f) for f in self.hom_iter(other)]
        pairs = [(f,g) for f in hom for g in hom if f<=g]
        return Poset(hom, pairs)


class Poset(PreOrder):
    def check(self):
        PreOrder.check(self)
        pairs = self.pairs
        for (a, b) in pairs:
            assert (b, a) not in pairs or a==b, (a, b)

    @classmethod
    def promote(cls, src, check=CHECK):
        if isinstance(src, Poset):
            return src.get_id()
        if not isinstance(src, PreOrder):
            raise TypeError

        # At this point we have a PreOrder, now we
        # identify elements to form a Poset.
        #print("promote:", src)
        ident = set()
        pairs = src.pairs
        for (a, b) in pairs:
            if a!=b and (b, a) in pairs:
                ident.add((a, b)) # identify these two
        relation = lambda a,b : (a,b) in ident
        els = list(set([a for (a,b) in pairs]+[b for (a,b) in pairs]))
        send = equ.quotient_rep(els, relation)
        #print("send:", send)
        get = send.get
        els = set(get(a, a) for a in els)
        pairs = set((get(a, a), get(b, b)) for (a, b) in pairs)
        tgt = Poset(els, pairs, check=check)
        hom = Hom(src, tgt, send)
        return hom

    @classmethod
    def generate(cls, pairs):
        P = PreOrder.generate(pairs)
        hom = cls.promote(P)
        return hom

    def add_pairs(self, pairs, check=CHECK):
        "add relations to form a new Poset"
        ups = self.ups
        for (a, b) in pairs:
            assert a in ups
            assert b in ups
        #pairs = self.items + pairs
        hom = PreOrder.add_pairs(self, pairs, check=check)
        P = hom.tgt
        hom = Poset.promote(P) * hom
        return hom

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

    def sup2(self, a, b): 
        "return colimit of a & b, if it exists, otherwise None"
        if (a, b) in self._sup:
            return self._sup[a, b]
        upa = set(self.ups[a])
        upb = set(self.ups[b])
        up = upa.intersection(upb)
        value = self.uniq_min(up)
        self._sup[a, b] = value
        return value

    def inf2(self, a, b): 
        "return limit of a & b, if it exists, otherwise None"
        if (a, b) in self._inf:
            return self._inf[a, b]
        dna = set(self.dns[a])
        dnb = set(self.dns[b])
        dn = dna.intersection(dnb)
        value = self.uniq_max(dn)
        self._inf[a, b] = value
        return value

    def sup(self, items):
        items = list(items)
        while len(items)>1:
            a = self.sup2(items.pop(), items.pop())
            items.append(a)
        if items:
            return items[0]
        return self.bot

    def inf(self, items):
        items = list(items)
        while len(items)>1:
            a = self.inf2(items.pop(), items.pop())
            items.append(a)
        if items:
            return items[0]
        return self.top

#    def tensor(L, R): # XXX
#        gen = [(l, r) for l in L for r in R if l!=L.bot and r!=R.bot]
#        pairs = []
#        for l in L:
#            for (r1, r2) in R.pairs:
#                pairs.append(((l,r1), (l,r2)))
#        for r in R:
#            for (l1, l2) in L.pairs:
#                pairs.append(((l1, r), (l2, r)))
#
#        F = Poset(pairs, els)
#        assert F.bot == (L.bot, R.bot)
#        for (l1, l2) in L.pairs:
#          for r in R:
#            lsup = L.sup2(l1, l2)
#            assert F.sup2((l1, r), (l2, r)) == (lsup, r)
#        for l in L:
#          for (r1, r2) in R.pairs:
#            rsup = R.sup2(r1, r2)
#            assert F.sup2((l, r1), (l, r2)) == (l, rsup)
#        return F
#    __matmul__ = tensor


class SupPoset(Poset):

    def __init__(self, els, pairs, check=CHECK):
        Poset.__init__(self, els, pairs, check=check)
        if check:
            self.check()

    def check(self):
        Poset.check(self)
        assert self.bot is not None
        for a in self.els:
          for b in self.els:
            if self.sup2(a, b) is None:
                print("SupPoset.check:")
                print(a)
                for up in self.ups[a]:
                    print("\t", up)
                print(b)
                for up in self.ups[b]:
                    print("\t", up)
                raise AssertionError

    @classmethod
    def free(cls, gens, check=CHECK):
        items = list(all_subsets(gens))
        pairs = [(a, b) for a in items for b in items if set(a).issubset(set(b))]
        try:
            # try to simplify this...
            assert "+" not in gens
            pairs = [('+'.join(a), '+'.join(b)) for (a,b) in pairs]
        except:
            pass
        return cls(None, pairs, check=check)

    @classmethod
    def promote(cls, item, check=CHECK):
        if isinstance(item, SupPoset):
            return item
        if not isinstance(item, PreOrder):
            raise TypeError(item)
        hom = Poset.promote(item, check)
        tgt = hom.tgt
        P = cls(tgt.els, tgt.pairs, check=check)
        hom = Hom(hom.src, P, hom.send, check=check)
        return hom

    def add_pairs(self, rels, check=True):
        "add relations to form a new SupPoset"
        ups = self.ups
        for (a, b) in rels:
            assert a in ups
            assert b in ups

        pairs = self.pairs
        els = []
        for x in self.els:
            for (lhs, rhs) in rels: # lhs <= rhs
                if (rhs,x) in pairs and not (lhs,x) in pairs:
                    break
            else:
                els.append(x)
        els = set(els)
        pairs = [(a, b) for (a, b) in self.pairs if a in els and b in els]
        P = SupPoset(els, pairs)
        # inclusion is a right adjoint
        radj = Hom(P, self, dict((a, a) for a in els))
        ladj = radj.get_ladj()
        return ladj

    def hom(self, other, check=CHECK):
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
                c = self.sup2(a, b)
                assert c is not None, "sup2(%s, %s) == None" % (a, b)
                fa, fb = f[a], f[b]
                d = other.sup2(fa, fb)
                if d!=f[c]:
                    #print("\t fail 2", a, b, c, fa, fb, d)
                    break
            else:
                #print("\tappend")
                hom.append(f)
        hom = [Hom(self, other, f, check=check) for f in hom]
        pairs = [(f,g) for f in hom for g in hom if f<=g]
        return SupPoset(hom, pairs, check=check)


class Hom(object):
    "Morphism between Poset's"
    def __init__(self, src, tgt, send, check=CHECK):
        assert isinstance(src, PreOrder), src
        assert isinstance(tgt, PreOrder), tgt
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
        s = ', '.join("%s:%s"%k for k in self.items)
        return "Hom(%s)"%(s,)
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
        pairs = self.tgt.pairs
        for a in self.src.els:
            lhs = self.send[a]
            rhs = other.send[a]
            if (lhs, rhs) not in pairs:
                return False
        return True
    
    def __lt__(self, other):
        assert self.src == other.src
        assert self.tgt == other.tgt
        return not self==other and self<=other
    
    def is_ladj(self):
        "is a left adjoint functor"
        src = self.src
        tgt = self.tgt
        send = self.send
        for a in src.els:
          for b in src.els:
            c = src.sup2(a, b)
            if send[c] != tgt.sup2(send[a], send[b]):
                return False
        return True

    def is_radj(self):
        "is a right adjoint functor"
        src = self.src
        tgt = self.tgt
        send = self.send
        for a in src.els:
          for b in src.els:
            c = src.inf2(a, b)
            if send[c] != tgt.inf2(send[a], send[b]):
                return False
        return True

    def get_ladj(self, check=True):
        src = self.src
        tgt = self.tgt
        radj = self.send
        # XXX check that src has all inf's ?
        ladj = {}
        for x in tgt.els:
            ladj[x] = tgt.inf([y for y in src.els if (x,radj[y]) in tgt.pairs])
        ladj = Hom(tgt, src, ladj)
        if check:
            #print(src, "-->")
            #print(tgt)
            #print(self)
            assert self.is_radj(), self
            assert ladj.is_ladj(), ladj
        return ladj


def main():

    P = PreOrder.generate('01 12 20'.split())
    assert len(P) == 3

    hom = P.add_pairs(['03'])
    assert hom.src is P
    assert len(hom.tgt) == 4

    make_poset = lambda desc:Poset.generate(desc.split()).tgt

    I = Poset.generate(['01']).tgt # 0 <= 1
    assert len(I) == 2

    hom = I.add_pairs(['10'])
    assert hom.src is I
    assert len(hom.tgt) == 1

    P = Poset.generate(['ab', 'ba']).tgt
    assert len(P) == 1, P

    P = Poset.generate(['ab', 'bc', 'ca']).tgt
    assert len(P) == 1, P

    P = make_poset('aa bb')
    assert P.sup2('a', 'b') == None

    P = make_poset('a1 b1 a0 b0')
    assert P.sup2('a', 'b') == None

    P = make_poset('a1 b1 a0 b0 12 02')
    assert P.sup2('a', 'b') == None

    P = make_poset('a1 b1 a0 b0 a2 b2 20 21')
    assert P.sup2('a', 'b') == '2'

    P = make_poset('0a 0b 0c a1 b1 c1')
    assert len(P) == 5

    assert P.sup2('a', 'b') == '1'
    assert P.sup2('a', '0') == 'a'
    assert P.sup2('a', '1') == '1'

    i = P.get_id()
    assert i*i == i

    F = SupPoset.free([])
    assert len(F)==1
    assert F.bot is not None

    F = SupPoset.free('ABC')
    assert len(F) == 8
    assert F.sup2("A", "A") == "A"
    assert F.sup2("", "A") == "A"
    assert F.sup2("", "A+B") == "A+B"
    assert ("A", "A+B") in F.pairs
    assert F.sup2("A", "A+B") == "A+B"

    sup = F.sup
    AB = sup("AB")
    AC = sup("AC")
    BC = sup("BC")
    ABC = sup("ABC")
    hom = F.add_pairs([("A", BC), ("B", AC), ("C", AB),], check=True)
    G = hom.tgt
    assert hom.is_ladj()
    assert hom["A"] == "A"
    assert hom["B"] == "B"
    assert hom["C"] == "C"
    assert hom[AB] == hom[ABC] == G.sup("ABC")
    assert hom[AC] == G.sup("ABC")
    assert hom[BC] == G.sup("ABC")
    assert len(G) == 5

    assert len(list(f for f in I.hom_iter(I))) == 3

    P = make_poset('0a 0b 0c a1 b1 c1')
    for f in P.hom_iter(I):
        Hom(P, I, f).check()

    Q = P.hom(I)
    for f in Q:
        f.check()
    assert len(Q) == 10
    assert len(list(P.hom_iter(P))) == 178

    P = SupPoset.promote(P).tgt
    I = SupPoset.promote(I).tgt

    Q = P.hom(I)
    for f in Q:
        f.check()
    assert len(Q) == 5

    R = SupPoset.generate('0a 0b a1 b1'.split()).tgt
    RR = R.hom(R)
    RR = SupPoset.promote(RR)

    PP = P.hom(P, check=False)
    assert isinstance(PP, SupPoset)
    #PP.check()
    #assert len(P.hom(P))==50 # slow...

#    print("II =", len(I@I))
#
#    PP = P@P
#    print(PP.els)
#    print(len(P@P))

    F = SupPoset.free("ABC", check=True)
    assert len(F.hom(I)) == len(F)

    print("OK")


if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")
    else:
        main()



