#!/usr/bin/env python3

import os
from random import choice

from bruhat.species import all_subsets
from bruhat.action import all_functions
from bruhat import equ
from bruhat.argv import argv

#CHECK = True # much slower...
CHECK = False

OPEN_COMMAND = "gvfs-open" # used to display pdf files

def closure(pairs, els):
    homs = set((a,b) for (a,b) in pairs) # set of (a,b) where a<=b
    for a in els:
        homs.add((a, a))
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
    def generate(cls, pairs, els=[], check=CHECK):
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

    def send_iter(self, other):
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
        hom = [Hom(self, other, f) for f in self.send_iter(other)]
        pairs = [(f,g) for f in hom for g in hom if f<=g]
        return Poset(hom, pairs)

    def is_iso(self, other):
        for hom in self.hom(other):
            if hom.is_iso():
                return True
        return False

    def get_op(self, check=True):
        els = self.els
        pairs = [(b, a) for (a, b) in self.pairs]
        return self.__class__(els, pairs, check=check)

    def get_skel(self):
        remain = set((a,b) for (a,b) in self.pairs if a!=b)
        for (a, b) in self.pairs:
            if a==b:
                continue
            for (c, d) in self.pairs:
                if c!=b or c==d:
                    continue
                if (a,d) in remain:
                    remain.remove((a,d))
        return remain

    def get_dot(self, filename=None, labels=True):
        pairs = self.get_skel()
        lines = ["digraph {"]
        if labels:
            lookup = dict((a, str(a)) for a in self.els)
        else:
            lookup = dict((a, idx+1) for (idx,a) in enumerate(self.els))
        for a,b in pairs:
            a = lookup[a]
            b = lookup[b]
            a = a or "_"
            line = '  "%s" -> "%s"'%(b, a)
            #line = line.replace("+", "")
            lines.append(line)
        lines.append("}")
        s = "\n".join(lines)
        if filename is None:
            return s
        f = open(filename, 'w')
        print(s, file=f)
        f.close()

    def show(self, labels=True):
        i = str(hash(self)).replace('-', 'm')
        name = "tmp.%s.dot"%i
        self.get_dot(name, labels=labels)
        os.system("dot -Tpdf %s > %s.pdf"% (name,name))
        os.system("%s %s.pdf" % (OPEN_COMMAND, name,))


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
        send = equ.quotient_rep(els, relation) # XXX random choice here XXX
        #print("send:", send)
        get = send.get
        els = set(get(a, a) for a in els)
        pairs = set((get(a, a), get(b, b)) for (a, b) in pairs)
        tgt = Poset(els, pairs, check=check)
        hom = Hom(src, tgt, send)
        return hom

    @classmethod
    def generate(cls, pairs, check=CHECK):
        P = PreOrder.generate(pairs, check=check)
        hom = cls.promote(P, check=check)
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

    def min(self, items): # cache ?
        #orig = items
        items = set(items)
        if not items:
            return items
        for (a, b) in self.pairs:
            if a!=b and a in items and b in items:
                items.remove(b)
        #assert items == self.slow_min(orig)
        return items

    def max(self, items): # cache ?
        items = set(items)
        if not items:
            return items
        for (a, b) in self.pairs:
            if a!=b and a in items and b in items:
                items.remove(a)
        return items

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
        if a>b:
            a, b = b, a
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
        if a>b:
            a, b = b, a
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
 
#    def sup(self, items):
#        if not items:
#            return self.bot
#        items = list(items)
#        a = items[0]
#        if len(items) == 1:
#            return a
#        #print("sup", items)
#        up = set(self.ups[a])
#        for i in range(1, len(items)):
#            a = items[i]
#            up = up.intersection(self.ups[a])
#        value = self.uniq_min(up)
#        #print("sup", items, value)
#        return value

    def inf(self, items):
        items = list(items)
        while len(items)>1:
            a = self.inf2(items.pop(), items.pop())
            items.append(a)
        if items:
            return items[0]
        return self.top

#    def inf(self, items):
#        if not items:
#            return self.top
#        items = list(items)
#        a = items[0]
#        if len(items) == 1:
#            return a
#        #print("inf", items)
#        dn = set(self.dns[a])
#        for i in range(1, len(items)):
#            a = items[i]
#            dn = dn.intersection(self.dns[a])
#        value = self.uniq_max(dn)
#        #print("inf", items, value)
#        return value


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

#    @classmethod
#    def free_on_poset(cls, P, check=CHECK):
        

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

    def get_gen(self):
        gen = set(self.els)
        gen.remove(self.bot)
        pairs = self.pairs
        for a in self.els:
            for b in self.els:
                if b==a:
                    break
                if (a, b) in pairs or (b, a) in pairs:
                    continue
                c = self.sup2(a, b)
                assert c!=a and c!=b
                if c in gen:
                    #print("%s v %s = %s" % (a, b, c))
                    gen.remove(c)
        return gen

    def add_pairs(self, rels, check=True):
        "add relations to form a new SupPoset"
        ups = self.ups
        for (a, b) in rels:
            assert a in ups, (a, list(ups.keys()))
            assert b in ups, (b, list(ups.keys()))

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

    def slow_hom(self, other, check=CHECK):
        "poset of hom's that preserve sup's, including empty sup == bot"
        els = self.els
        pairs = [(a, b) for a in els for b in els]
        assert self.bot is not None, "self is not sup-complete"
        assert other.bot is not None, "other is not sup-complete"
        hom = []
        for f in self.send_iter(other):
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

    def get_rels(self): # SLOOOW
        pairs = self.pairs
        gen = list(self.get_gen())
        gen.sort()
        rels = set()
        for rhs in all_subsets(gen):
            x = self.sup(rhs)
            for g in gen:
                if g in rhs:
                    continue # tautology
                if (g, x) in pairs:
                    rels.add((g, rhs))
        # now find a minimal set of relations
        pairs = set()
        for a in rels:
          for b in rels:
            if a[0]==b[0] and set(a[1]).issubset(b[1]):
                pairs.add((a, b)) # b implies a
        R = Poset(rels, pairs)
        rels = R.min(R.els) # woop
        return rels

    def send_iter(self, other, check=True):
        assert isinstance(other, SupPoset)
        gen = self.get_gen()
        #print("gen:", gen)
        els = other.els
        for send in all_functions(gen, els):
            for items in all_subsets(gen):
                lhs = self.sup(items)
                rhs = other.sup([send[g] for g in items])
                x = send.get(lhs)
                if x==rhs:
                    continue
                #print(lhs, "-->", rhs, x)
                if x is not None:
                    break
                send[lhs] = rhs
            else:
                #print(send)
                yield send

    def hom(self, other, check=True):
        hom = [Hom(self, other, send, check=check) 
            for send in self.send_iter(other, check)]
        pairs = [(f,g) for f in hom for g in hom if f<=g]
        return SupPoset(hom, pairs, check=check)

    def is_iso(self, other): # SLOOOOOOOOOOOOOW
        assert isinstance(other, SupPoset)
        gen = self.get_gen()
        #print("gen:", gen)
        els = other.els
        #els = other.get_gen()
        for send in all_functions(gen, els):
            send = dict(send)
            vals = list(send.values())
            if len(vals) < len(gen):
                continue
            if other.bot in vals:
                continue
            #print()
            #print(send)
            for items in all_subsets(gen):
                lhs = self.sup(items)
                rhs = other.sup([send[g] for g in items])
                assert rhs is not None
                x = send.get(lhs)
                if x==rhs:
                    continue
                #print(lhs, "-->", rhs, x)
                if x is not None:
                    break
                send[lhs] = rhs
            else:
                #print(send)
                f = Hom(self, other, send)
                if f.is_iso():
                    return True
        return False

    def tensor(L, R):
        lgen = L.get_gen()
        rgen = R.get_gen()
        pair = lambda a, b : (a, b)
        gen = [pair(l, r) for l in lgen for r in rgen]
        #print("gen:", gen)
        P = SupPoset.free(gen)
        send = dict((g, (g,)) for g in gen)
        pair = lambda a, b: send[(a, b)]
        lrels = L.get_rels()
        rrels = R.get_rels()
        rels = []
        for l in lgen:
            for (r, items) in rrels:
                assert pair(l, r) in P, pair(l, r)
                rel = (pair(l, r), tuple(pair(l, r1) for r1 in items))
                rels.append(rel)
        for r in rgen:
            for (l, items) in lrels:
                assert pair(l, r) in P, pair(l, r)
                rel = (pair(l, r), tuple(pair(l1, r) for l1 in items))
                rels.append(rel)

        #print("rels:", len(rels))
        pairs = [(a, P.sup(items)) for (a, items) in rels]
        #print("pairs")
        hom = P.add_pairs(pairs)
        P = hom.tgt

        return P

    __matmul__ = tensor

    def get_clean(P, check=True):
        gen = P.get_gen()
        rels = P.get_rels()
        Q = SupPoset.free(gen)
        pairs = [(a, Q.sup(items)) for (a, items) in rels]
        hom = Q.add_pairs(pairs)
        R = hom.tgt
        return R 

    def sum(L, R):
        els = [(l, r) for l in L.els for r in R.els]
        pairs = [((l1, r1), (l2, r2)) 
            for (l1,r1) in els for (l2,r2) in els
            if (l1,l2) in L.pairs and (r1,r2) in R.pairs]
        return SupPoset(els, pairs)
    __add__ = sum



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
        if other is None:
            return False
        assert self.src == other.src
        assert self.tgt == other.tgt
        return self.send == other.send

    def __ne__(self, other):
        if other is None:
            return True
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

    def is_iso(self):
        src = self.src
        tgt = self.tgt
        send = self.send
        if len(src) != len(tgt):
            return False
        vals = set(send.values())
        if len(vals) != len(tgt):
            return False
        # check inverse is a Hom
        inv = dict((b, a) for (a, b) in send.items())
        for (a, b) in tgt.pairs:
            if (inv[a], inv[b]) not in src.pairs:
                return False
        return True
    
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

    def get_ladj(self, check=False):
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

    P = PreOrder.generate([], list('012'))
    assert len(P) == 3

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

    assert F.get_gen() == set("ABC")

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
    assert G.get_gen() == set("ABC")
    G.get_dot("out.dot")

    P = SupPoset.promote(make_poset('0a 0b ac bc cd')).tgt
    assert P.get_gen() == set('abd')

    Q = P.get_op()
    assert Q.bot == 'd'
    assert Q.get_gen() == set('abc')
    assert Q.get_op() == P # equality on the nose.

    # hom's -----------------------------------------

    assert len(list(f for f in I.send_iter(I))) == 3

    P = make_poset('0a 0b 0c a1 b1 c1')
    for f in P.send_iter(I):
        Hom(P, I, f).check()

    Q = P.hom(I)
    for f in Q:
        f.check()
    assert len(Q) == 10
    assert len(list(P.send_iter(P))) == 178

    f = P.get_id()
    assert f.is_iso()

    P = SupPoset.promote(P).tgt
    I = SupPoset.promote(I).tgt

    Q = P.hom(I)
    for f in Q:
        f.check()
    assert len(Q) == 5

    # --------------------------
    # SupPoset's

    make_supposet = lambda desc : SupPoset.generate(desc.split(), check=True).tgt

    R = make_supposet('0a 0b a1 b1')
    RR = R.hom(R)
    RR = SupPoset.promote(RR)

    PP = P.hom(P, check=False)
    assert isinstance(PP, SupPoset)
    PP.check()
    assert len(PP) == 50

    F = SupPoset.free("ABC", check=True)
    assert len(F.hom(I)) == len(F)
    #F.get_dot("out.dot")

    P = SupPoset.free('abd')
    hom = P.add_pairs('bd'.split())
    a, b, d = [hom[i] for i in 'abd']
    Q = hom.tgt
    assert Q.sup2(b, d) == d
    c = Q.sup2(a, b)
    e = Q.sup2(c, d)
    assert set([Q.bot, a, b, c, d, e]) == set(Q.els)

    P = SupPoset.free('abcd')
    hom = P.add_pairs('ac bd'.split())
    Q = hom.tgt
    #Q.get_dot("out.dot")
    gen = set(hom[x] for x in 'abcd')
    assert Q.get_gen() == gen

    # ----------------------------------------

    P = SupPoset.free(range(8))
    assert P.bot == ()
    assert len(P) == 2**8

    # ----------------------------------------

    P = make_supposet('0a 0b ac ae be bd cf ef eg dg fh gh')
    R = P.get_clean()
    assert R.is_iso(P)

    # ----------------------------------------

    Q = SupPoset.free('DEFGH')
    rels = [('D', ('E', 'F')), ('F', ('G', 'H'))]
    pairs = [(a, Q.sup(items)) for (a, items) in rels]
    hom = Q.add_pairs(pairs)
    #R = hom.tgt
    #R.get_dot("out.dot")
    #print(len(R))

    Q = SupPoset.free('ABC')
    rels = [('B', ('A', 'C')),]
    pairs = [(a, Q.sup(items)) for (a, items) in rels]
    pairs.append(("A", "B"))
    hom = Q.add_pairs(pairs)
    #R = hom.tgt
    #R.get_dot("out.R.dot")
    #print(len(R))

    Q = SupPoset.free('ABC')
    rels = [('B', ('A', 'C')),('A', ('B', 'C')),  ('C', ('A', 'B')), ]
    pairs = [(a, Q.sup(items)) for (a, items) in rels]
    hom = Q.add_pairs(pairs)
    #R = hom.tgt
    #R.get_dot("out.P.dot")
    #print(len(R))

    #return

    # ----------------------------------------
    # tensor
    
    II = I@I
    assert len(II) == 2

    F2 = make_supposet('0a 0b a1 b1')
    assert F2.get_gen() == set('ab')
    assert F2.get_rels() == set()

    FF = F2@F2
    assert len(FF) == 2**4

    P = make_supposet('0a 0b 0c a1 b1 c1')
    assert P.get_gen() == set('abc')
    assert P.get_rels() == {('c',('a','b')), ('a',('b','c')), ('b',('a','c'))}
    #PP = P@P # slow.....
    #assert len(PP) == 50, len(PP)

    PI = P+I
    PI.check()
    assert len(PI) == 10
    assert len(PI.get_gen()) == 4
    #PI.show()

    M = make_supposet('0a ab b1 0c c1')
    Mop = M.hom(I)
    assert M.is_iso(Mop)

    V = make_supposet('0a ab ac b1 c1')
    Vop = make_supposet('0a 0b ac bc c1')
    assert not V.is_iso(Vop)
    assert Vop.is_iso(V.hom(I))
    assert V.is_iso(Vop.hom(I))

    if 1:
        #V_Vop = V.hom(Vop)
        #print(len(V_Vop)) # == 48
    
        V_V = V.hom(V)
        assert len(V_V) == 50
    
        VVop = V@Vop
        assert len(VVop) == 50

        #assert len(VVop) == 50
        #print(V_V.is_iso(VVop))

        #assert V_V.is_iso(VVop)
        #V_V.show(labels=False)
        #VVop.show(labels=False)
    
        #Vop_V = Vop.hom(V) # == 48
        #print(V_Vop.is_iso(Vop_V))

    F2op = F2.hom(I)
    #assert F2op.is_iso(F2)
    lhs = Vop.hom(F2)
    rhs = V @ F2
    #assert lhs.is_iso(rhs)

    #lhs.show(labels=False)
    #rhs.show(labels=False)

    print("OK")


if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")
    else:
        main()


