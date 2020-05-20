#!/usr/bin/env python3

import os
from random import choice

from bruhat.species import all_subsets
from bruhat.action import all_functions
from bruhat import equ, isomorph
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
        self.set_els = set(els)
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
        #return el in self.ups
        return el in self.set_els

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

    def slow_is_iso(self, other):
        for hom in self.hom(other):
            if hom.is_iso():
                return True
        return False

    def get_graph(self):
        "encode self into a graph structure for isomorph testing"
        points = []
        lookup = {}
        for idx, x in enumerate(self.els):
            p = isomorph.Point('', idx)
            points.append(p)
            lookup[x] = idx
        pairs = self.get_skel()
        graph = isomorph.Graph(points)
        for (a, b) in pairs:
            i = lookup[a]
            j = lookup[b]
            #graph.join(i, j)
            graph.add_directed(points[i], points[j])
        #print(graph)
        #print(lookup)
        return graph

    def is_iso(self, other, check=False):
        #if len(self.els) != len(other.els) or len(self.pairs) != len(other.pairs):
        if len(self.els) != len(other.els):
            return False

        #print("\n=========== is_iso ====================")
        lhs = self.get_graph()
        rhs = other.get_graph()

        for f in isomorph.search(lhs, rhs):
            send = dict((el, other.els[f[idx]]) for (idx, el) in enumerate(self.els))
            hom = Hom(self, other, send, check=check)
            if check:
                assert hom.is_iso()
            return True
        return False # no isomorph's found

    def get_op(self, check=False):
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

    def get_graph(self):
        "encode self into a graph structure for isomorph testing"
        points = []
        lookup = {}
        # decorate the tops and the bots is enough
        # for undirected graph to detect isomorphism of SupPoset's.
        bots = self.min(self.els)
        tops = self.max(self.els)
        names = {}
        for x in self.els:
            if x in bots:
                names[x] = 'b'
            elif x in tops:
                names[x] = 't'
        for idx, x in enumerate(self.els):
            desc = names.get(x, '')
            p = isomorph.Point(desc, idx)
            points.append(p)
            lookup[x] = idx
        pairs = self.get_skel()
        graph = isomorph.Graph(points)
        for (a, b) in pairs:
            i = lookup[a]
            j = lookup[b]
            graph.join(i, j)
            #graph.add_directed(points[i], points[j])
        #print(graph)
        #print(lookup)
        return graph


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
        assert len(gens) < 16, "too big: %d"%(len(gens),)
        items = list(all_subsets(gens))
        lookup = dict((a, set(a)) for a in items)
        pairs = [(a, b) for a in items for b in items 
            if lookup[a].issubset(lookup[b])]
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

    def add_pairs(self, rels, check=False):
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

    def get_rel(self): # SLOOOW
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

    def send_iter(self, other, check=False):
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

    def hom(self, other, check=False):
        hom = [Hom(self, other, send, check=check) 
            for send in self.send_iter(other, check)]
        pairs = [(f,g) for f in hom for g in hom if f<=g]
        return SupPoset(hom, pairs, check=check)

    def slow_tensor(L, R):
        lgen = L.get_gen()
        rgen = R.get_gen()
        #pair = lambda a, b : (a, b)
        gen = [(l, r) for l in lgen for r in rgen]
        #print("gen:", len(gen))
        #print("P:", len(P))
        send = dict((g, (g,)) for g in gen)
        pair = lambda a, b: send[(a, b)]
        lrels = L.get_rel()
        rrels = R.get_rel()
        rels = []
        for l in lgen:
            for (r, items) in rrels:
                #assert pair(l, r) in P, pair(l, r)
                rel = (pair(l, r), tuple(pair(l, r1) for r1 in items))
                rels.append(rel)
        for r in rgen:
            for (l, items) in lrels:
                #assert pair(l, r) in P, pair(l, r)
                rel = (pair(l, r), tuple(pair(l1, r) for l1 in items))
                rels.append(rel)

        # Only works for small examples
        P = SupPoset.free(gen)
        pairs = [(a, P.sup(items)) for (a, items) in rels]
        hom = P.add_pairs(pairs)
        P = hom.tgt
        return P

    def tensor(L, R, verbose=False):
        lgen = L.get_gen()
        rgen = R.get_gen()
        #pair = lambda a, b : (a, b)
        gen = [(l, r) for l in lgen for r in rgen]
        #print("gen:", len(gen))
        #print("P:", len(P))
        send = dict((g, (g,)) for g in gen)
        pair = lambda a, b: send[(a, b)]
        lrels = L.get_rel()
        rrels = R.get_rel()
        rels = []
        for l in lgen:
            for (r, items) in rrels:
                #assert pair(l, r) in P, pair(l, r)
                rel = (pair(l, r), tuple(pair(l, r1) for r1 in items))
                rels.append(rel)
        for r in rgen:
            for (l, items) in lrels:
                #assert pair(l, r) in P, pair(l, r)
                rel = (pair(l, r), tuple(pair(l1, r) for l1 in items))
                rels.append(rel)

        if verbose:
            print("tensor")
            print(gen)
            print(len(rels), rels)
        P = ConcretePoset(gen)
        pairs = [(a, P.sup(items)) for (a, items) in rels]
        if verbose:
            print(pairs)
        P = P.add_pairs(pairs) # does not construct hom ...
        return P
    __matmul__ = tensor

    def get_clean(P, check=False):
        gen = P.get_gen()
        rels = P.get_rel()
        Q = SupPoset.free(gen)
        pairs = [(a, Q.sup(items)) for (a, items) in rels]
        hom = Q.add_pairs(pairs)
        R = hom.tgt
        return R 

    def direct_sum(L, R):
        els = [(l, r) for l in L.els for r in R.els]
        pairs = [((l1, r1), (l2, r2)) 
            for (l1,r1) in els for (l2,r2) in els
            if (l1,l2) in L.pairs and (r1,r2) in R.pairs]
        return SupPoset(els, pairs)
    __add__ = direct_sum
    __mul__ = direct_sum


#class Set(object):
#    def __init__(self, items):
#        assert isinstance(items, tuple)
#        self.items = items

    def is_modular(self):
        pairs = self.pairs
        sup2, inf2 = self.sup2, self.inf2
        for (c, a) in pairs:
            for b in self.els:
                lhs = inf2(a, sup2(b, c))
                rhs = sup2(inf2(a, b), c)
                if lhs != rhs:
                    return False
        return True

    def is_distributive(self):
        els = self.els
        sup2, inf2 = self.sup2, self.inf2
        for a in els:
         for b in els:
          for c in els:
            lhs = inf2(a, sup2(b, c))
            rhs = sup2(inf2(a, b), inf2(a, c))
            if lhs != rhs:
                #print("is_distributive", a, b, c, lhs, rhs)
                return False
        return True


class ConcretePoset(Poset):
    def __init__(self, gen, check=True):
        els = list(all_subsets(gen))
        els.sort(key = lambda i:(len(i), i)) # shortest first
        set_els = set(els)
        gen = [(a,) for a in gen]
        for a in gen:
            assert a in set_els
        self.gen = gen
        assert len(set(gen)) == len(self.gen) # uniq
        assert len(els) == 2**len(gen)
        assert len(els[0])==0
        assert len(els[-1])==len(gen)
        self._bot = els[0]
        self._top = els[-1]
        self.els = els
        self.set_els = set_els

    def min(self, items):
        if items is self.els: # HACK
            return [self._bot]
        assert 0, "TODO"

    def max(self, items):
        if items is self.els: # HACK
            return [self._top]
        assert 0, "TODO"

    def get_skel(self):
        els = self.els
        pairs = []
        for a in els:
          n = len(a)
          sa = set(a)
          for b in els:
            if len(b) != n+1:
                continue
            sb = set(b)
            if not sb.issuperset(sa):
                continue
            if len(sb.difference(sa)) == 1:
                pairs.append((a, b)) # a <= b
        #print("ConcretePoset.get_skel")
        #print(pairs)
        return pairs

    def sup2(self, a, b): 
        top = self.top
        assert a in self.set_els
        assert b in self.set_els
        c = tuple(i for i in top if i in b or i in a) # union
        assert c in self.set_els
        return c

    def inf2(self, a, b): 
        assert a in self.set_els
        assert b in self.set_els
        c = tuple(i for i in a if i in b) # intersection
        assert c in self.set_els
        return c

    def add_pairs(self, rels, check=False):
        #if not rels:
        #    return self # XX convert to SupPoset

        "add relations to form a new SupPoset"
        set_els = self.set_els
        for (a, b) in rels:
            assert a in set_els, (a, list(set_els))
            assert b in set_els, (b, list(set_els))

        le = lambda a,b : set(a).issubset(b)
        #le = lambda a,b:self.sup2(a,b)==b # a bit slower
        els = []
        for x in self.els:
            for (lhs, rhs) in rels: # lhs <= rhs
                if le(rhs,x) and not le(lhs,x):
                    break
            else:
                els.append(x)
        els = set(els)
        #pairs = [(a, b) for (a, b) in self.pairs if a in els and b in els]
        pairs = [(a, b) for a in els for b in els if le(a, b)]
        P = SupPoset(els, pairs)
        ## inclusion is a right adjoint
        #radj = Hom(P, self, dict((a, a) for a in els))
        #ladj = radj.get_ladj()
        #return ladj
        return P


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
            if p not in tgt.pairs:
                print("not order preserving: %s --> %s"%(
                    (a,b), (send[a],send[b])))
                print(self)
                raise AssertionError

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
        els = src.els
        pairs = tgt.pairs
        for x in tgt.els:
            items = []
            up = tgt.ups[x]
            for y in els:
                if radj[y] in up:
                    items.append(y)
            ladj[x] = tgt.inf(items)
            # much slower:
            #ladj[x] = tgt.inf([y for y in src.els if (x,radj[y]) in tgt.pairs])
        ladj = Hom(tgt, src, ladj)
        if check:
            #print(src, "-->")
            #print(tgt)
            #print(self)
            assert self.is_radj(), self
            assert ladj.is_ladj(), ladj
        return ladj


def main():

    # ---------------------------------------
    # PreOrder

    make_preorder = lambda desc:PreOrder.generate(desc.split(), check=True)

    P = PreOrder.generate([], list('012'))
    assert len(P) == 3

    P = make_preorder('01 12 20')
    assert len(P) == 3

    hom = P.add_pairs(['03'])
    assert hom.src is P
    Q = hom.tgt
    assert len(Q) == 4

    assert not P.is_iso(Q)

    Q = make_preorder('01 10 02')
    assert len(P)==len(Q)
    assert not P.is_iso(Q)

    Q = make_preorder('12 23 31')
    assert P.is_iso(Q)

    # ---------------------------------------
    # Poset

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

    #P = make_poset('aa bb cc')

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

    # ---------------------------------------
    # SupPoset

    F = SupPoset.free([])
    assert len(F)==1
    assert F.bot is not None

    F3 = SupPoset.free('ABC')
    assert len(F3) == 8
    assert F3.sup2("A", "A") == "A"
    assert F3.sup2("", "A") == "A"
    assert F3.sup2("", "A+B") == "A+B"
    assert ("A", "A+B") in F3.pairs
    assert F3.sup2("A", "A+B") == "A+B"

    assert F3.get_gen() == set("ABC")

    sup = F3.sup
    AB = sup("AB")
    AC = sup("AC")
    BC = sup("BC")
    ABC = sup("ABC")
    hom = F3.add_pairs([("A", BC), ("B", AC), ("C", AB),], check=True)
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

    PP = P.hom(P)
    assert isinstance(PP, SupPoset)
    PP.check()
    assert len(PP) == 50

    F3 = SupPoset.free("ABC", check=True)
    assert len(F3.hom(I)) == len(F3)
    #F3.get_dot("out.dot")

    P = SupPoset.free('abd')
    hom = P.add_pairs('bd'.split())
    a, b, d = [hom[i] for i in 'abd']
    Q = hom.tgt
    assert Q.sup2(b, d) == d
    c = Q.sup2(a, b)
    e = Q.sup2(c, d)
    assert set([Q.bot, a, b, c, d, e]) == set(Q.els)

    # 3x3 trellis
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

    P = make_supposet('0a 0b a1 b1')
    assert P.is_iso(P)

    # ----------------------------------------

    # 3x3 trellis
    P = make_supposet('0a 0b ac ae be bd cf ef eg dg fh gh')
    assert P.hom(I).is_iso(P)
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
    # 

    P = ConcretePoset('abcd')
    assert len(P) == 2**4
    els = P.els
    for a in els:
      for b in els:
        c = P.sup2(a, b)
        assert set(c).issuperset(a)
        assert set(c).issuperset(b), (a, b, c)

    # ----------------------------------------
    # 
    
    II = I@I
    assert len(II) == 2

    F2 = make_supposet('0a 0b a1 b1')
    assert F2.get_gen() == set('ab')
    assert F2.get_rel() == set()

    T = make_supposet('0a ab') # a chain
    #print(T.get_gen())
    #print(T.get_rel())

    TT = T@T
    assert len(TT.get_gen()) == 4

    TT1 = T.slow_tensor(T)
    
    assert TT.is_iso(TT1)

    # -------------------------

    P = F2*F2*I
    assert len(P) == 32

    assert F2.is_modular()
    assert F2.is_distributive()
    assert (I@F2).is_iso(F2)
    assert (F2@F2).is_iso(SupPoset.free('ABCD'))
    
    assert len(F3) == 8
    assert F3.is_modular()
    assert F3.is_distributive()
    
    M3 = make_supposet('0a 0b 0c a1 b1 c1')
    assert M3.get_gen() == set('abc')
    assert M3.get_rel() == {('c',('a','b')), ('a',('b','c')), ('b',('a','c'))}
    assert M3.is_modular()
    assert not M3.is_distributive()

    #PP = M3@M3
    PP = M3.slow_tensor(M3)
    #print(PP)
    #print(PP.els)
    assert len(PP) == 50, len(PP)

    M3I = M3*I
    M3I.check()
    assert len(M3I) == 10
    assert len(M3I.get_gen()) == 4
    assert M3I.is_modular()
    assert not M3I.is_distributive()

    # smallest non-modular lattice
    N5 = make_supposet('0a ab b1 0c c1')
    N5op = N5.hom(I)
    assert N5.is_iso(N5op)
    assert not N5.is_modular()
    assert not N5.is_distributive()

    lhs = I.tensor(N5, verbose=False)
    #print(len(lhs))
    assert (I@N5).is_iso(N5)

    # tadpole's
    V = make_supposet('0a ab ac b1 c1')
    Vop = make_supposet('0a 0b ac bc c1')
    assert not V.is_iso(Vop)
    assert Vop.is_iso(V.hom(I))
    assert V.is_iso(Vop.hom(I))
    assert V.is_modular()
    assert V.is_distributive()
    assert Vop.is_modular()
    assert Vop.is_distributive()

    P = V@V
    #print(len(P))
    assert P.is_iso(V.slow_tensor(V))

    # 3x2 trellis
    B = make_supposet('0a 0b bc ac bd d1 c1')
    assert B.is_modular()
    assert B.is_distributive()

    L7 = make_supposet('0a 0b ac ad ae bd c1 d1 e1')
    #N6 = make_supposet('0a ab bc c1 0d d1')
    N6 = make_supposet('0a ab b1 0c cd d1')
    assert not N6.is_modular()
    assert not N6.is_distributive()

    def test_duals(A, B):
        Aop = A.get_op()
        Bop = B.get_op()
        lhs = A@B
        rhs = Aop.hom(B)
        a = lhs.is_iso(rhs)
        rhs = (Aop @ Bop).get_op()
        b = lhs.is_iso(rhs)
        return [a, b]

    # test that the distributive lattices have dual objects
    names = "I F2 T V B M3 N5".split()
    for sA in names:
      A = eval(sA)
      assert A.get_op().is_iso(A.hom(I))
      for sB in names:
        B = eval(sB)
        result = test_duals(A, B)
        result += [
            A.is_modular(), A.is_distributive(), 
            B.is_modular(), B.is_distributive()]
        result = [int(i) for i in result]
        #print(sA, sB, '\t', result)
        assert result[0] == result[1]
        if result[3] and result[5]:
            assert result[0] # distributive implies duals

    print("OK")


if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")
    else:
        main()



