#!/usr/bin/env python3


def closure(pairs):
    homs = set((a,b) for (a,b) in pairs) # set of (a,b) where a<=b
    els = set([a for (a,b) in pairs]+[b for (a,b) in pairs])
    for a in els:
        homs.add((a,a))
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
            if a!=b:
                assert a!=d, "loop found: %s = %s" % (a, b)
    ups = {}
    dns = {}
    for (a, b) in homs:
        up = ups.get(a, [])
        up.append(b)
        ups[a] = up
        dn = dns.get(b, [])
        dn.append(a)
        dns[b] = dn

    return homs, els, ups, dns


class Poset(object):
    def __init__(self, pairs):
        pairs, els, ups, dns = closure(pairs)
        items = list(pairs)
        items.sort() # canonical
        els = list(els)
        els.sort() # canonical
        self.items = items
        self.pairs = pairs
        self.els = els
        self.ups = ups
        self.dns = dns

    def __str__(self):
        return "Poset(%s, %s)"%(self.els, self.pairs)
    __repr__ = __str__

    def __eq__(self, other):
        return self.items == other.items

    def __ne__(self, other):
        return self.items != other.items

    def __hash__(self):
        return hash(tuple(self.items))

    def __getitem__(self, i):
        return self.els[i]

    def __len__(self):
        return len(self.els)

    def sup(self, a, b):
        pairs = self.pairs
        upa = set(self.ups[a])
        upb = set(self.ups[b])
        up = upa.intersection(upb)
        while len(up) > 1:
            items = list(up)
            n = len(up)
            for a in items:
              for b in items:
                if a==b:
                    continue
                if (a,b) in pairs and b in up:
                    up.remove(b)
            if len(up) == n:
                return None
        if len(up):
            return list(up)[0]


def main():

    I = Poset(['01']) # 0 <= 1
    #print(I)

    ok = False
    try:
        P = Poset(['ab', 'bc', 'ca'])
    except AssertionError:
        ok = True
    assert ok

    P = Poset('aa bb'.split())
    assert P.sup('a', 'b') == None

    P = Poset('0a 0b 0c a1 b1 c1'.split())
    assert len(P) == 5

    assert P.sup('a', 'b') == '1'
    assert P.sup('a', '0') == 'a'
    assert P.sup('a', '1') == '1'



if __name__ == "__main__":

    main()

