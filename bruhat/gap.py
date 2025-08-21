#!/usr/bin/env python3

"""
parse character table output from gap.
"""

import sys, os
from time import sleep
import select
from subprocess import Popen, PIPE
import re

from bruhat.argv import argv


PROMPT = "gap> "

# https://superuser.com/a/1657976
def escape_ansi(line):
    re1 = re.compile(r'\x1b\[[\x30-\x3f]*[\x20-\x2f]*[\x40-\x7e]')
    re2 = re.compile(r'\x1b[PX^_].*?\x1b\\')
    re3 = re.compile(r'\x1b\][^\a]*(?:\a|\x1b\\)')
    re4 = re.compile(r'\x1b[\[\]A-Z\\^_@]')
    # re5: zero-width ASCII characters
    # see https://superuser.com/a/1388860
    re5 = re.compile(r'[\x00-\x1f\x7f-\x9f\xad]+')

    for r in [re1, re2, re3, re4, re5]:
        line = r.sub('', line)

    return line


class Ref:
    idx = 0
    def __init__(self, ref=None):
        if ref is None:
            ref = "v%d"%Ref.idx
            Ref.idx += 1
        self.ref = ref
    def __str__(self):
        return self.ref
    gapstr = __str__

    def to_group(self):
        gap = Gap.the_gap
        
        


class Gap:
    the_gap = None # Singleton
    DEBUG = False

    def __new__(cls):
        if cls.the_gap is not None:
            return cls.the_gap
        gap = object.__new__(cls)
        gap.init()
        cls.the_gap = gap
        return gap

    def init(self):
        self.proc = Popen("gap", bufsize=0, stdin=PIPE, stdout=PIPE)
        sleep(0.1)
        #print(dir(proc.stdout))
        #print(proc.stdout.read(20))
        self.buf = ""
        self.pos = 0
        self.debug = open("gap_debug.dump", "w")
        self.expect()

    def read_nb(self, raw=False):
        proc = self.proc
        data = bytes()
        while select.select([proc.stdout],[],[],0)[0]!=[]:   
            data += proc.stdout.read(1)
        data = data.decode("utf-8")
        self.debug.write(data)
        if not raw:
            data = escape_ansi(data)
        #assert "Error" not in data # XXX read stderr XXX
        self.buf += data

    def send(self, data):
        proc = self.proc
        if self.DEBUG:
            print("Gap.send(%r)"%data)
        data = data + "\n"
        self.debug.write(data)
        data = data.encode("utf-8")
        proc.stdin.write(data)

    def expect(self, s=PROMPT, raw=False):
        while s not in self.buf[self.pos:]:
            self.read_nb(raw)
            sleep(0.001)

        if self.DEBUG:
            print("Gap.expect: buf=%r"%self.buf[self.pos:])
        pos = self.pos + self.buf[self.pos:].index(s)
        data = self.buf[self.pos : pos]
        self.pos = pos + len(s)
        if self.DEBUG:
            print("Gap.expect: return %r"%data)
        return data

    def define(self, item, ref=None):
        item = to_gap(item)
        ref = Ref(ref)
        assert not item.endswith(";")
        s = "%s := %s;;"%(ref, item)
        self.send(s)
        s = self.expect()
        return ref

    def __getattr__(self, name):
        return Method(self, name)

    def length(self, item):
        value = self.Length(item, get=True)
        return int(value)

    def iterate(self, item):
        item = to_gap(item)
        n = self.length(item)
        for i in range(n):
            s = "%s[%d];"%(item, i+1)
            self.send(s)
            yield self.expect()

    def tom(self, ref):
        self.SizeScreen([1000,1000])
        ref = self.TableOfMarks(ref)
        data = self.Display(ref, get=True, raw=True)
        from bruhat.tom import parse_tom
        tom = parse_tom(data)
        tom.gapref = ref
        return tom

    def to_group(gap, N, ref):
        from bruhat.gset import Group, Perm
        if N is None:
            assert 0, "arggghh"
            orbit = gap.OrbitPerms(ref, 1)
            N = gap.Maximum(orbit, get=True)
            N = int(N)
            print("N =", N)
        gens = gap.GeneratorsOfGroup(ref)
        perms = []
        for gapstr in gap.iterate(gens):
            perm = Perm.from_gapstr(N, gapstr)
            perms.append(perm)
        #G = Group.generate(perms, verbose=True)
        G = Group(None, perms)
        return G



class Method:
    def __init__(self, gap, name):
        self.gap = gap
        self.name = name
    def __call__(self, *args, ref=None, get=False, raw=False):
        args = [to_gap(arg) for arg in args]
        gapstr = "%s(%s)"%(self.name, ','.join(args))
        if get:
            gapstr = "%s;"%gapstr
        else:
            ref = Ref(ref)
            gapstr = "%s := %s;;"%(ref, gapstr)
        gap = self.gap
        gap.send(gapstr)
        s = gap.expect(raw=raw)
        if get:
            return s
        else:
            return ref


def to_gap(item):
    if type(item) is str:
        return item
    if type(item) is int:
        return str(item)
    if hasattr(item, "gapstr"):
        return item.gapstr()
    if hasattr(item, "gapref"):
        return str(item.gapref)
    if isinstance(item, list):
        return "[%s]"%(','.join(to_gap(i) for i in item)) # <-- recurse
    assert 0, "what's this: %r ?"%item


example_Irr = """\
[ Character( CharacterTable( <pc group of size 28 with 3 generators> ),
  [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ] ), Character( CharacterTable( <pc group of size 
    28 with 3 generators> ), [ 1, -1, 1, 1, -1, 1, 1, 1, 1, 1 ] ), 
  Character( CharacterTable( <pc group of size 28 with 3 generators> ),
  [ 1, -E(4), -1, 1, E(4), -1, 1, -1, 1, -1 ] ), 
  Character( CharacterTable( <pc group of size 28 with 3 generators> ),
  [ 1, E(4), -1, 1, -E(4), -1, 1, -1, 1, -1 ] ), 
  Character( CharacterTable( <pc group of size 28 with 3 generators> ),
  [ 2, 0, -2, E(7)^3+E(7)^4, 0, -E(7)^3-E(7)^4, E(7)+E(7)^6, -E(7)-E(7)^6, E(7)^2+E(7)^5, 
      -E(7)^2-E(7)^5 ] ), Character( CharacterTable( <pc group of size 28 with 
    3 generators> ),
  [ 2, 0, -2, E(7)^2+E(7)^5, 0, -E(7)^2-E(7)^5, E(7)^3+E(7)^4, -E(7)^3-E(7)^4, 
      E(7)+E(7)^6, -E(7)-E(7)^6 ] ), Character( CharacterTable( <pc group of size 28 with 
    3 generators> ),
  [ 2, 0, -2, E(7)+E(7)^6, 0, -E(7)-E(7)^6, E(7)^2+E(7)^5, -E(7)^2-E(7)^5, E(7)^3+E(7)^4, 
      -E(7)^3-E(7)^4 ] ), Character( CharacterTable( <pc group of size 28 with 
    3 generators> ),
  [ 2, 0, 2, E(7)^3+E(7)^4, 0, E(7)^3+E(7)^4, E(7)+E(7)^6, E(7)+E(7)^6, E(7)^2+E(7)^5, 
      E(7)^2+E(7)^5 ] ), Character( CharacterTable( <pc group of size 28 with 
    3 generators> ),
  [ 2, 0, 2, E(7)^2+E(7)^5, 0, E(7)^2+E(7)^5, E(7)^3+E(7)^4, E(7)^3+E(7)^4, E(7)+E(7)^6, 
      E(7)+E(7)^6 ] ), Character( CharacterTable( <pc group of size 28 with 
    3 generators> ),
  [ 2, 0, 2, E(7)+E(7)^6, 0, E(7)+E(7)^6, E(7)^2+E(7)^5, E(7)^2+E(7)^5, E(7)^3+E(7)^4, 
      E(7)^3+E(7)^4 ] ) ]
"""

def parse_Irr(s):
    print(s)
    s = s.replace("\n", " ")
    s = s.strip()
    assert s[:2] == "[ "
    assert s[-2:] == " ]"
    s = s[2:-2]

    s = s.replace("E(4)", "i")
    s = s.replace("E(8)-E(8)^3", r"\sqrt{2}")
    s = s.replace("-E(8)+E(8)^3", r"-\sqrt{2}")
    s = s.replace("-E(12)^7+E(12)^11", r"\sqrt{3}")
    s = s.replace("E(12)^7-E(12)^11", r"-\sqrt{3}")

    i = 3
    while "E(" in s:
        s = s.replace("E(%d)"%i, r"\zeta_{%d}"%i)
        i += 1

    # hack it:
    output = []
    items = s.split("Character( ")
    items.pop(0)
    #print(items)
    for item in items:
        # find matching close brace:
        item = item.split("CharacterTable(")[1]
        #print("scan:", repr(item))
        count = 1
        idx = 0
        while idx < len(item) and count:
            if item[idx] == "(":
                count += 1
            elif item[idx] == ")":
                count -= 1
            idx += 1
        assert count == 0, repr(item)
        item = item[idx:] # chop
        #print("chop:", repr(item))

        # look for characters:
        i = item.index("[")
        j = item.index("]")
        row = item[i+1:j]
        row = row.strip()
        output.append(row)
    return output


def latex_nosep(rows, desc):
    n = len(rows[0])
    lines = [("$$")]
    lines.append(r"\begin{array}{%s}"%(desc))
    for i, row in enumerate(rows):
        if type(row)==list:
            line = " & ".join(str(fld) for fld in row) + r" \\"
        else:
            line = str(row)
        lines.append(line)
    lines.append(r"\end{array}")
    lines.append("$$")
    s = "\n".join(lines)
    return s



def gen_latex(output, perm=None):
    rows = [r"\hline"]
    n = None
    for i, row in enumerate(output):
        flds = row.split(", ")
        flds = [fld.strip() for fld in flds]
        assert n is None or n==len(flds)
        n = len(flds)
        if perm is not None:
            assert len(perm)==n
            assert len(set(perm))==n
            flds = [flds[j] for j in perm]
        #print(flds)
        flds.insert(0, r"\rho_{%d}"%(i+1))
        rows.append(flds)
    row = [r"\mathrm{class}"] + ["1"]*n
    rows.insert(0, row)
    s = latex_nosep(rows, "c|" + "c"*n)
    print(s)


if __name__ == "__main__":
    cmd = argv.cmd
    if cmd:
        
        gap = Gap()
        
        #gap.send("G:=SmallGroup(28,1);")
    
        if cmd[-1] != ";":
            cmd += ";"
        gap.send(cmd)
        
        s = gap.expect()
        
        gap.send("Irr(G);")
        
        s = gap.expect()
        
        output = parse_Irr(s)
    
        perm = argv.perm
        gen_latex(output, perm)
    
    else:
    
        print("USAGE:")
        print('./gap.py  cmd="G:=SmallGroup(28,1);"')


