#!/usr/bin/env python3

"""
parse character table output from gap.
"""

import sys, os
from time import sleep
import select
from subprocess import Popen, PIPE

from argv import argv

class Gap(object):
    def __init__(self):
        self.proc = Popen("gap", bufsize=0, stdin=PIPE, stdout=PIPE)
        sleep(0.1)
        #print(dir(proc.stdout))
        #print(proc.stdout.read(20))
        self.buf = ""
        self.pos = 0

    def read_nb(self):
        proc = self.proc
        data = bytes()
        while select.select([proc.stdout],[],[],0)[0]!=[]:   
            data += proc.stdout.read(1)
        data = data.decode("utf-8")
        self.buf += data

    def send(self, data):
        proc = self.proc
        data = data + "\n"
        data = data.encode("utf-8")
        proc.stdin.write(data)

    def expect(self, s):
        while s not in self.buf[self.pos:]:
            self.read_nb()
            sleep(0.1)

        pos = self.pos + self.buf[self.pos:].index(s)
        data = self.buf[self.pos : pos]
        self.pos = pos + len(s)
        return data

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


cmd = argv.cmd
if cmd:
    
    gap = Gap()
    PROMPT = "gap> "
    gap.expect(PROMPT)
    
    #gap.send("G:=SmallGroup(28,1);")

    if cmd[-1] != ";":
        cmd += ";"
    gap.send(cmd)
    
    s = gap.expect(PROMPT)
    
    gap.send("Irr(G);")
    
    s = gap.expect(PROMPT)
    
    output = parse_Irr(s)

    perm = argv.perm
    gen_latex(output, perm)

else:

    print("USAGE:")
    print('./gap.py  cmd="G:=SmallGroup(28,1);"')


