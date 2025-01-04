#!/usr/bin/env python

"""
Extract table of marks from gap & do computations in the Burnside ring.

See also: relation.py where we explicitly build TOM's.
"""

import sys
from os import popen
from string import ascii_lowercase, ascii_uppercase
letters = ascii_uppercase+ascii_lowercase

import numpy
numpy.set_printoptions(threshold=sys.maxsize)

from bruhat.argv import argv
from bruhat.smap import SMap

all_names = [
    "A5", "A6", "A7", "A8", "A9", "A10", "A11", "J1", "J2",
    "J2.2", "J3", "J3.2", "L2(11)", "L2(19)", "L2(17)", "L2(23)",
    "L2(7)", "L2(8)", "L2(13)", "L2(16)", "L2(25)", "L2(27)",
    "L2(29)", "L2(31)", "L2(37)", "L2(32)", "L2(41)", "L2(43)",
    "L2(47)", "L2(49)", "L2(53)", "L2(59)", "L2(61)", "L2(67)",
    "L2(71)", "L2(73)", "L2(79)", "L2(64)", "L2(81)", "L2(83)",
    "L2(89)", "L2(97)", "L2(101)", "L2(103)", "L2(107)",
    "L2(109)", "L2(113)", "L2(121)", "L2(125)", "2^3:7:3",
    "2xA5", "19:6", "11:10", "D6xD10", "7:6", "3.PGL2(9)",
    "2^(1+4)-:A5", "2^(2+4):(3xS3)", "A4xA5", "A5xD10", "5^2:D12",
    "3.A6.2^2", "2^(1+4)-.S5", "2^(2+4):(S3xS3)", "(A4xA5):2",
    "(A5xD10).2", "L3(2).2x2", "5^2:(4xS3)", "2^4:(3xA5)",
    "(3xA6):2_2", "3^2.3^(1+2):8", "2^4:(3xA5).2", "L2(17)x2",
    "(3xM10):2", "3^2.3^(1+2):8.2", "19:18", "U3(3)", "U4(3)",
    "U3(5)", "U3(4)", "U3(11)", "U3(7)", "U3(8)", "U3(9)",
    "U4(2)", "U5(2)", "L2(7).2", "L2(16).2", "L2(16).4",
    "L2(11).2", "L2(8).3", "L2(13).2", "L2(25).2_2", "L2(32).5",
    "S5", "S8", "S6", "S7", "S9", "S10", "S4", "M11", "M12",
    "M12.2", "M22", "M22.2", "M23", "M24", "A6.2_3", "A6.2^2",
    "2.A8", "A6.2_2", "2.A5", "2.A6", "2.A7", "M9:2", "M8:S3",
    "4^2:D12.2", "S4xS3", "M8.S4", "4^2:D12", "A4xS3", "M9:S3",
    "2xS5", "(2^2xA5):2", "M8.S4.2", "3^(1+2)+:D8", "2^3:L3(2)",
    "2^3:L3(2)x2", "23:11", "L3(4)", "L3(3)", "L3(5)", "L3(7)",
    "L3(9)", "L4(3)", "L3(8)", "L3(11)", "2^4:A6", "2^4:A7",
    "2^4:A8", "2^4:A5", "2^4:A5`", "2^4:S5`", "2^4:S6", "2^5:S5",
    "2^4.S6", "2^4:S5", "L3(4).2_2", "L3(4).3.2_2", "L3(4).2_1",
    "L3(4).2^2", "L3(3).2", "L3(4).3", "L3(4).2_3", "L3(4).3.2_3",
    "L3(4).6", "L3(4).D12", "2^6:3.S6", "2^6:(L3(2)xS3)",
    "U3(5).2", "U3(3).2", "U4(3).2_3", "U3(4).2", "U3(5).3",
    "U3(5).S3", "U4(3).2_1", "U4(2).2", "4^3:L3(2)", "4.2^4:S5",
    "5:4xA5", "3^4:(M10x2)", "3^4:M10", "5^(1+2)+:3:8", "3^(1+4)+:4S5",
    "M11x2", "5^(1+2)+:3:8.2", "3^(1+4)+:2S5", "HS", "McL",
    "McL.2", "2.S8", "2.S5", "2.S6", "2.S5'", "2.S7", "He",
    "2^(1+6)+.L3(2)", "7^2:2L2(7)", "3.S7", "7^(1+2)+:(S3x3)",
    "S4xL3(2)", "7:3xL3(2)", "5^2:4A4", "S4(4).2", "2^2.L3(4).3.2_2",
    "(A7x3):2", "(A5xA5):4", "(A6xA4):2", "(A8x3):2", "(A7xA4):2",
    "(A6xA5):2", "(A6x3):2", "3^3:S4", "3^2:2A4", "2^4:(S3xS3)",
    "(A5x3):2", "S11", "S8x2", "S7xS3", "(S5xS5):2", "S6xS4",
    "2^5:S5`", "S7x2", "S6xS3", "S5xS4", "3^3:(2xS4)", "(S4xS4):2",
    "S6x2", "2^4:S4", "109:54", "(5x11).2", "27:2^2", "113:56",
    "(3x19).2", "28.2^2", "11^2:60", "D122", "30.2^2", "5^3:62",
    "63:2", "D62x2", "5^2:GL2(5)", "4^2:S3", "31:3", "7^2:2.L2(7):2",
    "(3xA4):2", "3^2:Q8", "19:3", "3^4:GL2(9)", "8^2:S3",
    "91:3", "31:5", "3^3:L3(3)", "3^4:2(A4xA4).2", "(4xA6):2",
    "S4xS4", "2^6:(7xL2(8))", "7^2:S3", "11^2:(5x2L2(11).2)",
    "133:3", "L5(2)", "S9x2", "S8xS3", "S7xS4", "S6xS5",
    "3.L3(4)", "3.L3(4).2_3", "3^2:2S4x2", "S3xS5", "7:6xS3",
    "S4(4)", "S4(5)", "S6(2)", "2^6:(3xA5)", "(A5xA5):2",
    "5^(1+2)+:4A5", "5^3:(2xA5).2", "2.(A5xA5).2", "2^6:L3(2)",
    "2^6:(3xA5):2", "(A5xA5):2^2", "2^5:S6", "2.[2^6]:(S3xS3)",
    "2.S6(2)", "2.(2^5:S6)", "2.2^6.L3(2)", "2.(S3xS6)",
    "2^2.[2^6].(S3xS3)", "2.U4(2).2", "2.U4(2)", "Sz(8)",
    "2F4(2)'", "Sz(32)", "2^(3+3):7", "13:4", "5:4", "D14",
    "2.[2^8]:5:4", "2^2.[2^8]:S3", "25:4", "D62", "2^(5+5):31",
    "41:4", "S3xL2(8)", "(7xL2(7)):2", "3^(1+2)+.2S4", "7^2:2A4",
    "2^(2+4):15", "5xA5", "5^2:S3", "13:3", "5^(1+2)+:8",
    "11^(1+2)+:40", "2(L2(11)x2).2", "(4^2x3):S3", "37:3",
    "7^(1+2):48", "2(L2(7)x4).2", "43:3", "2^(3+6):21", "3xL2(8)",
    "(9x3).S3", "3^(2+4):80", "5x2.A6.2_2", "10^2:S3", "73:3",
    "3^(1+2)+:2A4", "3^3:S4`", "2.(A4xA4).2", "2^(1+6)-:3^(1+2)+:2A4",
    "3xU4(2)", "2^(4+4):(3xA5)", "3^4:S5", "S3x3^(1+2)+:2A4",
    "3^4:A6", "3^(1+4)+.2S4", "2(A4xA4).4", "2^(2+4):(3xD10)",
    "13:6", "5^(1+2)+:24", "3x2S5", "6^2:S3", "7:3x3", "5^(1+2)+:8:2",
    "2S5.2", "(3x2S5).2", "6^2:D12", "(7:3x3):2", "3^(1+2)+:2S4",
    "2.(A4xA4).2.2", "3^4:(2xA6)", "3^(1+4)+.4S4", "U3(3)x2",
    "4(A4xA4).4", "3^(1+4)+.(2S4x2)", "2(A4xA4).4.2", "M10x2",
    "(4^2x2)S4", "G2(3)", "(3^(1+2)+x3^2):2S4", "2^3.L3(2)`",
    "2^(1+4)+:3^2.2", "A5xA5", "G2(3).2", "3^2.(3x3^(1+2)+):D8",
    "L2(8):3x2", "2^3.L3(2):2", "2^(1+4)+.(S3xS3)", "2.(2^4:S5)",
    "2.(S6x2)", "2.(3^(1+2)+:2S4)", "2.(3^3:(S4x2))", "2.(2.(A4xA4).2.2)",
    "2^2.L3(4).2_2", "2^2.L3(4)", "2^2.L3(4).3", "2^2.(2^4S5)",
    "2^2.S6", "2^2.(L2(7)x2)", "2^2.(3^2:Q8.2)", "2^2.(2^4:A5)",
    "2^2xA6", "2^2xL2(7)", "2^2.(3^2:Q8)", "2^2.(3^2:2A4)",
    "A4x7:3", "2^2.(3^2:2S4)", "2^2.(7:3xS3)", "U4(3).2^2_133",
    "3^(1+4)+.2^(1+4)-.S3", "U3(3):2x2", "2(A4xA4).4.2^2",
    "2xA6.2^2", "(4^2x2)(2xS4)", "(A9x3):2", "(A7xA5):2",
    "3^4:2^3.S4", "(A6xA6):2^2", "(A8xA4):2", "2^6:3^3:S4",
    "A12", "G2(4)", "3D4(2)", "2^(1+8)+:L2(8)", "Co3", "2^4`A8",
    "2xM12", "S3xL2(8):3", "A4xS5", "O8-(2)", "O8+(2)", "2^6:U4(2)",
    "(S3xS3xA5):2", "(3xU4(2)):2", "3^4:2^3.S4(a)", "2^6:A8",
    "S12", "A13", "HS.2", "He.2", "S13", "Sz(8).3", "2F4(2)" ]


def load_tom(name="M11", gapstr=None):
    assert name in all_names or gapstr

    if gapstr is None:
        cmd = """
        SizeScreen([1000,1000]);
        tom := TableOfMarks("%s");
        Display(tom);
        quit;
        """%name
    else:
        cmd = """
        SizeScreen([1000,1000]);
        G := %s;
        tom := TableOfMarks(G);
        Display(tom);
        quit;
        """%gapstr
    
    f = open("/tmp/tmp.gap", "w")
    print(cmd, file=f)
    f.close()
    
    #lines = popen("gap --norepl /tmp/tmp.gap").readlines()
    data = popen("gap --norepl /tmp/tmp.gap").read()
    #print(data)
    open("gapdump.out", "w").write(data)
    lines = data.split("\n")
    #print(data)
    #print(lines)
    
    rows = []
    start = False
    end = False
    for line in lines:
        line = line.strip()
        if line.startswith("1:"):
            start = True
        if not start:
            continue
        if not line:
            end = True
            continue
        #assert not end, "fix fix fix: see gapdump.out to parse"
        idx, data = line.split(":")
        idx = int(idx)-1
        assert 0 <= idx <= len(rows)
        data = data.replace(".", "0")
        data = data.split()
        flds = [int(i) for i in data]
        #print(idx, flds)
        if idx == len(rows):
            rows.append(flds)
        else:
            rows[idx] += flds
    
    debug = [len(item) for item in rows], len(rows)
    #print(rows)
    assert len(rows[-1]) == len(rows), debug
    
    N = len(rows)
    for item in rows:
        item += [0]*(N-len(item))
    
    lrev = lambda items:list(reversed(items))
    
    rows = [lrev(item) for item in reversed(rows)]
    rows = numpy.array(rows)
    return Tom(rows)
    

def get_letter(i, post=0):
    if i >= len(letters):
        i -= len(letters)
        post += 1
        return get_letter(i, post)
    else:
        return letters[i]+("_"+str(post) if post else "")


class Tom:
    def __init__(self, rows, names=None):
        self.rows = rows
        N = len(rows)
        if names is None:
            names = [get_letter(i) for i in range(N)]
        self.names = list(names)

    def str(self):
        lines = []
        for item in self.rows:
            line = " ".join((str(i) if i else ".") for i in item)
            lines.append(line)
        return '\n'.join(lines)

    def __getitem__(self, key):
        i = self.names.index(key)
        return self.rows[i]

    def __len__(self):
        return len(self.rows)

    def __str__(self):
        rows = self.rows
        names = self.names
        N = len(rows)
        smap = SMap()
        for j in range(N):
            smap[0,3*j] = "%3s"%names[j]
            smap[j+1,3*N+2] = names[j]
        for i,row in enumerate(rows):
            for j,c in enumerate(row):
                if j<i:
                    continue
                smap[i+1,3*j] = "%3s"%(c if c else ".")
        return str(smap)

    def get_desc(self, items):
        rows = self.rows
        names = self.names
        N = len(rows)
        assert len(items) == N
        items = items.copy()
        desc = []
        for i in range(N):
            assert items.min() >= 0
            if items[i] == 0:
                continue
            #print(items, "i=", i)
            n = items[i]
            row = rows[i]
            assert n % row[i] == 0
            j = n//row[i]
            items -= j*row
            s = "%d*%s"%(j, names[i]) if j>1 else names[i]
            desc.append(s)
        return "+".join(desc)

    def dump_maximal(self):
        rows = self.rows
        N = len(rows)
        print("maximal subgroups / minimal structures:")
        for i,row in enumerate(rows):
            col = rows[:,i]
            found = [j for j in col if j]
            if len(found)==2:
                print(self.names[i], row[N-1])

    def dump_burnside(self):
        N = len(self)
        names = self.names
        for a in names:
          for b in names:
            A = self[a]
            B = self[b]
            AB = A*B
            #print(a,"*", b)
            #print(AB)
            name = self.get_desc(AB)
            print("%6s"%name, end=" ")
          print()

    
    def dump_desc(self, row, desc):
        N = len(self)
        names = desc.split("+")
        div = row[N-1]
        total = 0
        for name in names:
            mul = 1
            if "*" in name:
                mul, name = name.split("*")
                mul = int(mul)
            row = self[name]
            count = row[N-1]
            print(name.ljust(6), mul, "\t*", count//div, "=", mul*count//div)
            total += mul * count // div
            assert count%div == 0
        print()
        print(total)


def test():

    tom = load_tom("A5")
    B = tom["B"]
    C = tom["C"]

    BC = B*C
    #assert tom.get_desc(BC) == "H"

    BG = tom["B"]*tom["G"]
    assert tom.get_desc(BG) == "2*G+I"

    tom.dump_burnside()
    #return

    tom = load_tom("A8")
    N = len(tom)

    rows = tom.rows
    #print(rows[:, N-1])
    tom.dump_maximal()



def test_M24():
    tom = load_tom("M24")
    N = len(tom)

    tom.dump_maximal()

    rows = tom.rows

    octad = tom.rows[4]
    desc = tom.get_desc(octad*octad)
    print("octad*octad =", desc)
    tom.dump_desc(octad, desc)

    line = tom["k"]
    desc = tom.get_desc(line*line)
    print("line*line =", desc)
    tom.dump_desc(line, desc)

    octern = tom["s_8"]
    desc = tom.get_desc(octern*octern)
    print(desc)
    tom.dump_desc(octern, desc)

    tom = load_tom(gapstr="GL(3,2)")
    print("tom:", len(tom))
    tom.dump_maximal()


def main():

    #test()

    #tom = load_tom("L2(7)")
    #tom = load_tom("M12")
    tom = load_tom("U3(3)");
    print(tom)

    E = tom["E"]
    F = tom["F"]
    print("E*E", tom.get_desc(E*E))
    print("F*F", tom.get_desc(F*F))

    tom.dump_maximal()

    #tom = load_tom("Co3")
    #print(len(tom))
    #tom.dump_maximal()
    #return
    #tom = load_tom("He")
    #print(len(tom))
    #return







if __name__ == "__main__":
    from time import time
    start_time = time()

    fn = argv.next() or "main"

    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)
    else:
        fn = eval(fn)
        fn()

    print("finished in %.3f seconds.\n"%(time() - start_time))

    

    
