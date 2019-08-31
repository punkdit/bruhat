#!/usr/bin/env python3

import numpy

from qupy.ldpc.gallagher import classical_distance
from bruhat.solve import array2, shortstr, find_kernel


# See: http://www.st.hirosaki-u.ac.jp/~betsumi/triply-even/


def get_all():
    f = open("matrixform.txt")
    
    line = f.readline()
    
    #items = []
    count = 0
    dim = None
    idx = None
    for line in f:
        line = line.strip()
        #print(line)
    
        if line in "[ ], ] ];".split():
            continue
        if line.startswith("<"):
            rows = []
            flds = line[1:].split(", ")
            dim = int(flds[0])
            idx = int(flds[1])
            count += 1
            continue
    
        assert line.startswith("["), repr(line)
        row = line[1:49]
        assert len(row) == 48
        rows.append([ord(c)-48 for c in row])
        if ">" in line:
            G = array2(rows)
            #items.append(G)
            yield (dim, idx, G)
            rows = None
            dim = None
            idx = None
        
items = list(get_all())


def get(dim, idx):
    for (_dim, _idx, G) in items:
        if dim==_dim and idx==_idx:
            return G

if 0:
    for (dim, idx, G) in get_all():
        print(dim, idx)
        print(shortstr(G))
        print()


maximal = [
    (15, 1, "Three copies of RM(1,4)"),
    (14, 1, "Direct sum of RM(1,4) and doubling of d+16."),
    (13, 1, "Doubling of g24 of order of Aut is 1002795171840."),
    (13, 2, "Doubling of (d10e27)+ of order of Aut is 443925135360."),
    (13, 3, "Doubling of d+24 of order of Aut is 12054469961318400."),
    (13, 4, "Doubling of d2+12 of order of Aut is 4348654387200."),
    (13, 5, "Doubling of d6+4 of order of Aut is 36238786560."),
    (13, 6, "Doubling of d4+6 of order of Aut is 32614907904."),
    (13, 7, "Doubling of d3+8 of order of Aut is 173946175488."),
    (9, 1, "The code induced from the triangular graph."),
]


for (dim, idx, cmt) in maximal:
    G = get(dim, idx)
    print(cmt)
    print(shortstr(G))
    print()

print("The Miyamoto's moonshine code, which contains the MacDonald [48, 6, 24] code:")
G = get(7, 144)
print(shortstr(G))



def get_dw_ge_4():
    # where dual code has weight >= 4
    yield (7, 144)
    dim = 8 
    idxs = [ 129, 130, 131, 132, 133, 1078, 1079 ]
    for idx in idxs:
        yield (dim, idx)
    dim = 9 
    idxs = [ 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 782, 783, 784, 785, 786,
          787, 788, 789, 790, 791, 798, 799, 800, 801, 805, 806, 840, 841, 878,
          879, 1109, 1711, 1712, 1713, 1714, 1715, 1716, 1960 ]
    for idx in idxs:
        yield (dim, idx)
    dim = 10 
    idxs = [ 16, 17, 18, 19, 20, 21, 22, 276, 277, 278, 279, 280, 281, 282, 287,
           288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 301, 302, 303, 304,
           336, 337, 338, 340, 342, 343, 345, 346, 350, 351, 358, 381, 382, 383,
           384, 454, 455, 456, 457, 464, 465, 548, 549, 550, 551, 552, 553, 554,
           615, 645, 939, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993,
           994, 995, 996, 997, 998, 999, 1000, 1001, 1002, 1003, 1240, 1241, 1242,
           1243, 1244, 1245, 1246, 1247 ]
    for idx in idxs:
        yield (dim, idx)
    dim = 11 
    idxs = [ 6, 7, 52, 54, 55, 56, 57, 58, 59, 61, 62, 63, 73, 74, 75, 76, 77, 80,
           81, 88, 89, 99, 100, 101, 102, 124, 125, 126, 130, 131, 132, 133, 150,
          151, 152, 153, 154, 160, 161, 162, 163, 186, 188, 189, 190, 208, 210, 211,
          319, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358,
          359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 508, 509, 510, 511, 512,
          513, 514, 515, 516, 517, 518, 519, 520 ]
    for idx in idxs:
        yield (dim, idx)
    dim = 12 
    idxs = [ 3, 7, 8, 10, 11, 13, 14, 18, 19, 23, 25, 26, 29, 31, 32, 36, 37, 42, 43,
          45, 71, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 149, 150, 151, 152,
          153, 154, 155, 156, 157, 158, 159 ]
    for idx in idxs:
        yield (dim, idx)
    dim = 13 
    idxs = [ 1, 2, 3, 4, 5, 6, 7, 11, 13, 14, 15, 16, 32, 33, 34, 35, 36, 37, 38, 39 ]
    for idx in idxs:
        yield (dim, idx)
    dim = 14 
    idxs = [ 1, 4, 5, 6 ]
    for idx in idxs:
        yield (dim, idx)
    yield (15, 1)


if 0:
    for (dim, idx) in get_dw_ge_4():
        G = get(dim, idx)
        H = list(find_kernel(G))
        H = array2(H)
        H = H.astype(numpy.int32)
        print("[%d, %d, %d]" % (G.shape[1], G.shape[0], classical_distance(H)))




