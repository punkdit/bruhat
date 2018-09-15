#!/usr/bin/env python

# Gaussian elimination over natural numbers

from __future__ import print_function

import sys, os

import numpy

#from solve import shortstr
shortstr = str

from argv import argv


def _parse(fld):
    "find the (sum of) coefficient(s) in front of a (the) symbol(s)"
    if len(fld)==1:
        return 1

    if "+" in fld:
        flds = fld.split("+")
        flds = [_parse(fld) for fld in flds]
        return sum(flds)

    #print("_parse", fld)
    i = 0
    while i<len(fld) and fld[i] in "0123456789":
        i += 1
    #assert i+1==len(fld)
    if i==0:
        return 1
    n = int(fld[:i])
    return n


def parse(table, perm=None):
    lines = table.split("\n")
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    global header # Lame ! haha..
    header = lines[0][2:].split()
    print("header:", header)
    lines = lines[2:]
    lines = [line[3:] for line in lines if line]
    rows = []
    for line in lines:
        flds = line.split()
        flds = [_parse(fld) for fld in flds]
        rows.append(flds)
    A = numpy.array(rows)
    choose = argv.choose
    if choose is not None:
        items = choose.split(",")
        idxs = [header.index(col) for col in items]
        print("idxs:", idxs)
        A = A[:, idxs]
        A = A[idxs, :]
        header = items
    if perm:
        # permute the cols
        idxs = [header.index(col) for col in perm]
        A = A[:, idxs]
    return A


def latex_table(table, rows, cols, upper=False):
    lines = [] 
    m, n = len(rows), len(cols)
    lines.append(r"\begin{array}{%s}"%('r'*n))
    lines.append(r" %s \\" % (' & '.join(cols)))
    lines.append(r"\hline")
    for i in range(m):
        row = rows[i]
        line = [str(table[i, j]) for j in range(n)]
        line = r"%s \\" % (' & '.join(line))
        lines.append(line)
    lines.append(r"\end{array}")
    s = '\n'.join(lines)
    return s




def zelim(A, verbose=False):

    A = A.copy()
    n = len(A)
    assert A.shape == (n, n)
    
    row = 0
    
    while row+1 < n:
    
        #print A
        #print "row = %d, col = %d" % (row, col)
        #assert A[row, :col].sum() == 0
        while row<n and A[row].sum() == 0:
            row += 1
    
        if row==n:
            break

        col = 0
        while A[row, col] == 0:
            col += 1

        val0 = A[row, col]
        assert val0
        for row1 in range(row + 1, n):
            val1 = A[row1, col]
            if val1 == 0:
                continue
            assert val1 % val0 == 0
            r = val1 // val0
            A[row1] -= r*A[row]
            assert A.min() >= 0
        #print
        #print shortstr(A)
    
        row += 1
    
    
    A = numpy.array([row for row in A if row.sum()])
    m, n = A.shape
    assert n==len(header)

#    if verbose:
#        print header
#        print A
#        
#        print "$$"
#        print latex_table(A, ' '*m, header)
#        print "$$"

    assert A.min() >= 0

    return A


def rank(A):
    A = zelim(A)
    idx = 0
    while idx<len(A):
        if A[idx].sum() == 0:
            break
        idx += 1
    return idx


tables = {}

tables["D_8"] = """
  | A B  C  D  E    F    G  H  
--+----------------------------
A | A B  C  D  E    F    G  H  
B | B 2B G  G  2E   H    2G 2H 
C | C G  2C G  H    2F   2G 2H 
D | D G  G  2D H    H    2G 2H 
E | E 2E H  H  2E+H 2H   2H 4H 
F | F H  2F H  2H   2F+H 2H 4H 
G | G 2G 2G 2G 2H   2H   4G 4H 
H | H 2H 2H 2H 4H   4H   4H 8H 
""" # D_8

tables["D_12"] = """
  | A B  C  D  E   F  G     H     I  J   
--+--------------------------------------
A | A B  C  D  E   F  G     H     I  J   
B | B 2B F  F  G   2F 2G    J     J  2J  
C | C F  2C F  H   2F J     2H    J  2J  
D | D F  F  2D I   2F J     J     2I 2J  
E | E G  H  I  E+I J  G+J   H+J   3I 3J  
F | F 2F 2F 2F J   4F 2J    2J    2J 4J  
G | G 2G J  J  G+J 2J 2G+2J 3J    3J 6J  
H | H J  2H J  H+J 2J 3J    2H+2J 3J 6J  
I | I J  J  2I 3I  2J 3J    3J    6I 6J  
J | J 2J 2J 2J 3J  4J 6J    6J    6J 12J 
"""

tables["D_16"] = """
  | A B  C  D  E    F    G  H  I     J     K   
--+--------------------------------------------
A | A B  C  D  E    F    G  H  I     J     K   
B | B 2B G  G  2E   H    2G 2H 2I    K     2K  
C | C G  2C G  H    2F   2G 2H K     2J    2K  
D | D G  G  2D H    H    2G 2H K     K     2K  
E | E 2E H  H  2E+H 2H   2H 4H 2I+K  2K    4K  
F | F H  2F H  2H   2F+H 2H 4H 2K    2J+K  4K  
G | G 2G 2G 2G 2H   2H   4G 4H 2K    2K    4K  
H | H 2H 2H 2H 4H   4H   4H 8H 4K    4K    8K  
I | I 2I K  K  2I+K 2K   2K 4K 2I+3K 4K    8K  
J | J K  2J K  2K   2J+K 2K 4K 4K    2J+3K 8K  
K | K 2K 2K 2K 4K   4K   4K 8K 8K    8K    16K 
"""


tables["Q_8"] = """
  | A B  C  D  E  F  
--+------------------
A | A B  C  D  E  F  
B | B 2B E  E  2E 2F 
C | C E  2C E  2E 2F 
D | D E  E  2D 2E 2F 
E | E 2E 2E 2E 4E 4F 
F | F 2F 2F 2F 4F 8F 
"""

tables["S_3"] = """
  | A B  C   D  
--+-------------
A | A B  C   D  
B | B 2B D   2D 
C | C D  C+D 3D 
D | D 2D 3D  6D 
""" # S_3

tables["S_4"] = """
  | A B  C   D    E     F  G     H     I     J     K   
--+----------------------------------------------------
A | A B  C   D    E     F  G     H     I     J     K   
B | B 2B F   H    J     2F J     2H    K     2J    2K  
C | C F  C+F I    E+J   3F G+J   K     I+K   3J    3K  
D | D H  I   D+I  2I    K  K     H+K   2I+K  2K    4K  
E | E J  E+J 2I   2E+K  3J J+K   2K    2I+2K 2J+2K 6K  
F | F 2F 3F  K    3J    6F 3J    2K    3K    6J    6K  
G | G J  G+J K    J+K   3J 2G+K  2K    3K    2J+2K 6K  
H | H 2H K   H+K  2K    2K 2K    2H+2K 4K    4K    8K  
I | I K  I+K 2I+K 2I+2K 3K 3K    4K    2I+5K 6K    12K 
J | J 2J 3J  2K   2J+2K 6J 2J+2K 4K    6K    4J+4K 12K 
K | K 2K 3K  4K   6K    6K 6K    8K    12K   12K   24K 
"""

tables["PS_4"] = """
  | A  D    E     I     K    
--+-------------------------
A | A  D    E     I     K    
D | D  D+I  2I    2I+K  4K   
E | E  2I   2E+K  2I+2K 6K   
I | I  2I+K 2I+2K 2I+5K 12K  
K | K  4K   6K    12K   24K  
""" # parabolics of S_4

tables["PS_5"] = """
  | A B    C       D       E        F      G    
--+---------------------------------------------
A | A B    C       D       E        F      G    
B | B B+D  D+E     2D+F    E+2F     3F+G   5G   
C | C D+E  C+E+F   D+3F    2E+2F+G  4F+3G  10G  
D | D 2D+F D+3F    2D+4F+G 6F+2G    6F+7G  20G  
E | E E+2F 2E+2F+G 6F+2G   2E+4F+5G 6F+12G 30G  
F | F 3F+G 4F+3G   6F+7G   6F+12G   6F+27G 60G  
G | G 5G   10G     20G     30G      60G    120G 
""" # parabolics of S_5


tables["D_4"] = open("D_4.table").read()

tables["PD_4"] = """
  | A B     C     D     E           F        G        H        I        J      K    
--+---------------------------------------------------------------------------------
A | A B     C     D     E           F        G        H        I        J      K    
B | B 2B+G  2F    2F    2G+J        2F+2J    4G+K     4J       4J       4J+2K  8K   
C | C 2F    2C+H  2F    2H+J        2F+2J    4J       4H+K     4J       4J+2K  8K   
D | D 2F    2F    2D+I  2I+J        2F+2J    4J       4J       4I+K     4J+2K  8K   
E | E 2G+J  2H+J  2I+J  2E+G+H+I+2K 6J+K     4G+2J+4K 4H+2J+4K 4I+2J+4K 6J+9K  24K  
F | F 2F+2J 2F+2J 2F+2J 6J+K        2F+6J+2K 8J+4K    8J+4K    8J+4K    8J+12K 32K  
G | G 4G+K  4J    4J    4G+2J+4K    8J+4K    8G+10K   8J+8K    8J+8K    8J+20K 48K  
H | H 4J    4H+K  4J    4H+2J+4K    8J+4K    8J+8K    8H+10K   8J+8K    8J+20K 48K  
I | I 4J    4J    4I+K  4I+2J+4K    8J+4K    8J+8K    8J+8K    8I+10K   8J+20K 48K  
J | J 4J+2K 4J+2K 4J+2K 6J+9K       8J+12K   8J+20K   8J+20K   8J+20K   8J+44K 96K  
K | K 8K    8K    8K    24K         32K      48K      48K      48K      96K    192K 
""" # parabolics of Weyl D_4




tables["A_5"] = """
  | A B    C     D     E     F     G     H      I   
--+-------------------------------------------------
A | A B    C     D     E     F     G     H      I   
B | B B+G  H     G+H   I     F+I   2G+I  H+2I   5I  
C | C H    C+H   2H    E+I   3H    2I    2H+2I  6I  
D | D G+H  2H    D+H+I 2I    3H+I  G+3I  2H+4I  10I 
E | E I    E+I   2I    2E+2I 3I    4I    6I     12I 
F | F F+I  3H    3H+I  3I    3F+3I 5I    3H+6I  15I 
G | G 2G+I 2I    G+3I  4I    5I    2G+6I 10I    20I 
H | H H+2I 2H+2I 2H+4I 6I    3H+6I 10I   2H+14I 30I 
I | I 5I   6I    10I   12I   15I   20I   30I    60I 
""" # alternating n=5

tables["GL32"] = """
  | A B     C     D     E      F      G          H     I       J        K       L        M      N      P    
--+--------------------------------------------------------------------------------------------------------
A | A B     C     D     E      F      G          H     I       J        K       L        M      N      P    
B | B B+J   G+I   M     J+M    F+N    G+J+N      P     I+2N    3J+P     K+N+P   L+3N     M+2P   3N+2P  7P   
C | C G+I   C+L   M     E+N    L+M    G+L+N      P     I+2N    J+3N     K+N+P   3L+P     M+2P   3N+2P  7P   
D | D M     M     D+M   2M     2M     P          H+P   M+P     2P       2P      2P       2M+2P  4P     8P   
E | E J+M   E+N   2M    2E+P   2M+N   J+N+P      2P    M+2N+P  2J+3P    N+3P    3N+2P    2M+4P  2N+6P  14P  
F | F F+N   L+M   2M    2M+N   2F+P   L+N+P      2P    M+2N+P  3N+2P    N+3P    2L+3P    2M+4P  2N+6P  14P  
G | G G+J+N G+L+N P     J+N+P  L+N+P  G+J+L+2N+P 3P    5N+P    3J+3N+3P K+2N+4P 3L+3N+3P 7P     5N+8P  21P  
H | H P     P     H+P   2P     2P     3P         3H+3P 4P      6P       6P      6P       8P     12P    24P  
I | I I+2N  I+2N  M+P   M+2N+P M+2N+P 5N+P       4P    I+3N+3P 6N+4P    2N+6P   6N+4P    M+9P   4N+12P 28P  
J | J 3J+P  J+3N  2P    2J+3P  3N+2P  3J+3N+3P   6P    6N+4P   6J+9P    3N+9P   9N+6P    14P    6N+18P 42P  
K | K K+N+P K+N+P 2P    N+3P   N+3P   K+2N+4P    6P    2N+6P   3N+9P    2K+10P  3N+9P    14P    2N+20P 42P  
L | L L+3N  3L+P  2P    3N+2P  2L+3P  3L+3N+3P   6P    6N+4P   9N+6P    3N+9P   6L+9P    14P    6N+18P 42P  
M | M M+2P  M+2P  2M+2P 2M+4P  2M+4P  7P         8P    M+9P    14P      14P     14P      2M+18P 28P    56P  
N | N 3N+2P 3N+2P 4P    2N+6P  2N+6P  5N+8P      12P   4N+12P  6N+18P   2N+20P  6N+18P   28P    4N+40P 84P  
P | P 7P    7P    8P    14P    14P    21P        24P   28P     42P      42P     42P      56P    84P    168P 
"""

tables["B_2"] = """
  | A B  C  D  E    F    G  H  
--+----------------------------
A | A B  C  D  E    F    G  H  
B | B 2B G  G  2E   H    2G 2H 
C | C G  2C G  H    2F   2G 2H 
D | D G  G  2D H    H    2G 2H 
E | E 2E H  H  2E+H 2H   2H 4H 
F | F H  2F H  2H   2F+H 2H 4H 
G | G 2G 2G 2G 2H   2H   4G 4H 
H | H 2H 2H 2H 4H   4H   4H 8H 
"""

tables["PB_2"] = """
  | A B    C    D  
--+----------------
A | A B    C    D  
B | B 2B+D 2D   4D 
C | C 2D   2C+D 4D 
D | D 4D   4D   8D 
"""


tables["B_3"] = """
  | A B  C  D  E   F    G  H  I     J     K     L     M     N     P     Q     R     S     T     U      V     W   X     Y     Z     a       b     c     d      e     f   g      h   
--+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A | A B  C  D  E   F    G  H  I     J     K     L     M     N     P     Q     R     S     T     U      V     W   X     Y     Z     a       b     c     d      e     f   g      h   
B | B 2B G  G  H   R    2G 2H S     V     V     W     S     W     b     b     2R    2S    e     f      2V    2W  e     e     e     c       2b    2c    h      2e    2f  h      2h  
C | C G  2C G  N   Q    2G W  Z     Z     T     W     T     2N    b     2Q    b     e     2T    g      e     2W  e     e     2Z    g       2b    h     h      2e    h   2g     2h  
D | D G  G  2D L   P    2G W  X     Y     X     2L    Y     W     2P    b     b     e     e     d      e     2W  2X    2Y    e     d       2b    h     2d     2e    h   h      2h  
E | E H  N  L  E+H U    W  3H I+S   J+V   K+V   L+W   M+S   N+W   d     g     f     3S    T+e   U+f    3V    3W  X+e   Y+e   Z+e   a+c     h     3c    d+h    3e    3f  g+h    3h  
F | F R  Q  P  U   F+U  b  f  f     d     g     d     2U    g     P+d   Q+g   R+f   2f    2g    2U+f   h     h   h     2d    h     d+g     b+h   2h    2d+h   2h    4f  2g+h   4h  
G | G 2G 2G 2G W   b    4G 2W e     e     e     2W    e     2W    2b    2b    2b    2e    2e    h      2e    4W  2e    2e    2e    h       4b    2h    2h     4e    2h  2h     4h  
H | H 2H W  W  3H  f    2W 6H 3S    3V    3V    3W    3S    3W    h     h     2f    6S    3e    3f     6V    6W  3e    3e    3e    3c      2h    6c    3h     6e    6f  3h     6h  
I | I S  Z  X  I+S f    e  3S 2I+f  Z+c   X+c   X+e   S+f   Z+e   h     h     2f    2S+2f e+h   3f     2c+e  3e  2X+h  e+h   2Z+h  c+h     2h    2c+2h 3h     2e+2h 6f  3h     6h  
J | J V  Z  Y  J+V d    e  3V Z+c   2J+c  V+c   Y+e   Y+c   Z+e   2d    h     h     2c+e  e+h   d+h    2V+2c 3e  e+h   2Y+h  2Z+h  2c+d    2h    4c+h  2d+2h  2e+2h 3h  3h     6h  
K | K V  T  X  K+V g    e  3V X+c   V+c   2K+c  X+e   T+c   T+e   h     2g    h     2c+e  2T+h  g+h    2V+2c 3e  2X+h  e+h   e+h   2c+g    2h    4c+h  3h     2e+2h 3h  2g+2h  6h  
L | L W  W  2L L+W d    2W 3W X+e   Y+e   X+e   2L+2W Y+e   3W    2d    h     h     3e    3e    d+h    3e    6W  2X+2e 2Y+2e 3e    d+h     2h    3h    2d+2h  6e    3h  3h     6h  
M | M S  T  Y  M+S 2U   e  3S S+f   Y+c   T+c   Y+e   2M+f  T+e   2d    2g    2f    2S+2f 2T+h  2U+2f  2c+e  3e  e+h   2Y+h  e+h   2a+h    2h    2c+2h 2d+2h  2e+2h 6f  2g+2h  6h  
N | N W  2N W  N+W g    2W 3W Z+e   Z+e   T+e   3W    T+e   2N+2W h     2g    h     3e    2T+2e g+h    3e    6W  3e    3e    2Z+2e g+h     2h    3h    3h     6e    3h  2g+2h  6h  
P | P b  b  2P d   P+d  2b h  h     2d    h     2d    2d    h     2P+2d b+h   b+h   2h    2h    2d+h   2h    2h  2h    4d    2h    2d+h    2b+2h 4h    4d+2h  4h    4h  4h     8h  
Q | Q b  2Q b  g   Q+g  2b h  h     h     2g    h     2g    2g    b+h   2Q+2g b+h   2h    4g    2g+h   2h    2h  2h    2h    2h    2g+h    2b+2h 4h    4h     4h    4h  4g+2h  8h  
R | R 2R b  b  f   R+f  2b 2f 2f    h     h     h     2f    h     b+h   b+h   2R+2f 4f    2h    4f     2h    2h  2h    2h    2h    2h      2b+2h 4h    4h     4h    8f  4h     8h  
S | S 2S e  e  3S  2f   2e 6S 2S+2f 2c+e  2c+e  3e    2S+2f 3e    2h    2h    4f    4S+4f 2e+2h 6f     4c+2e 6e  2e+2h 2e+2h 2e+2h 2c+2h   4h    4c+4h 6h     4e+4h 12f 6h     12h 
T | T e  2T e  T+e 2g   2e 3e e+h   e+h   2T+h  3e    2T+h  2T+2e 2h    4g    2h    2e+2h 4T+2h 2g+2h  2e+2h 6e  2e+2h 2e+2h 2e+2h 2g+2h   4h    6h    6h     4e+4h 6h  4g+4h  12h 
U | U f  g  d  U+f 2U+f h  3f 3f    d+h   g+h   d+h   2U+2f g+h   2d+h  2g+h  4f    6f    2g+2h 2U+5f  3h    3h  3h    2d+2h 3h    d+g+2h  4h    6h    2d+5h  6h    12f 2g+5h  12h 
V | V 2V e  e  3V  h    2e 6V 2c+e  2V+2c 2V+2c 3e    2c+e  3e    2h    2h    2h    4c+2e 2e+2h 3h     4V+4c 6e  2e+2h 2e+2h 2e+2h 4c+h    4h    8c+2h 6h     4e+4h 6h  6h     12h 
W | W 2W 2W 2W 3W  h    4W 6W 3e    3e    3e    6W    3e    6W    2h    2h    2h    6e    6e    3h     6e    12W 6e    6e    6e    3h      4h    6h    6h     12e   6h  6h     12h 
X | X e  e  2X X+e h    2e 3e 2X+h  e+h   2X+h  2X+2e e+h   3e    2h    2h    2h    2e+2h 2e+2h 3h     2e+2h 6e  4X+2h 2e+2h 2e+2h 3h      4h    6h    6h     4e+4h 6h  6h     12h 
Y | Y e  e  2Y Y+e 2d   2e 3e e+h   2Y+h  e+h   2Y+2e 2Y+h  3e    4d    2h    2h    2e+2h 2e+2h 2d+2h  2e+2h 6e  2e+2h 4Y+2h 2e+2h 2d+2h   4h    6h    4d+4h  4e+4h 6h  6h     12h 
Z | Z e  2Z e  Z+e h    2e 3e 2Z+h  2Z+h  e+h   3e    e+h   2Z+2e 2h    2h    2h    2e+2h 2e+2h 3h     2e+2h 6e  2e+2h 2e+2h 4Z+2h 3h      4h    6h    6h     4e+4h 6h  6h     12h 
a | a c  g  d  a+c d+g  h  3c c+h   2c+d  2c+g  d+h   2a+h  g+h   2d+h  2g+h  2h    2c+2h 2g+2h d+g+2h 4c+h  3h  3h    2d+2h 3h    2a+c+2h 4h    4c+4h 2d+5h  6h    6h  2g+5h  12h 
b | b 2b 2b 2b h   b+h  4b 2h 2h    2h    2h    2h    2h    2h    2b+2h 2b+2h 2b+2h 4h    4h    4h     4h    4h  4h    4h    4h    4h      4b+4h 8h    8h     8h    8h  8h     16h 
c | c 2c h  h  3c  2h   2h 6c 2c+2h 4c+h  4c+h  3h    2c+2h 3h    4h    4h    4h    4c+4h 6h    6h     8c+2h 6h  6h    6h    6h    4c+4h   8h    8c+8h 12h    12h   12h 12h    24h 
d | d h  h  2d d+h 2d+h 2h 3h 3h    2d+2h 3h    2d+2h 2d+2h 3h    4d+2h 4h    4h    6h    6h    2d+5h  6h    6h  6h    4d+4h 6h    2d+5h   8h    12h   4d+10h 12h   12h 12h    24h 
e | e 2e 2e 2e 3e  2h   4e 6e 2e+2h 2e+2h 2e+2h 6e    2e+2h 6e    4h    4h    4h    4e+4h 4e+4h 6h     4e+4h 12e 4e+4h 4e+4h 4e+4h 6h      8h    12h   12h    8e+8h 12h 12h    24h 
f | f 2f h  h  3f  4f   2h 6f 6f    3h    3h    3h    6f    3h    4h    4h    8f    12f   6h    12f    6h    6h  6h    6h    6h    6h      8h    12h   12h    12h   24f 12h    24h 
g | g h  2g h  g+h 2g+h 2h 3h 3h    3h    2g+2h 3h    2g+2h 2g+2h 4h    4g+2h 4h    6h    4g+4h 2g+5h  6h    6h  6h    6h    6h    2g+5h   8h    12h   12h    12h   12h 4g+10h 24h 
h | h 2h 2h 2h 3h  4h   4h 6h 6h    6h    6h    6h    6h    6h    8h    8h    8h    12h   12h   12h    12h   12h 12h   12h   12h   12h     16h   24h   24h    24h   24h 24h    48h 
"""

tables["PB_3"] = """
  | A B     C     D       E     F      G   
--+----------------------------------------
A | A B     C     D       E     F      G   
B | B 2B+E  2F    2E+F    4E+G  2F+2G  6G  
C | C 2F    2C+2F 2F+G    4G    4F+2G  8G  
D | D 2E+F  2F+G  2D+E+2G 4E+4G 2F+5G  12G 
E | E 4E+G  4G    4E+4G   8E+8G 12G    24G 
F | F 2F+2G 4F+2G 2F+5G   12G   4F+10G 24G 
G | G 6G    8G    12G     24G   24G    48G 
"""

# binary dicyclic group
tables["Dic3"] = """
  | A B  C   D  E  F   
--+--------------------
A | A B  C   D  E  F   
B | B 2B E   2D 2E 2F  
C | C E  C+E F  3E 3F  
D | D 2D F   4D 2F 4F  
E | E 2E 3E  2F 6E 6F  
F | F 2F 3F  4F 6F 12F 
"""

# binary dicyclic group
tables["Dic4"] = """
  | A B  C  D  E    F    G  H  I   
--+--------------------------------
A | A B  C  D  E    F    G  H  I   
B | B 2B G  G  2E   H    2G 2H 2I  
C | C G  2C G  H    2F   2G 2H 2I  
D | D G  G  2D H    H    2G 2H 2I  
E | E 2E H  H  2E+H 2H   2H 4H 4I  
F | F H  2F H  2H   2F+H 2H 4H 4I  
G | G 2G 2G 2G 2H   2H   4G 4H 4I  
H | H 2H 2H 2H 4H   4H   4H 8H 8I  
I | I 2I 2I 2I 4I   4I   4I 8I 16I 
"""

# binary tetrahedral group, aka SL(2, 3)
tables["2T"] = """
  | A B  C   D     E     F   G   
--+------------------------------
A | A B  C   D     E     F   G   
B | B 3B F   3D    G     3F  3G  
C | C F  C+F 2F    E+G   4F  4G  
D | D 3D 2F  2D+2F 2G    6F  6G  
E | E G  E+G 2G    2E+2G 4G  8G  
F | F 3F 4F  6F    4G    12F 12G 
G | G 3G 4G  6G    8G    12G 24G 
"""

# binary octahedral group
tables["2O"] = """
  | A B  C   D    E     F     G  H     I     J     K     L     M     N   P      Q   
--+---------------------------------------------------------------------------------
A | A B  C   D    E     F     G  H     I     J     K     L     M     N   P      Q   
B | B 2B G   H    L     L     2G 2H    M     M     N     2L    2M    2N  Q      2Q  
C | C G  C+G K    E+L   F+L   3G N     P     P     K+N   3L    Q     3N  P+Q    3Q  
D | D H  K   D+K  2K    N     N  H+N   I+P   J+P   2K+N  2N    M+Q   4N  2P+Q   4Q  
E | E L  E+L 2K   2E+N  L+N   3L 2N    2P    2P    2K+2N 2L+2N 2Q    6N  2P+2Q  6Q  
F | F L  F+L N    L+N   2F+N  3L 2N    Q     Q     3N    2L+2N 2Q    6N  3Q     6Q  
G | G 2G 3G  N    3L    3L    6G 2N    Q     Q     3N    6L    2Q    6N  3Q     6Q  
H | H 2H N   H+N  2N    2N    2N 2H+2N M+Q   M+Q   4N    4N    2M+2Q 8N  4Q     8Q  
I | I M  P   I+P  2P    Q     Q  M+Q   2I+Q  M+2P  2P+Q  2Q    2M+2Q 4Q  2P+3Q  8Q  
J | J M  P   J+P  2P    Q     Q  M+Q   M+2P  2J+Q  2P+Q  2Q    2M+2Q 4Q  2P+3Q  8Q  
K | K N  K+N 2K+N 2K+2N 3N    3N 4N    2P+Q  2P+Q  2K+5N 6N    4Q    12N 2P+5Q  12Q 
L | L 2L 3L  2N   2L+2N 2L+2N 6L 4N    2Q    2Q    6N    4L+4N 4Q    12N 6Q     12Q 
M | M 2M Q   M+Q  2Q    2Q    2Q 2M+2Q 2M+2Q 2M+2Q 4Q    4Q    4M+4Q 8Q  8Q     16Q 
N | N 2N 3N  4N   6N    6N    6N 8N    4Q    4Q    12N   12N   8Q    24N 12Q    24Q 
P | P Q  P+Q 2P+Q 2P+2Q 3Q    3Q 4Q    2P+3Q 2P+3Q 2P+5Q 6Q    8Q    12Q 2P+11Q 24Q 
Q | Q 2Q 3Q  4Q   6Q    6Q    6Q 8Q    8Q    8Q    12Q   12Q   16Q   24Q 24Q    48Q 
"""

# binary icosahedral group
tables["2I"] = """
  | A B    C     D     E     F     G     H     I      J      K   L    
--+-------------------------------------------------------------------
A | A B    C     D     E     F     G     H     I      J      K   L    
B | B B+G  I     G+I   K     F+K   2G+K  L     I+2K   2J+L   5K  5L   
C | C I    C+I   2I    E+K   3I    2K    H+L   2I+2K  2L     6K  6L   
D | D G+I  2I    D+I+K 2K    3I+K  G+3K  2L    2I+4K  J+3L   10K 10L  
E | E K    E+K   2K    2E+2K 3K    4K    2H+2L 6K     4L     12K 12L  
F | F F+K  3I    3I+K  3K    3F+3K 5K    3L    3I+6K  5L     15K 15L  
G | G 2G+K 2K    G+3K  4K    5K    2G+6K 4L    10K    2J+6L  20K 20L  
H | H L    H+L   2L    2H+2L 3L    4L    4H+4L 6L     8L     12L 24L  
I | I I+2K 2I+2K 2I+4K 6K    3I+6K 10K   6L    2I+14K 10L    30K 30L  
J | J 2J+L 2L    J+3L  4L    5L    2J+6L 8L    10L    4J+12L 20L 40L  
K | K 5K   6K    10K   12K   15K   20K   12L   30K    20L    60K 60L  
L | L 5L   6L    10L   12L   15L   20L   24L   30L    40L    60L 120L 
"""

def process(table):

    A = parse(table, perm=argv.perm)
    print(A)
    if A.shape[0]==A.shape[1]:
        print("det:", numpy.linalg.det(A))
    print()

    order = A[-1, -1]
    A = zelim(A)
    print(shortstr(A))

    reps = A[:, -1]
    print("order:", order)
    print("dimension of reps:", reps)
    print("sum of squares:", (reps*reps).sum())
    


if __name__ == "__main__":
    table = tables[argv.next()]
    process(table)


