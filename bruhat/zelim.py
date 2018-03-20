#!/usr/bin/env python

# Gaussian elimination over integers 

import sys, os

import numpy

from solve import shortstr
from argv import argv


def _parse(fld):
    if len(fld)==1:
        return 1
    if "+" in fld:
        flds = fld.split("+")
        flds = [_parse(fld) for fld in flds]
        return sum(flds)

    i = 0
    while i<len(fld) and fld[i] in "0123456789":
        i += 1
    assert i+1==len(fld)
    n = int(fld[:i])
    return n


def parse(table, perm=None):
    lines = table.split("\n")
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    global header # Lame ! haha..
    header = lines[0][2:].split()
    lines = lines[2:]
    lines = [line[3:] for line in lines if line]
    rows = []
    for line in lines:
        flds = line.split()
        flds = [_parse(fld) for fld in flds]
        rows.append(flds)
    A = numpy.array(rows)
    if perm:
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
    n = len(A)
    assert A.shape == (n, n)
    
    row = 0
    col = 0
    
    while row+1 < n:
    
        #print A
        #print "row = %d, col = %d" % (row, col)
        assert A[row, :col].sum() == 0
    
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
        print
        print shortstr(A)
    
        row += 1
        col += 1
    
    
    A = numpy.array([row for row in A if row.sum()])
    m, n = A.shape
    assert n==len(header)

    if verbose:
        print header
        print A
        
        print "$$"
        print latex_table(A, ' '*m, header)
        print "$$"

    assert A.min() >= 0

    return A



table = """
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
A = parse(table, perm=argv.get("perm"))
print shortstr(A)
A = zelim(A)
print
print shortstr(A)


table = """
  | A B  C   D  
--+-------------
A | A B  C   D  
B | B 2B D   2D 
C | C D  C+D 3D 
D | D 2D 3D  6D 
""" # S_3

A = parse(table)
A = zelim(A)

A = parse("""
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
""") # S_4
A = zelim(A)

A = parse("""
  | A  D    E     I     K    
--+-------------------------
A | A  D    E     I     K    
D | D  D+I  2I    2I+K  4K   
E | E  2I   2E+K  2I+2K 6K   
I | I  2I+K 2I+2K 2I+5K 12K  
K | K  4K   6K    12K   24K  
""") # parabolics of S_4
A = zelim(A)

A = parse("""
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
""") # alternating n=5
A = zelim(A)

A = parse("""
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
""")
A = zelim(A)

