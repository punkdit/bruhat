#!/usr/bin/env python3
"""
q-deformed Pascal triangles
"""

from bruhat.element import Z
from bruhat.poly import Poly
from bruhat.argv import argv

ring = Z
zero = Poly({}, ring)
one = Poly({():1}, ring)
q = Poly("q", ring)


def sl_pascal(row, col):
    assert 0 <= row
    assert 0 <= col <= row

    if col == 0 or col == row:
        return one

    left = sl_pascal(row-1, col-1)
    right = sl_pascal(row-1, col)
    p = q**(row - col)*left + right
    #print("sl_pascal: "
    return p
        
    
def sp_pascal(row, col):
    assert 0 <= row
    assert 0 <= col <= row

    if col == 0:
        return one

    left = sp_pascal(row-1, col-1)
    if col == row:
        right = 0
    else:
        right = sp_pascal(row-1, col)

    p = q**(2*col)*right + (q**col + 1)*left
    return p


def main():

    fn = argv.next()
    if fn is not None:
        fn = eval(fn)
    else:
        fn = sl_pascal

    value = argv.get("q")
    for row in range(7):
      for col in range(row+1):
        p = fn(row, col)
        if value is not None:
            p = p.substitute((("q", value),))
            s = str(p).replace(" + ", "+")
        else:
            s = p.qstr("q")
            #s = p.flatstr().replace(" + ", "+")
        print(s, end=" ")
      print()



if __name__ == "__main__":

    main()

