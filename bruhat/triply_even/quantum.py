#!/usr/bin/env python3


from bruhat.solve import shortstr, find_kernel, array2, dot2
from bruhat.triply_even import build

for (d, idx, Hx) in build.items:
    if d!=10:
        continue
    print(shortstr(Hx))
    print()
    Hz = array2(list(find_kernel(Hx)))
    print(shortstr(Hz))
    assert dot2(Hx, Hz.transpose()).sum() == 0

    #print(shortstr(dot2(Hz, Hz.transpose())))

    break


