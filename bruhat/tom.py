#!/usr/bin/env python

from os import popen

import numpy


def load_tom(name="M11"):
    cmd = """
    SizeScreen([20000,20000]);
    tom := TableOfMarks("M11");
    Display(tom);
    quit;
    """
    
    f = open("/tmp/tmp.gap", "w")
    print(cmd, file=f)
    f.close()
    
    lines = popen("gap --norepl /tmp/tmp.gap").readlines()
    
    tom = []
    start = False
    for line in lines:
        line = line.strip()
        if line.startswith("1:"):
            start = True
        if not start:
            continue
        if not line:
            continue
        flds = line.split(":")[1]
        flds = flds.replace(".", "0")
        flds = flds.split()
        flds = [int(i) for i in flds]
        tom.append(flds)
    
    print([len(item) for item in tom], len(tom))
    assert len(tom[-1]) == len(tom)
    
    N = len(tom)
    for item in tom:
        item += [0]*(N-len(item))
    
    lrev = lambda items:list(reversed(items))
    
    tom = [lrev(item) for item in reversed(tom)]
    
    for item in tom:
        s = " ".join((str(i) if i else ".") for i in item)
        print(s)
    





