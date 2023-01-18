#!/usr/bin/env python

import sys

data = open(sys.argv[1]).read()

data = data.replace(".", "")
data = data.replace("3", "")
data = data.replace(" ", "")
data = data.replace("\n\n", "\n")
data = data.replace("[[", "[")
data = data.replace("]]", "]\n")
data = data.replace("[", "")
data = data.replace("]", "")

#data = data.replace("\n\n", ".")
#data = data.replace("\n", "")
#data = data.replace(".", "\n")

f = open('oeqc.txt', 'w')
f.write(data)
f.close()


from bruhat.solve import parse, span

items = data.split("\n\n")
print(len(items))

found = set()
for item in items:
    item = item.strip()
    if not item or item[0] not in "01":
        print(repr(item))
        continue
    try:
        A = parse(item)
    except ValueError:
        print(item)
        raise
    n = A.shape[1]
    dist = [0 for i in range(n)]
    for v in span(A):
        dist[v.sum()] += 1
    dist = tuple(dist)
    if dist not in found:
        print(dist)
        found.add(dist)
        


print("done.\n")


