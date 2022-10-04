#!/usr/bin/env python

from time import sleep

from bruhat.gap import Gap, PROMPT



def parse(s, a, b, c, d=None):
    s = s.replace("^-1", "")
    s = s.replace("*", "+")
    s = s.replace("^", "*")
    s = s.replace(" ", "")
    return eval(s)



def make_gap(gap, orders, index=600):
    #print("make_gap", orders)
    n = len(orders)+1
    labels = "abcdef"[:n]
    lines = []
    s = ','.join('"%s"'%label for label in labels)
    lines.append('F := FreeGroup(%s);;' % s)
    lines.append('AssignGeneratorVariables(F);;')
    rels = []
    for i in range(n):
        rels.append('%s^2'%labels[i])
        if i+1<n and orders[i] is not None:
            rels.append('(%s*%s)^%s'%(labels[i], labels[i+1], orders[i]))
        for j in range(i+2, n):
            rels.append('(%s*%s)^2'%(labels[i], labels[j]))
    rels = ','.join(rels)
    #print(rels)
    lines.append("G := F/[%s];;"%rels)
    lines.append("gr := LowIndexNormalSubgroupsSearchForAll(G, %s);;"%index)
    lines.append("L := List(gr);;")

    for line in lines:
        #print(line)
        gap.send(line)
        s = gap.expect(PROMPT)
        #print(s)

    gap.send("Length(L);")
    count = gap.expect(PROMPT)
    #print("count:", count)
    count = int(count.strip())

    items = []
    for i in range(count):
        gap.send("S := Grp(L[%d]);"%(i+1))
        s = gap.expect(PROMPT)
        #print(s)
        gap.send("GeneratorsOfGroup(S);")
        item = gap.expect(PROMPT)
        item = item.replace(" ", "")
        item = item.replace("\n", "")
        items.append(item)

    return items
    



def load_db():
    gap = Gap()
    gap.expect(PROMPT)
    gap.send('LoadPackage("recog");')
    gap.expect(PROMPT)
    gap.send('LoadPackage("LINS");')
    gap.expect(PROMPT)

    db = {}
    keys = [(3, None), (5, 4), (3, 4, 4), (3, 3, 6)]
    #keys = [(5, 4), (5, 5)]
    for key in keys:
        items = make_gap(gap, key)
        #print("items:", key, len(items))
        db[key] = items
    return db



if __name__ == "__main__":

    db = load_db()

    print(db)


