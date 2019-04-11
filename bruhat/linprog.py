#!/usr/bin/env python3

from pulp import LpProblem, LpMinimize, LpVariable, LpStatus, value




system = LpProblem("system", LpMinimize)

vs = {}
def W(i, j):
    v = vs.get((i, j))
    if v is None:
        name = "W_%d_%d"%(i, j)
        v = LpVariable(name)
        vs[i,j] = v
    return v


for expr in [
    W(2,1) + W(3,1) == W(1,1) + 2*W(1,1),
#    2*W(1,1) + W(4,1) == W(2,1) + W(3,2),
#    W(3,1) + W(4,1) == W(1,1) + W(3,2),
#    W(3,2) + W(5,1) == W(2,1) + 3*W(1,1),
#    W(4,1) + W(5,1) == W(1,1) + 2*W(2,1),
#    W(3,2) + W(5,2) == 2*W(1,1) + W(4,3),
#    3*W(1,1) + W(6,1) == W(3,1) + W(4,3),
#    2*W(2,1) + W(6,1) == W(2,1) + W(4,3),
#    W(5,1) + W(6,1) == W(1,1) + W(5,2),
#    3*W(1,1) + 2*W(3,1) == W(3,2) + W(5,3),
#    2*W(2,1) + 2*W(3,1) == 2*W(1,1) + 4*W(1,1),
#    W(4,3) + W(7,1) == W(3,1) + 4*W(1,1),
#    W(5,2) + W(7,1) == W(2,1) + W(5,3),
#    W(6,1) + W(7,1) == W(1,1) + 2*W(3,1),
#    W(4,3) + W(7,2) == W(3,2) + W(5,4),
#    4*W(1,1) + W(8,1) == W(4,1) + W(5,4),
#    W(5,2) + W(7,2) == 2*W(1,1) + W(5,4),
#    W(5,3) + W(8,1) == W(3,1) + W(5,4),
    2*W(3,1) + W(8,1) == W(2,1) + 3*W(2,1),
#    W(7,1) + W(8,1) == W(1,1) + W(7,2),
    ]:

    system += expr

for v in vs.values():
    system += v >= 1.0

system += W(3,1) >= 2*W(1,1)

status = system.solve()
s = LpStatus[status]
print("pulp:", s)

#for v in vs.values():
#    print(v, "=", value(v))

if 0:
    names = {}
    vs = {}
    for (i, j) in keys:
        name = "W_%d_%d"%(i, j)
        names[i, j] = name
        var = LpVariable(name)
        vs[i, j] = var 
        system += var >= 1.0 

    for i, row in enumerate(rows):
        expr = 0 
        for j, key in enumerate(keys):
            value = row.get(key, 0)
            expr += value * vs[key]
        #print(expr)
        system += expr == 0.0 

    #system += vs[3, 1] >= vs[2, 2] # Infeasible at N=9
    #system += vs[5, 1] >= vs[4, 2] # Infeasible at N=25
    system += vs[4, 2] >= vs[3, 3]  # Infeasible at N=64

    status = system.solve()
    s = LpStatus[status]
    print("pulp:", s)
    if s == "Infeasible":
        return

    if 0:
        for key in keys:
            v = vs[key]
            print(v, "=", value(v))
     


