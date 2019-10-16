#!/usr/bin/env python2

from __future__ import print_function

from bruhat.dev.sage_gates import *

C1 = mulclose([X, Z, wI], mul)
print("|C1| =", len(C1)) # == 16

C2 = mulclose([X, S, H], mul)
print("|C2| =", len(C2)) # == 192


def is_cliff(A):
    Ai = inv(A)
    for g in [X, Z]:
        h = A*g*Ai
        h.set_immutable()
        if h not in C1:
            return False
    return True

#for A in C2:
#    assert is_cliff(A)
assert S not in C1
assert is_cliff(S)

def is_third_level(A):
    Ai = inv(A)
    for g in [X, Z]:
        h = A*g*Ai
        h.set_immutable()
        if not is_cliff(h):
            return False
    return True

assert not is_cliff(T)
assert is_third_level(T)

if 0:
    C3 = set(C2)
    for a in C2:
        aT = a*T
        for b in C2:
            A = aT*b
            A.set_immutable()
            C3.add(A)
    print("|C3| =", len(C3))
    #for A in C3:
    #    assert is_third_level(A)

    def src(a):
        return set([b for b in C3 if mul(a,b) in C3])

    def tgt(a):
        return set([b for b in C3 if mul(b,a) in C3])

    srcs = []
    for b in C3:
        src_b = src(b)
        if src_b not in srcs:
            print("|src_b| = ", len(src_b))
            srcs.append(src_b)
        if d==2 and len(srcs)==4: # there is only 4 of these to be found
            break
        if d==3 and len(srcs)==25: # conjecture that this is all we need
            break

    obs = list(srcs)
    tgts = []
    for b in C3:
        tgt_b = tgt(b)
        if tgt_b not in obs:
            obs.append(tgt_b)
        if tgt_b not in tgts:
            print("|tgt_b| = ", len(tgt_b))
            tgts.append(tgt_b)
        if d==2 and len(tgts)==4: # there is only 4 of these to be found
            break
        if d==3 and len(tgts)==25: # conjecture that this is all we need
            break

    done = False
    while not done:
        done = True
        print("obs:", len(obs))

        obs.sort(key = len, reverse=True)
        obs1 = list(obs)
        for s in obs:
          for t in obs:
            st = s.intersection(t)
            a = ' '
            if st not in obs1:
                a = '*'
                obs1.append(st)
                done = False
            print("%4d"%(len(st)), end=a)
          print()
        obs = obs1
        obs.sort(key = len, reverse=True)

    print("obs:", len(obs), [len(ob) for ob in obs])
    print(set(C2) == obs[-1])
    for ob in obs:
        #print(len(ob), end=" ")
        iob = set([inv(a) for a in ob])
        print(iob in obs, end=" ")
        #iob = [ia for a in iob if inv(a) in ob]
        #print((len(ob), len(iob)), end=" ")

    print()




#A = MS.matrix([1, 0, 0, 0, -1, 0, 0, 0, -1])
#print(is_cliff(A)) # nope
#print(is_third_level(A)) # nope

# --------------- two qudits -------------- #

XI = X.tensor_product(I)
IX = I.tensor_product(X)
ZI = Z.tensor_product(I)
IZ = I.tensor_product(Z)
wII = w2*II

gen = [XI, IX, ZI, IZ, wII]
for op in gen:
    op.set_immutable()
    #print()
    #print(op)
    assert op*inv(op) == II

C1_2 = mulclose(gen, mul)
print("|C1(2)| =", len(C1_2))

def is_cliff_2(A):
    Ai = inv(A)
    for g in [XI, IX, ZI, IZ]:
        h = A*g*Ai
        h.set_immutable()
        if h not in C1_2:
            return False
    return True

def is_third_level_2(A):
    Ai = inv(A)
    for g in [XI, IX, ZI, IZ]:
        h = A*g*Ai
        h.set_immutable()
        if not is_cliff_2(h):
            return False
    return True


def block_diag(blocks):
    A = copy(MS_2.zero_matrix())
    
    assert len(blocks)==d
    for i in range(d):
        A[d*i:d*(i+1), d*i:d*(i+1)] = blocks[i]
    
    A.set_immutable()
    return A

if 0:
    S_2 = block_diag([Z**(i**2) for i in range(d)])
    
    #print(S_2)
    print(S_2 in C1_2)
    print(is_cliff_2(S_2))
    print(is_third_level_2(S_2))



