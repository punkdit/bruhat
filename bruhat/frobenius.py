#!/usr/bin/env python

from bruhat import element
from bruhat.lin import Space, Lin
from bruhat.argv import argv


def test():
    ring = element.Q
    one = ring.one
    N = Space(ring, 0, 0, 'N')
    I = Space(ring, 1, 0, 'I')
    V = Space(ring, 2, 0, 'V')

    # direct sum Q+Q
    # green spider
    unit = Lin(V, I, [[1], [1]])
    mul = Lin(V, V@V, [[1, 0, 0, 0], [0, 0, 0, 1]])
    counit = Lin(I, V, [[1, 1]])
    cup = Lin(V@V, I, [[1],[0],[0],[1]])
    test_frobenius(mul, unit, counit, cup)

    # complex numbers as a Q-algebra
    unit = Lin(V, I, [[1], [0]])
    mul = Lin(V, V@V, [[1, 0, 0, -1], [0, 1, 1, 0]])
    counit = Lin(I, V, [[2, 0]])
    cup = Lin(V@V, I, [[one/2],[0],[0],[-one/2]])
    test_frobenius(mul, unit, counit, cup)

    # Q[sqrt(-5)]
    unit = Lin(V, I, [[1], [0]])
    mul = Lin(V, V@V, [[1, 0, 0, -5], [0, 1, 1, 0]])
    counit = Lin(I, V, [[2, 0]])
    cup = Lin(V@V, I, [[one/2],[0],[0],[-one/10]])
    test_frobenius(mul, unit, counit, cup)


def test_cyclotomic():
    ring = element.Q
    one = ring.one
    N = Space(ring, 0, 0, 'N')
    I = Space(ring, 1, 0, 'I')
    n = 4
    V = Space(ring, n, 0, 'V')
    #unit = Lin(V, I, [[1], [0]])
    
    from bruhat.slow_element import CyclotomicRing
    from bruhat import elim
    R = CyclotomicRing(8)
    w = R.x
    assert (w**4 + 1 == 0)

    unit = elim.zeros(ring, n, 1)
    unit[0, 0] = 1
    unit = Lin(V, I, unit)
    mul = elim.zeros(ring, n, n*n)
    basis = [w**i for i in range(n)]
    for i,a in enumerate(basis):
      for j,b in enumerate(basis):
        c = a*b
        if c in basis:
            k = basis.index(c)
            mul[k, i + n*j] = 1
        elif -c in basis:
            k = basis.index(-c)
            mul[k, i + n*j] = -1
        else:
            assert 0
        #print("%6s"%str(a*b), end=" ")
      #print()

    mul = Lin(V, V@V, mul)
    #print(mul)
    counit = elim.zeros(ring, 1, n)
    counit[0, 0] = n
    counit = Lin(I, V, counit)
    
    cap = counit * mul
    #print(cap)

    cup = elim.zeros(ring, n*n, 1)
    for i in range(n*n):
        v = cap[0, i]
        if v != 0:
            cup[i, 0] = 1//cap[0, i]
    cup = Lin(V@V, I, cup)
    #print(cup)
    #print(cap*cup)

    test_frobenius(mul, unit, counit, cup)



def test_frobenius(mul, unit, counit, cup):
    V, I = unit.hom

    i_iii = Lin(I, I@I@I, [[1]])
    i_ii = Lin(I, I@I, [[1]])
    iii_i = Lin(I@I@I, I, [[1]])
    ii_i = Lin(I@I, I, [[1]])
    i_i = Lin(I, I, [[1]])

    v_vi = (V@I).unitor()
    v_iv = (I@V).unitor()
    vi_v = (V@I).unitor(True)
    iv_v = (I@V).unitor(True)

    iV = V.identity()
    swap = (V@V).get_swap([1, 0])

    cap = counit * mul

    # _assoc
    lhs = mul * (iV @ mul)
    rhs = mul * (mul @ iV)
    assert lhs == rhs

    # unit
    assert iV == mul * (iV @ unit) * vi_v
    assert iV == mul * (unit @ iV) * iv_v

    # comm
    assert mul == mul * swap

    # V is self-adjoint 
    assert v_iv * (cap @ iV) * (iV @ cup) * vi_v == iV
    assert v_vi * (iV @ cap) * (cup @ iV) * iv_v == iV

    # build comul using mul, cup, etc.
    comul = (iV @ mul @ iV) * (cup @ cup)
    comul = (iV @ swap) * comul
    comul = (iV @ iV @ cap) * (comul @ iV)
    comul = comul.tgt.unitor() * comul * comul.src.unitor(True)

    print("comul")
    print(comul)

    lhs = comul * mul
    rhs = (mul @ iV) * (iV @ comul)
    assert lhs == rhs
    
    # special
    assert mul * comul == iV

    print("extra: counit*unit == %s" %( (counit * unit)[0,0]))

    if V.n == 2:
        # this holds for 2-dim frobenius 
        op = (i_iii * (counit @ counit @ counit)
            - i_ii * (counit @ cap)
            - i_ii * (counit @ cap) * (swap @ iV) 
            - i_ii * (cap @ counit)
            + 2 * counit * mul * (mul @ iV)
        )
    
        assert op.is_zero()


def main():
    test()
    test_cyclotomic()



if __name__ == "__main__":

    from time import time
    start_time = time()

    _seed = argv.get("seed")
    if _seed is not None:
        print("seed(%s)"%_seed)
        seed(_seed)

    profile = argv.profile
    fn = argv.next() or "main"

    print("%s()"%fn)

    if profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()

    print("\nOK: finished in %.3f seconds"%(time() - start_time))
    print()


