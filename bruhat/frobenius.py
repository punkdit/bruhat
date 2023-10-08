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

    unit = Lin(V, I, [[1], [0]])

    # direct sum Q+Q
    counit = Lin(I, V, [[1, 0]])
    mul = Lin(V, V@V, [[1, 0, 0, 0], [0, 0, 0, 1]])
    cup = Lin(V@V, I, [[1],[0],[0],[1]])
    test_frobenius(mul, unit, counit, cup)

    # complex numbers as a Q-algebra
    counit = Lin(I, V, [[2, 0]])
    mul = Lin(V, V@V, [[1, 0, 0, -1], [0, 1, 1, 0]])
    cup = Lin(V@V, I, [[one/2],[0],[0],[-one/2]])
    test_frobenius(mul, unit, counit, cup)

    # Q[sqrt(-5)]
    counit = Lin(I, V, [[2, 0]])
    mul = Lin(V, V@V, [[1, 0, 0, -5], [0, 1, 1, 0]])
    cup = Lin(V@V, I, [[one/2],[0],[0],[-one/10]])
    test_frobenius(mul, unit, counit, cup)


def test_frobenius(mul, unit, counit, cup):
    V, I = unit.hom
    iV = V.identity()
    swap = (V@V).get_swap([1, 0])

    cap = counit * mul
    #print( (cap*cup)[0,0] )
    #print(cap)
    #print(cup)

    # assoc
    lhs = mul * (iV @ mul)
    rhs = mul * (mul @ iV)
    assert lhs == rhs

    # unit
    assert iV == mul * (iV @ unit)
    assert iV == mul * (unit @ iV)

    # comm
    assert mul == mul * swap

    # V is self-adjoint 
    assert (cap @ iV) * (iV @ cup) == iV
    assert (iV @ cap) * (cup @ iV) == iV

    # build comul using mul, cup, etc.
    comul = (iV @ mul @ iV) * (cup @ cup)
    comul = (iV @ swap) * comul
    comul = (iV @ iV @ cap) * (comul @ iV)
    comul = comul.tgt.unitor() * comul
    
    # special
    assert mul * comul == iV

    print("extra:", (counit * unit)[0,0]==1)

    if V.n == 2:
        # this holds for 2-dim frobenius 
        u3 = Lin(I, I@I@I, [[1]])
        u2 = Lin(I, I@I, [[1]])
        op = (u3 * (counit @ counit @ counit)
            - u2 * (counit @ cap)
            - u2 * (counit @ cap) * (swap @ iV) 
            - u2 * (cap @ counit)
            + 2 * counit * mul * (mul @ iV)
        )
    
        assert op.is_zero()





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


