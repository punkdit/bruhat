
cdef extern from *:
    int __builtin_sadd_overflow (int a, int b, int *res)
    int __builtin_ssub_overflow (int a, int b, int *res)
    int __builtin_smul_overflow (int a, int b, int *res)

    int __builtin_add_overflow_p (int a, int b, int c)
    int __builtin_sub_overflow_p (int a, int b, int c)
    int __builtin_mul_overflow_p (int a, int b, int c)


cdef int gcd(int a, int b):
    cdef int c
    while b != 0:
        c = b
        b = a%b
        a = c
    return a


cdef class Type(object):
    pass


class Keyed(object): # mixin, used for Type's

    def __init__(self, key):
        if not type(key) is tuple:
            key = (key,)
        self.key = (self.__class__,) + key

    def __hash__(self):
        return hash(self.key)

    def __eq__(self, other):
        return self.__class__ == other.__class__ and self.key == other.key

    def __ne__(self, other):
        return self.__class__ != other.__class__ or self.key != other.key

    def __str__(self):
        return str(self.key)
    __repr__ = __str__


class Element(Type):
    def __init__(self, tp):
        assert isinstance(tp, Type)
        self.tp = tp
        #self.promote = tp.promote # ??


cdef class Ring(Type):
    pass

cdef class Fraction:

    cdef readonly int top
    cdef readonly int bot

    def __init__(self, int top, int bot):
        cdef int factor
        assert bot != 0
        factor = gcd(top, bot)
        top //= factor
        bot //= factor
        self.top = top
        self.bot = bot
        assert bot != 0

    @staticmethod
    cdef Fraction promote(object other):
        if type(other) is Fraction:
            return other
        if type(other) is long:
            return Fraction(other, 1)
        return None

    @property
    def tp(self):
        return Q

    def __hash__(self):
        return hash((self.top, self.bot))

    def __nonzero__(self):
        return self.top != 0

    def __eq__(self, object _other):
        cdef Fraction other
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        #return self.top * other.bot == other.top * self.bot
        cdef int topbot, bottop
        if __builtin_smul_overflow(self.top, other.bot, &topbot):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.top, &bottop):
            raise OverflowError()
        return topbot == bottop

    def __ne__(self, object _other):
        cdef Fraction other
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        #return self.top * other.bot != other.top * self.bot
        cdef int topbot, bottop
        if __builtin_smul_overflow(self.top, other.bot, &topbot):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.top, &bottop):
            raise OverflowError()
        return topbot != bottop

    def __lt__(self, object _other):
        cdef Fraction other
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        #return self.top * other.bot < other.top * self.bot
        cdef int topbot, bottop
        if __builtin_smul_overflow(self.top, other.bot, &topbot):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.top, &bottop):
            raise OverflowError()
        return topbot < bottop

    def __gt__(self, object _other):
        cdef Fraction other
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        #return self.top * other.bot > other.top * self.bot
        cdef int topbot, bottop
        if __builtin_smul_overflow(self.top, other.bot, &topbot):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.top, &bottop):
            raise OverflowError()
        return topbot > bottop

    def __le__(self, object _other):
        cdef Fraction other
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        #return self.top * other.bot <= other.top * self.bot
        cdef int topbot, bottop
        if __builtin_smul_overflow(self.top, other.bot, &topbot):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.top, &bottop):
            raise OverflowError()
        return topbot <= bottop

    def __ge__(self, object _other):
        cdef Fraction other
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        #return self.top * other.bot >= other.top * self.bot
        cdef int topbot, bottop
        if __builtin_smul_overflow(self.top, other.bot, &topbot):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.top, &bottop):
            raise OverflowError()
        return topbot >= bottop

    def __str__(self):
        top = str(self.top)
        if self.bot == 1:
            return top
        bot = str(self.bot)
        return "%s/%s"%(top, bot)

    def __repr__(self):
        return "Fraction(%s, %s)"%(self.top, self.bot)

    def __add__(_self, _other):
        cdef Fraction self, other
        cdef int top, bot
        self = Fraction.promote(_self)
        if self is None:
            return NotImplemented
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        #top = self.top * other.bot + other.top * self.bot
        #bot = self.bot * other.bot
        cdef int topbot, bottop, botbot
        if __builtin_smul_overflow(self.top, other.bot, &topbot):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.top, &bottop):
            raise OverflowError()
        if __builtin_sadd_overflow(topbot, bottop, &top):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.bot, &bot):
            raise OverflowError()
        return Fraction(top, bot)

    def __sub__(_self, _other):
        cdef Fraction self, other
        cdef int top, bot
        self = Fraction.promote(_self)
        if self is None:
            return NotImplemented
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        #top = self.top * other.bot - other.top * self.bot
        #bot = self.bot * other.bot
        cdef int topbot, bottop, botbot
        if __builtin_smul_overflow(self.top, other.bot, &topbot):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.top, &bottop):
            raise OverflowError()
        if __builtin_ssub_overflow(topbot, bottop, &top):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.bot, &bot):
            raise OverflowError()
        return Fraction(top, bot)

    def __pos__(self):
        return self

    def __neg__(self):
        return Fraction(-self.top, self.bot)

    def __mul__(_self, _other):
        cdef Fraction self, other
        cdef int top, bot
        self = Fraction.promote(_self)
        if self is None:
            return NotImplemented
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        #top = self.top * other.top
        #bot = self.bot * other.bot
        if __builtin_smul_overflow(self.top, other.top, &top):
            raise OverflowError()
        if __builtin_smul_overflow(self.bot, other.bot, &bot):
            raise OverflowError()
        return Fraction(top, bot)

    def __pow__(self, n, modulo):
        cdef int top, bot
        if modulo is not None:
            return NotImplemented
        top = self.top ** n
        bot = self.bot ** n
        return Fraction(top, bot)

    def __floordiv__(_self, _other):
        cdef Fraction self, other
        cdef int top, bot
        self = Fraction.promote(_self)
        if self is None:
            return NotImplemented
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        top = self.top * other.bot
        bot = self.bot * other.top
        return Fraction(top, bot)

    def __truediv__(_self, _other):
        cdef Fraction self, other
        cdef int top, bot
        self = Fraction.promote(_self)
        if self is None:
            return NotImplemented
        other = Fraction.promote(_other)
        if other is None:
            return NotImplemented
        top = self.top * other.bot
        bot = self.bot * other.top
        return Fraction(top, bot)


class Rational(Ring):
    def __init__(self):
        self.one = Fraction(1, 1)
        self.zero = Fraction(0, 1)

    def promote(self, value):
        return Fraction.promote(value)

    def __str__(self):
        return "Q"

Q = Rational()



cdef class FFElement(object):
    cdef readonly int value
    cdef readonly FiniteField ring
    def __init__(self, ring, value):
        self.ring = ring
        self.value = value % ring.p

    def __hash__(self):
        return hash((self.ring, self.value,))

    def __nonzero__(self):
        return self.value != 0

    def __eq__(self, object _other):
        cdef FFElement other
        other = self.ring.promote(_other)
        if other is None:
            return NotImplemented
        return self.value == other.value

    def __ne__(self, object _other):
        cdef FFElement other
        other = self.ring.promote(_other)
        if other is None:
            return NotImplemented
        return self.value != other.value

    def __lt__(self, object _other):
        cdef FFElement other
        other = self.ring.promote(_other)
        if other is None:
            return NotImplemented
        return self.value < other.value

    def __gt__(self, object _other):
        cdef FFElement other
        other = self.ring.promote(_other)
        if other is None:
            return NotImplemented
        return self.value > other.value

    def __le__(self, object _other):
        cdef FFElement other
        other = self.ring.promote(_other)
        if other is None:
            return NotImplemented
        return self.value <= other.value

    def __ge__(self, object _other):
        cdef FFElement other
        other = self.ring.promote(_other)
        if other is None:
            return NotImplemented
        return self.value >= other.value

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return "FFElement(%s)"%(self.value,)

    def __add__(_self, _other):
        cdef FFElement self, other
        if isinstance(_other, FFElement):
            other = _other
            self = other.ring.promote(_self)
        else:
            self = _self
            other = self.ring.promote(_other)
            if other is None:
                return NotImplemented
        return FFElement(self.ring, self.value + other.value)

    def __sub__(_self, _other):
        cdef FFElement self, other
        if isinstance(_other, FFElement):
            other = _other
            self = other.ring.promote(_self)
        else:
            self = _self
            other = self.ring.promote(_other)
            if other is None:
                return NotImplemented
        return FFElement(self.ring, self.value - other.value)

    def __mul__(_self, _other):
        cdef FFElement self, other
        if isinstance(_other, FFElement):
            other = _other
            self = other.ring.promote(_self)
        else:
            self = _self
            other = self.ring.promote(_other)
            if other is None:
                return NotImplemented
        return FFElement(self.ring, self.value * other.value)

    def __floordiv__(_self, _other):
        cdef FFElement self, other
        if isinstance(_other, FFElement):
            other = _other
            self = other.ring.promote(_self)
        else:
            self = _self
            other = self.ring.promote(_other)
            if other is None:
                return NotImplemented
        return self.ring.div(self, other)

    def __truediv__(_self, _other):
        cdef FFElement self, other
        if isinstance(_other, FFElement):
            other = _other
            self = other.ring.promote(_self)
        else:
            self = _self
            other = self.ring.promote(_other)
            if other is None:
                return NotImplemented
        return self.ring.div(self, other)

    def __pos__(self):
        return self

    def __neg__(self):
        return FFElement(self.ring, -self.value)

    def __pow__(self, int n, modulo):
        if modulo is not None:
            return NotImplemented
        return FFElement(self.ring, self.value ** n)



cdef class FiniteField(Ring):
    cdef readonly int p
    cdef readonly object zero
    cdef readonly object one
    cdef readonly object inv

    def __init__(self, int p):
        self.p = p
        one = FFElement(self, 1)
        mul = [FFElement(self, i) for i in range(1, p)]
        inv = {}
        for a in mul:
          for b in mul:
            if a*b == one:
                inv[a] = b
                inv[b] = a
        self.zero = FFElement(self, 0)
        self.one = one
        self.inv = inv

    cpdef FFElement promote(self, object other):
        if type(other) is FFElement:
            return other
        if type(other) is long:
            return FFElement(self, other)
        return None

    cdef FFElement div(self, FFElement left, FFElement right):
        if right == self.zero:
            raise ZeroDivisionError()
        return left * self.inv[right]


        

class IntegerRing(Type):
    one = 1
    zero = 0
    @classmethod
    def promote(cls, a):
        return a
Z = IntegerRing()



