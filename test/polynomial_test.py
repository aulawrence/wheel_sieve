import random
import unittest
from wheel_sieve.polynomial import Polynomial, inv


def mul(a, b, n):
    da = len(a)
    db = len(b)
    c = [0 for _ in range(da + db - 1)]
    for i in range(da):
        for j in range(db):
            c[i + j] = (c[i + j] + a[i] * b[j]) % n
    return c


def divmod_poly(a, b, n):
    da = len(a)
    db = len(b)
    if da < db:
        return [0], a
    bn_inv = inv(b[-1], n)
    quo = []
    rem = [pi for pi in a]
    for i in range(da - db + 1):
        quo_i = rem[da - i - 1] * bn_inv % n
        quo.append(quo_i)
        for j in range(db):
            rem[da - i - db + j] -= quo_i * b[j]
    for i in range(len(rem)):
        rem[i] %= n
    while len(rem) > 1 and rem[-1] == 0:
        rem.pop()
    return list(reversed(quo)), rem


def recip(a, n):
    d = len(a) - 1
    x_2d = [0 for _ in range(2 * d)]
    x_2d.append(1)
    return divmod_poly(x_2d, a, n)[0]


class TestPolynomialMethods(unittest.TestCase):
    def test_add(self):
        n = 4
        a = [0]
        b = [3]
        c = [3]
        self.assertEqual(Polynomial(a, n) + Polynomial(b, n), Polynomial(c, n))
        n = 4
        a = [0, 1, 2, 3, 3]
        b = [3, 3, 2, 2, 1]
        c = [3, 0, 0, 1]
        self.assertEqual(Polynomial(a, n) + Polynomial(b, n), Polynomial(c, n))
        n = 4
        a = [0, 1, 2, 3]
        b = [3, 3, 2, 2, 1]
        c = [3, 0, 0, 1, 1]
        self.assertEqual(Polynomial(a, n) + Polynomial(b, n), Polynomial(c, n))

    def test_sub(self):
        n = 4
        a = [0]
        b = [3]
        c = [1]
        self.assertEqual(Polynomial(a, n) - Polynomial(b, n), Polynomial(c, n))
        n = 4
        a = [0, 1, 2, 3, 0, 1]
        b = [3, 3, 2, 2, 1, 1]
        c = [1, 2, 0, 1, 3]
        self.assertEqual(Polynomial(a, n) - Polynomial(b, n), Polynomial(c, n))
        n = 5
        a = [0, 1, 2, 3]
        b = [3, 3, 2, 2, 1, 1]
        c = [2, 3, 0, 1, 4, 4]
        self.assertEqual(Polynomial(a, n) - Polynomial(b, n), Polynomial(c, n))

    def test_mul(self):
        random.seed(2)
        n = 5
        a = [3]
        b = [4]
        c = mul(a, b, n)
        self.assertEqual(Polynomial(a, n) * Polynomial(b, n), Polynomial(c, n))
        n = 5
        a = [4 for _ in range(10)]
        b = [4 for _ in range(100)]
        c = mul(a, b, n)
        self.assertEqual(Polynomial(a, n) * Polynomial(b, n), Polynomial(c, n))
        n = 97
        a = [random.randint(0, n - 1) for _ in range(100)]
        b = [random.randint(0, n - 1) for _ in range(1000)]
        c = mul(a, b, n)
        self.assertEqual(Polynomial(a, n) * Polynomial(b, n), Polynomial(c, n))

    def test_divmod(self):
        random.seed(2)
        n = 5
        a = [3]
        b = [4]
        q, r = divmod(Polynomial(a, n), Polynomial(b, n))
        self.assertEqual(len(r.coeff), 1)
        self.assertEqual(len(q.coeff), 1)
        self.assertEqual(q * Polynomial(b, n) + r, Polynomial(a, n))
        n = 5
        a = [4 for _ in range(100)]
        b = [4 for _ in range(10)]
        q, r = divmod(Polynomial(a, n), Polynomial(b, n))
        self.assertLess(len(r.coeff), len(b))
        self.assertEqual(len(q.coeff), len(a) - len(b) + 1)
        self.assertEqual(q * Polynomial(b, n) + r, Polynomial(a, n))
        n = 97
        a = [random.randint(0, n - 1) for _ in range(1000)]
        b = [random.randint(0, n - 1) for _ in range(100)]
        q, r = divmod(Polynomial(a, n), Polynomial(b, n))
        self.assertLess(len(r.coeff), len(b))
        self.assertEqual(len(q.coeff), len(a) - len(b) + 1)
        self.assertEqual(q * Polynomial(b, n) + r, Polynomial(a, n))

    def test_recip(self):
        random.seed(2)
        n = 5
        a = [3]
        self.assertEqual(Polynomial(recip(a, n), n), Polynomial(a, n).recip())
        n = 5
        a = [4 for _ in range(100)]
        self.assertEqual(Polynomial(recip(a, n), n), Polynomial(a, n).recip())
        n = 97
        a = [random.randint(0, n - 1) for _ in range(1000)]
        self.assertEqual(Polynomial(recip(a, n), n), Polynomial(a, n).recip())
        n = 97
        a = [random.randint(0, n - 1) for _ in range(63)]
        self.assertEqual(Polynomial(recip(a, n), n), Polynomial(a, n).recip())
        n = 97
        a = [random.randint(0, n - 1) for _ in range(64)]
        self.assertEqual(Polynomial(recip(a, n), n), Polynomial(a, n).recip())
        n = 97
        a = [random.randint(0, n - 1) for _ in range(65)]
        self.assertEqual(Polynomial(recip(a, n), n), Polynomial(a, n).recip())

    def test_mod_with_recip(self):
        n = 5
        a1 = Polynomial([3], n)
        a2 = Polynomial([2], n)
        b = Polynomial([4], n)
        b_recip = b.recip()
        self.assertEqual(divmod(a1, b)[1], a1.mod_with_recip(b, b_recip))
        self.assertEqual(divmod(a2, b)[1], a2.mod_with_recip(b, b_recip))
        n = 5
        a1 = Polynomial(list(range(5)) * 20, n)
        a2 = Polynomial(list(range(5)) * 10, n)
        b = Polynomial([1, 2, 3, 1], n)
        b_recip = b.recip()
        self.assertEqual(divmod(a1, b)[1], a1.mod_with_recip(b, b_recip))
        self.assertEqual(divmod(a2, b)[1], a2.mod_with_recip(b, b_recip))


if __name__ == "__main__":
    unittest.main()
