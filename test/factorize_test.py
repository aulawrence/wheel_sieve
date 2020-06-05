import unittest
from wheel_sieve.factorize import factorize


class TestFactorize(unittest.TestCase):
    def test_factorize(self):
        prime_dict, factor_dict = factorize(2**3 * 3**5 * 5**7 * 7**11 * 997)
        self.assertDictEqual(prime_dict, {2: 3, 3: 5, 5: 7, 7: 11, 997: 1})
        self.assertDictEqual(factor_dict, {})
        prime_dict, factor_dict = factorize(2 ** 64 - 1)
        num = 1
        for p, d in prime_dict.items():
            num *= p ** d
        for f, d in factor_dict.items():
            num *= f ** d
        self.assertEqual(num, 2 ** 64 - 1)


if __name__ == "__main__":
    unittest.main()
