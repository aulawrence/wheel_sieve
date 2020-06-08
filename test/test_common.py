import unittest
from wheel_sieve.common import inv, inv_multi, inv_power, gcd, InverseNotFound


class TestECMCommon(unittest.TestCase):
    def test_inv_multi_normal(self):
        n = 65537 * 65539
        element_list = range(1, 65537)
        target = {element: inv(element, n) for element in element_list}
        actual = inv_multi(element_list, n)
        self.assertEqual(target, actual)

    def test_inv_multi_error(self):
        n = 65537 * 65539
        element_list = range(60000, 65538)
        with self.assertRaises(InverseNotFound) as cm_target:
            _target = {element: inv(element, n) for element in element_list}
        with self.assertRaises(InverseNotFound) as cm_actual:
            _actual = inv_multi(element_list, n)
        self.assertEqual(gcd(cm_target.exception.x, n), 65537)
        self.assertEqual(gcd(cm_actual.exception.x, n), 65537)

    def test_inv_power(self):
        self.assertIsNone(inv_power(63, 3))
        self.assertEqual(inv_power(64, 3), 4)
        self.assertIsNone(inv_power(65, 3))
        self.assertEqual(inv_power(12 ** 6, 2), 1728)
        self.assertEqual(inv_power(2 ** 20000, 20000), 2)


if __name__ == "__main__":
    unittest.main()
