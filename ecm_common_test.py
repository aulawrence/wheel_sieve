import unittest
from ecm_common import inv, inv_multi, gcd, InverseNotFound


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
            target = {element: inv(element, n) for element in element_list}
        with self.assertRaises(InverseNotFound) as cm_actual:
            actual = inv_multi(element_list, n)
        self.assertEqual(gcd(cm_target.exception.x, n), 65537)
        self.assertEqual(gcd(cm_actual.exception.x, n), 65537)


if __name__ == "__main__":
    unittest.main()
