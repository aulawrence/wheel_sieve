import unittest
from wheel_sieve.ecm.ecm_weierstrass import get_curve, mul_pt_exn, mul_pt_multi, InverseNotFound, gcd


class TestECMWeierstrass(unittest.TestCase):
    def test_mul_pt_multi_normal(self):
        n = 65537 * 65539
        pt = (1, 1)
        curve = get_curve(pt, 133, n)
        target = [mul_pt_exn(pt, curve, k) for k in range(-1000, 1001)]
        actual = mul_pt_multi(pt, curve, range(-1000, 1001))
        self.assertEqual(target, actual)

    def test_mul_pt_multi_exception(self):
        n = 65537 * 65539
        pt = (1, 1)
        curve = get_curve(pt, 133, n)
        with self.assertRaises(InverseNotFound) as cm_target:
            target = [mul_pt_exn(pt, curve, k) for k in range(8000, 9000)]
        with self.assertRaises(InverseNotFound) as cm_actual:
            actual = mul_pt_multi(pt, curve, range(8000, 9000))
        self.assertEqual(gcd(cm_target.exception.x, n), 65539)
        self.assertEqual(gcd(cm_actual.exception.x, n), 65539)


if __name__ == "__main__":
    unittest.main()
