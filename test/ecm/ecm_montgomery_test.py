import unittest
from wheel_sieve.ecm.ecm_montgomery import (
    get_curve_a,
    get_curve_suyama,
    add_pt,
    dbl_pt,
    mul_pt_exn,
)


class TestECMMontgomery(unittest.TestCase):
    def test_2x_pt(self):
        n = 65537
        pt1, curve1 = get_curve_a(1, 3, n)
        pt2, curve2 = get_curve_suyama(7, n)
        for pt, curve in [(pt1, curve1), (pt2, curve2)]:
            pt_2x_a = dbl_pt(pt, curve)
            pt_2x_b = mul_pt_exn(pt, curve, 2)
            self.assertEqual(pt_2x_a, pt_2x_b)

    def test_3x_pt(self):
        n = 65537
        pt1, curve1 = get_curve_a(1, 3, n)
        pt2, curve2 = get_curve_suyama(7, n)
        for pt, curve in [(pt1, curve1), (pt2, curve2)]:
            pt_3x_a = add_pt(dbl_pt(pt, curve), pt, pt, curve)
            pt_3x_b = mul_pt_exn(pt, curve, 3)
            self.assertEqual(pt_3x_a, pt_3x_b)


if __name__ == "__main__":
    unittest.main()
