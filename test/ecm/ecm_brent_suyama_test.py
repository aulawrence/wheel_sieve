import unittest
from wheel_sieve.ecm.ecm_brent_suyama import apply_polynomial, get_difference_seq


def step_seq(seq):
    for i in range(len(seq) - 1):
        seq[i] += seq[i + 1]


class TestECMBrentSuyama(unittest.TestCase):
    def test_get_difference_seq(self):
        polynomial = (0, 1, 0, 2, 0, 3)
        a0 = 5
        an = 100
        d = 7
        diff_seq = get_difference_seq(polynomial, a0, d)
        for i in range(a0, an, 7):
            target = apply_polynomial(polynomial, i)
            self.assertEqual(diff_seq[0], target)
            step_seq(diff_seq)


if __name__ == "__main__":
    unittest.main()
