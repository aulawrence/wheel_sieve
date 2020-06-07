import unittest
from wheel_sieve.wheel_sieve_byte import wheel_sieve, wheel_sieve_count


class TestWheelSieveByte(unittest.TestCase):
    def test_wheel_sieve(self):
        self.assertEqual(list(wheel_sieve(11, 31)), [11, 13, 17, 19, 23, 29])
        self.assertEqual(list(wheel_sieve(11, 32)), [11, 13, 17, 19, 23, 29, 31])

    def test_wheel_sieve_count(self):
        self.assertEqual(wheel_sieve_count(1, 101), 25)
        self.assertEqual(wheel_sieve_count(1, 102), 26)


if __name__ == "__main__":
    unittest.main()
