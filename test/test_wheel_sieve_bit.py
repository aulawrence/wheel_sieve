import unittest
from wheel_sieve.wheel_sieve_bit import wheel_sieve_count


class TestWheelSieveBit(unittest.TestCase):
    def test_wheel_sieve_count(self):
        self.assertEqual(wheel_sieve_count(1, 101), 25)
        self.assertEqual(wheel_sieve_count(1, 102), 26)


if __name__ == "__main__":
    unittest.main()
