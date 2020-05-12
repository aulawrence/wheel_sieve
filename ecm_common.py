from math import gcd
import numpy as np
from wheel_sieve_byte import PRIME_GEN, wheel_sieve


class InverseNotFound(Exception):
    """Inverse Not Found Exception. Cannot find inverse of x (mod n).

    Args:
        x (int): Number that could not be inverted.
        n (int): Modulus.
    """

    def __init__(self, x, n):
        super(InverseNotFound, self).__init__("Inverse of {0:d} (mod {1:d}) not found.".format(x, n))
        self.x = x
        self.n = n


class CurveInitFail(Exception):
    """Curve Init Fail Exception. Cannot init curve with given parameters.
    """


def inv(x, n):
    """Compute inverse of x (mod n).

    Args:
        x (int): Number to be inverted.
        n (int): Modulus.

    Raises:
        InverseNotFound: Thrown when x cannot be inverted.

    Returns:
        int: Inverse of x.
    """
    a = x % n
    b = n
    tb = 0
    ta = 1
    while a > 0:
        (q, a), b = divmod(b, a), a
        ta, tb = tb - q * ta, ta
    assert (tb * x - b) % n == 0
    if b != 1:
        raise InverseNotFound(x % n, n)
    return tb % n


def init_wheel(b1, b2, wheel):
    """Initialize Wheel. Generates:
    j_list: list of j, where 1 <= j < wheel // 2 and j coprime to wheel.
    prime_array: bitarray of whether each number in range [b1, b2) is a prime.
        The entries are stored in 8-bit little-endian format.
        With c1 = b1 // wheel, if either of n1 = c * wheel + j or n2 = c * wheel - j is prime,
        the following bit is set to 1:
            axis 0: c - c1
            axis 1: (index of j in j_list) // 8
            bit # : (index of j in j_list) % 8

    Args:
        b1 (int): Lower bound of prime range.
        b2 (int): Upper bound of prime range.
        wheel (int): Wheel. Typically primorial numbers like 30, 210, 2310.

    Returns:
        tuple(list of int, np.array): (j_list, prime_array).
    """
    j_list = [j for j in range(1, wheel // 2) if gcd(j, wheel) == 1]
    j_index = {j: i for i, j in enumerate(j_list)}
    c1 = b1 // wheel
    c2 = b2 // wheel + 2
    prime_array = np.zeros((c2 - c1, (len(j_list) - 1) // 8 + 1), dtype=np.uint8)
    for p in wheel_sieve(b1, b2):
        c = p // wheel
        j = p % wheel
        if j > wheel // 2:
            c += 1
            j = wheel - j
        jq, jr = divmod(j_index[j], 8)
        prime_array[c - c1, jq] |= 1 << jr
    return j_list, prime_array
