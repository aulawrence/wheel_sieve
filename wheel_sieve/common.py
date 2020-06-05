from math import gcd
from gmpy2 import invert, iroot
import numpy as np
from wheel_sieve.wheel_sieve_byte import PRIME_GEN, wheel_sieve


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
    try:
        return invert(x, n)
    except ZeroDivisionError:
        raise InverseNotFound(x % n, n)


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


def inv_multi(element_list, n):
    """Compute inverse (mod n) of multiple elements.
    Uses Montgomery's trick so that for a list of length k, it only takes 1 modular inverse and O(k) modular multiplications
    instead of k modular inverses.

    Args:
        element_list (list of int): List of elements to be inverted (mod n).
        n (int): Modulus.

    Returns:
        dict (int, int): Dictionary mapping elements to their inverse.
    """
    inv_dict = dict()
    d = len(element_list)
    if d == 1:
        inv_dict[element_list[0]] = inv(element_list[0], n)
        return inv_dict
    k = 1 << (d - 1).bit_length()
    element_tree = [None] * (k - 1) + [1] * k
    element_tree[k - 1:k + d - 1] = element_list
    while k > 1:
        for i in range(k // 2 - 1, k - 1):
            element_tree[i] = element_tree[2 * i + 1] * element_tree[2 * i + 2] % n
        k //= 2
    element_tree[0] = inv(element_tree[0], n)
    k = 2
    while k < d:
        for i in range(k // 2 - 1, k - 1):
            inv1 = element_tree[i] * element_tree[2 * i + 2] % n
            inv2 = element_tree[i] * element_tree[2 * i + 1] % n
            element_tree[2 * i + 1] = inv1
            element_tree[2 * i + 2] = inv2
        k *= 2
    for i in range(k - 1, k + d - 1):
        element = element_tree[i]
        if i % 2 == 0:
            inv_element = element_tree[(i - 1) // 2] * element_tree[i - 1] % n
        else:
            inv_element = element_tree[(i - 1) // 2] * element_tree[i + 1] % n
        inv_dict[element] = inv_element
    return inv_dict


def inv_power(x, d):
    """Computes the d-th root of x if it is an integer.
    x must be >= 0 and d must be >= 1.

    Args:
        x (int): x.
        d (int): d.

    Raises:
        ValueError: Thrown when x < 0 or d < 1.

    Returns:
        int: The d-th root of x if it is an integer, otherwise None.
    """
    if x < 0 or d < 1:
        raise ValueError
    if x == 0:
        return 0
    x_d, is_perfect = iroot(x, d)
    if is_perfect:
        return x_d
    return None
