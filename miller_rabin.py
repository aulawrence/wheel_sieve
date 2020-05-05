import numpy as np
from wheel_sieve_byte import PRIME_GEN
import pdb


def powmod(x, r, n):
    """Computes (x ** r) % n

    Args:
        x (int): Base
        r (int): Power
        n (int): Modulo

    Returns:
        int: (x ** r) % n
    """
    y = 1
    x1 = x % n
    while r >= 1:
        if r % 2 == 1:
            y = (y * x1) % n
        r //= 2
        x1 = (x1 ** 2) % n
    return y


def witness_uniform(n, k):
    """Choose k distinct numbers from the range [2, n-2] with equal weights.

    Args:
        n (int): Number determining the range, n >= 4
        k (int): Number of elements to be selected, k <= n - 3

    Returns:
        list of int: List of chosen numbers
    """
    assert n >= 4
    assert k <= n - 3
    return np.random.choice(n - 3, k, replace=False) + 2


def witness_prime(ubound):
    """Returns the list of primes in [1, ubound).

    Args:
        ubound (int): Upper bound, not inclusive

    Returns:
        list of int: List of primes
    """
    return list(PRIME_GEN(ubound))


def miller_rabin(n, witness_list):
    """Miller-Rabin Test on n with a given list of witnesses.

    Args:
        n (int): Number to be tested for primality.
        witness_list (list of int): List of witnesses.

    Returns:
        bool: Whether n passes the test:
              True indicates possible prime
              False indicates composite
    """
    if (n == 2 or n == 3):
        return True
    elif (n < 5 or n % 2 == 0):
        return False
    r = 0
    d = n - 1
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in witness_list:
        if 2 <= a <= n - 2:
            x = powmod(a, d, n)
            if x not in (1, n - 1):
                done = True
                for _ in range(r):
                    x = (x * x) % n
                    if x == n - 1:
                        done = False
                        break
                if done:
                    return False
    return True


def probable_primes(n, d, ubound, witness_list):
    """Generates a list of probable primes in the range [n, n+d).
       1. Sieve the list of multiples of primes in [1, ubound).
       2. Apply the Miller-Rabin Test with witness_list on the remaining values and yield the numbers passing the test.

    Args:
        n (int): Lower bound of range.
        d (int): Size of the range.
        ubound (int): Upper bound of primes used in the sieve in step 1.
        witness_list (list of int): list of witnesses used in the Miller-Rabin Test in step 2.

    Yields:
        int: Probable prime.
    """
    b = ~np.zeros((d,), dtype=bool)
    for p in PRIME_GEN(ubound):
        i = ((n - 1) // p + 1) * p - n
        b[i::p] = False
    for i in np.nonzero(b)[0]:
        if miller_rabin(n + int(i), witness_list):
            yield n + i


if __name__ == "__main__":
    np.random.seed(2)
    n = 3 * 2 ** 2046
    for i, v in enumerate(np.random.randint(2**32, size=(2046 // 32))):
        n += v * 2 ** (i * 32)
    n += np.random.randint(2**(2046 % 32)) * 2 ** (2046 // 32 * 32)
    d = 10000

    for v in probable_primes(n - d, d, 3000, witness_prime(1000)):
        print(v)
