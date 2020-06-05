import numpy as np
from wheel_sieve.wheel_sieve_byte import PRIME_GEN

# Memory (bytes) used in sieving numpy array
MEM = 8_000_000


def wheel_sieve_count(lbound, ubound, p_list=(2, 3, 5)):
    """Count primes in [lbound, ubound) using wheel sieve.

    Args:
        lbound (int): lower bound of range.
        ubound (int): upper bound of range.
        p_list (tuple, optional): wheel primes. Defaults to (2, 3, 5).

    Returns:
        int: number of primes in the range.
    """
    assert isinstance(lbound, int) and isinstance(ubound, int) and 1 <= lbound <= ubound
    assert all(isinstance(x, int) and x > 1 for x in p_list)
    PRIME_GEN.set_n(int(np.sqrt(ubound)) + 1)
    # wheel : product of primes in p_list
    # k_list : numbers from 1 to wheel coprime to wheel
    # kl_list : multiplication table of k_list mod wheel
    wheel = 1
    for wheel_prime in p_list:
        wheel *= wheel_prime
    k_list = [k for k in range(1, wheel) if all(k % wheel_prime != 0 for wheel_prime in p_list)]
    kl_list = {k: [k_list.index(k * l % wheel) for l in k_list] for k in k_list}
    # Count primes between lbound and max(p_list)
    count = 0
    if lbound <= max(p_list):
        for wheel_prime in PRIME_GEN(max(p_list) + 1):
            if lbound <= wheel_prime and wheel % wheel_prime == 0:
                count += 1

    # Split into smaller partitions to reduce memory usage
    step_size = max(1, MEM // len(k_list)) * wheel * 8
    prev = lbound
    while prev < ubound:
        curr = ((prev + step_size - 1) // (wheel * 8) + 1) * wheel * 8
        curr = min(curr, ubound)
        res = _wheel_sieve_bit(prev, curr, wheel, k_list, kl_list)
        count += np.sum(np.bincount(np.ndarray.flatten(res), minlength=256) * BIT_COUNT)
        prev = curr
    return count


CLEAR_BIT = (
    np.uint8(~(1 << 7)),
    np.uint8(~(1 << 6)),
    np.uint8(~(1 << 5)),
    np.uint8(~(1 << 4)),
    np.uint8(~(1 << 3)),
    np.uint8(~(1 << 2)),
    np.uint8(~(1 << 1)),
    np.uint8(~(1 << 0)),
)

CLEAR_BIT_FROM = (
    np.uint8(-1 << 8),
    np.uint8(-1 << 7),
    np.uint8(-1 << 6),
    np.uint8(-1 << 5),
    np.uint8(-1 << 4),
    np.uint8(-1 << 3),
    np.uint8(-1 << 2),
    np.uint8(-1 << 1),
)

# [sum(c == "1" for c in "{:>08}".format(bin(np.uint8(x))[2:])) for x in range(256)]

BIT_COUNT = np.array([
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
])


def _wheel_sieve_bit(lbound, ubound, wheel, k_list, kl_list):
    # Ignore numbers outside of range [lbound, ubound)
    am = lbound // wheel
    bm = (ubound - 2) // wheel + 1
    q, r = divmod(bm - am, 8)
    if r != 0:
        q += 1
    c_array = ~np.zeros((q, len(k_list)), dtype=np.uint8)
    for i, k in enumerate(k_list):
        if am * wheel + k < lbound or am * wheel + k == 1:
            c_array[0, i] &= CLEAR_BIT[0]
        if (bm - 1) * wheel + k >= ubound:
            c_array[-1, i] &= CLEAR_BIT_FROM[r - 1]
        elif r != 0:
            c_array[-1, i] &= CLEAR_BIT_FROM[r]
    # Sieve
    for prime in PRIME_GEN(int(np.sqrt(ubound)) + 1):
        if wheel % prime == 0:
            continue
        km = max(prime // wheel, am // prime)
        kp_list = kl_list[prime % wheel]
        for k, idx in zip(k_list, kp_list):
            nm = prime * km + prime * k // wheel - am
            if nm < 0 or (km == 0 and k == 1):
                nm += prime
            for _ in range(8):
                qi, ri = divmod(nm, 8)
                if qi > q:
                    break
                c_array[qi::prime, idx] &= CLEAR_BIT[ri]
                nm += prime
    return c_array


if __name__ == "__main__":
    import time
    START = time.time()
    print(wheel_sieve_count(1, 1_000_000_000))
    # print(wheel_sieve_count(1,4_000_000_000))
    # print(wheel_sieve_count(1, 1+60*10_000_000, (2,3,5,7,11)))
    # print(wheel_sieve_count(2**45, 2**45+2**30))
    END = time.time()
    print(END - START)
