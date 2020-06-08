"""Prime sieve with wheel factorization.
"""
import numpy as np

# Memory (bytes) used in sieving numpy array
MEM = 8_000_000


class Worker(object):
    """Prime generating worker. Keeps a list of generated primes.
    """

    def __init__(self, ubound=None):
        super(Worker, self).__init__()
        self.ubound = ubound
        self.primes = [2, 3, 5]
        self.lim = 6

    def set_n(self, ubound):
        """Set upper bound for worker. Worker will not generate primes >= upper bound.

        Args:
            ubound (int): upper bound.
        """
        if self.ubound is None or ubound > self.ubound:
            self.ubound = ubound

    def gen(self):
        """Generate next batch of primes.
        """
        prime_tail = self.primes[-1]
        new_lim = prime_tail * prime_tail - 1
        if self.ubound is not None:
            new_lim = min(new_lim, self.ubound + 1)
        if new_lim > self.lim:
            for prime in wheel_sieve(self.lim, new_lim):
                self.primes.append(prime)
            self.lim = new_lim

    def __call__(self, ubound):
        """Returns a generator for primes in [1, ubound).

        Args:
            ubound (int): upper bound.

        Yields:
            int: primes in ascending order.
        """
        self.set_n(ubound)
        i = 0
        while True:
            if i >= len(self.primes):
                self.gen()
            if i >= len(self.primes) or self.primes[i] >= ubound:
                break
            else:
                yield self.primes[i]
                i += 1


PRIME_GEN = Worker()


def wheel_sieve_count(lbound, ubound, p_list=(2, 3, 5)):
    """Count primes in [lbound, ubound) using wheel sieve.

    Args:
        lbound (int): lower bound of range.
        ubound (int): upper bound of range.
        p_list (tuple, optional): wheel primes. Defaults to (2, 3, 5).

    Raises:
        ValueError: Thrown when input is bad.

    Returns:
        int: number of primes in the range.
    """
    if (
        not isinstance(lbound, int)
        or not isinstance(ubound, int)
        or not 1 <= lbound <= ubound
    ):
        raise ValueError
    if not all(isinstance(x, int) and x > 1 for x in p_list):
        raise ValueError
    PRIME_GEN.set_n(int(np.sqrt(ubound)) + 1)
    # wheel : product of primes in p_list
    # k_list : numbers from 1 to wheel coprime to wheel
    # kl_list : multiplication table of k_list mod wheel
    wheel = 1
    for wheel_prime in p_list:
        wheel *= wheel_prime
    k_list = [
        k
        for k in range(1, wheel)
        if all(k % wheel_prime != 0 for wheel_prime in p_list)
    ]
    kl_list = {k: [k_list.index(k * l % wheel) for l in k_list] for k in k_list}
    # Count primes between lbound and max(p_list)
    count = 0
    if lbound <= max(p_list):
        for wheel_prime in PRIME_GEN(max(p_list) + 1):
            if lbound <= wheel_prime and wheel % wheel_prime == 0:
                count += 1
    # Split into smaller partitions to reduce memory usage
    step_size = max(1, MEM // len(k_list)) * wheel
    prev = lbound
    while prev < ubound:
        curr = ((prev + step_size - 1) // wheel + 1) * wheel
        curr = min(curr, ubound)
        res = _wheel_sieve(prev, curr, wheel, k_list, kl_list)
        count += np.count_nonzero(res)
        prev = curr
    return count


def wheel_sieve(lbound, ubound, p_list=(2, 3, 5)):
    """Generate primes in [lbound, ubound) using wheel sieve.

    Args:
        lbound (int): lower bound of range.
        ubound (int): upper bound of range.
        p_list (tuple, optional): wheel primes. Defaults to (2, 3, 5).

    Raises:
        ValueError: Thrown when input is bad.

    Yields:
        int: primes in ascending order.
    """
    if (
        not isinstance(lbound, int)
        or not isinstance(ubound, int)
        or not 1 <= lbound <= ubound
    ):
        raise ValueError
    if not all(isinstance(x, int) and x > 1 for x in p_list):
        raise ValueError
    PRIME_GEN.set_n(int(np.sqrt(ubound)) + 1)
    # wheel : product of primes in p_list
    # k_list : numbers from 1 to wheel coprime to wheel
    # kl_list : multiplication table of k_list mod wheel
    wheel = 1
    for wheel_prime in p_list:
        wheel *= wheel_prime
    k_list = [
        k
        for k in range(1, wheel)
        if all(k % wheel_prime != 0 for wheel_prime in p_list)
    ]
    kl_list = {k: [k_list.index(k * l % wheel) for l in k_list] for k in k_list}
    # Yield primes between lbound and max(p_list)
    if lbound <= max(p_list):
        for prime in PRIME_GEN(max(p_list) + 1):
            if lbound <= prime and wheel % prime == 0:
                yield prime
    # Split into smaller partitions to reduce memory usage
    step_size = max(1, MEM // len(k_list)) * wheel
    prev = lbound
    while prev < ubound:
        curr = prev + step_size
        curr = min(curr, ubound)
        res = _wheel_sieve(prev, curr, wheel, k_list, kl_list)
        for km, k in zip(*np.where(res)):
            yield int((prev // wheel + km) * wheel + k_list[k])
        prev = curr


def _wheel_sieve(lbound, ubound, wheel, k_list, kl_list):
    # Ignore numbers outside of range [lbound, ubound)
    am = lbound // wheel
    bm = (ubound - 2) // wheel + 1
    c_array = np.ones((bm - am, len(k_list)), dtype=np.bool)
    for i, k in enumerate(k_list):
        if am * wheel + k < lbound or am * wheel + k == 1:
            c_array[0, i] = False
        if (bm - 1) * wheel + k >= ubound:
            c_array[-1, i] = False
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
            c_array[nm::prime, idx] = False
    return c_array


if __name__ == "__main__":
    print(wheel_sieve_count(1, 1_000_000_000))
