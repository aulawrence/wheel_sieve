from collections import defaultdict
import random
from wheel_sieve.miller_rabin import miller_rabin, witness_prime
from wheel_sieve.common import PRIME_GEN, inv_power
from wheel_sieve.ecm.ecm_polyeval import ecm


def factor_small_primes(n, ubound):
    """Factor n with small primes into n = \\prod{i} {p_i}^{d_{p_i}} * x.

    Args:
        n (int): Number to be factored.
        ubound (int): Upper bound for small primes, not inclusive.

    Returns:
        tuple(dict, int): (prime_factors, x), where
            prime_factors: dict(int, int) mapping p_i to d_{p_i},
            x: (int) x.

    """
    prime_factors = dict()
    x = n
    for prime in PRIME_GEN(ubound):
        power = 0
        while x % prime == 0:
            power += 1
            x //= prime
        if power > 0:
            prime_factors[prime] = power
    return prime_factors, x


def factor_power(n, ubound):
    """Find a factor of n if n is a perfect k-th power for some k < ubound.
    i.e. n = d ^ k for some integers d, k, where k >= 2.

    Args:
        n (int): Number to be factored.
        ubound (int): Upper bound for k, not inclusive.

    Returns:
        tuple(int, int): (d, k), or (None, None) if n is not a perfect power.
    """
    for power in PRIME_GEN(ubound):
        factor = inv_power(n, power)
        if factor is not None:
            return factor, power
    return None, None


def factor_ecm(n, ecm_kwargs_list, seed=None):
    """Find a factor of a number n using ECM.

    Args:
        n (int): Number to be factored.
        ecm_kwargs_list (list of dict): List of dicts containing keyword arguments to be passed to ecm call.
        seed (int, optional): Random seed to be set every ecm call. Defaults to None.

    Returns:
        int: Factor, or None if not found.
    """
    for ecm_kwargs in ecm_kwargs_list:
        if seed is not None:
            random.seed(seed)
        factor = ecm(n, **ecm_kwargs)
        if factor is not None:
            return factor
    return None


def factorize(n, witness=witness_prime(100)):
    """Factorize a number n, where n >= 2, with ECM into n = \\prod{i} {p_i}^{d_{p_i}} * \\prod{j} {f_j}^{d_{f_j}}.
    Each p_i passes the Miller Rabin Primality Test and is (probably) prime.
    Each f_j are known composite that cannot be factored because we try a fixed number of curves.

    Args:
        n (int): Integer to factorize.
        witness (list of int, optional): Witness to be used in Miller Rabin Primality Test. Defaults to witness_prime(100).

    Raises:
        ValueError: Thrown when n < 2

    Returns:
        tuple(dict, dict): (prime_factors, remaining_factors), where
            prime_factors: dict(int, int) mapping p_i to d_{p_i},
            remaining_factors: dict(int, int) mapping f_j to d_{f_j}.
    """
    if n < 2:
        raise ValueError
    prime_factors, factor = factor_small_primes(n, 1033)
    if factor == 1:
        return prime_factors, dict()
    remaining_factors = defaultdict(int)
    working_dict = defaultdict(int)
    working_dict[factor] = 1
    while working_dict:
        factor_i, power_i = working_dict.popitem()
        if miller_rabin(factor_i, witness):
            for dt in [working_dict, remaining_factors]:
                for factor_j in list(dt.keys()):
                    if factor_j % factor_i == 0:
                        power_j = dt.pop(factor_j)
                        while factor_j % factor_i == 0:
                            factor_j //= factor_i
                            power_i += power_j
                        if factor_j > 1:
                            working_dict[factor_j] += power_j
            prime_factors[factor_i] = power_i
        else:
            factor, power = factor_power(factor_i, factor_i.bit_length() // 10 + 1)
            if factor is not None:
                working_dict[factor] += power_i * power
                continue
            ecm_kwargs_list = [
                {
                    "rounds": 10,
                    "b1": 2_000,
                    "b2": 50_000,
                    "wheel": 210,
                    "output": False,
                },
                {
                    "rounds": 40,
                    "b1": 11_000,
                    "b2": 600_000,
                    "wheel": 2310,
                    "output": False,
                },
                {
                    "rounds": 100,
                    "b1": 50_000,
                    "b2": 4_000_000,
                    "wheel": 2310,
                    "output": False,
                },
                {
                    "rounds": 200,
                    "b1": 250_000,
                    "b2": 40_000_000,
                    "wheel": 2310,
                    "output": False,
                },
            ]
            factor = factor_ecm(factor_i, ecm_kwargs_list, seed=2)
            if factor is not None:
                working_dict[factor_i // factor] += power_i
                working_dict[factor] += power_i
                continue
            remaining_factors[factor_i] += power_i
    return prime_factors, dict(remaining_factors)


if __name__ == "__main__":
    print(factorize((2 ** 256 - 1) * (2 ** 64 - 1) ** 3))
