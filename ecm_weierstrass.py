import random
import numpy as np
from math import gcd
from wheel_sieve_byte import PRIME_GEN


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
    assert tb * x % n == b
    if b != 1:
        raise InverseNotFound(x % n, n)
    return tb % n


def get_curve(pt0, a, n):
    """Given a, n, and point pt0, returns the Elliptic Curve that pt0 is on.

    Args:
        pt0 (tuple(int, int)): Point (x0, y0).
        a (int): Parameter a of the Elliptic Curve.
        n (int): Modulus.

    Returns:
        tuple(int, int, int): (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).
    """
    x0, y0 = pt0
    b = (y0**2 - x0**3 - a * x0) % n
    return (a, b, n)


def get_delta(curve):
    """Computes the discriminant (4 * a**3 + 27 * b**2) % n of the Elliptic Curve.

    Args:
        curve (tuple(int, int, int)): (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).

    Returns:
        int: The discriminant.
    """
    a, b, n = curve
    delta = (4 * a**3 + 27 * b**2) % n
    return gcd(delta, n) % n


def neg_pt(pt, curve):
    """Negate a point.

    Args:
        pt (tuple(int, int)): Point (x, y). Use (None, None) for point at infinity.
        curve (tuple(int, int, int)): (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).

    Returns:
        tuple(int, int): Point -pt.
    """
    if pt == (None, None):
        return (None, None)
    x, y = pt
    a, b, n = curve
    return (x, n - y)


def add_pt_exn(pt1, pt2, curve):
    """Adds two points pt1 and pt2 on curve.

    Args:
        pt1 (tuple(int, int)): Point (x1, y1). Use (None, None) for point at infinity.
        pt2 (tuple(int, int)): Point (x2, y2). Use (None, None) for point at infinity.
        curve (tuple(int, int, int)): (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).

    Raises:
        InverseNotFound: Throws InverseNotFound when the sum is the point at infinity.

    Returns:
        tuple(int, int): Point pt1 + pt2.
    """
    x1, y1 = pt1
    x2, y2 = pt2
    a, b, n = curve
    if pt1 == (None, None):
        return pt2
    elif pt2 == (None, None):
        return pt1
    elif pt1 == pt2:
        s = ((3 * x1 * x1 + a) * inv(2 * y1, n)) % n
    else:
        s = (y2 - y1) * inv(x2 - x1, n) % n
    xr = (s * s - x1 - x2) % n
    yr = (s * (x1 - xr) - y1) % n
    return (xr, yr)


def mul_pt_exn(point, curve, k):
    """Multiplies point by k times on curve.

    Args:
        point (tuple(int, int)): Point (x, y). Use (None, None) for point at infinity.
        curve (tuple(int, int, int)): (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).
        k (int): Multiplier.

    Raises:
        InverseNotFound: Thrown when a number cannot be inverted during the calculation.

    Returns:
        tuple(int, int): Point k * point.
    """
    if k < 0:
        return mul_pt_exn(neg_pt(point, curve), curve, -k)
    res = (None, None)
    while k >= 1:
        if k % 2 == 1:
            res = add_pt_exn(res, point, curve)
        k //= 2
        point = add_pt_exn(point, point, curve)
    return res


def ecm(n, rounds, b1, b2):
    """Elliptic Curve Factorization Method.
    For each round:
        0. Generate random point and curve.
        1. Repeatedly multiply the current point by small primes raised to some power, determined by b1.
        2. Repeatedly try to multiply the point from step 1 by possible primes (with wheel of 210) between b1 and b2.
    Returns when a non-trivial factor is found.

    Args:
        n (int): Number to be factorized.
        rounds (int): Number of random curves to try.
        b1 (int): Bound for primes used in step 1.
        b2 (int): Bound for primes searched for in step 2. b1 < b2.

    Returns:
        int: Non-trivial factor if found, otherwise returns None.
    """
    k_ls = []
    for p in PRIME_GEN(b1):
        k_ls.append(p ** int(round(np.log(b1) / np.log(p))))
    for roundi in range(rounds):
        print("Round {}...".format(roundi))
        count = 0
        delta = 0
        while delta == 0 and count < 20:
            count += 1
            x0 = random.randint(1, n - 1)
            y0 = random.randint(1, n - 1)
            a = random.randint(1, n - 1)
            pt = (x0, y0)
            curve = get_curve(pt, a, n)
            delta = get_delta(curve)
        if delta == 0:
            break
        if 1 < delta < n:
            return delta
        try:
            # Step 1
            print(" - Step 1")
            for k in k_ls:
                pt = mul_pt_exn(pt, curve, k)
            # Step 2
            print(" - Step 2")
            q = pt
            wheel = 210
            mq = mul_pt_exn(q, curve, wheel)
            jq_list = []
            for j in [k for k in range(1, wheel // 2) if gcd(k, wheel) == 1]:
                jq = mul_pt_exn(q, curve, j)
                jq_list.append(jq)
                res = gcd(jq[1], n)
                if 1 < res < n:
                    return res
            c = (b1 // wheel) * wheel
            cq = mul_pt_exn(q, curve, c)
            while c < b2 + wheel:
                s = cq[1] if cq[1] != 0 else 1
                for jq in jq_list:
                    if cq[0] != jq[0]:
                        s = s * (cq[0] - jq[0]) % n
                res = gcd(s, n)
                if 1 < res < n:
                    return res
                elif res == n:
                    res = gcd(cq[1], n)
                    if 1 < res < n:
                        return res
                    for jq in jq_list:
                        res = gcd(cq[0] - jq[0], n)
                        if 1 < res < n:
                            return res
                    # s is a multiple of n while each of cq[1] and {(cq[0] - jq[0]) % n} is not.
                    # There must be at least 2 non-trivial factors. The function should have returned.
                    assert False
                c += wheel
                cq = add_pt_exn(cq, mq, curve)
        except InverseNotFound as e:
            res = gcd(e.x, n)
            if 1 < res < n:
                return res
    return None


if __name__ == "__main__":
    random.seed(2)
    n = 310739457793333465418548557523014289  # (413198756866051421 * 752033864163021509)
    print(ecm(n, 100, 10000, 800000))
