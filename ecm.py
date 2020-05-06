import random
import numpy as np
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


def gcd(a, b):
    a = abs(a)
    b = abs(b)
    while a > 0:
        a, b = b % a, a
    return b


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
    sb = 1
    tb = 0
    sa = 0
    ta = 1
    while a > 0:
        (q, a), b = divmod(b, a), a
        sa, sb = sb - q * sa, sa
        ta, tb = tb - q * ta, ta
        assert sa * n + ta * (x % n) == a and sb * n + tb * (x % n) == b
    if b != 1:
        raise InverseNotFound(x % n, n)
    return (tb % n + n) % n


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


def add_pt(pt1, pt2, curve):
    """Adds two points pt1 and pt2 on curve.

    Args:
        pt1 (tuple(int, int)): Point (x1, y1). Use (None, None) for point at infinity.
        pt2 (tuple(int, int)): Point (x2, y2). Use (None, None) for point at infinity.
        curve (tuple(int, int, int)): (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).

    Returns:
        tuple(int, int): Point pt1 + pt2.
    """
    try:
        return add_pt_exn(pt1, pt2, curve)
    except InverseNotFound:
        return (None, None)


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
    res = (None, None)
    while k >= 1:
        if k % 2 == 1:
            res = add_pt_exn(res, point, curve)
        k //= 2
        point = add_pt_exn(point, point, curve)
    return res


def mul_pt(point, curve, k):
    """Multiplies point by k times on curve.

    Args:
        point (tuple(int, int)): Point (x, y). Use (None, None) for point at infinity.
        curve (tuple(int, int, int)): (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).
        k (int): Multiplier.

    Returns:
        tuple(int, int): Point k * point. Returns (None, None) for point at infinity.
    """
    res = (None, None)
    while k >= 1:
        if k % 2 == 1:
            res = add_pt(res, point, curve)
        k //= 2
        point = add_pt(point, point, curve)
    return res


def ecm(n, rounds, b):
    """Elliptic Curve Factorization Method.
    For each round:
        1. Generate random point and curve.
        2. Repeatedly multiply the current point by small primes raised to some power.
    Returns when a non-trivial factor is found.

    Args:
        n (int): Number to be factorized.
        rounds (int): Number of random curves to try.
        b (int): Bound for primes used in step 2.

    Returns:
        int: Non-trivial factor if found, otherwise returns None.
    """
    k_ls = []
    for p in PRIME_GEN(b):
        k_ls.append(p ** int(np.log(b) / np.log(p)))
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
            for k in k_ls:
                pt = mul_pt_exn(pt, curve, k)
        except InverseNotFound as e:
            res = gcd(e.x, n)
            if 1 < res < n:
                return res
    return None


if __name__ == "__main__":
    random.seed(2)
    n = 10648244288842058842742264007469181  # (103190330403778789 * 103190330403788729)
    print(ecm(n, 100, 10000))
