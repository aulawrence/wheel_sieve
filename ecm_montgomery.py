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


def get_curve_suyama(sigma, n):
    """Given parameter sigma, generate an Elliptic Curve (mod n) and a point on it using Suyama's parametrization.
       The constructed cruve's group order is a multiple of 12, compared to 4 guaranteed for Montgomery Curves.

    Args:
        sigma (int): The sigma parameter.
        n (int): Modulus.

    Raises:
        CurveInitFail: Thrown when the curve generated by the given parameters fails the necessary conditions.

    Returns:
        tuple(tuple(int, int), tuple(int, int, int)): (Point, Curve), where
            Point = (x0, z0) in projective coordinates ignoring y.
            Curve = (A, s, n), representing B * (y/z) ** 2 == (x/z) ** 3 + A * (x/z) ** 2 + (x/z) (mod n),
                    ignoring B and y.
                    s = (A+2)/4 % n is precomputed for point doubling.
    """
    if sigma % n in (n - 5, n - 3, n - 1, 0, 1, 3, 5) or sigma * 3 % n in (n - 5, 5):
        raise CurveInitFail()
    u = sigma ** 2 - 5 % n
    v = 4 * sigma % n
    x0 = u ** 3 % n
    z0 = v ** 3 % n
    A = ((v - u) ** 3 * (3 * u + v) * inv(4 * u**3 * v, n) - 2) % n
    if A in (n - 2, 2):
        raise CurveInitFail()
    s = (A + 2) * inv(4, n) % n
    # For completeness...
    # B = u * inv(z0, n) % n
    # y = (sigma ** 2 - 1) * (sigma ** 2 - 25) * (sigma ** 4 - 25) % n
    # x0_norm = (x0 * inv(z0, n)) % n
    # y0_norm = (y * inv(z0, n)) % n
    # assert B * y0_norm ** 2 % n == (x0_norm ** 3 + A * x0_norm ** 2 + x0_norm) % n
    return (x0, z0), (A, s, n)


def get_curve_a(x, A, n):
    """Given parameters x and A, generate an Elliptic Curve (mod n) and a point on it.

    Args:
        x (int): Desited x coordinate of the point.
        A (int): Parameter A of Montgomery Curve.
        n (int): Modulus.

    Raises:
        CurveInitFail: Thrown when the curve generated by the given parameters fails the necessary conditions.

    Returns:
        tuple(tuple(int, int), tuple(int, int, int)): (Point, Curve), where
            Point = (x0, z0) in projective coordinates ignoring y.
            Curve = (A, s, n), representing B * (y/z) ** 2 == (x/z) ** 3 + A * (x/z) ** 2 + (x/z) (mod n),
                    ignoring B and y.
                    s = (A+2)/4 % n is precomputed for point doubling.
    """
    if A % n in (n - 2, 2):
        raise CurveInitFail()
    x0 = x % n
    z0 = 1
    s = (A + 2) * inv(4, n) % n
    # For completeness...
    # x0_norm = x0
    # y0_norm = 2
    # B = (x0_norm ** 3 + A * x0_norm ** 2 + x0_norm) * inv(y0_norm ** 2, n) % n
    # assert B * y0_norm ** 2 % n == (x0_norm ** 3 + A * x0_norm ** 2 + x0_norm) % n
    return (x0, z0), (A, s, n)


def add_pt(ptp, ptq, pt_, curve):
    """Computes point P+Q given points P, Q and P-Q, and curve.
       Does not return correct result when P == Q, use dbl_pt instead.

    Args:
        ptp (tuple(int, int)): Point P.
        ptq (tuple(int, int)): Point Q.
        pt_ (tuple(int, int)): Point P-Q.
        curve (tuple(int, int, int)): Curve.

    Returns:
        (tuple(int, int)): Point P+Q.
    """
    xp, zp = ptp
    xq, zq = ptq
    x_, z_ = pt_
    A, s, n = curve
    u = (xp - zp) * (xq + zq) % n
    v = (xp + zp) * (xq - zq) % n
    xr = z_ * (u + v) ** 2 % n
    zr = x_ * (u - v) ** 2 % n
    return (xr, zr)


def dbl_pt(pt, curve):
    """Computes point 2P given point P and curve.

    Args:
        pt (tuple(int, int)): Point P.
        curve (tuple(int, int, int)): Curve.

    Returns:
        (tuple(int, int)): Point 2P.
    """
    x, z = pt
    A, s, n = curve
    a = (x + z) ** 2 % n
    b = (x - z) ** 2 % n
    t = a - b
    xr = a * b % n
    zr = t * (b + s * t) % n
    return (xr, zr)


def mul_pt(pt, curve, k):
    """Computes point kP given point P, curve and k using Montgomery Ladder.

    Args:
        pt (tuple(int, int)): Point P.
        curve (tuple(int, int, int)): Curve.
        k (int): Multiplier.

    Returns:
        (tuple(int, int)): Point kP.
    """
    if k <= 1:
        if k < 0:
            # x and z coordinates are the same for P and -P.
            return mul_pt(pt, curve, -k)
        if k == 0:
            return (0, 0)
        return pt
    res0 = pt
    res1 = dbl_pt(pt, curve)
    j = k.bit_length() - 2
    while j >= 0:
        if (k >> j) % 2 == 1:
            res0 = add_pt(res1, res0, pt, curve)
            res1 = dbl_pt(res1, curve)
        else:
            res1 = add_pt(res1, res0, pt, curve)
            res0 = dbl_pt(res0, curve)
        j -= 1
    return res0


def check(pt, curve):
    """Given point P (x, z), check that P is not the point at infinity, i.e. gcd(z, n) == 1, and return P.
       If gcd(z, n) > 1, throws InverseNotFound.

    Args:
        pt (tuple(int, int)): Point P.
        curve (tuple(int, int, int)): Curve.

    Raises:
        InverseNotFound: Thrown when point P is the point at infinity.

    Returns:
        tuple(int, int): Point P.
    """
    x, z = pt
    A, s, n = curve
    if gcd(z, n) > 1:
        raise InverseNotFound(z, n)
    return pt


def ecm(n, rounds, b1, b2):
    """Elliptic Curve Factorization Method.
    For each round:
        0. Generate random point and curve.
        1. Repeatedly multiply the current point by small primes raised to some power, determined by b1.
        2. Repeatedly try to multiply the point from step 1 by possible primes (with wheel of 2310) between b1 and b2.
    Returns when a non-trivial factor is found.

    Args:
        n (int): Number to be factorized. n >= 12.
        rounds (int): Number of random curves to try.
        b1 (int): Bound for primes used in step 1.
        b2 (int): Bound for primes searched for in step 2. b1 < b2.

    Returns:
        int: Non-trivial factor if found, otherwise returns None.
    """
    assert n >= 12
    k_ls = []
    for p in PRIME_GEN(b1):
        k_ls.append(p ** int(round(np.log(b1) / np.log(p))))
    for roundi in range(rounds):
        print("Round {}...".format(roundi))
        count = 0
        success = False
        while not success and count < 20:
            try:
                count += 1
                sigma = random.randint(6, n - 6)
                pt, curve = get_curve_suyama(sigma, n)
                success = True
            except InverseNotFound as e:
                res = gcd(e.x, n)
                if 1 < res < n:
                    return res
            except CurveInitFail:
                pass
        if not success:
            print(" - Curve Init Failed.")
            break
        try:
            # Step 1
            print(" - Step 1")
            for k in k_ls:
                pt = check(mul_pt(pt, curve, k), curve)
            # Step 2
            print(" - Step 2")
            q = pt
            wheel = 2310
            mq = check(mul_pt(q, curve, wheel), curve)
            xj_list = []
            for j in [k for k in range(1, wheel // 2) if gcd(k, wheel) == 1]:
                xj, zj = check(mul_pt(q, curve, j), curve)
                xj_list.append(xj * inv(zj, n) % n)
            c = (b1 // wheel) * wheel
            cq = check(mul_pt(q, curve, c), curve)
            cq_ = check(mul_pt(q, curve, c - wheel), curve)
            while c < b2 + wheel:
                s = 1
                for xj in xj_list:
                    t = (xj * cq[1] - cq[0]) % n
                    if t != 0:
                        s = s * t % n
                res = gcd(s, n)
                if 1 < res < n:
                    return res
                elif res == n:
                    for xj in xj_list:
                        res = gcd(xj * cq[1] - cq[0], n)
                        if 1 < res < n:
                            return res
                    # s is a multiple of n while each of {(xj *  cq[1] - cq[0]) % n} is not.
                    # There must be at least 2 non-trivial factors. The function should have returned.
                    assert False
                c += wheel
                cq, cq_ = check(add_pt(cq, mq, cq_, curve), curve), cq
        except InverseNotFound as e:
            res = gcd(e.x, n)
            if 1 < res < n:
                return res
    return None


if __name__ == "__main__":
    random.seed(2)
    n = 32795254512039893688982510946402035163613484930081  # (4110718401119136339711691 * 7977986160061812908308291)
    print(ecm(n, 300, 50_000, 5_000_000))
