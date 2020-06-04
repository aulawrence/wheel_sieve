import random
import time
from math import gcd
import numpy as np
from ecm_common import PRIME_GEN, InverseNotFound, CurveInitFail, inv, init_wheel


def get_curve_suyama(sigma, n):
    """Given parameter sigma, generate an Elliptic Curve (mod n) and a point on it using Suyama's parametrization.
       The constructed curve's group order is a multiple of 12, compared to 4 guaranteed for Montgomery Curves.

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
        x (int): Desired x coordinate of the point.
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
    xr = z_ * ((u + v) ** 2 % n) % n
    zr = x_ * ((u - v) ** 2 % n) % n
    return (xr, zr)


def to_weierstrass(pt, curve):
    """Given a point P and an Montgomery Curve it is on, computes the equivalent point and curve in weierstrass form.
    Note: Multiple calls for same curve with different P will produce different output curves. This
        is due to y-coordinates being omitted in the representation. Without the ability to square-root
        y (mod n) by fixing B, the natural thing to do is to fix y and calculate B. So different point P
        produces different B.

    Args:
        pt (tuple(int, int)): Point P in XZ form.
        curve (tuple(int, int, int)): Curve in Montgomery form.

    Returns:
        tuple(tuple(int, int), tuple(int, int, int)): (Point, Curve), where
            Point = (t, v) in XY form.
            Curve = (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).
    """
    x, z = pt
    A, s, n = curve
    y_norm = 1
    x_norm = x * inv(z, n)
    B = (x_norm ** 3 + A * x_norm ** 2 + x_norm) % n
    assert B * y_norm ** 2 % n == (x_norm ** 3 + A * x_norm ** 2 + x_norm) % n
    B_inv = inv(B, n)
    three_inv = inv(3, n)
    t = (x_norm * B_inv + A * three_inv * B_inv) % n
    v = (y_norm * B_inv) % n
    a = (3 - A ** 2) * three_inv * B_inv * B_inv % n
    b = (2 * A ** 3 - 9 * A) * (three_inv * B_inv % n) ** 3 % n
    assert v ** 2 % n == (t ** 3 + a * t + b) % n
    return (t, v), (a, b, n)


def add_pt_exn(ptp, ptq, pt_, curve):
    """Computes point P+Q given points P, Q and P-Q, and curve.
       Does not return correct result when P == Q, use dbl_pt instead.

    Args:
        ptp (tuple(int, int)): Point P.
        ptq (tuple(int, int)): Point Q.
        pt_ (tuple(int, int)): Point P-Q.
        curve (tuple(int, int, int)): Curve.

    Raises:
        InverseNotFound: Thrown when point P+Q is the point at infinity.

    Returns:
        (tuple(int, int)): Point P+Q.
    """
    return check(add_pt(ptp, ptq, pt_, curve), curve)


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
    zr = t * ((b + s * t) % n) % n
    return (xr, zr)


def mul_pt_exn(pt, curve, k):
    """Computes point kP given point P, curve and k using Montgomery Ladder.

    Args:
        pt (tuple(int, int)): Point P.
        curve (tuple(int, int, int)): Curve.
        k (int): Multiplier.

    Raises:
        InverseNotFound: Thrown when point kP is the point at infinity.

    Returns:
        (tuple(int, int)): Point kP.
    """
    if k <= 2:
        if k < 0:
            # x and z coordinates are the same for P and -P.
            return mul_pt_exn(pt, curve, -k)
        if k == 0:
            # InverseNotFound will be thrown
            return check((0, 0), curve)
        if k == 1:
            return check(pt, curve)
        return check(dbl_pt(pt, curve), curve)
    res0 = pt
    res1 = dbl_pt(pt, curve)
    j = k.bit_length() - 2
    while j >= 1:
        if (k >> j) % 2 == 1:
            res0 = add_pt(res1, res0, pt, curve)
            res1 = dbl_pt(res1, curve)
        else:
            res1 = add_pt(res1, res0, pt, curve)
            res0 = dbl_pt(res0, curve)
        j -= 1
    if k % 2 == 1:
        res0 = add_pt(res1, res0, pt, curve)
    else:
        res0 = dbl_pt(res0, curve)
    return check(res0, curve)


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
        2. Repeatedly try to multiply the point from step 1 by primes (with wheel of 2310) between b1 and b2.
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
    wheel = 2310
    st = time.time()
    j_list, prime_array = init_wheel(b1, b2, wheel)
    print("Init time: {:.2f}".format(time.time() - st))
    for round_i in range(rounds):
        st = time.time()
        print("Round {}...".format(round_i))
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
            print("{:>5.2f}: Step 1".format(time.time() - st))
            for p in PRIME_GEN(b1):
                for _ in range(int(np.log(b1) / np.log(p))):
                    pt = mul_pt_exn(pt, curve, p)
            # Step 2
            print("{:>5.2f}: Step 2".format(time.time() - st))
            q = pt
            mq = mul_pt_exn(q, curve, wheel)
            xj_list = []
            for j in j_list:
                xj, zj = mul_pt_exn(q, curve, j)
                xj_list.append(xj * inv(zj, n) % n)
            c1 = b1 // wheel
            c2 = b2 // wheel + 2
            c = 0
            cq = mul_pt_exn(q, curve, c1 * wheel)
            cq_ = mul_pt_exn(q, curve, (c1 - 1) * wheel)
            while c < c2 - c1:
                s = 1
                for xj, is_prime in zip(xj_list, np.unpackbits(prime_array[c, :], bitorder="little")):
                    if is_prime:
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
                c += 1
                cq, cq_ = add_pt_exn(cq, mq, cq_, curve), cq
            print("{:>5.2f}: End".format(time.time() - st))
        except InverseNotFound as e:
            res = gcd(e.x, n)
            if 1 < res < n:
                return res
    return None


if __name__ == "__main__":
    random.seed(2)
    n = 294636370796972331405770334382449402989049465216208991677129  # (406724252548875212358759885439 * 724413085648406196306771670711)
    print(ecm(n, 430, 250_000, 40_000_000))
