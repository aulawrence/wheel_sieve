import random
from math import gcd
import numpy as np
from ecm_common import PRIME_GEN, InverseNotFound, inv


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
    a, b, n = curve
    gen = add_pt_gen(pt1, pt2, curve)
    inv_req = gen.send(None)
    if inv_req is None:
        return gen.send(None)
    else:
        return gen.send(inv(inv_req, n))


def add_pt_gen(pt1, pt2, curve):
    """Adds two points pt1 and pt2 on curve. Returns a generator that generates exactly 2 elements.
    Step 0:
        Send None to start generator.
    Step 1: 
        Generator yields the number to be inverted (mod n), or None if no inversion is required.
        Send the inverted number to generator, or None if no inversion is required.
    Step 2:
        Generator yields the point pt1 + pt2.

    Args:
        pt1 (tuple(int, int)): Point (x1, y1). Use (None, None) for point at infinity.
        pt2 (tuple(int, int)): Point (x2, y2). Use (None, None) for point at infinity.
        curve (tuple(int, int, int)): (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).

    Yields:
        int: The number to be inverted (mod n), or None if no inversion is required.
        tuple(int, int): Point pt1+pt2.
    """
    x1, y1 = pt1
    x2, y2 = pt2
    a, b, n = curve
    if pt1 == (None, None):
        yield None
        yield pt2
        return
    elif pt2 == (None, None):
        yield None
        yield pt1
        return
    elif pt1 == pt2:
        res = yield 2 * y1
        s = (3 * x1 * x1 + a) * res % n
    else:
        res = yield x2 - x1
        s = (y2 - y1) * res % n
    xr = (s * s - x1 - x2) % n
    yr = (s * (x1 - xr) - y1) % n
    yield (xr, yr)


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


def mul_pt_gen(point, curve, k):
    """Multiplies point by k times on curve. Returns a generator that generates multiple elements.
    Step 0:
        Send None to start generator.
    Step 1: 
        Generator yields the number to be inverted (mod n), or None if no inversion is required.
        Send the inverted number to generator, or None if no inversion is required.
        If the generator did not yield None, Repeat step 1.
    Step 2:
        Generator yields the point k * point.

    Args:
        point (tuple(int, int)): Point (x, y). Use (None, None) for point at infinity.
        curve (tuple(int, int, int)): (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).
        k (int): Multiplier.

    Yields:
        int: The number to be inverted (mod n), or None if no inversion is required.
        tuple(int, int): Point k * point.
    """
    if k < 0:
        yield from mul_pt_gen(neg_pt(point, curve), curve, -k)
        return
    res = (None, None)
    while k >= 1:
        if k % 2 == 1:
            gen = add_pt_gen(res, point, curve)
            inv_req = gen.send(None)
            if inv_req is None:
                res = gen.send(None)
            else:
                inv_res = yield inv_req
                res = gen.send(inv_res)
        k //= 2
        gen = add_pt_gen(point, point, curve)
        inv_req = gen.send(None)
        if inv_req is None:
            point = gen.send(None)
        else:
            inv_res = yield inv_req
            point = gen.send(inv_res)
    yield None
    yield res


def mul_pt_multi(point, curve, k_ls):
    """Multiply point by k times on curve, for k in k_ls.
    Uses Montgomery's trick to reduce number of modular inversions.

    Args:
        point (tuple(int, int)): Point (x, y). Use (None, None) for point at infinity.
        curve (tuple(int, int, int)): (a, b, n) representing the Elliptic Curve y**2 = x**3 + a*x + b (mod n).
        k_ls (list of int): List of multiplier.

    Raises:
        InverseNotFound: Thrown when a number cannot be inverted during the calculation.

    Returns:
        list of tuple(int, int): List of points [k * point for k in k_ls].
    """
    a, b, n = curve
    gen_list = [mul_pt_gen(point, curve, k) for k in k_ls]
    denom_list = [gen.send(None) for gen in gen_list]
    working_set = set()
    for i in range(len(k_ls)):
        if denom_list[i] is not None:
            working_set.add(i)
    denom_set = set()
    denom_inv_dict = dict()
    while working_set:
        prod = 1
        for i in working_set:
            if denom_list[i] is not None and denom_list[i] not in denom_set:
                denom_set.add(denom_list[i])
                prod = prod * denom_list[i]
        prod_inv = inv(prod, n)
        for denom in denom_set:
            denom_inv_dict[denom] = prod_inv * (prod // denom) % n
        remove_set = set()
        for i in working_set:
            while denom_list[i] is not None and denom_list[i] in denom_set:
                denom_list[i] = gen_list[i].send(denom_inv_dict[denom_list[i]])
            if denom_list[i] is None:
                remove_set.add(i)
        working_set -= remove_set
        denom_set.clear()
        denom_inv_dict.clear()
    return [gen.send(None) for gen in gen_list]


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
        k_ls.append(p ** int(np.log(b1) / np.log(p)))
    for round_i in range(rounds):
        print("Round {}...".format(round_i))
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
            print(" - Curve Init Failed.")
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
            print(" - End")
        except InverseNotFound as e:
            res = gcd(e.x, n)
            if 1 < res < n:
                return res
    return None


if __name__ == "__main__":
    random.seed(2)
    n = 310739457793333465418548557523014289  # (413198756866051421 * 752033864163021509)
    print(ecm(n, 100, 10000, 800000))
