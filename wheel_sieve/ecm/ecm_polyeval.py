import random
import time
from math import gcd
import numpy as np
from wheel_sieve.common import PRIME_GEN, InverseNotFound, CurveInitFail
from wheel_sieve.ecm.ecm_brent_suyama import (
    apply_polynomial,
    get_difference_seq,
    step_difference_seq_exn,
)
from wheel_sieve.polynomial import Polynomial
import wheel_sieve.ecm.ecm_montgomery as mnt
import wheel_sieve.ecm.ecm_weierstrass as wst


def product_tree(poly_list, n):
    """Product Tree Algorithm. Multiply a list of polynomials, poly_list.
    The leaf nodes are populated by polynomials in poly_list. If the length of poly_list is not a power of 2,
    the remaining leaf nodes are takes the value of the polynomial f(x) = 1.
    The value of each internal node is the product of its children.

    Args:
        poly_list (list of Polynomial): List of polynomials to multiply.
        n (int): Modulus.

    Returns:
        list of Polynomial: The Product Tree, a complete binary tree in list form.
            The root node is at position 0 of the list. The children of node i are node 2*i+1 and node 2*i+2.
    """
    poly_num = len(poly_list)
    k = 1
    while k < poly_num:
        k *= 2
    res = [None for _ in range(k * 2 - 1)]
    for i in range(poly_num):
        res[k + i - 1] = poly_list[i]
    for i in range(poly_num, k):
        res[k + i - 1] = Polynomial([1], n)
    while k > 1:
        for i in range(k // 2 - 1, k - 1):
            res[i] = res[2 * i + 1] * res[2 * i + 2]
        k //= 2
    return res


def recip_tree(prod_tree):
    """Recip Tree Algorithm. Compute the reciprocal polynomail of each element in the given product tree.

    Args:
        prod_tree (list of Polynomial): Product Tree.

    Returns:
        list of Polynomial: The Recip Tree, a complete binary tree in list form.
            The root node is at position 0 of the list. The children of node i are node 2*i+1 and node 2*i+2.
            r_tree[i] = prod_tree[i].recip()
    """
    r_tree = [prod_tree[0].recip()]
    for i in range(len(prod_tree) // 2):
        gi_recip = r_tree[i]
        g1 = prod_tree[2 * i + 1]
        g2 = prod_tree[2 * i + 2]
        d1 = len(g1.coeff) - 1
        d2 = len(g2.coeff) - 1
        g1_recip = (gi_recip[d2:] * g2)[d2:]
        r_tree.append(g1_recip)
        g2_recip = (gi_recip[d1:] * g1)[d1:]
        r_tree.append(g2_recip)
    return r_tree


def remainder_tree(f, g_tree, g_recip_tree, n):
    """Remainder Tree Algorithm. Given polynomials f, g_1, g_2, ..., g_m where g_i(x) = x - x_i,
    use a product tree to compute :math:\prod_{i=0}^{m}(f \mod g_i) \mod n, which is :math:\prod_{i=0}^{m} f(x_i) \mod n.
    If one f(x_i) is zero, the product will be zero. To prevent this, a slight change is added to omit those terms from the product.
    The leaf nodes are populated by f mod g_i, which is f(x_i). If the length of poly_list is not a power of 2,
    the remaining leaf nodes takes the value of 1.
    The value of each internal node is the product of its children mod n.

    Args:
        f (Polynomial): Polynomial f.
        g_tree (list of Polynomial): Product Tree of polynomial [g_1, g_2, ..., g_m]. Each of g_i is assumed to be of degree 1.
        g_recip_tree (list of Polynomial): Recip Tree of g_tree.
        n (int): Modolus.

    Returns:
        list of Polynomial: The Remainder Tree, a complete binary tree in list form.
            The root node is at position 0 of the list. The children of node i are node 2*i+1 and node 2*i+2.
    """
    f_mod_g_tree = []
    for i in range(len(g_tree)):
        gi = g_tree[i]
        gi_recip = g_recip_tree[i]
        if i == 0:
            f_mod_g = f
        else:
            f_mod_g = f_mod_g_tree[(i - 1) // 2]
        while True:
            di = len(gi.coeff) - 1
            hi = (f_mod_g[di:] * gi_recip)[di:]
            f_mod_gi = f_mod_g - gi * hi
            if len(f_mod_gi.coeff) < len(gi.coeff) or (
                len(f_mod_gi.coeff) == 1 and f_mod_gi.coeff[0] == 0
            ):
                break
            f_mod_g = f_mod_gi
        f_mod_g_tree.append(f_mod_gi)
    # assert len(g_recip_tree) == len(g_tree)
    # assert all(g_recip == g.recip() for g_recip, g in zip(g_recip_tree, g_tree))
    k = len(f_mod_g_tree) + 1
    for i in range(k // 2 - 1, k - 1):
        f_mod_g_tree[i] = f_mod_g_tree[i].coeff[0]
        if f_mod_g_tree[i] == 0:
            f_mod_g_tree[i] = 1
    k //= 2
    while k > 1:
        for i in range(k // 2 - 1, k - 1):
            f_mod_g_tree[i] = f_mod_g_tree[2 * i + 1] * f_mod_g_tree[2 * i + 2] % n
        k //= 2
    return f_mod_g_tree


def ecm(n, rounds, b1, b2, wheel=2310, output=True):
    """Elliptic Curve Factorization Method.
    For each round:
        0. Generate random point and curve.
        1. Repeatedly multiply the current point by small primes raised to some power, determined by b1.
        2. Standard continuation from b1 to b2 with Brent-Suyama's Extension and Polyeval.
    Returns when a non-trivial factor is found.

    Args:
        n (int): Number to be factorized. n >= 12.
        rounds (int): Number of random curves to try.
        b1 (int): Bound for primes used in step 1.
        b2 (int): Bound for primes searched for in step 2. b1 < b2.
        wheel (int, optional): Wheel, where only numbers coprime to wheel will be considered in step 2. Defaults to 2310.
        output (bool, optional): Whether to print progress to stdout. Defaults to True.

    Returns:
        int: Non-trivial factor if found, otherwise returns None.
    """
    assert n >= 12
    j_list = [j for j in range(1, wheel // 2) if gcd(j, wheel) == 1]
    block_size = 1 << (len(j_list) - 1).bit_length() - 1
    for round_i in range(rounds):
        if output:
            st = time.time()
            print("Round {}...".format(round_i))
        count = 0
        success = False
        while not success and count < 20:
            try:
                count += 1
                sigma = random.randint(6, n - 6)
                mnt_pt, mnt_curve = mnt.get_curve_suyama(sigma, n)
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
            if output:
                print("{:>5.2f}: Step 1".format(time.time() - st))
            for p in PRIME_GEN(b1):
                for _ in range(int(np.log(b1) / np.log(p))):
                    mnt_pt = mnt.mul_pt_exn(mnt_pt, mnt_curve, p)
            # Step 2
            if output:
                print("{:>5.2f}: Step 2".format(time.time() - st))
            polynomial = (2, 0, 9, 0, 6, 0, 1)  # f(x) = x^6 + 6x^4 + 9x^2 + 2
            q, wst_curve = mnt.to_weierstrass(mnt_pt, mnt_curve)
            c1 = b1 // wheel
            c2 = b2 // wheel + 2
            c = 0
            k_ls = [
                apply_polynomial(polynomial, j) for j in j_list
            ] + get_difference_seq(polynomial, c1 * wheel, wheel)
            mul_res = wst.mul_pt_multi(q, wst_curve, k_ls)
            xj_list = []
            for i in range(len(j_list)):
                xj_list.append(mul_res[i][0])
            cq_list = mul_res[len(j_list) :]
            f_tree = product_tree([Polynomial([n - xj, 1], n) for xj in xj_list], n)
            f_recip_tree = recip_tree(f_tree)
            H = Polynomial([1], n)
            g_poly_list = []
            while c < c2 - c1:
                for _ in range(min(block_size, c2 - c1 - c)):
                    g_poly_list.append(Polynomial([n - cq_list[0][0], 1], n))
                    step_difference_seq_exn(cq_list, wst_curve)
                    c += 1
                G = product_tree(g_poly_list, n)[0]
                H = (H * G).mod_with_recip(f_tree[0], f_recip_tree[0])
                g_poly_list.clear()
            rem_tree = remainder_tree(H, f_tree, f_recip_tree, n)
            res = gcd(rem_tree[0], n)
            if 1 < res < n:
                return res
            elif res == n:
                for rem in rem_tree[len(rem_tree) // 2 :]:
                    res = gcd(rem, n)
                    if 1 < res < n:
                        return res
                assert False
            if output:
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
