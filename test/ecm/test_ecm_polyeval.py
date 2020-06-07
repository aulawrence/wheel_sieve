import random
import unittest
from wheel_sieve.ecm.ecm_polyeval import (
    Polynomial,
    product_tree,
    recip_tree,
    remainder_tree,
)


def product_loop(f_list):
    prod = f_list[0]
    for f in f_list[1:]:
        prod *= f
    return prod


def remainder_loop(s_list, t_list, n):
    prod = 1
    for s in s_list:
        prod_t = 1
        for t in t_list:
            prod_t = prod_t * (s - t) % n
        if prod_t != 0:
            prod = prod * prod_t % n
    return prod


class TestProductTree(unittest.TestCase):
    def test1(self):
        n = 257
        f_list = [Polynomial([2], n)]
        target = product_loop(f_list)
        p_tree = product_tree(f_list, n)
        self.assertEqual(target, p_tree[0])

    def test2(self):
        random.seed(2)
        n = 310739457793333465418548557523014289
        f_list = [Polynomial([random.randint(0, n - 1), 1], n) for _ in range(100)]
        target = product_loop(f_list)
        p_tree = product_tree(f_list, n)
        self.assertEqual(target, p_tree[0])

    def test3(self):
        random.seed(2)
        n = 310739457793333465418548557523014289
        f_list = [
            Polynomial(
                [random.randint(0, n - 1) for _ in range(random.randint(1, 5))], n
            )
            for _ in range(100)
        ]
        target = product_loop(f_list)
        p_tree = product_tree(f_list, n)
        self.assertEqual(target, p_tree[0])


class TestRecipTree(unittest.TestCase):
    def test1(self):
        n = 257
        f_list = [Polynomial([2], n)]
        p_tree = product_tree(f_list, n)
        r_tree = recip_tree(p_tree)
        self.assertEqual(r_tree, [p.recip() for p in p_tree])

    def test2(self):
        random.seed(2)
        n = 310739457793333465418548557523014289
        f_list = [Polynomial([random.randint(0, n - 1), 1], n) for _ in range(100)]
        p_tree = product_tree(f_list, n)
        r_tree = recip_tree(p_tree)
        self.assertEqual(r_tree, [p.recip() for p in p_tree])

    def test3(self):
        random.seed(2)
        n = 310739457793333465418548557523014289
        f_list = [
            Polynomial(
                [random.randint(0, n - 1) for _ in range(random.randint(1, 5))], n
            )
            for _ in range(100)
        ]
        target = product_loop(f_list)
        p_tree = product_tree(f_list, n)
        r_tree = recip_tree(p_tree)
        self.assertEqual(r_tree, [p.recip() for p in p_tree])


class TestRemainderTree(unittest.TestCase):
    def test1(self):
        random.seed(2)
        n = 257
        s_list = [random.randint(1, n - 1) for _ in range(256)]
        t_list = [random.randint(1, n - 1) for _ in range(256)]
        target = remainder_loop(s_list, t_list, n)
        Fx = product_tree([Polynomial([n - i, 1], n) for i in t_list], n)[0]
        g_tree = product_tree([Polynomial([n - i, 1], n) for i in s_list], n)
        g_recip_tree = recip_tree(g_tree)
        r_tree = remainder_tree(Fx, g_tree, g_recip_tree, n)
        self.assertEqual(target, r_tree[0])

    def test2(self):
        random.seed(2)
        n = 310739457793333465418548557523014289
        s_list = [random.randint(1, n - 1) for _ in range(5)]
        t_list = [random.randint(1, n - 1) for _ in range(501)]
        target = remainder_loop(s_list, t_list, n)
        Fx = product_tree([Polynomial([n - i, 1], n) for i in t_list], n)[0]
        g_tree = product_tree([Polynomial([n - i, 1], n) for i in s_list], n)
        g_recip_tree = recip_tree(g_tree)
        r_tree = remainder_tree(Fx, g_tree, g_recip_tree, n)
        self.assertEqual(target, r_tree[0])

    def test3(self):
        random.seed(2)
        n = 310739457793333465418548557523014289
        s_list = [random.randint(1, n - 1) for _ in range(501)]
        t_list = [random.randint(1, n - 1) for _ in range(3)]
        target = remainder_loop(s_list, t_list, n)
        Fx = product_tree([Polynomial([n - i, 1], n) for i in t_list], n)[0]
        g_tree = product_tree([Polynomial([n - i, 1], n) for i in s_list], n)
        g_recip_tree = recip_tree(g_tree)
        r_tree = remainder_tree(Fx, g_tree, g_recip_tree, n)
        self.assertEqual(target, r_tree[0])

    def test4(self):
        random.seed(2)
        n = 257
        s_list = [2, 3]
        t_list = [4]
        target = remainder_loop(s_list, t_list, n)
        Fx = product_tree([Polynomial([n - i, 1], n) for i in t_list], n)[0]
        g_tree = product_tree([Polynomial([n - i, 1], n) for i in s_list], n)
        g_recip_tree = recip_tree(g_tree)
        r_tree = remainder_tree(Fx, g_tree, g_recip_tree, n)
        self.assertEqual(target, r_tree[0])


if __name__ == "__main__":
    unittest.main()
