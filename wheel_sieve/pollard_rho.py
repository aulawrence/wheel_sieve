"""Pollard's Rho Algorithm on finding a factor of an integer.
"""
from math import gcd


def pollard(x0, g, n):
    """Pollard's Rho algorithm.
    With Pollard and Brent's improvement of computing gcd every 100 multiples.

    Args:
        x0 (int): Initial x.
        g (int->int): A polynomial modulo n.
        n (int): Number to be factored.

    Returns:
        int: A nontrivial factor if found, otherwise None is returned.
    """
    x = x0
    y = x0
    d = 1
    done = False
    while d == 1 and not done:
        for _ in range(100):
            x = g(x)
            y = g(g(y))
            e = d
            d = (d * (x - y)) % n
            if d == 0:
                d = e
                done = True
                break
        d = gcd(d, n)
    if d == n or d == 1:
        return None
    return d


def pollard2(x0, g, n):
    """Pollard's Rho algorithm.
    With Brent's cycle finding method.
    With Pollard and Brent's improvement of computing gcd every 100 multiples.

    Args:
        x0 (int): Initial x.
        g (int->int): A function representing a polynomial modulo n.
        n (int): Number to be factored.

    Returns:
        int: A nontrivial factor if found, otherwise None is returned.
    """
    x = x0
    y = g(x0)
    power = 1
    lam = 1
    d = 1
    done = False
    while d == 1 and not done:
        for _ in range(100):
            if power == lam:
                x = y
                power *= 2
                lam = 0
            y = g(y)
            lam += 1
            e = d
            d = (d * (x - y)) % n
            if d == 0:
                d = e
                done = True
                break
        d = gcd(d, n)
    if d == n or d == 1:
        return None
    return d


if __name__ == "__main__":
    num = 3176919178044924099992905369  # (41610855623299 * 76348326186931)
    x_init = 2
    print(pollard2(x_init, lambda x: (x * x + 3) % num, num))
