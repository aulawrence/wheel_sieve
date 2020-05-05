def gcd(a, b):
    a = abs(a)
    b = abs(b)
    while a > 0:
        a = b % a
        b = a
    return b


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
    while d == 1:
        for _ in range(100):
            x = g(x)
            y = g(g(y))
            d = (d * (x - y)) % n
        d = gcd(d, n)
    if d == n:
        return None
    return d


def pollard2(x0, g, n):
    """Pollard's Rho algorithm.
       With Brent's cycle finding method.
       With Pollard and Brent's improvement of computing gcd every 100 multiples.

    Args:
        x0 (int): Initial x.
        g (int->int): A polynomial modulo n.
        n (int): Number to be factored.

    Returns:
        int: A nontrivial factor if found, otherwise None is returned.
    """
    x = x0
    y = g(x0)
    power = 1
    lam = 1
    d = 1
    while d == 1:
        for _ in range(100):
            if power == lam:
                x = y
                power *= 2
                lam = 0
            y = g(y)
            lam += 1
            d = (d * (x - y)) % n
        d = gcd(d, n)
    if d == n:
        return None
    return d


if __name__ == "__main__":
    n = 10648244288842058842742264007469181  # (103190330403778789 * 103190330403788729)
    x = 2

    def g(x):
        return (x * x + 3) % n
    print(pollard2(x, g, n))
