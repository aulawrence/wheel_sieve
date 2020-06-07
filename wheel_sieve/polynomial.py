"""Polynomial arithmetic (mod n).
"""
from wheel_sieve.common import inv


class Polynomial(object):
    """Polynomial. f(x) = a0 + a1*x + a2*x**2 + ... + ad*x**d (mod n).

    Args:
        coeff (list(int)): Coefficient list [a0, a1, a2, ..., ad].
        n (int): Modulus.
        copy (bool, optional): Make a copy of the coefficient list. Defaults to False.
    """

    def __init__(self, coeff, n, copy=False):
        self.coeff = coeff if not copy else coeff.copy()
        self.n = n

    def __eq__(self, other):
        if not isinstance(other, Polynomial) or other.n != self.n:
            return False
        return self.coeff == other.coeff

    def __getitem__(self, slice_index):
        """Slice the coefficient list by slice_index, then return the Polynomial with the resultant coefficient list.

        Args:
            slice_index (int or slice): Index to access, or slice to apply on the original coefficients list.

        Returns:
            Polynomial: Polynomial with coefficient list as indexed/ sliced.
        """
        return Polynomial(self.coeff[slice_index], self.n)

    def __str__(self):
        return "Polynomial({}, mod {})".format(self.coeff, self.n)

    def __repr__(self):
        return "Polynomial({}, {})".format(self.coeff, self.n)

    def __add__(self, other):
        if not isinstance(other, Polynomial) or other.n != self.n:
            raise ValueError
        res = []
        i = 0
        while i < len(self.coeff) and i < len(other.coeff):
            res.append((self.coeff[i] + other.coeff[i]) % self.n)
            i += 1
        while i < len(self.coeff):
            res.append(self.coeff[i])
            i += 1
        while i < len(other.coeff):
            res.append(other.coeff[i])
            i += 1
        while len(res) > 1 and res[-1] == 0:
            res.pop()
        return Polynomial(res, self.n)

    def __sub__(self, other):
        if not isinstance(other, Polynomial) or other.n != self.n:
            raise ValueError
        res = []
        i = 0
        while i < len(self.coeff) and i < len(other.coeff):
            res.append((self.coeff[i] - other.coeff[i]) % self.n)
            i += 1
        while i < len(self.coeff):
            res.append(self.coeff[i])
            i += 1
        while i < len(other.coeff):
            res.append(self.n - other.coeff[i])
            i += 1
        while len(res) > 1 and res[-1] == 0:
            res.pop()
        return Polynomial(res, self.n)

    def __mul__(self, other):
        """Multiplies two polynomials self and other.
        Relies on fast integer multiplication in Python implementation.

        Args:
            other (Polynomial): Multiplier.

        Raises:
            ValueError: Thrown when other is not Polynomial or is incompatible.

        Returns:
            Polynomial: self * other.
        """
        if not isinstance(other, Polynomial) or other.n != self.n:
            raise ValueError
        d = max(len(self.coeff), len(other.coeff))
        k = (d * self.n ** 2 + 1).bit_length()
        k_8 = (k - 1) // 8 + 1
        k = k_8 * 8
        bt_self = bytes.join(
            b"", (ai.to_bytes(k_8, byteorder="little") for ai in self.coeff)
        )
        t_self = int.from_bytes(bt_self, byteorder="little")
        if self == other:
            t_other = t_self
        else:
            bt_other = bytes.join(
                b"", (ai.to_bytes(k_8, byteorder="little") for ai in other.coeff)
            )
            t_other = int.from_bytes(bt_other, byteorder="little")
        t_res = t_self * t_other
        res = []
        bt_res = t_res.to_bytes((t_res.bit_length() - 1) // 8 + 1, byteorder="little")
        i = 0
        while i < len(bt_res):
            res.append(int.from_bytes(bt_res[i : i + k_8], byteorder="little") % self.n)
            i += k_8
        return Polynomial(res, self.n)

    def __divmod__(self, other):
        """Divides polynomials self by other, return quotient and remainder.

        Args:
            other (Polynomial): Divisor.

        Raises:
            ValueError: Thrown when other is not Polynomial or is incompatible.

        Returns:
            tuple(Polynomial, Polynomial): (quotient, remainder).
        """
        if not isinstance(other, Polynomial) or other.n != self.n:
            raise ValueError
        pn = len(self.coeff)
        qn = len(other.coeff)
        if pn < qn:
            poly_quo = Polynomial([0], self.n)
            poly_rem = Polynomial(self.coeff, self.n, copy=True)
        else:
            other_recip = other.recip()
            d = qn - 1
            dividend = self
            poly_quo = Polynomial([0], self.n)
            while True:
                quo = (dividend[d:] * other_recip)[d:]
                poly_quo += quo
                rem = dividend - quo * other
                if len(rem.coeff) < len(other.coeff) or (
                    len(rem.coeff) == 1 and rem.coeff[0] == 0
                ):
                    break
                dividend = rem
            poly_rem = rem
        return poly_quo, poly_rem

    def mod_with_recip(self, other, other_recip):
        """Compute polynomial remainder self % other given reciprocal polynomial other_recip.

        Args:
            other (Polynomial): Modulus.
            other_recip (Polynomial): The reciprocal polynomial of the modulus.

        Raises:
            ValueError: Thrown when other is not Polynomial or is incompatible.

        Returns:
            Polynomial: Remainder.
        """
        if not isinstance(other, Polynomial) or other.n != self.n:
            raise ValueError
        pn = len(self.coeff)
        qn = len(other.coeff)
        if pn < qn:
            poly_rem = Polynomial(self.coeff, self.n, copy=True)
        else:
            d = qn - 1
            dividend = self
            while True:
                quo = (dividend[d:] * other_recip)[d:]
                rem = dividend - quo * other
                if len(rem.coeff) < len(other.coeff) or (
                    len(rem.coeff) == 1 and rem.coeff[0] == 0
                ):
                    break
                dividend = rem
            poly_rem = rem
        return poly_rem

    def recip(self):
        """Get the reciprocal polynomial. For f(x) of degree n, recip(f)(x) = :math:`\\lfloor{\\frac{x^{2n}}{f(x)}}\\rfloor`.
        Uses Montgomery's RECIP algorithm.

        Returns:
            Polynomial: The reciprocal polynomial.
        """
        d = len(self.coeff) - 1
        inv_fn = inv(self.coeff[d], self.n)
        R_curr = Polynomial([inv_fn], self.n)
        e_curr = -self.coeff[d - 1] * inv_fn % self.n
        k = 2
        while k < d * 2:
            R_prev = R_curr
            H = (R_prev * R_prev) * Polynomial(
                [self.coeff[d - k + j + 1] for j in range(k)], self.n
            )
            R_curr_coeff = [0 for _ in range(k // 2)]
            for ai in R_prev.coeff:
                R_curr_coeff.append(2 * ai)
            for j in range(k):
                R_curr_coeff[j] = (R_curr_coeff[j] - H.coeff[j + k - 2]) % self.n
            R_curr = Polynomial(R_curr_coeff, self.n)
            e_prev = e_curr
            if k == 2:
                e_curr = (e_prev * e_prev - self.coeff[d - k] * inv_fn) % self.n
            elif k <= d:
                e_curr = (
                    e_prev * e_prev
                    - H.coeff[k - 3] * self.coeff[d]
                    - self.coeff[d - k] * inv_fn
                ) % self.n
            k *= 2
        res = R_curr
        if k == d * 2:
            # Only needed when k is a power of 2.
            res.coeff.insert(0, e_curr * inv_fn % self.n)
        res.coeff = res.coeff[-d - 1 :]
        return res
