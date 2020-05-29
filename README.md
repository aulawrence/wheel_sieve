# Wheel Sieve and More

## Files
|File|Description|Remark|
|--|--|--|
|[ecm_brent_suyama.py](ecm_brent_suyama.py)|ECM: Stage 1 in Montgomery Form, Stage 2 Standard Continuation with Brent-Suyama's Extension||
|[ecm_montgomery.py](ecm_montgomery.py)|Lenstra Elliptic Curve Factorization in Montgomery Form and XZ coordinates||
|[ecm_polyeval.py](ecm_polyeval.py)|ECM: Stage 1 in Montgomery Form, Stage 2 Standard Continuation with Brent-Suyama's Extension and Polyeval||
|[ecm_weierstrass.py](ecm_weierstrass.py)|Lenstra Elliptic Curve Factorization in Weierstrass Form and XY coordinates|Slower than Montgomery Curve due to high cost of inverse.|
|[miller_rabin.py](miller_rabin.py)|Miller-Rabin Prime Test||
|[pollard_rho.py](pollard_rho.py)|Pollard's Rho Algorithm||
|[wheel_sieve_bit.py](wheel_sieve_bit.py)|Prime sieving with wheel factorization (bit array)|Implementation is slower than the byte version with the same memory constraint.|
|[wheel_sieve_byte.py](wheel_sieve_byte.py)|Prime sieving with wheel factorization||

## Branches
|Branch|Description|
|--|--|
|master|Master branch|
|gmp|Uses gmpy2 to speed up integer arithmetic|

## Dependencies
Python 3.5+

As listed in [requirements.txt](requirements.txt):
```
gmpy2
numpy
```

## Referenced Links

### Elliptic Curve Method
- https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization
- https://www.rieselprime.de/ziki/Elliptic_curve_method
- https://en.wikipedia.org/wiki/Montgomery_curve
- https://www.alpertron.com.ar/ECM.HTM

### Others
- https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test
- https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
- https://en.wikipedia.org/wiki/Wheel_factorization

## Referenced Papers

### Elliptic Curve Method
 - Zimmermann, P., & Dodson, B. (2006, July). 20 years of ECM. In International Algorithmic Number Theory Symposium (pp. 525-542). Springer, Berlin, Heidelberg.
   - [Link](https://hal.inria.fr/inria-00070192v2)
 - Brent, R., Kruppa, A., & Zimmermann, P. (2017). FFT extension for algebraic-group factorization algorithms.
   - [Link](https://hal.inria.fr/hal-01630907/document)
 - Montgomery, P. L. (1992). An FFT extension of the elliptic curve method of factorization (Doctoral dissertation, UCLA).
   - [Link](http://cr.yp.to/bib/1992/montgomery.ps)
