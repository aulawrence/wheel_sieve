# Wheel Sieve and More

## Files
|File|Description|Remark|
|--|--|--|
|ecm_montgomery.py|Lenstra Elliptic Curve Factorization in Montgomery Form||
|ecm_weierstrass.py|Lenstra Elliptic Curve Factorization in Weierstrass Form|Slower than Montgomery Curve due to high cost of inverse.|
|miller_rabin.py|Miller-Rabin Prime Test||
|pollard_rho.py|Pollar's Rho Algorithm||
|wheel_sieve_bit.py|Prime sieving with wheel factorization (bit array)|Implementation is slower than the byte version with the same memory constraint.|
|wheel_sieve_byte.py|Prime sieving with wheel factorization||

## Dependencies
Python 3.5+
```
numpy
```

## Referenced Links

### Elliptic Curve Method
- https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization
- https://www.rieselprime.de/ziki/Elliptic_curve_method
- https://en.wikipedia.org/wiki/Montgomery_curve
- https://www.alpertron.com.ar/ECM.HTM (Click Help!)

### Others
- https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test
- https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
- https://en.wikipedia.org/wiki/Wheel_factorization

## Referenced Papers

### Elliptic Curve Method
 - Paul Zimmermann, Bruce Dodson. 20 years of ECM. 7th Algorithmic Number Theory Symposium (ANTS VII), 2006, Berlin/Germany, Germany. pp.525–542. inria-00070192v2
   - [Link](https://hal.inria.fr/inria-00070192v2)
