# Lenstra's Elliptic Curve Factorization

This project implements **Lenstra's Elliptic Curve Factorization Algorithm** (ECM) in Python. ECM is a number-theoretic algorithm used for integer factorization, specifically designed to find smaller prime factors of large composite numbers. It uses properties of elliptic curves over finite fields and is especially efficient in finding smaller factors (20–25 digits) of large numbers.

### Features

- Implementation of Lenstra’s ECM to factor large composite numbers.
- Ability to handle multiple elliptic curves to increase the chances of successful factorization.
- Modular arithmetic implemented for operations on elliptic curves.

### Example

```python
N = (3209622181 * 6727426213 * 2810645183)
print(lenstraFactorial(N))
>>> 2810645183 #Takes ~ (0.5822016000020085 seconds)
```
