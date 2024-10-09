"""Functions in support of EC class and Lenstra's EC factorization."""

from Point import identity

def extGCD(n, m):
    """Return the GCD and Bezout coefficients of two integers.

    Parameters
    ----------
    n, m : int : Non-zero integers, order does not matter.

    Returns
    -------
    res : tuple : Integers a, b, c such that n*a + m*b = gcd(n, m) = c.

    Example(s)
    ----------
    >>> extGCD(15, 2)
    >>> (1, -7, 1)

    >>> extGCD(12, 14)
    >>> (-1, 1, 2)
    
    """
    prevX, x = 0, 1
    prevY, y = 1, 0

    while n != 0:
        (q, n), m = divmod(m, n), n
        prevY, y = y, prevY - q * y
        prevX, x = x, prevX - q * x

    res = (prevX, prevY, m)
    return res

def modInv(n, m):
    """Attempt to find the multiplicative inverse of an element in Z/mZ.

    Parameters
    ----------
    n : int : Number to be inverted (not required to be principal value).
    m : int : Modulus being considered.

    Returns
    -------
    exists : bool : If the element is invertible.
    res    : int  : If invertible the inverse else the obstruction (=gcd(n,m)).

    Example(s)
    ----------
    >>> modInv(2, 6)
    >>> (False, 2)

    >>> modInv(3, 8)
    >>> (True, 3)
    
    """
    n %= m
    n, _, gcf = extGCD(n, m)
    if gcf != 1:
        exists = False
        res = gcf
    else:
        exists = True
        res = n % m
        
    return exists, res

def factorial(n):
    """Return n!."""
    return 1 if (n in {0, 1}) else n*factorial(n-1)

def isSmooth(a, b, m):
    """Return if (y**2 = x**3 + ax + b) for a, b in Z/mZ is smooth.

    Parameters
    ----------
    a, b : int : Coefficients used to define: (y**2 = x**3 + ax + b).
    m    : int : Defines the ring of scalars a and b live in.

    Returns
    -------
    - : bool : If the defined curve is smooth.
    
    """
    negDisc = 64 * a**3 + 432 * b**2
    return (negDisc % m != 0)

def modularSlope(P, Q):
    """Return the slope of the line connecting 2 distinct points in (Z/mZ)^2.

    Parameters
    ----------
    P, Q : Point : Unique points in (Z/mZ)^2.

    Returns
    -------
    exists : bool : If the slope is well-defined.
    res    : int  : If well-defined the slope else the obstruction (=gcd(n,m)).

    Notes
    -----
    In (Z/mZ)^2, the slope of a line connecting points P=(a,b) and Q=(x,y) is
    given by (b-y) * (x-a)**(-1). This is analogous to the classic
    interpretation of the slope, with the added nuisance that (x-a) may not
    be invertible in Z/mZ.

    Example(s)
    ----------
    >>> P = Point(2, 3, 5)
    >>> Q = Point(3, 1, 5)
    >>> modularSlope(P, Q)
    >>> (True, 3)

    """
    dy = Q.y - P.y
    dx = Q.x - P.x
    modulus = P.modulus
    
    exists, dxInv = modInv(dx, modulus)
    if exists:
        res = (dy * dxInv) % modulus
    else:
        res = dxInv
        
    return exists, res

def isCurveInd(P, Q=identity, k=1):
    """Return if [k]P + Q is curve-independent.

    Parameters
    ----------
    P, Q : Point : Unique points in (Z/mZ)^2.
    k    : int   : Multiple of P. 

    Returns
    -------
    - : bool : If [k]P and Q can be added without curve specific information.

    """
    return (P.ID or Q.ID or P == Q.inverse) and (k in {-1, 0, 1})

def resCurveInd(P, Q=identity, k=1):
    """Return [k]P + Q given the sum is curve-independent.

    Parameters
    ----------
    P, Q : Point : Unique points in (Z/mZ)^2.
    k    : int   : Multiple of P. 

    Returns
    -------
    res : Point : Is equal to [k]P + Q.

    Notes
    -----
    The logic in this function is based off of backend use cases, namely
    it is never the case: (k != 1) and (Q.ID == False).

    """
    if Q.ID and k == 1:
        res = P
    elif (P == Q.inverse) or (k == 0):
        res = identity
    elif k == -1:
        res = P.inverse
    if P.ID:
        res = Q

    return res

