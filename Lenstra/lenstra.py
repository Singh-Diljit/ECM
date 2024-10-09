"""Two Implementations of Lenstra's elliptic-curve factorization."""

from random import randint
from Point import Point
from EllipticCurve import EC
from helperFuncs import factorial, isSmooth

def generateCurvePoint(m, effort=10**5):
    """Attempt to return a random point on an elliptic curve over Z/mZ.

    Parameters
    ----------
    m      : int : Modulus being considered, defines ring of scalars Z/mZ.
    effort : int : Number of attempts to generate a point and smooth curve.

    Returns
    -------
    P     : Point : Random point in (Z/mZ)^2.
    curve : EC    : Random elliptic curve over Z/mZ such that P is on it.

    Raises
    ------
    Exception if no elliptic curve is found - none of the generated Weierstrass
    curves were non-singular.

    """
    haveCurve = False
    for _ in range(effort):
        #Generate a random point in (Z/mZ)^2
        s, t = [randint(0, m-1) for _ in range(2)]

        #Given (s, t) is on the Weierstrass curve there is 1 degree of freedom
        a = randint(0, m-1)
        b = (t**2 - s**3 - a*s) % m
        
        if isSmooth(a, b, m):
            haveCurve = True
            P = Point(s, t, m)
            curve = EC(a, b, m)
            break

    if not haveCurve:
        raise Exception('No EC found, try increasing effort parameter.')

    return P, curve

def lenstraFactorial(N, bound=500, effortObs=500, effortGen=10**5):
    """Attempt to return a factor of N.

    Parameters
    ----------
    N         : int : Integer whose factor is to be found.
    bound     : int : [bound!]P for some point P on an EC is computed.
    effortObs : int : Maximum times EC mult is attempted.
    effortGen : int : Number of attempts to generate a point and smooth curve.

    Returns
    -------
    res : int : A factor of N (res=N, if factor can't be found).

    Example(s)
    ----------
    >>> N = (3209622181 * 6727426213 * 2810645183)
    >>> lenstraFactorial(N)
    >>> 2810645183 (0.5822016000020085 seconds)

    """
    num = factorial(bound)
    searching = True

    #Search for local obstructions to elliptic curve multiplication
    for _ in range(effortObs):
        if not searching:
            break
        P, curve = generateCurvePoint(N, effortGen)
        searching, res = curve.mult(P, num)

    return N if searching else res
