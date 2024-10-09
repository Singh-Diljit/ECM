"""Implement 'EC' class, for elliptic curves over Z/mZ."""

from Point import Point, identity
from helperFuncs import isCurveInd, resCurveInd, modularSlope, modInv

class EC:
    """Implements elliptic curves over Z/mZ."""

    def __init__(self, a, b, modulus):
        """Initialize an elliptic curve.

        Parameters
        ----------
        a, b    : int  : Coefficients in: (y**2 = x**3 + ax + b).
        modulus : int  : (x, y) is in (Z/mZ)^2 for m = modulus.

        Initializes
        -----------
        self.a       : int : Coefficient of linear term in associated cubic.
        self.b       : int : Constant term in associated cubic.
        self.modulus : int : The modulus.
        
        """
        self.a = a % modulus
        self.b = b % modulus
        self.modulus = modulus

    def tangent(self, P):
        """Return slope of the line tangent to self and a given point.

        Parameters
        ----------
        P : Point : A point on the curve.

        Returns
        -------
        exists : bool : If the slope is well-defined.
        res    : int  : If well-defined the slope else the obstruction.
        
        """
        dx = (3 * (P.x ** 2) + self.a) % self.modulus
        dy = (2 * P.y)
        
        exists, dyInv = modInv(dy, self.modulus)        
        if exists:
            res = (dx * dyInv) % self.modulus
        else:
            res = dyInv

        return exists, res

    def __add_sumExists(self, P, Q, slope):
        """Given the slope of the secant and that the sum exists, return P + Q.

        Parameters
        ----------
        P, Q  : Point : Unique points in (Z/mZ)^2.
        slope : int   : Modular slope the secant line.

        Returns
        -------
        - : Point : The sum of P and Q.
        
        """
        resX = slope**2 - P.x - Q.x
        resY = slope*(P.x - resX) - P.y
        
        return Point(resX, resY, self.modulus)
           
    def double(self, P):
        """Attempt to return [2]P.

        Parameters
        ----------
        P : Point : A point on the curve.

        Returns
        -------
        exists : bool       : If 2[P] is well-defined.
        res    : Point, int : If well-defined [2]P else the obstruction.
        
        """      
        if P == P.inverse:
            exists = True
            res = identity

        else:
            exists, tangentAttempt = self.tangent(P)
            if exists:
                res = self.__add_sumExists(P, P, tangentAttempt)
            else:
                res = tangentAttempt

        return exists, res

    def add(self, P, Q):
        """Attempt to return the sum of two points.

        Parameters
        ----------
        P, Q  : Point : Point on the curve.

        Returns
        -------
        exists : bool       : If P+Q is well-defined.
        res    : Point, int : If well-defined P+Q else the obstruction.
        
        """    
        if P == Q:
            exists, res = self.double(P)

        elif isCurveInd(P, Q):
            exists, res = True, resCurveInd(P, Q)
            
        else:
            exists, slopeAttempt = modularSlope(P, Q)
            if exists:
                res = self.__add_sumExists(P, Q, slopeAttempt)
            else:
                res = slopeAttempt

        return exists, res
                    
    def mult(self, P, k):
        """Attempt to return k[P].

        Parameters
        ----------
        P : Point : Point on the curve.
        k : int   : Multiple of P.  

        Returns
        -------
        exists : bool       : If [k]P is well-defined.
        res    : Point, int : If well-defined [k]P else the obstruction.
        
        """
        if k < 0:
            P, k = P.inverse, -k

        if isCurveInd(P, k=k):
            exists, res = True, resCurveInd(P, k=k) 

        else:
            exists, res = True, identity
            for bit in bin(k)[:1:-1]:
                if not exists:
                    break
                if bit == '1':
                    exists, res = self.add(P, res)
                    if not exists: break
                
                exists, P = self.double(P)
                
        return exists, res
