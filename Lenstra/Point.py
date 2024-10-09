"""Implement 'Point' class, for points on an EC with components in Z/mZ."""

class Point:
    """Implements points (on an EC) of the form: (x, y) for x, y in Z/mZ."""
                  
    def __init__(self, x=0, y=0, modulus=1, ID=False):
        """Initialize a point.

        Parameters
        ----------
        x, y    : int  : x and y coordinates of the point.
        modulus : int  : (x, y) is in (Z/mZ)^2 for m = modulus.
        ID      : bool : If representing the point at infinity.

        Initializes
        -----------
        self.ID      : bool : If point represents point at infinity. 
        self.x       : int  : x-coordinate of point.
        self.y       : int  : y-coordinate of point.
        self.modulus : int  : The modulus.
        
        """
        self.ID = ID
        self.modulus = modulus
            
        if self.ID:
            self.x, self.y = float('inf'), float('inf')
        else:
            self.x, self.y = x % modulus, y % modulus
            
    @property
    def inverse(self):
        """Return the additive inverse wrt the elliptic curve group law.

        Returns
        -------
        P : Point : The point, P, s.t. P + self == ID under the group law.

        """
        if self.ID:
            res = Point(ID=True)
        else:
            newY = self.modulus-self.y if (self.y != 0) else 0
            res = Point(self.x, newY, self.modulus)
        return res
        
identity = Point(ID=True)
