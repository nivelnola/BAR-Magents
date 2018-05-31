import math

# Coordinate grid management
# Set maximum coordinates
def setMaxCoords(xnew, ynew, znew):
    x_max = xnew
    y_max = ynew
    z_max = znew

class Coordinate:
    'Discrete coordinate points in 3-space. Effectively a 3-vector'

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def getCoord(self):
        return [self.x, self.y, self.z]

    def distance(self, end):
        return [end.x - self.x, end.y - self.y, end.z - self.z]

# Current-carrying segment
class CurrSegment:
    'Wire segment of a certain length, with current moving from beg-->end'

    def __init__(self, beg_new, end_new, I_new):
        self.beg = beg_new
        self.end = end_new

        self.curr = I_new

        vect = self.beg.distance(self.end)
        self.length = math.sqrt(vect[0]**2 + vect[1]**2 + vect[2]**2)
        self.mid = Coordinate(self.beg.x + vect[0]/2, self.beg.y + vect[1]/2, self.beg.z + vect[2]/2)

# Current-carrying dL
#class CurrSegDelta(CurrSegment):
#    'dL wire segment with current moving through'

#    def __init__(self, 
