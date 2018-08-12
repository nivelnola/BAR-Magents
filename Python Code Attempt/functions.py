import math

def cross(vectA, vectB):
    try:
        assert len(vectA)==3 and len(vectB)==3
        x =   vectA[1]*vectB[2] - vectA[2]*vectB[1]
        y = -(vectA[0]*vectB[2] - vectA[2]*vectB[2])
        z =   vectA[0]*vectB[1] - vectA[1]*vectB[0]
        return [x,y,z]
    except TypeError,e:
        print "Vector contents must be numbers."
    except:
        print "Vectors must be of length 3."

def displacement(start, end):
    try:
        assert len(start)==3 and len(end)==3
        x = end[0] - start[0]
        y = end[1] - start[1]
        z = end[2] - start[2]

        return [x, y, z]
    except TypeError,e:
        print "Point coordinates must be numeric."
    except:
        print "Points must have 3 coordinate values."

def unitVector(start, end):
    try:
        vect = displacement(start, end)
        norm = math.sqrt(vect[0]**2 + vect[1]**2 + vect[2]**2)
        vect[0] /= norm
        vect[1] /= norm
        vect[2] /= norm

        return vect
    except TypeError,e:
        raise
    except:
        raise

        
def midpoint(start, end):
    try:
        assert len(start)==3 and len(end)==3
        x = (end[0] + start[0])/2.0
        y = (end[1] + start[1])/2.0
        z = (end[2] + start[2])/2.0

        return [x, y, z]
    except TypeError,e:
        print "Point coordinates must be numeric."
    except:
        print "Points must have 3 coordinate values."

# def biotSavart(start, end, point, current)
