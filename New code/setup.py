class Coord:
    'Discrete coordinate points in 3-space'
    
    # Determines size of memory management & display size
    x_max = 100
    y_max = 100
    z_max = 100

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def getCoord(self):
        return [self.x, self.y, self.z]
