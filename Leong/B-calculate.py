'''
Calculates magnetic field lines for a set of wires 
using Biot-Savart.

Created on Dec 4, 2011

@author: Eric Leong
with guidance from Professor Wolf
'''

from numpy.ma.core import sqrt
import inspect
import math
import pickle
import time

# constants
I = 1
mu0 = 1e-7
k = mu0 / (4 * math.pi)

def bfield_segment(x, y, z, lx0, ly0, lz0, lx1, ly1, lz1):
    '''Calculates magnetic field using Biot-Savart'''
    # Midpoint of the line segment
    mx = (lx0 + lx1) / 2
    my = (ly0 + ly1) / 2
    mz = (lz0 + lz1) / 2
    
    # Distance from the midpoint of the line segment to the point
    r = sqrt((x - mx) ** 2 + (y - my) ** 2 + (z - mz) ** 2)
    
    # Check for divide by zero!
    if r == 0:
        return 0, 0, 0
    
    # Multiply by constants.
    c = k * I * r ** (-3)
    
    # Algebraic cross product
    dbx = c * ((ly1 - ly0) * (z - mz) - (y - my) * (lz1 - lz0))
    dby = -c * ((lx1 - lx0) * (z - mz) - (x - mx) * (lz1 - lz0))
    dbz = c * ((lx1 - lx0) * (y - my) - (x - mx) * (ly1 - ly0))
    
    return dbx, dby, dbz

def bfield(x, y, z, w):
    '''Calculates the b field at the given point from the wires.'''
    
    # Initial b field
    bx, by, bz = 0, 0, 0
    
    for wire in w:
        # Iterate using an index in order to access the next element
        for p in range(len(wire) - 1):
            dbx, dby, dbz = bfield_segment(x, y, z,
                                wire[p][0], wire[p][1], wire[p][2],
                                wire[p + 1][0], wire[p + 1][1], wire[p + 1][2])
            
            bx += dbx
            by += dby
            bz += dbz
                    
    return bx, by, bz

def bfield_mag(x, y, z, w):
    '''Magntitude of b field at the given point from the wires.'''
    bx, by, bz = bfield(x, y, z, w)
    return sqrt(bx ** 2 + by ** 2 + bz ** 2)

def fieldpoint(i, t, h, d, w):
    '''Calculates a point in the field line.
    Note that there is no time dependence.
    '''
    # Determine the next position.
    o = map(lambda i, d: [i[0] + d[0] * h, i[1] + d[1] * h], i, d)
    
    # Determine the magnetic field at the new position.
    bx, by, bz = bfield(o[0][0], o[1][0], o[2][0], w)
    b = sqrt(bx ** 2 + by ** 2 + bz ** 2)
    
    # Avoid divide by zero errors
    if b == 0:
        # It's probably [[0, 0], [0, 0], [0, 0]]
        # But due to floating point error, it might not actually be.
        return [[bx, 0], [by, 0], [bz, 0]]
    
    # Divide by magnetic field strength to arrive at the gradient.
    return [[bx / b, 0], [by / b, 0], [bz / b, 0]]

def integrate(o, t, h, w):
    '''Calculates the next position given initial position o.
    Uses a forth order Runge-Kutta implementation.
    Note that o is this format: [[x, vx], [y, vy], [z, vz]]
    Also note that there is no time dependence, so t and h don't depend on anything.
    '''
    k1 = fieldpoint(o, t, 0, [[0, 0]] * len(o), w)
    k2 = fieldpoint(o, t + h / 2., h / 2., k1, w)
    k3 = fieldpoint(o, t + h / 2., h / 2., k2, w)
    k4 = fieldpoint(o, t + h, h, k3, w)
    
    # Return it back into this form [[x, vx], [y, vy], [z, vz]]
    return map(lambda o, k1, k2, k3, k4: 
               [o[0] + 1 / 6. * h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
                o[1] + 1 / 6. * h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])],
               o, k1, k2, k3, k4)

def fieldline(x, y, z, w, sstart=0, send=100, sstep=.5):
    '''Creates one field line.'''
    # Store magnitude of b
    b = bfield_mag(x, y, z, w)

    # Initial position + strength of field line.
    line = [[x, y, z, b]]
    
    s = sstart
    # Initial position + velocity
    pos = [[x, 0], [y, 0], [z, 0]]
    
    while s <= send + sstep:
        # Determine next position.
        pos = integrate(pos, s, sstep, w)
        
        # Remove velocity from position
        pos = [[pos[0][0], 0], [pos[1][0], 0], [pos[2][0], 0]]

        # Store magnitude of b field
        b = bfield_mag(pos[0][0], pos[1][0], pos[2][0], w)

        # Add point to line
        line.append([pos[0][0], pos[1][0], pos[2][0], b])
        
        # if near the start point (again), don't go on
        dist_sq = (x - pos[0][0]) ** 2 + (y - pos[1][0]) ** 2 + (z - pos[2][0]) ** 2
        if s > sstart + sstep and dist_sq < sstep ** 2:
            # Integrate once more to 'cover' the hole in case the line hasn't been finished
            pos = integrate(pos, s, sstep, w)
            # Remove velocity from position
            pos = [[pos[0][0], 0], [pos[1][0], 0], [pos[2][0], 0]]
            # Store magntitude of b field
            b = bfield_mag(pos[0][0], pos[1][0], pos[2][0], w)
            # Add point to line
            line.append([pos[0][0], pos[1][0], pos[2][0], b])
            break
        
        s += sstep
        
    return line

def split(tstart, tend, tstep, f, offsets=[], nlines=5, valid=True, starts=[], at=0):
    '''Split the wire into segments.
    offsets - the numerical values to offset the the function by.
    nlines - number of field lines for this wire.
    valid - whether or not the wire is active at this time
    starts - exact positions that the magnetic field lines should start at
    '''
    
    # Check if the wire is valid.
    if hasattr(valid, '__call__'):
        if not valid(at):
            return []
    elif not valid:
        return []
    
    points = []

    # Number of segments between field lines
    numbetweenflines = math.ceil((tend - tstart) / (nlines * tstep))
    # Number of arguments to f
    numargs = len(inspect.getargspec(f)[0])

    n = 0
    t = tstart
    while t <= tend + tstep:
        n += 1

        point = None
        
        if t <= tend:
            # Determine correct number of arguments for f
            if numargs == 4:
                point = f(t, 0, at, False)
            elif numargs == 3:
                point = f(t, 0, at)
            elif numargs == 2:
                point = f(t, 0)
            else:
                point = f(t)
            # Add point to wire
            if point:
                points.append(point)
        
        # Check if it's time for a field line
        if (n - 1) % numbetweenflines == 0:
            
            # Store start points for field lines.
            # Iterate through requested offsets,
            for offset in offsets:
                
                s = None
                # Determine correct number of arguments for f
                if numargs == 4:
                    s = f(t, offset, at, True)
                elif numargs == 3:
                    s = f(t, offset, at)
                elif numargs == 2:
                    s = f(t, offset)

                if s:
                    starts.append(s)
            
        t += tstep
        
    return points

def calculate(wires, at=0):
    '''Calculate the wire segments and field line points.'''
    
    # Reset field line starting points.
    starts = []
    
    # Build the wires
    w = [split(*wire, starts=starts, at=at) for wire in wires]
    
    # Create the fieldlines
    fieldlines = [ fieldline(p[0], p[1], p[2], w) for p in starts ]
    
    return w, fieldlines

def animate(wires, atstart, atend, atstep, prepend='', n=0):
    '''Iterate through the animation and calculate wire segments 
    and field line points for different moments in time.'''
    at = atstart
    
    while at < atend:
        
        # Time how long it takes to calculate
        stime = time.time()
        w, l = calculate(wires, at)
        print n, time.time() - stime
        
        # Dump the data into the file
        f = open('%s%04d.txt' % (prepend, n), 'w')
        pickle.dump([w, l], f)
        f.close()
        
        # Increment counters
        at += atstep
        n += 1

if __name__ == "__main__":
    # animation
    n = 0
    atstep = .01 # rate
    atstart = 0 + n * atstep
    atend = 1
    prepend = 'line_'
    
    # Single short wire.
    def a(t, o, at, fl):
        if fl:
            return (t, o, 0)
        else:
            if abs(t) <= 1:
                return (t, o, 0)
            else:
                return None
    
    # wires
    
    # single wire
#    rs = -2
#    re = 2
#    ri = .025 # rate
#    wires = [
#             [rs, re, ri,
#              # parametric variable, offset, animation var, (whether or not offset is being computed)
#              a, [.25, .75], 9, lambda at: at >= 0],
#             ]

    # toroid
    from numpy.ma.core import cos, sin
    
    toroid_r = 1
    r = 1
    rs = 0
    re = 2 * math.pi
    ri = math.pi / 30 # rate
    offset = [-.25, -.5, -.75]
    start_toroid_r = 5
    nl = 2
    wires = [
         [rs, re, ri,
          lambda t, o, at: (cos(0 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), sin(0 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), (r + o) * sin(t) ), offset, nl, lambda at: at >= 0],
         [rs, re, ri,
          lambda t, o, at: (cos(1 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), sin(1 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), (r + o) * sin(t) ), offset, nl, lambda at: at >= 0],
         [rs, re, ri,
          lambda t, o, at: (cos(2 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), sin(2 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), (r + o) * sin(t) ), offset, nl, lambda at: at >= 0],
         [rs, re, ri,
          lambda t, o, at: (cos(3 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), sin(3 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), (r + o) * sin(t) ), offset, nl, lambda at: at >= 0],
         [rs, re, ri,
          lambda t, o, at: (cos(4 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), sin(4 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), (r + o) * sin(t) ), offset, nl, lambda at: at >= 0],
         [rs, re, ri,
          lambda t, o, at: (cos(5 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), sin(5 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), (r + o) * sin(t) ), offset, nl, lambda at: at >= 0],
         [rs, re, ri,
          lambda t, o, at: (cos(6 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), sin(6 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), (r + o) * sin(t) ), offset, nl, lambda at: at >= 0],
         [rs, re, ri,
          lambda t, o, at: (cos(7 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), sin(7 * 2 * math.pi / 8) * (toroid_r + (-start_toroid_r * (at - 1) if at <= 1 else 0) + (r + o) * cos(t)), (r + o) * sin(t) ), offset, nl, lambda at: at >= 0],
         ]

    animate(wires, atstart, atend, atstep, prepend, n)
