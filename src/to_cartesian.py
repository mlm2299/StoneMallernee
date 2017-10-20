'''
Created on Jun 1, 2015

@author: maggiemallernee
'''

import numpy as np
from scipy import constants
from math import acos, sin, asin, cos, degrees, radians

    
class O_Orbit(object):
    
    perspective = "orbital elements"
    '''
    Initialize orbit object with
    @param m1: mass of the central body
    @param m2: mass of the orbiting body
    @param a: semi-major axis (in AU)
    @param e: eccentricity
    @param i: inclination (radians)
    @param omega: longitude of ascending node (radians)
    @param w: argument of percentre (radians)
    @param f: true anomaly (radians)
    '''
    def __init__(self, m1, m2, a, e, i, omega, w, f):
        self.m1 = m1
        self.m2 = m2
        self.a = a
        self.e = e
        self.i = i
        self.omega = omega
        self.w = w
        self.f = f
        
    def get_mu(self):
        return (constants.G * constants.year**2 * 2 * 10**30 / constants.au**3) * (self.m1 + self.m2)
        
    def print_out(self):
        print "Central mass: " + str(self.m1)
        print "Orbiting mass: " + str(self.m2)
        print "a: " + str(self.a)
        print "e: " + str(self.e)
        print "i: " + str(self.i)
        print "omega: " + str(self.omega)
        print "w: " + str(self.w)
        print "f: " + str(self.f)
        

class C_Orbit(object):
    
    perspective = "cartesian coordinates"
    '''
    Initialize orbit object with
    @param m1: mass of the central body
    @param m2: mass of the orbiting body
    @param X, Y, Z: cartesian coordinates
    @param U, V, W: cartesian velocities
    '''
    def __init__(self, m1, m2, X, Y, Z, U, V, W):
        self.m1 = m1
        self.m2 = m2
        self.X = X 
        self.Y = Y
        self.Z = Z
        self.U = U 
        self.V = V 
        self.W = W   
    
    def get_mu(self):
        return (constants.G * constants.year**2 * 2 * 10**30 / constants.au**3) * (self.m1 + self.m2)
    
    def print_out(self):
        print "Central mass: " + str(self.m1)
        print "Orbiting mass: " + str(self.m2)
        print "X: " + str(self.X)
        print "Y: " + str(self.Y)
        print "Z: " + str(self.Z)
        print "U: " + str(self.U)
        print "V: " + str(self.V)
        print "W: " + str(self.W)
    
    #get magnitude of R
    def get_R(self):
        return (self.X ** 2 + self.Y ** 2 + self.Z ** 2) ** .5
    
    #get magnitude of V squared
    def get_V_sqrd(self):
        return self.U ** 2 + self.V ** 2 + self.W ** 2

    #get h as a vector in the form of a numpy array object
    def get_H(self):
        h1 = self.Y * self.W - self.Z * self.V
        h2 = self.Z * self.U - self.X * self.W
        h3 = self.X * self.V - self.Y * self.U
        H = np.array([h1, h2, h3], float)
        return H
        
    #get the magnitude of h
    def get_h(self):
        H = self.get_H()
        return (H[0] ** 2 + H[1] ** 2 + H[2] ** 2) ** .5
            
            
    
'''
function to convert an orbit from cartesian coordinates to orbital elements
    
@param orbit: a Cartesian_Orbit object
@print: the orbital elements
@return: a new Orbital_Orbit object 
'''
def convert_from_cartesian(orbit):
    #define cartesian elements from orbit
    h = orbit.get_h()
    H = orbit.get_H()
    R = orbit.get_R()
    V_sqrd = orbit.get_V_sqrd()
    V = V_sqrd ** .5
    mu = orbit.get_mu()
    X = orbit.X
    Z = orbit.Z
    
    #get the value of a
    term1 = 2 / R
    term2 = V_sqrd / mu
    a = (term1 - term2) ** -1
    
    #DEBUG print "h: %s" % (h) 
    
    #get the value of e
    term = h ** 2 / (mu * a)
    e = (1 - term) ** .5
    
    #get the value of i
    i = acos(H[2] / h)   ##radians vs degrees!
    
    #DEBUG
    #print "hx: %s" % (H[0])
    #print "hy: %s" % (H[1])
    #print "hz: %s" % (H[2])
    
    #get the value of omega -- not sure if both sin and cos would matter
    hz_pos = True if H[2] > 0 else False
    cos_omega = H[1] / (h * sin(i))
    sin_omega = -1* (H[0] / (h *sin(i)))
    if hz_pos:
        cos_omega = cos_omega * -1
        sin_omega = sin_omega * -1
    omega = degrees(acos(cos_omega))
    
    
    
    #DEBUG
    print "sin of omega: %s" % (sin_omega)
    print "inverse sin omega: %s" % (asin(sin_omega))
    print "asin omega in degrees: %s" % (degrees(asin(sin_omega)))
    
    print "cos of omega: %s" % (-1*H[1] / (h * sin(i)))
    print "inverse cos omega: %s" % (acos(-1*H[1] / (h * sin(i))))
    print "acos omega in degrees: %s" % (degrees(acos(-1*H[1] / (h * sin(i)))))
    
    #get the value of f
    cos_f = (1/e) * ((a * (1 - e**2) / R) - 1)
    f = degrees(acos(cos_f))
    
    '''#get the value of w
    sin_wf = Z / (R * sin(i)) 
    cos_wf = (1/cos(omega)) * ((X/R) + sin(omega)*sin_wf*cos(i))
    wf = degrees(acos(cos_wf))     ##domain issues!
    w = wf - f'''
    
    '''sin_f = a * (1 - e ** 2) * V / (h * e) 
    f = degrees(asin(sin_f))   ##deal with domain issues!
    '''
    #get the value of w from the expression for sin(w+f)
    sin_wf = Z / (R * sin(i)) 
    wf = degrees(asin(sin_wf))
    w = wf - f

    return O_Orbit(orbit.m1, orbit.m2, a, e, i, omega, w, f)

'''
Function to convert from orbital elements to cartesian coordinates
    
@param: an Orbital_Orbit object
@print: the cartesian position and velocity components
@return: a Cartesian_Orbit object
'''
def convert_from_orbital(orbit):
    #define orbital elements from orbit
    a = orbit.a
    e = orbit.e
    i = orbit.i
    omega = orbit.omega 
    w = orbit.w
    f = orbit.f
    mu = orbit.get_mu() 
    
    #define p
    p = a * (1 - e ** 2)
    
    #define r
    r = p / (1 + e * cos(f))
    
    #get X, Y, Z
    X = r * (cos(omega)*cos(w+f) - sin(omega)*sin(w+f)*cos(i))
    Y = r * (sin(omega)*cos(w+f) + cos(omega)*sin(w+f)*cos(i)) #Nick's paper says -
    Z = r * sin(w+f) * sin(i)
    
    #get U, V, W -- hardcoded version
    num = (mu / p) ** 0.5 
    
    #the likely wrong matrix transform with z not = 0
    #U = num * (-1*sin(f)*cos(omega)*cos(w) - cos(omega)*sin(w)*(e + cos(f)) + sin(omega)*cos(i)*sin(f)*sin(w) - sin(omega)*cos(i)*cos(w)*(e + cos(f)) + sin(omega)*sin(i))
    #V = num * (-1*sin(f)*sin(omega)*cos(w) - sin(omega)*sin(w)*(e + cos(f)) - cos(omega)*cos(i)*sin(f)*sin(w) + cos(omega)*cos(i)*cos(w)*(e + cos(f)) - cos(omega)*sin(i))
    #W = num * (-1*sin(i)*sin(f)*sin(w) + sin(i)*cos(w)*(e + cos(f)))
    
    #from Nick's paper
    U = -1 * num * (cos(omega)*(e*sin(w) + sin(w+f)) + cos(i)*sin(omega)*(e*cos(w)) + cos(w+f))
    V = -1 * num * (sin(omega)*(e*sin(w) + sin(w+f)) - cos(i)*cos(omega)*(e*cos(w)) + cos(w+f))
    W = num * (e*cos(w)*sin(i) + cos(w+f)*sin(i))
    
    #when I did the matrix transforms
    #U = num * (-1*sin(f)*cos(omega)*cos(w) - cos(omega)*sin(w)*(e + cos(f) + sin(omega)*cos(i)*sin(f)*sin(w) - sin(omega)*cos(i)*cos(w)*(e + cos(f))))
    #V = num * (-1*sin(f)*sin(omega)*cos(w) - sin(omega)*sin(w)*(e + cos(f) - cos(omega)*cos(i)*sin(f)*sin(w) + cos(omega)*cos(i)*cos(w)*(e + cos(f))))
    #W = num * (-1*sin(i)*sin(f)*sin(w) + sin(i)*cos(w)*(e + cos(f)))
    
    #get U, V, W -- numpy version
    '''p1 = np.array([[cos(w), -1*sin(w), 0], [sin(w), cos(w), 0], [0, 0, 1]], float)
    p2 = np.array([[1, 0, 0], [0, cos(i), -1*sin(i)], [0, sin(i), cos(i)]], float)
    p3 = np.array([[cos(omega), -1*sin(omega), 0], [sin(omega), cos(omega), 0], [0, 0, 1]], float)
    
    uvw = np.array([[sin(f), e + cos(f), 0], [0, 0, 0], [0, 0, 0]], float)
    
    step1 = p1.dot(uvw)
    return step1
    step2 = p2.dot(step1)
    step3 = p3.dot(step2)
    
    U = num * step3[0]
    V = num * step3[1]
    W = num * step3[2]'''
    
    
    return C_Orbit(orbit.m1, orbit.m2, X, Y, Z, U, V, W)

#MAIN TEST OPERATION BELOW

#full circle test
'''orbit_in = O_Orbit(100, .03, 2.517, .9896, 1.305, 2.09, .49, 2.8)
orbit_mid = convert_from_orbital(orbit_in)
orbit_out = convert_from_cartesian(orbit_mid)

print "Input:"
orbit_in.print_out()
print "Intermediate:"
orbit_mid.print_out()
print "Output:"
orbit_out.print_out()
'''

#test against mercury
orbit_o = O_Orbit(100, .0009548, 2.51728, .989626, radians(1.3047), radians(100.5), radians(66.05), radians(182.8))
orbit_c = C_Orbit(100, .0009548, 4.407, -.829, -.095, -.035, .015, .0007)

o_to_c = convert_from_orbital(orbit_o)
print "From orbital to cartesian:"
o_to_c.print_out()

c_to_o = convert_from_cartesian(orbit_c)
print "\nFrom cartesian to orbital:"
c_to_o.print_out()

#unit issues
#print constants.G
    
'''TODO

5. add a way to specify precision
6. specify units & make sure trig functions are using the correct units 
-- currently all using radians
7. fix expressions for U, V, W
'''    
    
    
    

    