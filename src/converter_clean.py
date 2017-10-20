'''
Created on Jun 1, 2015

@author: maggiemallernee
'''

import numpy as np
from scipy import constants
from math import acos, sin, cos, degrees, radians

    
class O_Orbit(object):
    
    '''
    Initialize orbit object with
    @param m1: mass of the central body (solar masses)
    @param m2: mass of the orbiting body (solar masses)
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
        return (constants.G * constants.day**2 * 1.98855 * 10**30 / constants.au**3) * (self.m1 + self.m2)
        
    def print_out_r(self):      #print out angles in radians
        print "Central mass: " + str(self.m1)
        print "Orbiting mass: " + str(self.m2)
        print "a: " + str(self.a)
        print "e: " + str(self.e)
        print "i: " + str(self.i)
        print "omega: " + str(self.omega)
        print "w: " + str(self.w)
        print "f: " + str(self.f)
        
    def print_out_d(self):      #print out angles in degrees
        print "Central mass: " + str(self.m1)
        print "Orbiting mass: " + str(self.m2)
        print "a: " + str(self.a)
        print "e: " + str(self.e)
        print "i: " + str(degrees(self.i))
        print "omega: " + str(degrees(self.omega))
        print "w: " + str(degrees(self.w))
        print "f: " + str(degrees(self.f))
        
    def element_val_list_r(self):   #with angles in radians
        return [self.m1, self.m2, self.a, self.e, self.i, self.omega, self.w, self.f]
        
    def element_val_list_d(self):     #with angles in degrees
        return [self.m1, self.m2, self.a, self.e, degrees(self.i), degrees(self.omega), degrees(self.w), degrees(self.f)]    
        
    def element_list(self):
        return ['m1', 'm2', 'a', 'e', 'i', 'omega', 'w', 'f']
        

class C_Orbit(object):
    
    '''
    Initialize orbit object with
    @param m1: mass of the central body (solar masses)
    @param m2: mass of the orbiting body (solar masses)
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
        return (constants.G * constants.day**2 * 1.98855 * 10**30 / constants.au**3) * (self.m1 + self.m2)
    
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
            
    #get R dot
    def get_R_dot(self):
        sign = (self.X * self.U + self.Y * self.V + self.Z * self.W) / abs(self.X * self.U + self.Y * self.V + self.Z * self.W)
        return sign * (self.get_V_sqrd() -  (self.get_h()**2 / self.get_R()**2)) ** 0.5
        
    def element_val_list(self):
        return [self.m1, self.m2, self.X, self.Y, self.Z, self.U, self.V, self.W]
        
    def element_list(self):
        return ['m1','m2','X', 'Y', 'Z', 'U', 'V', 'W']
    
'''
function to convert an orbit from cartesian coordinates to orbital elements
    
@param orbit: a Cartesian_Orbit object
@print: the orbital elements
@return: a new Orbital_Orbit object 
'''
def convert_from_cartesian(orbit):
    #define cartesian elements from input orbit
    h = orbit.get_h()
    H = orbit.get_H()
    R = orbit.get_R()
    V_sqrd = orbit.get_V_sqrd()
    R_dot = orbit.get_R_dot()
    mu = orbit.get_mu()
    X = orbit.X
    Z = orbit.Z
    
    #get the value of a
    term1 = 2 / R
    term2 = V_sqrd / mu
    a = (term1 - term2) ** -1
    
    #get the value of e
    term = h ** 2 / (mu * a)
    e = (1.0 - term) ** .5
    
    #get the value of i
    i = acos(H[2] / h)
    
    #get the value of omega
    hz_pos = True if H[2] > 0 else False
    cos_omega = H[1] / (h * sin(i))
    if hz_pos:
        cos_omega = cos_omega * -1
    omega = acos(cos_omega)
    
    #get the value of f
    sin_f = (a*(1 - e**2)*R_dot) / (h*e)
    cos_f = (1/e) * ((a * (1 - e**2) / R) - 1)
    
    f = np.arctan2(sin_f, cos_f)
    
    #get the value of w from the expression for w+f 
    sin_wf = Z / (R * sin(i)) 
    cos_wf = (1/cos(omega)) * ((X/R) + sin(omega)*sin_wf*cos(i))
    wf = np.arctan2(sin_wf, cos_wf)    
    w = wf - f
    
    #ensure that positive angles returned
    if omega < 0:
        omega = omega + 2*constants.pi        
    if f < 0:
        f = f + 2*constants.pi
    if w < 0:
        w = w + 2*constants.pi
        
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
    Y = r * (sin(omega)*cos(w+f) + cos(omega)*sin(w+f)*cos(i)) 
    Z = r * sin(w+f) * sin(i)
    
    #get U, V, W
    num = (mu / p) ** 0.5 
    
    U = -1 * num * (cos(omega)*(e*sin(w) + sin(w+f)) + cos(i)*sin(omega)*(e*cos(w) + cos(w+f)))
    V = -1 * num * (sin(omega)*(e*sin(w) + sin(w+f)) - cos(i)*cos(omega)*(e*cos(w) + cos(w+f)))
    W = num * (e*cos(w)*sin(i) + cos(w+f)*sin(i))

  
    return C_Orbit(orbit.m1, orbit.m2, X, Y, Z, U, V, W)
    
'''
Function to convert cartesian positions and velocities from the heliocentric reference frame to barycentric

@param: Orbital orbit object
@return: New orbital object for m2
'''
'''def convert_to_barycentric(orbit):
    m1 = orbit.m1
    m2 = orbit.m2
    a = orbit.a
    e = orbit.e
    i = orbit.i
    omega = orbit.omega
    w = orbit.w
    f = orbit.f
    
    #get center of mass --- what's the orbiting mass???
        
    #define new orbit for m1 - e, i, omega remain the same
    #a1 = (m2 / (m1 + m2)) * a
    #w1 = w + constants.pi
    #f1 = f + constants.pi

    #define changes for the orbit m2
    a2 = (m1 / (m1 + m2)) * a
    
        
    #return O_Orbit(m2, m1, a1, e, i, omega, w1, f1) #return m1 orbit
    return O_Orbit(m1, m2, a2, e, i, omega, w, f)   #return m2 orbit
'''

#MAIN TEST OPERATION BELOW

#orbit_J = C_Orbit(1.0, .001, 4.84143144246472090, -1.16032004402742839, -.103622044471123109, .00166007664274403694, .00769901118419740425,-.0000690460016972063023)
#orbital_J = convert_from_cartesian(orbit_J)
#orbital_J.print_out_d() 

#full circle test
'''
orbit_in = C_Orbit(10000000, 100000, 4.84143144246472090, -1.16032004402742839, -0.103622044471123109, 0.00166007664274403694, 0.00769901118419740425, -0.0000690460016972063023)
orbit_mid = convert_from_cartesian(orbit_in)
orbit_out = convert_from_orbital(orbit_mid)
'''
#orbit_in = O_Orbit(10000000, 100000, 2000, 0.999999896694, radians(1.30465899538), radians(100.47168534), radians(66.05626), radians(180.000000134))
#orbit_mid = convert_from_orbital(orbit_in)
#orbit_out = convert_from_cartesian(orbit_mid)

#print "Input:"
#orbit_in.print_out_d()
#print "Intermediate:"
#orbit_mid.print_out()
#print "Output:"
#orbit_out.print_out_d()


'''
#test against mercury
orbit_o = O_Orbit(2, .0009547919, 3.36860, .478659, radians(1.3047), radians(100.4748), radians(64.6367), radians(197.8617))
orbit_c = C_Orbit(2, .0009547919, 4.762297, 0.247190, -0.107677, -0.002640, 0.008095, .000026)

o_to_c = convert_from_orbital(orbit_o)
print "From orbital to cartesian:"
o_to_c.print_out()

c_to_o = convert_from_cartesian(orbit_c)
print "\nFrom cartesian to orbital:"
c_to_o.print_out()
'''
'''
#test Barycentric converter
helio_orbit_o = O_Orbit(2.00000, 1.0000000, 3.146817, 0.6522820, radians(1.30468025), radians(100.4746748), radians(65.3623465), radians(197.5708520))
print "\nThe original orbit:\n"
helio_orbit_o.print_out()
convert_from_orbital(helio_orbit_o).print_out()
print "\nWith the conversion to barycentric:\n"
bary_orbit_o = convert_to_barycentric(helio_orbit_o)
bary_orbit_o.print_out()
convert_from_orbital(bary_orbit_o).print_out()
'''
    

    
