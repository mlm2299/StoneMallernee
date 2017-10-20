'''
Created on June 11, 2015

@author: Maggie Mallernee

    1. file output project I:
        - create a function that takes the parameters of an orbiting body (name, mass, cartesian elements)
        - output a file with the same format as big.in that is readable by mercury
    2. file output project II:
        - modified version of I that converts to the frame of reference used by mercury
'''

import converter_clean as c
from scipy import constants

#G = constants.G * constants.day**2 * 1.98855 * 10**30 / constants.au**3
'''
Function to format big.in files for mercury (in the lab frame)

@param: name
@param: mass
@param: cartesian orbital elements      #too bad not just an orbit object? --- this does not require a value for m1

@return: big.in file to be input to mercury
'''
def format(name, mass, density, x, y, z, u, v, w):
    file = open("big.in", 'w')
    file.write(')O+_06 Big-body initial data  (WARNING: Do not delete this line!!) \n')
    file.write(") Lines beginning with `)' are ignored. \n")
    file.write(")--------------------------------------------------------------------- \n")
    file.write(" style (Cartesian, Asteroidal, Cometary) = Cartesian \n")
    file.write(" epoch (in days) = 2451000.5 \n")
    file.write(")--------------------------------------------------------------------- \n")
    file.write(' %s    m=%s r=3.d0 d=%s \n' % (name, mass, density))  #previously, d=1.33
    file.write("  %s %s %s \n" % (x, y, z))
    file.write("  %s %s %s \n" % (u, v, w))
    file.write("  0. 0. 0.\n")
    file.close()
    
    return file
    
#format("JUPITER", 0.00100000000000000000, 5.2027536, 0.00, 0.00, 
 #   0.00, 0.0075454, 0.00)

'''
Function to format big.in files in the correct reference frame for mercury
Initial conditions to set up a binary:
'''  
def frame_format(name, m_cen, m_big, a_big, r_start):
    G = constants.G * constants.day**2. * 1.98855 * 10.**30. / constants.au**3.
    mu = G * (m_cen + m_big)
    v = (mu / a_big) ** (0.5)
    
    #to get proper density
    r_star = 0.00464913 #AU
    r_t = r_star * (m_big) ** (1./3.)
    density = (3. / (4. * constants.pi)) * (m_big / (r_t ** 3.)) #might want to format/cut down
    #convert from solar masses / cubic AU to g / cm^3
    density = density * ((1.989 * 10.**30.) / constants.au**3) * (1000. / 1.) * (1. / 100.)**3. # (solar mass / cubic AU) (cubic AU / solar mass) (g / kg) (cm / m)^3

   # r_start = 30.0  #for run E -- choose start at peri
    v = (mu * ((2./r_start) - (1./a_big)))**0.5
    print v
    
    return format(name, m_big, density, r_start, 0.0, 0.0, 0.0, v, 0.0)
    
 
    
#frame_format("JUPITER", 10000000, 100000, 2000)
   
'''
 Function to format big.in files for mercury in the correct reference frame (mercury frame)
 Assumes input is heliocentric
 
 @param: name
 @param: cartesian orbit object
 
 @return: big.in file to be input to mercury

def frame_format(name, orbit):
     m1 = orbit.m1
     m2 = orbit.m2
     R = orbit.get_R()
     mu = orbit.get_mu()
     X = orbit.X
     Y = orbit.Y
     Z = orbit.Z

     #get new positions, now in barycentric frame
     num = m1 / (m1 + m2)
     Xb = num * X
     Yb = num * Y
     Zb = num * Z
    
     #get new velocities, adjusted to barycentric then mercury frame
    
     num2 = m2 / (m1 + m2)
     R1 = [-1*num2 * X, -1*num2*Y, -1*num2*Z]   #list of position of the more massive object
     R2 = [Xb, Yb, Zb]
    
     omegak = (mu / R**3) ** 0.5
    
     Ub = omegak * (R2[0] - R1[0])
     Vb = omegak * (R2[1] - R1[1])
     Wb = omegak * (R2[2] - R1[2])
    
     #create the file
     format(name, m2, Xb, Yb, Zb, Ub, Vb, Wb)
     

orbit = c.C_Orbit(2.000000000, 1.0000000000, 4.84143144246472090, -1.16032004402742839, -1.03622044471123109, 
    1.66007664274403694, 7.69901118419740425, -6.90460016972063023)
     
    
#frame_format("JUPITER", orbit)
'''
     
    
     
     
 
 
 
 
 
 
 
 
