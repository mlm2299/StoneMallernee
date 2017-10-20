'''
Created on June 30, 2015

@author: Maggie Mallernee

To sample stars with the hope of populating some p+q:p (interior) resonances with a larger body
To create a small.in file that works with Mercury
'''

import random
import anomaly_convert as ac
from scipy import constants as c

#@param: a_big the semi-major axis of the perturbing body
def sample_format(num_stars, m_cen, m_big, a_big, p, q, e, i):    #NEED PLACE TO DENOTE P & Q IN THE FILE
    file = open("small.in", 'w')
    file.write(')O+_06 Small-body initial data  (WARNING: Do not delete this line!!) \n')
    file.write(") Lines beginning with `)' are ignored. \n")
    file.write(")--------------------------------------------------------------------- \n")
    file.write(" style (Cartesian, Asteroidal, Cometary) = Ast \n")
    file.write(")--------------------------------------------------------------------- \n")
    
    a = (p / (p+q)) ** (2./3.) * a_big * (m_cen / (m_cen+m_big))**(1./3.)     #without the last term: assumes m1/(m1+m2) = 1
    
    for x in range(1, num_stars+1):
        name = 'S' + str(x)
        omega = random.uniform(0.0, 360.)   #DEGREES
        w = random.uniform(0.0, 360.)
        M = random.uniform(0.0, 360.)
        E = ac.newton(e, 1.0, M)      
        f = ac.get_f(e, E)
        file.write(' %s    ep=2450400.5 \n' % (name))
        file.write(" %s %s %s %s %s %s 0 0 0\n" % (a, e, i, omega, w, f))
    
    
    file.close()
    
    return file
    
#sample_format(25, 1.00, 0.1, 5.2027536, 1.00, 2.00, 0.0000, 0.0000)