'''
Created on June 21, 2015

@author: Maggie Mallernee

Converts mean anomaly (M) to true anomaly (f) via the eccentric anomaly (E)
'''
from math import sin, cos, atan, tan, degrees, isnan
#import numpy as np
#from scipy import constants as c

#Functions to get E from M with Newton's method

def f(e, E, M):
    if(isnan(E - e * sin(E) - M)):
        print "WARNING: F ISNAN"
    return E - e * sin(E) - M
    
def df(e, E):
    return 1 - e * cos(E)
    
def newton(e, x, M):
    counter1 = 0
    counter2 = 0
    done = False 
    #x_init = x
    while counter1 < 101 and not done:
        while abs(f(e, x, M)) > .00000001:
            if(isnan(x)):
                print "WARNING: X ISNAN"
            x = x - (f(e, x, M)/df(e, x))
            counter1 += 1
        if abs(f(e, x, M)) < .00000001:
            done = True
        counter2 += 1
    #print "The approximation of E given e %s, M %s, and initial guess %s is: %s" % (e, M, x_init, x)
    #print str(counter1) + ' ' + str(counter2)
    return x
    
#Functions to get f from E in DEGREES

def get_f(e, E):
    
    #need to get correct sign - see that both in the same quadrant (f and E)
    
    tan_f2 = ((1+e)/(1-e)) ** 0.5 * tan(E/2.)
    f = 2 * atan(tan_f2)
    #print "%s %s" % (f, E)
    #if E > 0.5*c.pi and E < 1.5*c.pi:
     #   f = f + c.pi
    #print f
    
    return degrees(f)
    
#TEST
#newton(.999, 1.0, 0.)
