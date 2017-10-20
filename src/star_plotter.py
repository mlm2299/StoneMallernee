'''
Created on June 18, 2015

@author: Maggie Mallernee

Script to plot sampled star orbits and orbital elements
'''
import numpy as np
import matplotlib.pyplot as plt

#Function that calculates the radius of influence (of a black hole), in parsecs
def get_r_infl(mass):
    return (mass / (10**6 * 1.98855 * 10**30)) ** 0.5

def FA(a, r_infl):
    return a**(5./4.)/r_infl**(5./4.)
    
def get_data(infile):

    with open(infile) as fi:   #assume only read (r) mode
        #create list of each non-empty line in a given file
        lines = []
        for line in fi:
            if not line.isspace():
                lines.append(line)
    
    #put star data in a list - each element should be a list of name(str) + orbital elements(float)
    stars = []
    
    for n in range(1, len(lines)):
        x = lines[n].split()
        for y in range(1, len(x)):    
            x[y] = float(x[y])
        stars.append(x)
        
    a = []
    e = []
    i = []
    omega = []
    w = []
    f = []
    
    for n in range(len(stars)):
        a.append(stars[n][1])
        e.append(stars[n][2])
        i.append(stars[n][3])
        omega.append(stars[n][4])
        w.append(stars[n][5])
        f.append(stars[n][6])    
    
    plt.figure(1)
    
    plt.subplot(231)
    #r_infl = get_r_infl(7.9564 * 10**36)    #hardcode
    #t = np.arange(0.0, r_infl, 0.0001)
    plt.plot(sorted(a), 'ro')
    #plt.plot(t, FA(t, r_infl), 'k')
    plt.title("Semi-Major Axis Distribution")
    
    plt.subplot(232)
    plt.plot(sorted(e), 'g^')
    plt.title("Eccentricity Distribution")
    
    plt.subplot(233)
    plt.plot(sorted(i), 'bs')
    plt.title("Inclination Distribution") 
      
    plt.subplot(234)
    plt.plot(sorted(omega), 'r+')
    plt.title("Ascending Node Distribution")
    
    plt.subplot(235)
    plt.plot(sorted(w), 'g--')
    plt.title("Argument of Pericentre Distribution")
    
    plt.subplot(236)
    plt.plot(sorted(f), 'b^')
    plt.title("True Anomaly Distribution")
    
    plt.show()    
    
get_data('stars.dat')
    
'''
    #separate into lists of elements
    a = np.array([0.0], float)
    e = np.array([0.0], float)
    i = np.array([0.0], float)
    omega = np.array([0.0], float)
    w = np.array([0.0], float)
    M = np.array([0.0], float)
    
    for n in range(len(stars)):
        a = np.append(a, float(stars[n][1]))
        e = np.append(e, float(stars[n][2]))
        i = np.append(i, float(stars[n][3]))
        omega = np.append(omega, float(stars[n][4]))
        w = np.append(w, float(stars[n][5]))
        M = np.append(M, float(stars[n][6]))                    #this could be MUCH more efficient, all in stars list
        
    a = a[1:]
    e = e[1:]
    i = i[1:]
    omega = omega[1:]
    w = w[1:]
    M = M[1:]
'''
    
    
