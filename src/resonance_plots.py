'''
Created on June 24, 2015

@author: Maggie Mallernee

Script to plot orbital elements over time in order to look for resonances
1. read in .aei file for each star
2. plot orbital elements vs. time
    - one plot per element
    - different color line per star
'''

import matplotlib.pyplot as plt
import random
import numpy as np
#from scipy import constants as c
import width_calculator as wid

#go ahead and create a star object? - could write as subclass of orbit

def get_data(aei_f):
    
    with open(aei_f) as f:   #assume only read (r) mode
        #create list of each non-empty line in a given file
        lines = []
        for line in f:
            if not line.isspace():
                lines.append(line)
    
    over_time = []
    
    for n in range(2, len(lines)):
        x = lines[n].split()    #a list of t, a, e, i, omega, w, f, mass, x, y, z, u, v, w, M from mercury
        for y in range(0, len(x)):  #might not work with the scientific notation?
            if x[y] == '***************':
                x[y] = '400'
            x[y] = float(x[y])
        over_time.append(x)     #a list of lists of elements, one per time interval
        
    t = []
    a = []
    e = []
    i = []
    omega = []
    w = []
    f = []
    M = []
    
    for n in range(len(over_time)):     #maybe make each element a set of key value pairs with time?
        t.append(over_time[n][0])
        a.append(over_time[n][1])
        e.append(over_time[n][2])
        i.append(over_time[n][3])
        omega.append(over_time[n][4])
        w.append(over_time[n][5])
        f.append(over_time[n][6])   
        M.append(over_time[n][14])
        
    return [t, a, e, i, omega, w, f, M]
'''    
#Function to get element for all of the stars, takes a list of file names, returns 3D array
def aggregate(list_files):
    stars = []
    
    for f in list_files:
        stars.append(get_data(f))
    
    return stars
'''
import resonance_ID_p1 as rcheck

def plot_res_ratio(file_name, big_file_name, p, q, m_big, m_cen):
    
    data = get_data(file_name)
    big_data = get_data(big_file_name)  #A OF BIG BODY??
    
    t = data[0]
    a = data[1]
    e = data[2]
    #i = data[3]
   # omega = data[4]
   # w = data[5]
   # f = data[6]
   
    a_big = big_data[1]
    
    res = []
    min21 = []
    max21 = []
    for x in range(len(a)):
        res.append( pow((abs(a[x]) / a_big[x]),(3./2.)) )
        min_max = rcheck.get_res_ratio_min_max(p, q, a[x], e[x], a_big[x], m_big, m_cen)   
        min21.append(min_max[0])
        max21.append(min_max[1])
    plt.figure("Elements over time (years): " + file_name)
    
    plt.scatter(t, res)
    plt.plot(t, min21, 'r-')
    plt.plot(t, max21, 'r-')
    plt.ylabel("Resonance p:q ratio")
    
    plt.show()

def plot_a_widths(file_name, big_file_name, p, q, m_big, m_cen):
    data = get_data(file_name)
    big_data = get_data(big_file_name)  #A OF BIG BODY??
    
    t = data[0]
    a = data[1]
    e = data[2]
    #i = data[3]
   # omega = data[4]
   # w = data[5]
   # f = data[6]
   
    a_big = big_data[1]
    
    
    wid_min_max = []
    for x in range(len(a)):
       # a = (p / (p+q)) ** (2./3.) * a_big
        wid_min_max.append(wid.get_width(m_cen, m_big, a_big[x], p, q, e[x]))

def plot_elements(file_name, big_file_name, m_cen):
    
    data = get_data(file_name)
    big_data = get_data(big_file_name)  #A OF BIG BODY??
    
    t = data[0]
    a = data[1]
    e = data[2]
    i = data[3]
   # omega = data[4]
   # w = data[5]
   # f = data[6]
   
    a_big = big_data[1]
    
    res = []
    for x in range(len(a)):
        res.append( (abs(a[x]) / a_big[x])**(3./2.) )
        
    plt.figure("Elements over time (years): " + file_name)
    plt.title("Elements over time for %s" % (file_name))
    
    
    ax1 = plt.subplot(411)
    plt.scatter(t, res)
    plt.ylabel("Resonance p:p+q ratio", fontsize = 8)
    plt.setp(ax1.get_xticklabels(), visible = False)
    plt.setp(ax1.get_yticklabels(), fontsize = 6)
   # ax1.set_ylim(ymin = 0)
    
    ax2 = plt.subplot(412, sharex = ax1)
    plt.scatter(t, a)
    plt.ylabel("Semi-Major Axis (AU)", fontsize = 8)       #ideally might also include TD radius
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), fontsize = 6)
    
    
    ax3 = plt.subplot(413, sharex = ax1)
    plt.scatter(t, e)
    plt.ylabel("Eccentricity", fontsize = 8)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), fontsize = 6)
   # ax3.set_ylim(ymin = 0)
    
    #to show upper limit of e (=1)
    e_max = []
    for y in range(len(t)):
        e_max.append(1.0)
    
    plt.plot(t, e_max, '-')
    
    ax4 = plt.subplot(414, sharex = ax1)
    plt.scatter(t, i)
    plt.ylabel("Inclination", fontsize = 8)
    plt.xlabel("Time (Years)")
    plt.setp(ax4.get_xticklabels(), fontsize = 10)
    plt.setp(ax4.get_yticklabels(), fontsize = 6)
   # ax4.set_ylim(ymin = 0)
    #C_pq = []
        
    #for t in range(len(data[0])):
        #C_pq.append( (m_cen * a[t])**0.5 * ((p+q) - p * (1 - e[t]**2)**0.5 * math.cos(math.radians(i[t]))) ) DEPENDS ON P AND Q
    
    #plt.subplot(515)
    #plt.scatter(t, t)
    #plt.ylabel("Adiabatic invariant in progress")
    
    
    
    #plt.subplot(234)
    #plt.scatter(t, omega)
    #plt.title("Ascending Node over Time")
    #
    #plt.subplot(235)
    #plt.scatter(t, w)
    #plt.title("Argument of Pericentre over Time")
    #
    #plt.subplot(236)
    #plt.scatter(t, f)
    #plt.title("True Anomaly over Time")
    plt.xlim(xmin = 0, xmax = max(t))
    plt.savefig('%s Elements' % file_name, format = 'eps', dpi = 1000)
    plt.show()
    
def plot_all(num_stars):
    
    file_name = 'S'
    
    temp_data = []
    
    t = []
    a = []
    e = []
    i = []
    omega = []
    w = []
    f = []
    
    for x in range(1, num_stars+1):
        file_name = 'S' + str(x) + '.aei'
        temp_data = get_data(file_name)
        t.append(temp_data[0])
        a.append(temp_data[1])
        e.append(temp_data[2])
        i.append(temp_data[3])
        omega.append(temp_data[4])
        w.append(temp_data[5])
        f.append(temp_data[6])
        
    plt.figure(1)
    
    plt.subplot(231)
    plt.scatter(t, a)
    plt.title("Semi-Major Axis over Time")
    
    plt.subplot(232)
    plt.scatter(t, e)
    plt.title("Eccentricity over Time")
    
    plt.subplot(233)
    plt.scatter(t, i)
    plt.title("Inclination over Time")
    
    plt.subplot(234)
    plt.scatter(t, omega)
    plt.title("Ascending Node over Time")
    
    plt.subplot(235)
    plt.scatter(t, w)
    plt.title("Argument of Pericentre over Time")
    
    plt.subplot(236)
    plt.scatter(t, f)
    plt.title("True Anomaly over Time")
    
    plt.show()

def plot_by_star(num_stars):
    
    file_name = ''
    data = []

    plt.figure(1)
    
    for x in range(1, num_stars+1):
        file_name = 'S' + str(x) + '.aei'
        data = get_data(file_name)
        color = [random.random(), random.random(), random.random()]
        plt.subplot(231)
        plt.plot(data[0], data[1], ':', c = color)
            
        plt.subplot(232)
        plt.plot(data[0], data[2], ':', c = color)
            
        plt.subplot(233)
        plt.plot(data[0], data[3], ':', c = color)
        
        plt.subplot(234)
        plt.plot(data[0], data[4], ':', c = color)
            
        plt.subplot(235)
        plt.plot(data[0], data[5], ':', c = color)
    
        plt.subplot(236)
        plt.plot(data[0], data[6], ':', c = color)
        
        
    plt.subplot(231)
    plt.title("Semi-Major Axis over Time")
    plt.xlabel("t (years)")
    plt.ylabel("a (AU)")
    plt.subplot(232)
    plt.title("Eccentricity over Time")
    plt.xlabel("t (years)")
    plt.ylabel("e")
    plt.subplot(233)
    plt.title("Inclination over Time")
    plt.xlabel("t (years)")
    plt.ylabel("i (degrees)")
    plt.subplot(234)
    plt.title("Ascending Node over Time")
    plt.xlabel("t (years)")
    plt.ylabel("Omega (degrees)")       #put in actual symbol with unicode
    plt.subplot(235)
    plt.title("Argument of Pericentre over Time")
    plt.xlabel("t (years)")
    plt.ylabel("w (degrees)")           #again, unicode
    plt.subplot(236)
    plt.title("True Anomaly over Time")
    plt.xlabel("t (years)")
    plt.ylabel("f (degrees)")
    
    plt.show()

    
    
#same as plot_by_star except creates a different figure (window) for each orbital element
#also includes (width) boundaries of 2:1 resonances
#@param time in years
def plot_big(num_stars, time):
    
    file_name = ''
    data = []
    
    plt.figure(1)
    plt.title("Semi-Major Axis over Time")
    plt.xlabel("t (years)")
    plt.ylabel("a (AU)")
    
    #put in resonance width boundaries
    '''t = np.arange(0, time, 1.0)
    #for 2:1 resonances
    a_min_max = wid.get_width_list(1.0, .001, 5.2027536, 1.0, 1.0, 0.1, time)
    plt.plot(t, a_min_max[0], 'r-', label='2:1')
    plt.plot(t, a_min_max[1], 'r-')
    #for 3:2 resonances
    a_min_max = wid.get_width_list(1.0, .001, 5.2027536, 2.0, 1.0, 0.1, time)
    plt.plot(t, a_min_max[0], 'b-', label='3:2')
    plt.plot(t, a_min_max[1], 'b-')
    #for 3:1 resonances
    a_min_max = wid.get_width_list(1.0, .001, 5.2027536, 1.0, 2.0, 0.1, time)
    plt.plot(t, a_min_max[0], 'g-', label='3:1')
    plt.plot(t, a_min_max[1], 'g-')
    
    plt.legend()
    '''
    
    plt.figure(2)
    plt.title("Eccentricity over Time")
    plt.xlabel("t (years)")
    plt.ylabel("e")
    plt.ylim(0.0, 1.0)
    
    plt.figure(3)
    plt.title("Inclination over Time")
    plt.xlabel("t (years)")
    plt.ylabel("i (degrees)")
    
    plt.figure(4)
    plt.title("Ascending Node over Time")
    plt.xlabel("t (years)")
    plt.ylabel("Omega (degrees)")       #put in actual symbol with unicode
    
    plt.figure(5)
    plt.title("Longitude of Pericentre over Time")
    plt.xlabel("t (years)")
    plt.ylabel("w (degrees)")           #again, unicode
    
    plt.figure(6)
    plt.title("True Anomaly over Time")
    plt.xlabel("t (years)")
    plt.ylabel("f (degrees)")
    
    plt.figure(7)
    plt.title("Mean Anomaly over Time")
    plt.xlabel("t (years)")
    plt.ylabel("M (degrees)")
    
    #astrofest_set = [1, 2, 7]   #sample set of 2:1 mJ = .001 stars with 1 librating 7 circulating and 2 somewhere in between
    for x in range(1, num_stars + 1):
    #for x in astrofest_set:
        file_name = 'S' + str(x) + '.aei'
        data = get_data(file_name)
        color = [random.random(), random.random(), random.random()]
        '''colors = ['magenta', 'blue', 'green']   #
        if x == 1:
            c = colors[0]
        elif x == 2:
            c = colors[1]
        else:
            c = colors[2] '''  #
        plt.figure(1)
        plt.plot(data[0], data[1], '-o', c = color)  #
        
        plt.figure(2)
        plt.plot(data[0], data[2], '-o', c = color)
            
        plt.figure(3)
        plt.plot(data[0], data[3], '-o', c = color)
        
        plt.figure(4)
        plt.plot(data[0], data[4], '-o', c = color)
            
        plt.figure(5)
        plt.plot(data[0], data[5], '-o', c = color)
    
        plt.figure(6)
        plt.plot(data[0], data[6], '-o', c = color)
        
        plt.figure(7)
        plt.plot(data[0], data[7], '-o', c = color)
        
        
    plt.figure(1)
    plt.show()
    plt.figure(2)
    plt.show()
    plt.figure(3)
    plt.show()
    plt.figure(4)
    plt.show()
    plt.figure(5)
    plt.show()
    plt.figure(6)
    plt.show()
    plt.figure(7)
    plt.show()
        
        
def plotI(num_stars):
    
    file_name = ''
    data = []
        
    plt.figure(1)
    plt.title("Inclination over Time")
    plt.xlabel("t (years)")
    plt.ylabel("i (degrees)")   
    
    for x in range(1, num_stars+1):
        file_name = 'S' + str(x) + '.aei'
        data = get_data(file_name)
        color = [random.random(), random.random(), random.random()]
        
        plt.figure(1)
        plt.plot(data[0], data[3], '-o', c = color)
        
    plt.figure(1)
    plt.show()

#TABLE 8.1 ROW 1
#@param: big_body a file name for the outer resonant body, here Jupiter    
def plot_res_arg1(num_stars, big_body):
    
    file_name = ''
    data = []
    
    big_data = get_data(big_body)
    #print big_data
    w_big = []
    M_big = []
    #print len(big_data[0])
    for x in range(len(big_data[0])):
        w_big.append(big_data[5][x])
        M_big.append(big_data[7][x])
    lambda_big = []
    for x in range(len(w_big)):
        lambda_big.append(w_big[x]+M_big[x])
    
    
    
    for x in range(1, num_stars+1):
        file_name = 'S' + str(x) + '.aei'
        data = get_data(file_name)
        lambda_little = []
        phi = []
        phi_mod = []
        w = data[5]
        M = data[7]
        color = [random.random(), random.random(), random.random()]
        for y in range(len(data[0])):
            lambda_little.append(w[y]+M[y])
            #print "%s %s %s" % (lambda_big[y], lambda_little[y], w[y])
            phi.append(2*lambda_big[y] - lambda_little[y] - w[y])     #from Table 8.1 ROW 1
            phi_mod.append(phi[y] % 360)                              #phi mod 360 degrees
            #phi_mod.append(2*(lambda_big[y]%360) - (lambda_little[y]%360) - (w[y]%360)) 
        plt.figure(x)
        plt.subplot(211)
        plt.title("Resonant Argument over Time")
        #plt.xlabel("t (years)")
        plt.ylabel("phi (degrees)")
        plt.plot(data[0], phi, '-o', c = color)
        
        plt.subplot(212)
        #plt.title("Resonant Argument over Time [0, 360]")
        plt.xlabel("t (years)")
        plt.ylabel("phi mod 360 (degrees)")
        #plt.ylim(0., 360.)
        plt.plot(data[0], phi_mod, '-o', c = color)
        
        plt.show()

    #for x in range(1, num_stars+1):        
    #    plt.figure(x)    
    #    plt.show()

#TABLE 8.1 ROW 2
#@param: big_body a file name for the outer resonant body, here Jupiter    
def plot_res_arg2(num_stars, big_body):
    
    file_name = ''
    data = []
    
    big_data = get_data(big_body)
    #print big_data
    w_big = []
    M_big = []
    #print len(big_data[0])
    for x in range(len(big_data[0])):
        w_big.append(big_data[5][x])
        M_big.append(big_data[7][x])
    lambda_big = []
    for x in range(len(w_big)):
        lambda_big.append(w_big[x]+M_big[x])
    
    
    
    for x in range(1, num_stars+1):
        file_name = 'S' + str(x) + '.aei'
        data = get_data(file_name)
        lambda_little = []
        phi = []
        phi_mod = []
        w = data[5]
        M = data[7]
        color = [random.random(), random.random(), random.random()]
        for y in range(len(data[0])):
            lambda_little.append(w[y]+M[y])
            #print "%s %s %s" % (lambda_big[y], lambda_little[y], w[y])
            phi.append(2*lambda_big[y] - lambda_little[y] - w_big[y])     #from Table 8.1 ROW 2
            phi_mod.append(phi[y] % 360)                              #phi mod 360 degrees
        plt.figure(x)
        plt.subplot(211)
        plt.title("Resonant Argument over Time")
        #plt.xlabel("t (years)")
        plt.ylabel("phi (degrees)")
        plt.plot(data[0], phi, '-o', c = color)
        
        plt.subplot(212)
        #plt.title("Resonant Argument over Time [0, 360]")
        plt.xlabel("t (years)")
        plt.ylabel("phi mod 360 (degrees)")
        plt.ylim(0., 360.)
        plt.plot(data[0], phi_mod, '-o', c = color)
        
        plt.show()
#MAIN OPERATION      
#plot_all(6)
#plot_by_star(25)
#    plot_big(25,3000)
#plot_res_arg1(25, 'Jupiter.aei')
#plot_res_arg2(25, 'Jupiter.aei')
#get_width21()
#plotI(6)
#plot_elements("S85.aei")