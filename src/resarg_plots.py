# -*- coding: utf-8 -*-
'''
Created July 13, 2015

@author: Maggie Mallernee
Function to plot resonant argument
'''

import matplotlib.pyplot as plt
import random
#import math
#import numpy as np
import resonance_plots as rp
#from scipy import constants as c
#import width_calculator as wid

#TABLE 8.1 ROW 1
#@param: big_body _ a file name for the outer resonant body, here Jupiter 
#@param: row _ denotes which formula to use to calculate phi, from Table 8.1 (Solar System Dynamics)    
def plot_res_arg1(num_stars, big_body, p, q, row):
    
    file_name = ''
    data = []
    
    big_data = rp.get_data(big_body)
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
    
    #astrofest_set = [1, 2, 7]       #sample set of 2:1 mJ = .001 stars with 1 librating 7 circulating and 2 somewhere in between
    #get data from each star, plot phi for each star
    for x in range(1, num_stars+1):
    #for x in astrofest_set:
        file_name = 'S' + str(x) + '.aei'
        data = rp.get_data(file_name)
        lambda_little = []
        phi = []
        phi_mod = []
        w = data[5]
        M = data[7]
        color = [random.random(), random.random(), random.random()]
        '''colors = ['magenta', 'blue', 'green']       #
        if x == 1:
            c = colors[0]
        elif x == 2:
            c = colors[1]
        else:
            c = colors[2]        '''       #
        for y in range(len(data[0])):
            lambda_little.append(w[y]+M[y])
            if row == 1:
                phi.append((p+q)*lambda_big[y] + (1-p-q)*lambda_little[y] - w[y] + 180)     #from Table 8.1 ROW 1, add 180 to get between -180 and 180
            elif row == 2:
                phi.append((p+q)*lambda_big[y] + (1-p-q)*lambda_little[y] - w_big[y] + 180)     #from Table 8.1 ROW 2
            elif row == 3:
                phi.append((p+q)*lambda_big[y] + (2-p-q)*lambda_little[y] - 2*w[y] + 180)
            elif row == 4:
                phi.append((p+q)*lambda_big[y] + (2-p-q)*lambda_little[y] - w_big[y] - w[y] +180)
            phi_mod.append((phi[y] % 360) - 180)                              #phi mod 360 degrees, minus 180 to get correct frame
            #phi_mod.append(2*(lambda_big[y]%360) - (lambda_little[y]%360) - (w[y]%360)) 
        plt.figure(x)
        #plt.subplot(211)
        #plt.title("Resonant Argument over Time")
        #plt.xlabel("t (years)")
        #plt.ylabel("phi (degrees)")
        #plt.plot(data[0], phi, '-o', c = color)
        
        #plt.subplot(212)
        #plt.title("Resonant Argument over Time [0, 360]")
        plt.xlabel("t (years)")
        plt.ylabel("phi mod 360 (degrees)")
        plt.ylim(-180., 180.)
        plt.plot(data[0], phi_mod, ':o', c = color)     #
        
        plt.show()

    #for x in range(1, num_stars+1):        
    #    plt.figure(x)    
    #    plt.show()
    
def plot_phi_all(num_stars, big_body, num):
    
    file_name = ''
    data = []
    
    big_data = rp.get_data(big_body)
    w_big = []
    M_big = []
    for x in range(len(big_data[0])):
        w_big.append(big_data[5][x])
        M_big.append(big_data[7][x])
    lambda_big = []
    for x in range(len(w_big)):
        lambda_big.append(w_big[x]+M_big[x])
    
    if num_stars % 5 == 0:
        counter = num_stars / 5
    else:
        counter = (num_stars / 5) + 1
    
    for x in range(1, counter+1):
        
        plt.figure(x)
        plt.title("Resonant Argument over Time [-180, 180]")
        plt.xlabel("t (years)")
        plt.ylabel("phi mod 360 (degrees)")
        plt.ylim(-180., 180.)
        
    #get data from each star, plot phi for each star
    for x in range(1, num_stars+1):
        file_name = 'S' + str(x) + '.aei'
        data = rp.get_data(file_name)
        lambda_little = []
        phi = []
        phi_mod = []
        w = data[5]
        M = data[7]
        color = [random.random(), random.random(), random.random()]
        for y in range(len(data[0])):
            lambda_little.append(w[y]+M[y])
            if num == 1:
                phi.append(2*lambda_big[y] - lambda_little[y] - w[y] + 180)     #from Table 8.1 ROW 1, add 180 to get between -180 and 180
            elif num == 2:
                phi.append(2*lambda_big[y] - lambda_little[y] - w_big[y] + 180)     #from Table 8.1 ROW 2
            phi_mod.append((phi[y] % 360) - 180)                              #phi mod 360 degrees, minus 180 to get correct frame
        num_plot = (x % counter) + 1
        plt.figure(num_plot)
        plt.plot(data[0], phi_mod, '-o', c = color)
        
    plt.show()
    
def plot_all_grid(num_stars, big_body, num):
    
    file_name = ''
    data = []
    
    big_data = rp.get_data(big_body)
    w_big = []
    M_big = []
    for x in range(len(big_data[0])):
        w_big.append(big_data[5][x])
        M_big.append(big_data[7][x])
    lambda_big = []
    for x in range(len(w_big)):
        lambda_big.append(w_big[x]+M_big[x])
    
    if num_stars % 9 == 0:
        counter = num_stars / 9
    else:
        counter = (num_stars / 9) + 1
    
    for x in range(1, counter+1):
        
        plt.figure(x)
        plt.title("Resonant Argument over Time [-180, 180]")   
    
        
    #get data from each star, plot phi for each star
    for x in range(1, num_stars+1):
        file_name = 'S' + str(x) + '.aei'
        data = rp.get_data(file_name)
        lambda_little = []
        phi = []
        phi_mod = []
        w = data[5]
        M = data[7]
        color = [random.random(), random.random(), random.random()]
        for y in range(len(data[0])):
            lambda_little.append(w[y]+M[y])
            if num == 1:
                phi.append(2*lambda_big[y] - lambda_little[y] - w[y] + 180)     #from Table 8.1 ROW 1, add 180 to get between -180 and 180
            elif num == 2:
                phi.append(2*lambda_big[y] - lambda_little[y] - w_big[y] + 180)     #from Table 8.1 ROW 2
            phi_mod.append((phi[y] % 360) - 180)                              #phi mod 360 degrees, minus 180 to get correct frame
        num_fig = (x % counter) + 1
        plt.figure(num_fig)
        num_plot = '33' + str((x % 9) + 1)
        plt.subplot(num_plot)
        plt.xlabel("t (years)")
        plt.ylabel("phi (degrees)")
        plt.ylim(-180., 180.)
        plt.plot(data[0], phi_mod, '-o', c = color)
        
    plt.show()

def plot_phi_phi_dot(num_stars, big_body, num):
    
    file_name = ''
    data = []
    
    big_data = rp.get_data(big_body)
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
    
    #get data from each star, plot phi for each star
    for x in range(1, num_stars+1):
        file_name = 'S' + str(x) + '.aei'
        data = rp.get_data(file_name)
        lambda_little = []
        phi = []
        phi_mod = []
        t = data[0]
        w = data[5]
        M = data[7]
        phi_dot = []
        #phi_dot_mod = []
        phi_prime = []
        phi_prime_mod = []
        t_prime = []
        color = [random.random(), random.random(), random.random()]
        for y in range(len(data[0])):
            lambda_little.append(w[y]+M[y])
            if num == 1:
                phi.append(2*lambda_big[y] - lambda_little[y] - w[y] + 180.)     #from Table 8.1 ROW 1
            elif num == 2:
                phi.append(2*lambda_big[y] - lambda_little[y] - w_big[y] + 180.)     #from Table 8.1 ROW 2
            phi_mod.append((phi[y] % 360.) - 180.)                              #phi mod 360 degrees
        for y in range(len(data[0]) - 1):
            phi_dot.append( (phi[y+1]-phi[y]) / (t[y+1]-t[y]) ) 
            phi_prime.append( (phi[y+1]+phi[y]) / 2. )
            t_prime.append( (t[y+1]+t[y]) / 2.)
            phi_prime_mod.append((phi_prime[y] % 360.) - 180.)
        plt.figure(x)
        plt.subplot(211)
        plt.title("Phi (mod) vs. Phi Dot")    #unicode
        plt.xlabel("phi (degrees)")
        #plt.xlim(-180., 180.)
        plt.ylabel("phi dot (degrees)")
        plt.plot(phi_prime_mod, phi_dot, 'o', c = color)
        plt.subplot(212)
        plt.title("Phi (mod) over Time")
        plt.xlabel("t (years)")
        plt.ylabel("phi mod 360 (degrees)")
        plt.plot(t_prime, phi_prime_mod, '-o', c = color)
        
        plt.show()
        

#MAIN




#plot_phi_all(25, 'Jupiter.aei', 1)
#plot_all_grid(25, 'Jupiter.aei', 1)
#plot_phi_phi_dot(25, 'Jupiter.aei', 1)
            

