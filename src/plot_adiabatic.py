'''
Created on November 30, 2015
A script to plot the adiabatic invariant over time, 
as given by Tremaine & Yu
'''

import resonance_plots as rp
import math
import random
import matplotlib.pyplot as plt

def plot_adiabatic(file_name, m_cen, p, q):     
    
    
        data = rp.get_data(file_name)
        a = data[1]
        e = data[2]
        i = data[3]
        
        
        color = [random.random(), random.random(), random.random()]
        
        C_pq = []
        
        for t in range(len(data[0])):
            C_pq.append( (m_cen * a[t])**0.5 * ((p+q) - p * (1 - e[t]**2)**0.5 * math.cos(math.radians(i[t]))) )
        
        plt.figure(1)
        
        plt.xlabel("t (years)")
        plt.ylabel("Adiabatic invariant C_pq")  
      
        plt.plot(data[0], C_pq, 'o', c = color)

        plt.show()
        
