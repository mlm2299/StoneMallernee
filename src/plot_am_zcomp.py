

import matplotlib.pyplot as plt
import random
import math
import resonance_plots as rp


#PLOT Z COMP. OF ANGULAR MOMENTUM   
def plot_z_ang_mom(num_stars):
    
    file_name = ''
    data = []
    
    #get data from each star, plot z component of angular momentum
    for x in range(1, num_stars+1):
        file_name = 'S' + str(x) + '.aei'
        data = rp.get_data(file_name)
        t = data[0]
        e = data[2]
        i = data[3]
        z_comp = []
        
        for n in range(len(t)):
            z_comp.append((1. - (e[n]**2)) ** 0.5 * math.cos(math.radians(i[n])))
        
        
        color = [random.random(), random.random(), random.random()]
        plt.figure(x)
        plt.plot(t, z_comp, ':o', c = color)
        plt.xlabel("time t (years)")
        plt.ylabel(" (1-e^2) ^ (1/2) * cos(i) ")
        plt.show()




