'''
Created on June 16, 2015

@author: Maggie Mallernee

To sample stellar orbits assuming isotropy and spherical distribution, with angles in DEGREES
'''

#from scipy import constants as c
import random
from math import cos, asin, degrees
import anomaly_convert as ac

#Function that calculates the radius of influence (of a black hole), in parsecs
#Need to verify units & empirics (r_infl ~ 0.76 ( M_bh / 10^6 M_s)^(0.58)
def get_r_infl(mass):
    return (mass / (10**6 * 1.98855 * 10**30)) ** 0.5
    
#def get_r_max
    
def star_sampler(mass_bh, num_stars):
    
    sample_array = []
    
    #name = ""
    
    for x in range(1, num_stars+1):
        
        #name = 'S' + str(x)
        
        #randomly sample omega, w, M (mean anomaly) **in DEGREES
        omega = random.uniform(0.0, 360.)
        w = random.uniform(0.0, 360.)
        M = random.uniform(0.0, 360.)
        
        #randomly sample i
        sin_i = random.uniform(0.0, 2.0) - 1.0
        i = degrees( asin(sin_i) )
        
        #randomly sample e using its cdf F(X) = e**2
        #thermal distribution of e
        is_ok = False
        while(is_ok != True):
            X = random.random()
            e = X ** 0.5
            if(e != (1.0 / cos(1.0) )):
                is_ok = True
        
        
        #get f from M
        E = ac.newton(e, 1.0, M)      #need method to choose the initial guess
        f = ac.get_f(e, E)
        
        #randomly sample a using the Bahcall-Wolf Cusp distribution
        X = random.random()
        #a = X ** (4./5.) * 1240     #4:2 outer bound?
        #for initial GR computational testing -- bound to the 3:1 resonance 
        #(48.0750 - for 100 AU binary separation) 
        #(14.4225 - for 30 AU binary separation)
        a = X ** (4./5.) * (48.0750)   
        
        #radius of pericentre -- for sorting purposes only
        r_p = a * (1 - e)
        
        sample_array.append([a, e, i, omega, w, f, r_p])
        
    return sample_array
  
#sorts by pericentre, ascending  
def sort_sample_by_peri(sample):
    
    sample.sort(key=lambda x: x[6])
    
    return sample
        

def star_sample_file_writer(file_name, mass_bh, num_stars_to_sample, num_to_include):
    
    sample = sort_sample_by_peri(star_sampler(mass_bh, num_stars_to_sample))
    
    get_file(file_name, sample, num_to_include)
    
#to_include ought to be a list of different numbers to sample
def sample_multiple(mass_bh, num_stars, to_include):
    
    sample = sort_sample_by_peri(star_sampler(mass_bh, num_stars))
    
    for n in range(len(to_include)):
        
        file_name = "small" + str(to_include[n]) + ".in"
        
        get_file(file_name, sample, to_include[n])
        

#returns a file with the smallest num_included number of stars with pericentres greater than twice the tidal radius
#r_t currently set to 1.0016247683 AU
def get_file(file_name, sample, num_included):
    
    with open(file_name, 'w') as file:
            
            file.write(')O+_06 Small-body initial data  (WARNING: Do not delete this line!!) \n')
            file.write(") Lines beginning with `)' are ignored. \n")
            file.write(")--------------------------------------------------------------------- \n")
            file.write(" style (Cartesian, Asteroidal, Cometary) = Ast \n")
            file.write(")--------------------------------------------------------------------- \n")
            
            name = ""
            
            r_t = 1.0016247683
            
            n = 0
            i = 1
            while i <= num_included: 
                if sample[n][6] > (2 * r_t):
                    name = "S" + str(i)
                    i += 1
                    
                    #write to the file
                    #TEST file.write("Test - radius of pericentre: %s \n" % (sample[n][6]))
                    file.write(' %s    ep=2450400.5 \n' % (name))
                    
                    if i > num_included:
                        file.write(" %s %s %s %s %s %s 0 0 0\n\n" % (sample[n][0], sample[n][1], sample[n][2], sample[n][3], sample[n][4], sample[n][5]))
                    else:
                        file.write(" %s %s %s %s %s %s 0 0 0\n" % (sample[n][0], sample[n][1], sample[n][2], sample[n][3], sample[n][4], sample[n][5]))

                n += 1   
            #name = "S" + str(num_included) 
            #file.write(' %s    ep=2450400.5 \n' % (name))
            #file.write(" %s %s %s %s %s %s 0 0 0  " % (sample[n][0], sample[n][1], sample[n][2], sample[n][3], sample[n][4], sample[n][5]))
                    
    return file
    
    
#MAIN

            
            
            
        