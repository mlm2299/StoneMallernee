'''
Created on September 17, 2015
@author: Maggie Mallernee

This file automates the process of resonance testing; given file of orbital elements, it assesses whether or not the body is in resonance
Given a .aei file
'''
import aei_reader as r
import width_calculator as wc
import get_resarg

#use width calculator to produce ranges -- need to know range of eccentricities & resonances to sample
#eliminate stars who do not stay in width


#returns a multidimensional array of stars by resonance by in/out at a given timestep (index)
class Res_tester(object):
    
    def __init__(self):
        all_stars = []
        
        res21 = []
        res32 = []
        res43 = []
        res31 = []
        res53 = []
            
        resonances = [res21, res32, res43, res31, res53]    #0 to 4
    
    def get_resonances(self):   return self.resonances
    def get_allstars(self):     return self.all_stars
    
    def width_tester(num_stars, m_cen, m_big, time):    #do per star
        
        all_stars = Res_tester.get_allstars()
        resonances = Res_tester.get_resonances()
        
        for x in range(1, num_stars+1):     #for every star ...
            file_name = 'S' + str(x) + '.aei'
            data = r.get_data(file_name)        #such a waste of space, list of lists
            a = r.get_a(data)               #list of a over time
            e = r.get_e(data)               #list of e over time
            
            #test for each resonance
            #function get_width(m_cen, m_big, a_big, p, q, e, time) - need to determine which of these should be parameters
            #^^also need to know which value of e to use
            #^^returns array of [min, max]
            
            
            #compute and compare instantaneous width with actual values of a
            for t in range(time):       #need to ensure that is the correct range
                bounds21 = wc.get_width(m_cen, m_big, 1., 1., e[t], time)     #array of [min, max]
                bounds32 = wc.get_width(m_cen, m_big, 2., 1., e[t], time)
                bounds43 = wc.get_width(m_cen, m_big, 3., 1., e[t], time)
                bounds31 = wc.get_width(m_cen, m_big, 1., 2., e[t], time)
                bounds53 = wc.get_width(m_cen, m_big, 3., 2., e[t], time)
                
                bounds = [bounds21, bounds32, bounds43, bounds31, bounds53]     #indexed 0 to 4
                
                for num, minmax in enumerate(bounds):
                    if a[t] >= minmax[0] and a[t] <= minmax[1]:
                        resonances[num].append(1)               #choosing 1s == yes in resonance (that way can test numerically)
                    else:
                        resonances[num].append(0)
            
            all_stars.append(resonances)
        
        return all_stars
        
    #function to eliminate stars not within the resonant width -- written to read the above produced array
    def in_width():
        all_stars = Res_tester.get_allstars()
        resonances = Res_tester.get_resonances()        #do these get the updated versions?
        
        strike = 0
        
        for star in all_stars:
            for res in resonances:
                if sorted(res)[0] == 0:
                    strike = strike + 1
            if strike == 5:                 #if the star is outside the width for every resonance at a given time -- seems too general? look instead for continuity?
                all_stars.remove(star)      #will lose the ability to track individual stars -- important?
            strike = 0
                 
    def circ_test(time):
        
        all_stars = Res_tester.get_allstars()
                
        
        
