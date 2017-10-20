'''
Created July 9, 2015

@author: Maggie Mallernee
Function that calculates widths of various resonances
'''

from scipy import constants as c

#09/17/15 -- separate the one for plots and just for values + switch to min/max

#@param: m_big mass of the perturbing body in solar masses
#@param: m_cen mass of central body
def get_width(m_cen, m_big, a_big, p, q, e):
    a = (p / (p+q)) ** (2./3.) * a_big
    G = c.G * c.day**2 * 1.98855 * 10**30 / c.au**3
    n = (a**3/(G))**-0.5              #assumes central body is one solar mass
    if q == 1.:                         #different values of alpha f sub d of alpha -- **could also make a dictionary
        if p == 1.:
            aFda = -0.749964
        elif p == 2.:
            aFda = -1.54553
        elif p == 3.:
            aFda = -2.34472
    elif q == 2.:
        if p == 1.:
            aFda = 0.287852
        elif p == 3.:
            aFda = 2.32892
    abs_Cr = abs((m_big/m_cen) * n * aFda)
    j2 = -1.*p
    
    max_del = a * (((16./3.)*(abs_Cr/n)*e)** 0.5 * (1. + (abs_Cr/(27.*j2**2.*e**3.)))**.5 - (2.*abs_Cr/(9.*j2*e*n)))    #assumes j4 = -1
    
    return [a - max_del, a + max_del]       #min, max
    
    
#returns a list of minimum/maximum a values, with as many values as time --> for plotting    
def get_width_list(m_cen, m_big, a_big, p, q, e, time):
    vals = get_width(m_cen, m_big, a_big, p, q, e, time)
    mini = vals[0]
    maxi = vals[1]
    
    a_max = [] 
    a_min = [] 
    for x in range(time):
        a_max.append(maxi)
        a_min.append(mini)
    
    return [a_min, a_max]