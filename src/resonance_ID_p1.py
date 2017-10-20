'''
Created on January 26, 2016
Part 1 of a two-part algorithm to identify and isolate resonances:
    For a given .aei file, goes through a series of resonances and eliminates impossible conclusions based on
    resonant width over a minimum time interval based on the libration period of the specific resonance.
    Should leave each star with at most a list of possible resonances and corresponding time intervals.
'''
from scipy import constants as c

#gets the value of a over time for a given aei file
def get_data_part1(file_name):
    
    with open(file_name) as f:   #assume only read (r) mode
        #create list of each non-empty line in a given file
        lines = []
        for line in f:
            if not line.isspace():
                lines.append(line)
    
    over_time = []
    
    for n in range(2, len(lines)):
        x = lines[n].split()    #a list of t, a from mercury
        for y in range(0, len(x)):  #might not work with the scientific notation?
            if x[y] == '***************':
                x[y] = '400'
            x[y] = float(x[y])
        over_time.append(x)     #a list of lists of elements, one per time interval, so [t0,a0], [t1,a1] ...
        
    t = []
    a = []
    e = []
    
    for n in range(len(over_time)): 
        t.append(over_time[n][0])
        a.append(over_time[n][1])
        e.append(over_time[n][2])
        
    return [t, a, e]

#reads info.out for the central mass, returns as a float
def get_central_mass():
    m_cen = 0.0
    with open("info.out") as f:
        words = []
        for line in f:
            if("Central mass:" in line):
                    #get the central mass / convert from sci notation
                    words = line.split()
                    m = words[2]
    
    m = m.split('E+')
    m_cen = float(m[0]) * pow(10.0, float(m[1]))
    return m_cen
    
#reads big.in for the mass of the secondary body
def get_big_mass():
    m_big = 0.0
    with open("big.in") as f:
        words = []
        for line in f:
            if("m=" in line):
                    words = line.split()
                    m = words[1]
    
    m = m.strip('m=')
    m_big = float(m)
    
    return m_big
    
    
#verified
#NOT resonance width, but rather error bounds for the ratio ( i.e. 2:1 ~ 0.5 but not quite )
def get_res_ratio_min_max(p, q, sub_res, a, e, a_big, m_big, m_cen):  
    ideal = 0.0
    minimum = 0.0
    maximum = 0.0
    
    reduce(p, q)
    
    ideal = (float) (p) / (p + q)
    
    if q == 1.:   
        if sub_res == 0.:
            j4 = -1.    #always = -1 for first order while we assume jupiter does not precess
        elif sub_res == 1.:              
            j4 = 0.
    elif q == 2.:
        if sub_res == 0.:
            j4 = -2.    #a bit of an oversimplification but fine for now
        elif sub_res == 1.:
            j4 = 0.     #ROW 6
    else:
        j4 = -1.
        
    aFda = get_aFda(p, q)
    
    n = (a**3/(G * m_cen))**-0.5  
    C_r = (m_big/m_cen) * n * aFda
    j2 = -1. * p    
    abs_Cr = abs(C_r)
    
    if q == 1.:
        dA = a * (((16./3.)*(abs_Cr/n)*e)** 0.5 * (1. + (abs_Cr/(27.*j2**2.*e**3.)))**.5 - (2.*abs_Cr/(9.*j2*e*n)))    #assumes j4 = -1
    elif q == 2.:
        dA = a * (((16./3.)*(abs_Cr/n)* pow(e, abs(j4)))** 0.5)
            
    dP = 3 * c.pi * pow(a/(G * m_cen), 0.5) * dA      
            
    p_big = 2 * c.pi * pow(a_big, 1.5) * pow(G*m_cen, -0.5)     #period of Jupiter
            
    minimum = ideal - (dP/p_big)
    maximum = ideal + (dP/p_big)
   # if p == 1. and q == 1.:
   #     print ideal             #returning 0.5 for the graph but 0.0 below
    return [minimum, maximum]
        
            
    #res_ratio_widths.append([p, q, minimum, maximum])   #stores a row in res_ratio_widths containing p, q, and min/max ratio values
#verified
def get_lib_period(p, q, sub_res, a, e, a_big, m_big, m_cen):
    reduce(p, q)
    
    n = (a**3/(G * m_cen))**-0.5  
    #print n  
    aFda = get_aFda(p, q)   
    if q == 1.:   
        if sub_res == 0.:
            j4 = -1.    #always = -1 for first order while we assume jupiter does not precess
        elif sub_res == 1.:              
            j4 = 0.
    elif q == 2.:
        if sub_res == 0.:
            j4 = -2.    #a bit of an oversimplification but fine for now
        elif sub_res == 1.:
            j4 = 0.     #ROW 6
    else:
        j4 = -1. #ASSUMED DEFAULT  
    C_r = (m_big/m_cen) * n * aFda
   # print C_r
    j2 = -1. * p    
    w_ = pow(abs(-3.* (j2 ** 2) * C_r * n * (e ** abs(j4))) , 0.5)     #ADDED ABS, for when q = 2, aFda > 0, so trying to take sqrt of a negative value  
    #print w_ #TEST
    if w_ == 0.0:
        w_ = 5000. 
    T_lib = 4.*c.pi / w_    #Multiplied by 2 twice to approximate the elliptic integral
    T_lib = T_lib / 365.0 #convert from days to years
    
    return T_lib
    
#tabulate p | q | t_start | t_stop | t1 line # | t2 line #
#for a minimum number of libration periods (maybe start with 4)

def get_res_timetable(star_file_name, big_file_name):     #NEED WAY TO READ OTHER MASSES AUTOMATICALLY - info.out
    star_data = get_data_part1(star_file_name)   #3D array [t, a, e]
    big_data = get_data_part1(big_file_name)
    m_cen = get_central_mass()
    m_big = get_big_mass()
    
    print "Central mass: %s" % m_cen
    print "Secondary mass: %s" % m_big
    print "Small body: %s" % star_file_name
    
    res_timetable = []
    
    t_start = 0
    t_end = 0
    
    line_start = 0
    line_end = 0
    
    
    for p in range(MAX_P_VAL):
        p = p+1
        for q in range(MAX_Q_VAL):
            q = q + 1
            reduce(p, q)
            for sub_res in range(NUM_SUB_RESONANCES[q-1]): 
                in_res = False
                meets_min = False
                # print "%s, %s" % (p, q)
                for t in range(len(star_data[0])):
                    a_big = big_data[1][t]
                    a_star = star_data[1][t]
                    if a_star < 0: 
                        print star_data[0][t]
                        a_star = abs(a_star)
                    e_star = star_data[2][t]
                    T_lib = get_lib_period(p, q, sub_res, a_star, e_star, a_big, m_big, m_cen)      #UPDATE THIS FUNCTION
                    # print T_lib
                    min_max = get_res_ratio_min_max(p, q, sub_res, a_star, e_star, a_big, m_big, m_cen)      #THIS ONE TOO
                    
                    T_ratio = pow((a_star / a_big), (3./2.))
                    
                    if T_ratio > min_max[0] and T_ratio < min_max[1]: 
                        #print in_res
                        if not in_res:
                            #print in_res
                            in_res = True
                            meets_min = False
                            t_start = star_data[0][t]
                            line_start = t
                        elif star_data[0][t] - t_start > MIN_NUM_LIB * T_lib:
                            meets_min = True
                            if t == len(star_data[0]) - 1:  #if at the end & meets min
                                t_end = star_data[0][t]
                                line_end = t
                                res_timetable.append([p, q, sub_res, t_start, t_end, line_start, line_end])         #UPDATE WHAT DEPENDS ON THIS
                                #print "p: %s, q: %s, t_start: %s, t_end: %s, line_start: %s, line_end: %s" % (p, q, t_start, t_end, line_start, line_end)
                    elif meets_min:
                        t_end = star_data[0][t - 1]  #make note of real t vs. line #
                        line_end = t - 1
                        res_timetable.append([p, q, sub_res, t_start, t_end, line_start, line_end])         #THIS TOO
                        #print "p: %s, q: %s, t_start: %s, t_end: %s, line_start: %s, line_end: %s" % (p, q, t_start, t_end, line_start, line_end)
                        in_res = False
                        meets_min = False
                    else:
                        meets_min = False
    
    return res_timetable
    
                    

#HOW TO DEAL WITH REPEATS (IE 2:1 == 4:2)  
def reduce(p, q):
    for c in range(MAX_P_VAL + 1):
        if( c>1 and (p+q)%c == 0 and p%c == 0 ):
            p_ = p / c
            q_ = (p_*c + q - p) / c
            p = p_
            q = q_
            #return true
            
def is_duplicate(p, q):
    for c in range(MAX_P_VAL + 1):
        if( c>1 and (p+q)%c == 0 and p%c == 0 ):
            return True;
    return False;
            
def get_aFda(p, q): #***NEED DIFFERENTIATION BETWEEN DIFFERENT SUBRESONANCES***#
    reduce(p, q)
    aFda = 0.0
    
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
    
    return aFda
    
#def plot
        
             
            
G = c.G * c.day**2 * 1.98855 * 10**30 / c.au**3  
#print G                                                    

#potential_resonances = []

MAX_P_VAL = 3
MAX_Q_VAL = 2

MIN_NUM_LIB = 4

NUM_SUB_RESONANCES = [2, 2]


