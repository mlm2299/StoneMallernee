'''
Created on September 17, 2015
@author: Maggie Mallernee

Reads .aei files and returns the data enclosed
'''

#function taken from resonance_plots, written June 24, 2015
#@param: aei_f the filename of a .aei file from Mercury
#assumes output format of t, a, e, i, omega, w, f, mass, x, y, z, u, v, w, M in that order
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
    #print aei_f
    #print over_time[0]
    t = []
    a = []
    e = []
    i = []
    omega = []
    w = []
    f = []
    M = []
    
    for n in range(len(over_time)-1):     #maybe make each element a set of key value pairs with time?
        t.append(over_time[n][0])
        a.append(over_time[n][1])
        e.append(over_time[n][2])
        i.append(over_time[n][3])
        omega.append(over_time[n][4])
        w.append(over_time[n][5])
        f.append(over_time[n][6])   
        M.append(over_time[n][14])
        
    return [t, a, e, i, omega, w, f, M]
    
#set of functions to get individual elements
#takes parameter data, the output of get_data
def get_t(data):
    return data[0]
    
def get_a(data):
    return data[1]
    
def get_e(data):
    return data[2]
    
def get_i(data):
    return data[3]
    
def get_omega(data):
    return data[4]
    
def get_w(data):
    return data[5]
    
def get_f(data):
    return data[6]
    
def get_M(data):
    return data[7]
