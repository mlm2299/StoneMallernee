'''
Created on June 10, 2015

@author: Maggie Mallernee

file input project:
    - read in .aei files
    - take first 6 (orbital) elements
    - convert them to cartesian
    - compare these to the last 6 (cartesian) elements from the file
    - return % error for each element
'''
import converter_clean as c
from math import radians

'''
Function to compute experimental error
@param: exp The experimental value
@param: acc The accepted value

@return the value of the error in percent
'''
def exp_error(exp, acc):
    return abs(exp - acc)/acc * 100
    
'''
TODO #1 function

@param: file containing all 12 orbital elements in the format of mercury .aei files
@param: mass of the more massive object in the binary
'''
def read(file_name, m1): 
    
    with open(file_name) as f:   #assume only read (r) mode
        #create list of each non-empty line in a given file
        lines = []
        for line in f:
            if not line.isspace():
                lines.append(line)
    
    #separate orbital elements -- assuming the name and variable names are on separate lines
    m = lines[2].split()     #should be a list of t, a, e, i, omega, w, f, mass, x, y, z, u, v, w from mercury
    
    #FROM ORBITAL TO CARTESIAN
    #create and print orbit object with orbital elements
    orbit_o = c.O_Orbit(float(m1), float(m[7]), float(m[1]), float(m[2]), radians(float(m[3])), radians(float(m[4])), radians(float(m[5])), radians(float(m[6])))
    print "\nThe orbital elements from the infile:\n"
    orbit_o.print_out_d()                                  
    
    #create and print orbit object with elements converted to cartesian
    orbit_c = c.convert_from_orbital(orbit_o)
    print "\nCartesian elements converted from orbital input:\n"
    orbit_c.print_out()
    print "\n"
    
    elem = orbit_c.element_list()
    x = orbit_c.element_val_list()  #list of 'experimental' values for cartesian elements
    
    #calculate % error between converted cartesian elements and those read from file
    for n in range(2, 8):
        print "Percent error of %s conversion: %s" % (elem[n], str(exp_error(x[n], float(m[n+6]))))
    
    #FROM CARTESIAN TO ORBITAL
    orbit_c = c.C_Orbit(float(m1), float(m[7]), float(m[8]), float(m[9]), float(m[10]), float(m[11]), float(m[12]), float(m[13]))
    print "\nThe cartesian elements from the infile: \n"
    orbit_c.print_out()
    
    orbit_o = c.convert_from_cartesian(orbit_c)
    print "\nThe orbital elements converted from cartesian input:\n"
    orbit_o.print_out_d()
    print "\n"
    
    elem = orbit_o.element_list()
    x = orbit_o.element_val_list_d()  #'experimental' values for orbital elements
    
    for n in range(2, 8):
        print "Percent error of %s conversion: %s" % (elem[n], str(exp_error(x[n], float(m[n-1]))))
    
read('JUPITER.txt', 2.00000)