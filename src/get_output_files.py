'''
Created on October 18th, 2016
@author: Maggie Mallernee

This file includes methods to create summarizing output files.
'''

import resonance_ID_p1 as algo1
import resonance_ID_p2 as algo2


#method to create an executive summary file
#for each star, outputs the times it was in which resonances
#i.e. the resonances that pass through both parts of the res detection algo
def get_sum_file(num_stars, big_file_name, fudge):
    
    file = open("summary.out", 'w')
    file.write("Summary of Stars in Resonance\n")
    for x in range(1, num_stars+1):
        
        name = 'S' + str(x) + ".aei"
        write_to_sum(file, name, big_file_name, fudge)
        
    return file
        
        
#writes to a given file
#only for librating resonances, [ p, q, sub_res, start time, end time ]
def write_to_sum(out_file, star_file_name, big_file_name, fudge):

    out_file.write(star_file_name + '\n')    
    res_table = algo1.get_res_timetable(star_file_name, big_file_name)
    phi_table = algo2.get_phi_table(star_file_name, big_file_name, res_table, fudge)
    #row { p, q, sub_res, time[], phi[], circ test result }
    for n in range(len(phi_table)):
        if(phi_table[n][5] == 0):
            out_file.write("     p: %s, q: %s, sub_res: %s, t_start: %s, t_end: %s\n" % (phi_table[n][0], phi_table[n][1], phi_table[n][2], phi_table[n][3][0], phi_table[n][3][len(phi_table[n][3]) - 1]))