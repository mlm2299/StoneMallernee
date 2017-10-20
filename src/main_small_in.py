'''
Main for sampling stars & generating small.in
@author: Maggie Mallernee
November 28, 2016
'''

import star_sampler1 as get_small_in


'''CREATING INPUT FILES'''
#for x in range(32):
    #name = "small.in" + str(x)
    
get_small_in.star_sample_file_writer("small.in", 10**7, 100, 30)     #m_cen, num_stars_to_sample, num_to_include (smallest by peri)