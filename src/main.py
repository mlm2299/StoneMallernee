'''
Created on December 9, 2015
@author: Maggie Mallernee
'''

import parse
#import plot_adiabatic as pa
#import star_sampler1 as get_small_in
#import plot_am_zcomp as L
#import file_o_projects as get_big_in  #THIS THE ONE THAT PRINTS
#import resarg_plots as resp
#import resonance_plots as rp
#import resonance_ID_p1 as algo1
#import resonance_ID_p2 as algo2
#import get_output_files as out


'''CREATING INPUT FILES'''
#get_small_in.star_sample_file_writer("small.in", 10**7, 100, 30)     #m_cen, num_stars_to_sample, num_to_include (smallest by peri)
#get_small_in.sample_multiple(10**7, 100, [1, 3, 10, 30])
#get_big_in.frame_format('JUPITER', 10000000.0, 100000.0, 100.0, 100.0)       #name, m_cen, m_big, a_big, r_start (equality of a_big and r_start means eccentricity of 0)

'''READING OUTPUT FILES'''
#parse.get_special("info.out")   #reads info.out
#tde_file_names = parse.get_tde_file_names()

'''PLOTTING DATA'''
#parse.plot_special("info.out")  #plots a, e, i of the special cases

#parse.plot_TDE("info.out")

#rp.plot_elements("S70.aei", "JUPITER.aei", 10000000.)

#pa.plot_adiabatic("S5.aei", 10000000., 1., 1.)   #num_stars, m_cen, p, q

#rp.plot_res_ratio("S70.aei", "JUPITER.aei", 1., 1., 100000., 10000000.)

#L.plot_z_ang_mom(25)  #num_stars

#resp.plot_res_arg1(95, 'Jupiter.aei', 1., 1., 1.)#num_stars, big_body, p, q, row

'''AUTOMATIC RESONANCE CHECKING'''
#star_file_name = "S1.aei"
#res_table = algo1.get_res_timetable(star_file_name, "JUPITER.aei")

#phi_table = algo2.get_phi_table(star_file_name, "JUPITER.aei", res_table, 10.)
#algo2.plot_all_phi(phi_table)

#for special bodies
#for star_file_name in tde_file_names:
#    res_table = algo1.get_res_timetable(star_file_name, "JUPITER.aei")
#
#    phi_table = algo2.get_phi_table(star_file_name, "JUPITER.aei", res_table)
#    algo2.plot_all_phi(phi_table)

'''OUTPUT FILE PRODUCTION'''
#out.get_sum_file(25, "JUPITER.aei", 1.)

#parse.get_special_file("info.outE")

#parse.get_info_single_batch("analysisTest1.out", 32, 30)
#parse.plot_tde_vs_t(32)
#parse.write_tde_times(32)
#parse.plot_tde_vs_e(32, 10000.)
   

