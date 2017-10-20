'''
Created on September 28, 2015
@author: Maggie Mallernee

The main class where all function calls are placed
'''

import resonance_infiles
import file_o_projects

file_o_projects.frame_format('JUPITER', 10000000.0, 10000.0, 2000.0)           #name, m_cen, m_big, a_big
resonance_infiles.sample_format(25, 10000000.0, 10000.0, 2000.0, 1.0, 1.0, 0.3, 10.0)                 #num_stars, m_cen, m_big, a_big, p, q, e, i