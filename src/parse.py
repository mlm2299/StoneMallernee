'''
Created on December 7, 2015
@author Maggie Mallernee

Designed to parse info.out and plot special (TDE/Ejection stars)
'''
#import resonance_plots as rp
import aei_reader
import math
import matplotlib as plt
import pylab as P

tde_bodies = []
sec_tde_bodies = []
eject_bodies = []

def get_tde_times(num_runs):
    tde_bodies = []
    sec_tde_bodies = []

    for x in range(num_runs):
        
        info_out_name = "run" + str(x) + "Info.out"
        
        with open(info_out_name) as f:   #assume only read (r) mode
            #create list of each non-empty line in a given file
            lines = []
                
            for line in f:
                if not line.isspace():
                    lines.append(line)
            
            for line in lines:
                if("collided with the central body at" in line):
                    if(not "JUPITER" in line):
                        tde_bodies.append(int(float(line[48:62].strip())))
                elif("was hit by" in line):
                    sec_tde_bodies.append(int(float(line[37:48].strip())))
    
    return [tde_bodies, sec_tde_bodies]     
    
def get_tde_times_names(num_runs):
    tde_bodies = []
    sec_tde_bodies = []

    for x in range(num_runs):
        
        info_out_name = "run" + str(x) + "Info.out"
        
        with open(info_out_name) as f:   #assume only read (r) mode
            #create list of each non-empty line in a given file
            lines = []
                
            for line in f:
                if not line.isspace():
                    lines.append(line)
            
            for line in lines:
                if("collided with the central body at" in line):
                    if(not "JUPITER" in line):
                        tde_bodies.append(["run" + str(x) + line[1:4], line[48:62]])
                elif("was hit by" in line):
                    sec_tde_bodies.append(["run" + str(x) + line[21:24], line[37:48]])
    
    return [tde_bodies, sec_tde_bodies]  
    
def write_tde_times(num_runs):
    tdes = get_tde_times(num_runs)
    tde_bodies = tdes[0]
    sec_tde_bodies = tdes[1]
    tde_bodies.sort()
    sec_tde_bodies.sort()
    
    
    file = open("TDE.out", 'w')
    file.write("Times of TDEs, primary then secondary\n")
    for x in range(len(tde_bodies)):
        file.write("%s\n" % tde_bodies[x])
    file.write("\n")
    for y in range(len(sec_tde_bodies)):
        file.write("%s\n" % sec_tde_bodies[y])
    return file
    
def plot_tde_vs_t(num_runs):
    tdes = get_tde_times(num_runs)
    tde_bodies = tdes[0]
    sec_tde_bodies = tdes[1]
                    
    P.figure()
    n, bins, patches = P.hist([tde_bodies, sec_tde_bodies], 2600, stacked=True)
    P.title("TDEs over Time")
    maxim = max(max(tde_bodies), max(sec_tde_bodies))
    P.xlim(xmin = maxim - 500., xmax = maxim)
    P.xlabel("t (years)")
    P.ylabel("# of TDEs")
    
    # FOR HEIGHT LABELS
    width = bins[1] - bins[0]
    for x in range(len(n[1])):
        height = n[1][x]   # total height including both primary and secondary TDEs
        P.text(bins[x] + width/2, height + 5, int(height), ha='center', va='bottom')
        
    P.show()

# Plot histogram of final eccentricities with TDEs in the last ~x years
def plot_tde_vs_e(num_runs, x):
    tdes = get_tde_times_names(num_runs)
    tde_bodies = tdes[0]
    final_ecc = []
    for tde in tde_bodies:
        file_name = tde[0].strip() + ".aei"
        t = tde[1]
        if t > (26000. - x):
            data = aei_reader.get_data(file_name)
            e = aei_reader.get_e(data)
            #print len(e)
            if len(e) > 0:
                final_ecc.append(e[len(e) - 1])
    
    P.figure()
    n, bins, patches = P.hist(final_ecc, 250)
    P.title("TDEs by Final Eccentricity")
    P.xlim(xmin = 0, xmax = max(e))
    P.xlabel("Final Eccentricity")
    P.ylabel("# of TDEs")
    
    # FOR HEIGHT LABELS
    #width = bins[1] - bins[0]
    #for x in range(len(n[1])):
    #    height = n[1][x]   # total height including both primary and secondary TDEs
    #    P.text(bins[x] + width/2, height + 5, int(height), ha='center', va='bottom')
        
    P.show()
        
def get_info_single_batch(info_file_name, num_runs, num_stars):
        
    #read info.outs
    tde_bodies = []
    sec_tde_bodies = []
    eject_bodies = []
    mergers = []
    num_stars = 0
             
    for x in range(num_runs):
        
        tde_bodies_this_run = []
        sec_tde_bodies_this_run = []
        ### PARSE FROM INFO.OUT FILES ###
        info_out_name = "run" + str(x) + "Info.out"
            
        with open(info_out_name) as f:   #assume only read (r) mode
            #create list of each non-empty line in a given file
            lines = []
                
            for line in f:
                if not line.isspace():
                    lines.append(line)
                
            head_line = ""
            
                
            for line in lines:
                if("collided with the central body at" in line):
                    if(not "JUPITER" in line):
                        #tde_bodies.append(["run" + str(x) + line[1:4], line[48:62]])
                        tde_bodies_this_run.append(["run" + str(x) + line[1:4], line[48:62]])
                    else:
                        mergers.append(["run" + str(x) + "JUPITER", line[48:62]])
                elif("was hit by" in line):
                    sec_tde_bodies_this_run.append(["run" + str(x) + line[21:24], line[37:48]])
                elif("ejected at" in line):
                    eject_bodies.append(["run" + str(x) + line[1:4], line[25:37]])
                elif("Number of Small bodies" in line):
                    #head_line = line
                    num_stars += int(float(line[30:]))
        survivors = num_stars - len(tde_bodies) - len(sec_tde_bodies) - len(eject_bodies)
        
        #add times to / from merger
        if x > 0:
            print "RUN: %s" % x
            jupiter_aei_data = aei_reader.get_data("run" + str(x) + "JUPITER.aei")
            #if x == 1:
            #    for i in range(len(jupiter_aei_data[0])):
            #        print i
            for n in range(len(tde_bodies_this_run)):
                star_aei_data = aei_reader.get_data((tde_bodies_this_run[n][0]).strip() + ".aei")
                #print len(star_aei_data[0])   # why not ???
                delta_t = 10.0239562 
                int_time_tde = int(round(float(tde_bodies_this_run[n][1]),-1)/delta_t)  #convert to an index
                #if x == 0:
                print "N: %s" % n
                print "Time of TDE: %s" % tde_bodies_this_run[n][1]
                print "Time of TDE index: %s" % int_time_tde
                print "Resulting time from aei: %s" % jupiter_aei_data[0][int_time_tde]
                print "TDE: %s" % tde_bodies_this_run[n][0]
                if int_time_tde > 1:
                    #if x == 0:
                    #    for i in range(len(jupiter_aei_data[0])):
                    #        print "Binary time vs index: %s" % jupiter_aei_data[0][i]
                    tde_bodies.append([tde_bodies_this_run[n][0],  #filename
                                    tde_bodies_this_run[n][1], #time of TDE
                                    int(float(mergers[len(mergers)-1][1])) - int_time_tde, #time of TDE wrt merge time
                                    jupiter_aei_data[1][int_time_tde], #binary separation --- assumes aei data exists for each year
                                    jupiter_aei_data[2][int_time_tde], #binary eccentricity
                                    star_aei_data[3][len(star_aei_data)-1]]) #final orbital inclination
                               
        #for n in range(len(sec_tde_bodies_this_run)):
        #    star_aei_data = aei_reader.get_data((sec_tde_bodies_this_run[n][0]).strip() + ".aei")
        #    int_time_tde = int(float(sec_tde_bodies_this_run[n][1]))/10
        #    sec_tde_bodies.append([sec_tde_bodies_this_run[n][0],  #filename
        #                       sec_tde_bodies_this_run[n][1], #time of TDE
        #                       int(float(mergers[len(mergers)-1][1])) - int_time_tde, #time of TDE wrt merge time
        #                       jupiter_aei_data[int_time_tde][1], #binary separation --- assumes aei data exists for each year
        #                       jupiter_aei_data[int_time_tde][2], #binary eccentricity
        #                       star_aei_data[len(star_aei_data)-1][3]]) #final orbital inclination
        
    file = open(info_file_name, 'w')
    file.write("Number of stars: %s\n" % (str(num_stars)))
    file.write(str(len(tde_bodies)) + " central TDE stars \n")
    file.write(str(len(sec_tde_bodies)) + " secondary TDE stars \n")
    file.write(str(len(eject_bodies)) + " Ejections \n")
    file.write(str(survivors) + " Survived\n\n")  
    
    #add percentage time since start & percentage until merge
    
    #also include physical time until end of inspiral
    #can prob eliminate physical time since start
    file.write(str(len(tde_bodies)) + " central TDE stars: \n")
    file.write("    " + "TDE" + "    " + "Time of Disruption" + "\n")
    for n in range(len(tde_bodies)):
        file.write("    " + tde_bodies[n][0] + "    " + tde_bodies[n][1] + "\n")
    file.write(str(len(sec_tde_bodies)) + " secondary TDE stars: \n")
    for n in range(len(sec_tde_bodies)):
        file.write("    " + sec_tde_bodies[n][0] + "    " + sec_tde_bodies[n][1] + "\n")
    file.write(str(len(eject_bodies)) + " ejections: \n")
    for n in range(len(eject_bodies)):
        file.write("    " + eject_bodies[n][0] + "    " + eject_bodies[n][1] + "\n")
    #ejection velocity -- convert from velocity as leave the grid to velocity at infinity (list both)
    file.write('\n' + str(survivors) + " Survived\n")   
    
    #to plot --- number of disruptions vs. percent of time since start histogram
    #second plot --- number of disruption vs. percent of time until merger histogram (might want log axis depending)
    return file
        
      
            
        
    #count of primary TDEs, secondary TDEs, ejections, survivors
    
    #include time of disruption for TDEs, time of merger
    
    #ejection velocity
    
    #last time-step orbital elements for survivors


def get_special_file(info_file_name):
    
    with open(info_file_name) as f:   #assume only read (r) mode
        #create list of each non-empty line in a given file
        lines = []
        
        for line in f:
            if not line.isspace():
                lines.append(line)
        
        del tde_bodies[:]   #clear the lists
        del sec_tde_bodies[:]
        del eject_bodies[:]
        
        head_line = ""
        num_stars = 0
        
        for line in lines:
            if("collided with the central body at" in line):
                if(not "JUPITER" in line):
                    tde_bodies.append(line[2:4])
            elif("was hit by" in line):
                sec_tde_bodies.append(line[22:24])
            elif("ejected at" in line):
                eject_bodies.append(line[2:4])
            elif("Number of Small bodies" in line):
                head_line = line
                num_stars = int(float(line[30:]))
    survivors = num_stars - len(tde_bodies) - len(sec_tde_bodies) - len(eject_bodies)
    file = open("basic_run_stats.dat", 'w')
    file.write(head_line + '\n')
    file.write(str(len(tde_bodies)) + " central TDE stars: " + ', '.join(tde_bodies) + '\n')
    file.write(str(len(sec_tde_bodies)) + " secondary TDE stars: " + ', '.join(sec_tde_bodies) + '\n')
    #add time of disruption & separate by TDE with the primary vs. TDE with the secondary
    file.write('\n' + str(len(eject_bodies)) + " Ejections: " + ', '.join(eject_bodies) + '\n')
    #ejection velocity -- convert from velocity as leave the grid to velocity at infinity
    file.write('\n' + str(survivors) + " Survived")
    #will want their orbital elements at the last time-step -- can apply results of Seto paper
    #a and e most important
    
    #write in such a way that can combine with 31 other similar files for same run on different
    #star sample
        
    return file
    
def get_special(info_file_name):
    
    with open(info_file_name) as f:   #assume only read (r) mode
        #create list of each non-empty line in a given file
        lines = []
        
        for line in f:
            if not line.isspace():
                lines.append(line)
        
        del tde_bodies[:]   #clear the lists
        del eject_bodies[:]
        
        for line in lines:
            if("collided with the central body at" in line):
                tde_bodies.append(line[0:4])
            elif("ejected at" in line):
                eject_bodies.append(line[0:4])
        
        print str(len(tde_bodies)) + " TDE stars: " + ', '.join(tde_bodies)
        print str(len(eject_bodies)) + " Ejections: " + ', '.join(eject_bodies)
        
def plot_TDE(info_file_name):
    file_name = ""
    
    for star in tde_bodies:
        file_name = star.strip() + ".aei"
        #rp.plot_elements(file_name, "JUPITER.aei", 10000000)
    
def plot_special(info_file_name):
    
    #get_special(info_file_name)
    
    file_name = ""
    
    for star in tde_bodies:
        file_name = star.strip() + ".aei"
        #rp.plot_elements(file_name, "JUPITER.aei", 10000000)
        
    for star in eject_bodies:
        file_name = star.strip() + ".aei"
        #rp.plot_elements(file_name, "JUPITER.aei", 10000000)
        
#returns a list of file_names        
def get_tde_file_names():
    
    names = []
    for star in tde_bodies:
        names.append(star.strip() + ".aei")
    return names
    
def get_eject_file_names():
    
    names = []
    for star in eject_bodies:
        names.append(star.strip() + ".aei")
    return names
    
    
        
    
        
