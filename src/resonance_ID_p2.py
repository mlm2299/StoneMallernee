'''
Resonance determination algorithm part 2:
    Takes information from resonance_ID_p1 (which gives a table of possible resonances in given timetables),
    then calculates the resonant argument for each of those resonances in those time periods,
    and makes a final determination of whether it is librating or circulating.
    This should hopefully narrow results to a single resonance (or no resonance at all).

@author Maggie Mallernee
Created February 26, 2016
'''
import resonance_plots as rp
import matplotlib.pyplot as plt
import numpy as np

#algo1_table should be the table from resonance_ID_p1, containing rows of the form [p, q, t_start, t_end, line_start, line_end]
#seeks to create a [ p, q, sub_res, time[], phi[] ] 
def get_phi_table(file_name, big_file_name, algo1_table, fudge): #automate ROW
    
    all_phi_table = []
    
    time_phi = []
    for x in range(len(algo1_table)):
        #print algo1_table[x]
        p = algo1_table[x][0]
        q = algo1_table[x][1]
        sub_res = algo1_table[x][2]
        #t_start = algo1_table[x][3]
        #t_end = algo1_table[x][4]
        line_start = algo1_table[x][5]
        line_end = algo1_table[x][6]
        time_phi = get_res_arg1(file_name, big_file_name, p, q, sub_res, line_start, line_end) #gives line #s to find time
        all_phi_table.append([p, q, sub_res, time_phi[0], time_phi[1], test_for_circ(time_phi, fudge)]) #adding a row { p, q, sub_res, time[], phi[], circ test result }
        print "The maximum phi for this set is: %s" % (get_phi_max(time_phi))
        print "The minimum phi for this set is: %s" % (get_phi_min(time_phi))
        print "Phi changed sign %s times during the run" % (get_flips(time_phi))
        print "Result of circulation test is: %s" % (test_for_circ(time_phi, .1))
    
    return all_phi_table
    

#needs to be mindful of time conversions -- algo 1 should give real time values, not row numbers
#@return an array of [time, phi-mod]
def get_res_arg1(file_name, big_file_name, p, q, sub_res, line_start, line_end):      #SHOULD get row from p/q?
    
    #get necessary data about secondary black hole
    big_data = rp.get_data(big_file_name)
    
    w_big = big_data[5]
    M_big = big_data[7]

    
    lambda_big = []
    for x in range(len(w_big)):
        lambda_big.append(w_big[x]+M_big[x])
    
    #get data about the star/test particle
    data = []
    data = rp.get_data(file_name)
    lambda_little = []
    phi = []
    phi_mod = []
    time = []  
    Om = data[4]     
    w = data[5]
    M = data[7]
    
    for x in range(len(w)):
        lambda_little.append(w[x]+M[x])

    #to set row from Table 8.1 ------CORRECT THIS --- add row 6 for 2nd order
    if q == 1.:
        row = 1    
    if q == 2.: 
        if sub_res == 0.:
            row = 3    #--- also check Row 6 --- plot both for 2nd order
        elif sub_res == 1.:
            row = 6
 
    #to calculate resonant argument over time
    for y in range(line_end - line_start + 1):  #gets the proper number of lines to be considered
        t = ((int) (line_start + y))            #shifts to get the right interval
        time.append(data[0][t]) #will have indices 0 - y_max
        if row == 1:
            phi.append((p+q)*lambda_big[t] + (1-p-q)*lambda_little[t] - w[t] + 180)     #from Table 8.1 ROW 1, add 180 to get between -180 and 180
        elif row == 2:
            phi.append((p+q)*lambda_big[t] + (1-p-q)*lambda_little[t] - w_big[t] + 180)     #from Table 8.1 ROW 2
        elif row == 3:
            phi.append((p+q)*lambda_big[t] + (2-p-q)*lambda_little[t] - 2*w[t] + 180)
        elif row == 4:
            phi.append((p+q)*lambda_big[t] + (2-p-q)*lambda_little[t] - w_big[t] - w[t] +180)
        elif row == 6:
            phi.append((p+q)*lambda_big[t] + (2-p-q)*lambda_little[t] - 2*Om[t] + 180)  #ADDED -- VERIFY
        phi_mod.append((phi[y] % 360) - 180)   #to go to correct index (y not t)          #phi mod 360 degrees, minus 180 to get correct frame
        #phi_mod.append(2*(lambda_big[t]%360) - (lambda_little[t]%360) - (w[t]%360)) 
    
    return [time, phi_mod]
   
def plot_all_phi(all_phi_table): 
    
    for x in range(len(all_phi_table)):
        plt.figure(x)
        p = all_phi_table[x][0]
        q = all_phi_table[x][1]
        sub_res = all_phi_table[x][2]
        time = all_phi_table[x][3]
        phi_mod = all_phi_table[x][4]
        
        plt.title("Phi over Time for %s : %s Resonance, Subresonance %s" % (p+q, p, sub_res))
        plt.xlabel("t (years)")
        plt.ylabel("phi mod 360 (degrees)")
        plt.ylim(-180., 180.)
        plt.plot(time, phi_mod, ':o', c = 'blue')     
    
        plt.show()


###### Take in phi vs time and do meat of algo pt 2 ######

#Calculate maximum phi in the table
def get_phi_max(phi_v_time):
    p_max = phi_v_time[1][0];
    
    for x in range(len(phi_v_time[0])): #in the range of years
        if(phi_v_time[1][x] > p_max):
            p_max = phi_v_time[1][x]
            
    return p_max
    
#Calculate the minimum phi in the table
def get_phi_min(phi_v_time):
    p_min = phi_v_time[1][0];
    
    for x in range(len(phi_v_time[0])): #in the range of years
        if(phi_v_time[1][x] < p_min):
            p_min = phi_v_time[1][x]
            
    return p_min
    
#Calculate the number of sign flips (includes crossing 0 and crossing the -180/+180 boundary)
#should ignore 0
def get_flips(phi_v_time):
    #ensure that non initialized to 0
    i = 0
    while(np.sign(phi_v_time[1][i]) == 0):
        i += 1
    curr_sign = np.sign(phi_v_time[1][i]);
    num_flips = 0
    for x in range(len(phi_v_time[0]) - i): #in the range of years
        x += i
        sign = np.sign(phi_v_time[1][x])
        if(sign != curr_sign and sign != 0):
            num_flips += 1
            curr_sign = sign
            
    return num_flips
    
#get our theoretical delta phi
def get_ideal_d_phi(phi_v_time):
    num_pts = len(phi_v_time[0])
    flips = get_flips(phi_v_time)
    
    #return 180. / (num_pts / flips)
    return 180. * (float(flips) / float(num_pts))
    
#do the actual test 
#return 1 if circulating
#return 0 otherwise
def test_for_circ(phi_v_time, fudge):
    bound = fudge * (get_ideal_d_phi(phi_v_time))
    min_p = abs(get_phi_min(phi_v_time))
    max_p = abs(get_phi_max(phi_v_time))
    
    print "Bound: %s" % bound
    
    if(min_p > max_p):
        real_max = min_p 
    else: 
        real_max = max_p
        
    if(180. - real_max < bound):
        return 1
    else:
        return 0
    

