'''
Created on September 21, 2015
@author: Maggie Mallernee

Script to calculate/get the resonant argument over time
'''

import resonance_plots as rp

#@return: 2D array of resonant arguments over time for each star
def get_resargs(num_stars, big_body, p, q, row):
    
    file_name = ''
    data = []
    
    big_data = rp.get_data(big_body)

    w_big = []
    M_big = []

    for x in range(len(big_data[0])):
        w_big.append(big_data[5][x])
        M_big.append(big_data[7][x])
    lambda_big = []
    for x in range(len(w_big)):
        lambda_big.append(w_big[x]+M_big[x])
    
    resargs = []

    for x in range(1, num_stars+1):
        file_name = 'S' + str(x) + '.aei'
        data = rp.get_data(file_name)
        lambda_little = []
        phi = []
        phi_mod = []
        w = data[5]
        M = data[7]
        
        for y in range(len(data[0])):
            lambda_little.append(w[y]+M[y])
            if row == 1:
                phi.append((p+q)*lambda_big[y] + (1-p-q)*lambda_little[y] - w[y] + 180)     #from Table 8.1 ROW 1, add 180 to get between -180 and 180
            elif row == 2:
                phi.append((p+q)*lambda_big[y] + (1-p-q)*lambda_little[y] - w_big[y] + 180)     #from Table 8.1 ROW 2
            elif row == 3:
                phi.append((p+q)*lambda_big[y] + (2-p-q)*lambda_little[y] - 2*w[y] + 180)
            elif row == 4:
                phi.append((p+q)*lambda_big[y] + (2-p-q)*lambda_little[y] - w_big[y] - w[y] +180)
            phi_mod.append((phi[y] % 360) - 180) 
        
        resargs.append(phi_mod)
    return resargs