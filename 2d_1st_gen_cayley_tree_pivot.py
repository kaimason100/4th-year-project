#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 16:10:13 2019

@author: kai
"""

 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 17:15:52 2019

@author: kai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as pt
import time
import pandas as pd


def SAW(max_algs, X,Y, X_branches, Y_branches):
    global neighbours,rejections, table, epsilon, points, branches 
    length = temp_length
    branch_length = b_length
    neighbours = 0
    epsilon = 0.9
    points = 4
    branches = 8
    
    #below are the 7 transformation matrices which can be applied
    matrices = np.array([[[0,-1],[1,0]],[[-1,0],[0,-1]],
                        [[0,1],[-1,0]], [[1,0],[0,-1]], [[-1,0], [0,1]], [[0,1], [1,0]], [[0,-1], [-1,0]]])
    global X2, Y2, X2_branches, Y2_branches
    X2, X2_branches = np.copy(X), np.copy(X_branches)
    Y2, Y2_branches = np.copy(Y), np.copy(Y_branches)             
    
    #loop for number of times pivot applied
    rejections = 0
    table = {}
    for i in range(max_algs):     
        table2 = table
        X3, X3_branches = np.copy(X2), np.copy(X2_branches)
        Y3, Y3_branches = np.copy(Y2), np.copy(Y2_branches)
        rand_matrix = np.random.randint(0,len(matrices)-1)
        trans_matrix = matrices[rand_matrix]
        
        
        #loop for applying pivot to end of walk
        arm_or_branch = np.random.randint(0,2)
        if arm_or_branch == 0:
            selector = np.random.randint(0,points)
            pivot = np.random.randint(1,length - 2)
            j = pivot + 1
            k = 0
            while j < length:
                [X2[selector][j], Y2[selector][j]] = trans_matrix.dot(([X2[selector][j] - X2[selector][pivot], Y2[selector][j] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]   
                j = j+1 
            while k < branch_length:
                [X2_branches[2*selector][k], Y2_branches[2*selector][k]] = trans_matrix.dot(([X2_branches[2*selector][k] - X2[selector][pivot], Y2_branches[2*selector][k] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]
                [X2_branches[2*selector+1][k], Y2_branches[2*selector+1][k]] = trans_matrix.dot(([X2_branches[2*selector+1][k] - X2[selector][pivot], Y2_branches[2*selector+1][k] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]
                k = k+1
            
        if arm_or_branch == 1:
            selector = np.random.randint(0, branches)
            pivot  = np.random.randint(1, branch_length - 2)
            j = pivot + 1
            while j < branch_length:
                [X2_branches[selector][j], Y2_branches[selector][j]] = trans_matrix.dot(([X2_branches[selector][j] - X2_branches[selector][pivot], Y2_branches[selector][j] - Y2_branches[selector][pivot]])) + [X2_branches[selector][pivot], Y2_branches[selector][pivot]]
                j = j+1
        table = {}
        overlap = False
        for i in range(points):
            for j in range(length):
                table[X2[i][j],Y2[i][j]] = True           
        for i in range(branches):
            for j in range(branch_length):
                table[X2_branches[i][j], Y2_branches[i][j]] = True              
        if len(table) < points*length + branches*branch_length - 11:
            overlap = True
            table = table2
 
        '''if overlap == False:
            old_neighbours = neighbours      
            neighbours = 0 
            for (i, j) in table2:
                if (i+1, j) in table2:
                    neighbours = neighbours + 1
                if (i-1, j) in table2:
                    neighbours = neighbours + 1
                if (i, j+1) in table2:
                    neighbours = neighbours + 1
                if (i, j-1) in table2:
                    neighbours = neighbours + 1
                    
            #neighbours = neighbours -2*points*end_length + 2*points
            #neighbours = neighbours/2
            
                    
            if neighbours < old_neighbours:
                acc = np.random.rand()
                if acc > np.exp((epsilon*(neighbours-old_neighbours))):
                    overlap = True
                    neighbours = old_neighbours'''
                    
        
                
        if overlap:
            X2, X2_branches = np.copy(X3), np.copy(X3_branches)
            Y2, Y2_branches = np.copy(Y3), np.copy(Y3_branches)
            rejections = rejections + 1        
        
        
        
        
#Defining function to calculate radius of gyration 
def Rg(x,y,x_b,y_b):
    global R_cm
    N = len(x)*len(x[0]) + len(x_b)*len(x_b[0])
    R = np.zeros((N,2))
    x2,y2 = [], []
    for i in range(len(x)):
        for j in range(len(x[0])):
            x2.append(x[i][j])
            y2.append(y[i][j])
            
            
    for k in range(len(x_b)):
        for l in range(len(x_b[0])):
            x2.append(x_b[k][l])
            y2.append(y_b[k][l])
            
            
    R[:,0] = x2
    R[:,1] = y2
    R_cm = [0,0]
    j = 0
    while j < N:
        R_cm = R_cm + R[j]
        j = j + 1

    R_cm = R_cm/N

    R_g = 0
    i = 0
    while i < N:
        R_temp = R[i] - R_cm
        R_g = R_g + R_temp[0]**2 + R_temp[1]**2
        i = i + 1
    R_g = np.sqrt(R_g/N)
    return R_g        
        

#Defining initial lengths of arms and branches
temp_length = 100
b_length = int(temp_length/5)

t0 = time.time()
zeros = np.zeros(temp_length)
line = np.arange(temp_length)
branch = np.arange(b_length)
b_loc1 = [temp_length-1]*b_length
b_loc2 = [-temp_length+1]*b_length

x_branches = [branch, -branch, b_loc1, b_loc1, branch, -branch, b_loc2, b_loc2]
y_branches = [b_loc1, b_loc1, branch, -branch, b_loc2, b_loc2, branch, -branch]


X0 = [zeros, line, zeros, -line] 
Y0 = [line, zeros, -line, zeros] 



'''SAW(3000, X0, Y0, x_branches, y_branches)



#Plotting initial configuration
pt.figure()
pt.title('1st Generation Cayley Tree', fontsize=15)
pt.xlabel('$x$',fontsize=11)
pt.ylabel('$y$', fontsize=11)
for i in range(points):
    pt.plot(X0[i],Y0[i])
    
for j in range(branches):
    pt.plot(x_branches[j], y_branches[j])
    
        

#Plotting chain after pivot steps
pt.figure()
pt.title('1st Generation Cayley Tree', fontsize=15)
pt.xlabel('$x$',fontsize=11)
pt.ylabel('$y$', fontsize=11)
for i in range(points):
    pt.plot(X2[i],Y2[i])
    
for j in range(branches):
    pt.plot(X2_branches[j], Y2_branches[j])
    '''
       
#Equilibration
overlap = False
rsq = []
x_axis = []
rg_in = Rg(X0,Y0, x_branches, y_branches)
print(rg_in)
x_axis.append(0)
rsq.append(rg_in)
print(R_cm)
for i in range(1,1000):
    if i%50 == 0:
        print(i)
    SAW(6*i, X0, Y0, x_branches, y_branches)
    x_axis.append(6*i)
    if overlap ==False:
        rsq.append(Rg(X2, Y2, X2_branches, Y2_branches))
    else:
        rsq.append(rsq[i-1])
        
print('finished')

#Plotting equilibration 
pt.figure()
pt.title('Equilibrium of 1st Generation Cayley Tree')
pt.xlabel('Number of pivots applied')   
pt.ylabel('$R_g$')     
pt.plot(x_axis,rsq)




