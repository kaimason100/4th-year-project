#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 10:15:06 2018

@author: kai
"""

import numpy as np
import matplotlib.pyplot as pt
from mpl_toolkits import mplot3d
import pandas as pd
ax = pt.axes(projection='3d') 

#Function to generate SAW
def SAW(max_algs, X,Y,Z):
    #Define length
    length = len(X)
    #Define global variables to be called outside the function
    global epsilon, neighbours, X2, Y2, Z2, overlap, table
    #define counter for interacting neighbours
    neighbours = 0
    #define epsilon value
    epsilon = 0
    #below are the 47 transformation matrices which can be applied
    matrices = np.array([[[1,0,0],[0,1,0],[0,0,-1]], [[1,0,0],[0,-1,0],[0,0,-1]],
                            [[-1,0,0],[0,-1,0],[0,0,-1]], [[-1,0,0],[0,-1,0],[0,0,1]], 
                            [[-1,0,0],[0,1,0],[0,0,1]], [[1,0,0],[0,-1,0],[0,0,1]], 
                            [[-1,0,0],[0,1,0],[0,0,-1]], [[0,1,0],[0,0,1],[1,0,0]], 
                            [[0,-1,0],[0,0,1],[1,0,0]], [[0,-1,0],[0,0,-1],[1,0,0]], 
                            [[0,-1,0],[0,0,-1],[-1,0,0]], [[0,1,0],[0,0,-1],[-1,0,0]], 
                            [[0,1,0],[0,0,1],[-1,0,0]], [[0,-1,0],[0,0,1],[-1,0,0]], 
                            [[0,1,0],[0,0,-1],[1,0,0]], [[0,0,1],[1,0,0],[0,1,0]], 
                            [[0,0,-1],[1,0,0],[0,1,0]], [[0,0,-1],[-1,0,0],[0,1,0]], 
                            [[0,0,-1],[-1,0,0],[0,-1,0]], [[0,0,1],[-1,0,0],[0,-1,0]], 
                            [[0,0,1],[1,0,0],[0,-1,0]], [[0,0,-1],[1,0,0],[0,-1,0]], 
                            [[0,0,1],[-1,0,0],[0,1,0]], [[1,0,0],[0,0,1],[0,1,0]], 
                            [[-1,0,0],[0,0,1],[0,1,0]], [[-1,0,0],[0,0,-1],[0,1,0]], 
                            [[-1,0,0],[0,0,-1],[0,-1,0]], [[1,0,0],[0,0,-1],[0,-1,0]], 
                            [[1,0,0],[0,0,1],[0,-1,0]], [[1,0,0],[0,0,-1],[0,1,0]], 
                            [[-1,0,0],[0,0,1],[0,-1,0]], [[0,0,1],[0,1,0],[1,0,0]], 
                            [[0,0,-1],[0,1,0],[1,0,0]], [[0,0,-1],[0,-1,0],[1,0,0]], 
                            [[0,0,-1],[0,-1,0],[-1,0,0]], [[0,0,1],[0,-1,0],[-1,0,0]], 
                            [[0,0,1],[0,1,0],[-1,0,0]], [[0,0,-1],[0,1,0],[-1,0,0]], 
                            [[0,0,1],[0,-1,0],[1,0,0]], [[0,1,0],[1,0,0],[0,0,1]], 
                            [[0,-1,0],[1,0,0],[0,0,1]], [[0,-1,0],[-1,0,0],[0,0,1]], 
                            [[0,-1,0],[-1,0,0],[0,0,-1]], [[0,1,0],[-1,0,0],[0,0,-1]], 
                            [[0,1,0],[1,0,0],[0,0,-1]], [[0,-1,0],[1,0,0],[0,0,-1]], 
                            [[0,1,0],[-1,0,0],[0,0,1]]])
    #Arrays to store coordinates of updated SAW
    X2 = np.copy(X)
    Y2 = np.copy(Y)
    Z2 = np.copy(Z)
        
    #loop for number of times pivot applied
    algs = 0
    table = {}
    while algs < max_algs:
        #Store value of neighbours before transformation is applied
        old_neighbours = neighbours 
        #Empty dictionary which will be used to check for self-avoidance and count interacting neighbours
        table2 = dict(table)
        #copies of coordinate arrays to revert back to if walk is rejected
        X3 = np.copy(X2)
        Y3 = np.copy(Y2)
        Z3 = np.copy(Z2)
        #random pivot point selected
        pivot = np.random.randint(1,length - 2)
        #random matrix selected 
        rand_matrix = np.random.randint(0,len(matrices)-1)
        trans_matrix = matrices[rand_matrix]
        #loop for applying transformation from pivot point to end of walk
        j = pivot + 1
        while j < length:
            [X2[j], Y2[j], Z2[j]] = trans_matrix.dot(([X2[j] - X2[pivot], Y2[j] - Y2[pivot], Z2[j] - Z2[pivot]])) + [X2[pivot], Y2[pivot], Z2[pivot]] 
            j = j + 1
            
        #check for self avoidance
        overlap = False
        table = {}
        for i in range(0, length):
            table[X2[i],Y2[i], Z2[i]] = True
        if len(table) < length:
            overlap = True
            
        #code to count number of interacting neighbours
        if overlap == False:
            neighbours = 0 
            
            for (i, j, k) in table:
                if (i+1, j, k) in table:
                    neighbours = neighbours + 1
                if (i-1, j, k) in table:
                    neighbours = neighbours + 1
                if (i, j+1, k) in table:
                    neighbours = neighbours + 1
                if (i, j-1, k) in table:
                    neighbours = neighbours + 1
                if (i, j, k+1) in table:
                    neighbours = neighbours + 1
                if (i, j, k-1) in table:
                    neighbours = neighbours + 1
                    
            neighbours = neighbours - 2*length + 2
            neighbours = neighbours/2
                    
            #metropolis monte carlo step to determine whether walk is accpeted
            if neighbours < old_neighbours:
                acc = np.random.rand()
                if acc > np.exp((epsilon*(neighbours-old_neighbours))):
                    overlap = True
                    

        #if not self avoiding then revert back to config at beginning of loop
        if overlap == True:
            neighbours = old_neighbours
            table = dict(table2)
            X2 = np.copy(X3)
            Y2 = np.copy(Y3)
            Z2 = np.copy(Z3) 
        algs = algs + 1
        
#Function to calculate radius of gyration
def Rg(x,y,z):
    N = len(x)
    R = np.zeros((N,3))
    
    R[:,0] = x
    R[:,1] = y
    R[:,2] = z
    
    R_cm = [0,0,0]
    j = 0
    while j < N:
        R_cm = R_cm + R[j]
        j = j + 1

    R_cm = R_cm/N

    R_g = 0
    i = 0
    while i < N:
        R_temp = R[i] - R_cm
        R_g = R_g + R_temp[0]**2 + R_temp[1]**2 + R_temp[2]**2
        i = i + 1
    R_g = np.sqrt(R_g/N)
    return R_g        
        
######################################################       
        
#Initialising data 
temp_length =1000
rsq = []
x = np.arange(0,500)
X0 = np.arange(temp_length)
Y0 = np.zeros(temp_length)
Z0 = np.zeros(temp_length)
# =============================================================================
# #Equilibration
# for i in range(0,500):
#     print(i)       
#     SAW(i, X0, Y0, Z0),
#     rsq.append(X2[temp_length - 1]**2 + Y2[temp_length - 1]**2)
# 
# #Plotting equilibration
# pt.figure()
# pt.title('Equilibration of a 100 Node Chain in 3-D')
# pt.xlabel('Number of pivot steps')
# pt.ylabel('$r^2$')
# pt.plot(x,rsq)
# =============================================================================

#Collecting data about radius of gyration
for j in range(10):
    rsq = []
    SAW(8000, X0, Y0, Z0)
    X4, Y4, Z4 = np.copy(X2), np.copy(Y2), np.copy(Z2)
    for i in range(0, 1000):
        SAW(5, X4, Y4, Z4)
        rsq.append(Rg(X2, Y2, Z2))
    print(np.mean(rsq))





#Plotting the random walk
pt.figure()
ax.set_xlabel('$x$', fontsize = 15)
ax.set_ylabel('$y$', fontsize = 15)
ax.set_zlabel('$z$', fontsize = 15)
ax.set_title('3-D Self Avoiding Random Walk')
ax.plot3D(X2, Y2 ,Z2)

