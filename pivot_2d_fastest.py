#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 11:43:58 2018

@author: kai
"""

import numpy as np
import matplotlib.pyplot as pt
import time
import pandas as pd

def SAW(max_algs, X,Y):
    global neighbours, rejections, epsilon, overlap, table, X2, Y2
    rejections = 0
    neighbours = 0
    epsilon = 2.0
    length = len(X)
    
    #below are the 7 transformation matrices which can be applied
    matrices = np.array([[[0,-1],[1,0]],[[-1,0],[0,-1]],
                        [[0,1],[-1,0]], [[1,0],[0,-1]], [[-1,0], [0,1]], [[0,1], [1,0]], [[0,-1], [-1,0]]])
    X2 = np.copy(X)
    Y2 = np.copy(Y)
        
        
    #loop for number of times pivot applied
    algs = 0
    table = {}
    while algs < max_algs:
        table2 = dict(table)
        X3 = np.copy(X2)
        Y3 = np.copy(Y2)
        pivot = np.random.randint(1,length - 2)
        rand_matrix = np.random.randint(0,len(matrices)-1)
        trans_matrix = matrices[rand_matrix]
        
        #loop for applying pivot to end of walk
        j = pivot + 1
        while j < length:
            [X2[j], Y2[j]] = trans_matrix.dot(([X2[j] - X2[pivot], Y2[j] - Y2[pivot]])) + [X2[pivot], Y2[pivot]]      
            j = j + 1  
        
        #Check for self avoidance
        overlap = False
        table = {}
        for i in range(length):
            table[X2[i],Y2[i]] = True
        if len(table) < length:
            overlap = True
            
        #Inclusion of attractive interaction
        if overlap == False:
            old_neighbours = neighbours      
            neighbours = 0 
            for (i, j) in table:
                if (i+1, j) in table:
                    neighbours = neighbours + 1
                if (i-1, j) in table:
                    neighbours = neighbours + 1
                if (i, j+1) in table:
                    neighbours = neighbours + 1
                if (i, j-1) in table:
                    neighbours = neighbours + 1
            #Correction for overcounting       
            neighbours = neighbours -2*length + 2
            neighbours = neighbours/2
            
            #monte carlo step using boltzmann factor       
            if neighbours < old_neighbours:
                acc = np.random.rand()
                if acc > np.exp((epsilon*(neighbours-old_neighbours))):
                    overlap = True
                    neighbours = old_neighbours

         #If self avoiding or energetically rejected then revert back to previous config       
        if overlap:
            table = dict(table2)
            X2 = np.copy(X3)
            Y2 = np.copy(Y3)
            rejections = rejections + 1
        algs = algs + 1
        
 
#Define function to calculate radius of gyration
def Rg(x,y):
    global R_cm
    N = len(x)
    R = np.zeros((N,2))
            
    R[:,0] = x
    R[:,1] = y
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
    
    

#Set initil lengths
temp_length = 100
zeros = np.zeros(temp_length)
line = np.arange(temp_length)
X0 = np.arange(temp_length)
Y0 = np.zeros(temp_length)


N = temp_length

#Equilibration

rsq = []
x_axis = []
rg_in = Rg(X0,Y0)
print(rg_in)
x_axis.append(0)
rsq.append(rg_in)
print(R_cm)
for i in range(1,1000):
    if i%50 == 0:
        print(i)
    SAW(10*i, X0, Y0)
    x_axis.append(10*i)
    if overlap ==False:
        rsq.append(Rg(X2, Y2))
    else:
        rsq.append(rsq[i-1])

print('finished')


#Plotting equilibration graph
pt.figure()
pt.plot(x_axis, rsq)
pt.xlabel('Number of pivot steps')
pt.ylabel('$R_g$')
pt.title('Equilibration of a 100 Node Chain')


#Plotting radius of gyration as function of number of nodes
df = pd.read_excel('2d_pivot_data.xlsx')
log_length = df['log(nodes)'].tolist()
log_Rg = df['log(Rg)'].tolist()

pt.figure()
pt.title('Relationship Between $N$ and $R_g$ for 2D Linear Chains')
pt.xlabel('$\log_{10}(N)$')
pt.ylabel('$\log_{10}(R_g)$')
pt.plot(log_length, log_Rg, 'x')
m, b = np.polyfit(log_length, log_Rg, 1)
print(m, b)
x_axis = np.linspace(1,3)
y_axis = m*x_axis + b
pt.plot(x_axis, y_axis)










        
    









