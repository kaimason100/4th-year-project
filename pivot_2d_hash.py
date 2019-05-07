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
    global neighbours, rejections, epsilon, overlap, table
    rejections = 0
    neighbours = 0
    epsilon = 2.0
    global table
    length = len(X)
    
    #below are the 47 transformation matrices which can be applied
    matrices = np.array([[[0,-1],[1,0]],[[-1,0],[0,-1]],
                        [[0,1],[-1,0]], [[1,0],[0,-1]], [[-1,0], [0,1]], [[0,1], [1,0]], [[0,-1], [-1,0]]])
    global X2,Y2
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
        
        overlap = False
        table = {}
        for i in range(length):
            table[X2[i],Y2[i]] = True
        if len(table) < length:
            overlap = True
            

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
                    
            neighbours = neighbours -2*length + 2
            neighbours = neighbours/2
            
                    
            if neighbours < old_neighbours:
                acc = np.random.rand()
                if acc > np.exp((epsilon*(neighbours-old_neighbours))):
                    overlap = True
                    neighbours = old_neighbours

                
        if overlap:
            table = dict(table2)
            X2 = np.copy(X3)
            Y2 = np.copy(Y3)
            rejections = rejections + 1
        algs = algs + 1
        
 
temp_length = 100
zeros = np.zeros(temp_length)
line = np.arange(temp_length)
X0 = np.arange(temp_length)
Y0 = np.zeros(temp_length)


N = temp_length

# =============================================================================
# #Data collection of radius of gyration
# SAW(250000, X0,Y0)
# print(rejections)
# 
# 
# for count1 in range(5):
#     rg = []
#     for count2 in range(1000):
#        SAW(10, X2, Y2)
#        r = np.sqrt(X2[temp_length-1]**2 + Y2[temp_length-1]**2)
#        rg.append(r)
#     print(np.mean(rg))
#     
# =============================================================================
 

       
#Equilibration
x_axis=[]
rsq = []
table = {}
for j in range(temp_length):
    table[X0[j],Y0[j]] = True
Rg = 0
R_cm = 0
for (i,j) in table:
    R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2))
for (l,m) in table:
    Rg = Rg + np.sqrt((1/N)*(np.sqrt(l**2+m**2) - R_cm)**2)


# =============================================================================
# for count in range(1,1000):
#     print(count)
#     SAW(count, X0,Y0)
#     x_axis.append(count)
#     if overlap == False:
#         Rg = 0
#         R_cm = 0
#         for (i,j) in table:
#             R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2))
#         for (l,m) in table:
#             Rg = Rg + np.sqrt((1/N)*(np.sqrt(l**2+m**2) - R_cm)**2)
#         rsq.append(Rg)
#     else: 
#         rsq.append(Rg)
# 
# =============================================================================
#Plotting equilibration graph
pt.figure()
pt.plot(x_axis, rsq)
pt.xlabel('Number of pivot steps')
pt.ylabel('$R_g$')
pt.title('Equilibration of a 100 Node Chain')




#Plotting radius of gyration as function of number of nodes
df = pd.read_excel('2d_pivot_data.xlsx')
log_length = df['log(length)'].tolist()
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








        
    









