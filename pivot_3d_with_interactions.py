#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 10:15:06 2018

@author: kai
"""

import numpy as np
import matplotlib.pyplot as pt
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit
from numpy.polynomial.polynomial import polyfit
import time
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
        
        
 ######################################################       
        
        
        

#Initialising data 
temp_length =1000
rsq = []
x = np.arange(0,500)
X0 = np.arange(temp_length)
Y0 = np.zeros(temp_length)
Z0 = np.zeros(temp_length)
#Equilibration
'''
for i in range(0,500):
    print(i)       
    SAW(i, X0, Y0, Z0),
    rsq.append(X2[temp_length - 1]**2 + Y2[temp_length - 1]**2)

pt.figure()
pt.title('Equilibration of a 100 Node Chain in 3-D')
pt.xlabel('Number of pivot steps')
pt.ylabel('$r^2$')
pt.plot(x,rsq)



SAW(100000, X0, Y0, Z0)
N = temp_length

#Collecting data
for j in range(5):
    rsq = []
    for i in range(0, 1000):
        SAW(10, X2, Y2, Z2)
        if overlap == False:
            Rg = 0
            R_cm = 0
            for (i,j,k) in table:
                R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2 + k**2))
            for (l,m,n) in table:
                Rg = Rg + np.sqrt((1/N)*(np.sqrt(l**2+m**2+n**2) - R_cm)**2)
            rsq.append(Rg)
        else:
            rsq.append(Rg)
    
    print(np.mean(rsq))
print('length: ', temp_length)
print('neighbours: ', neighbours)

pt.figure()
ax.set_xlabel('$x$', fontsize = 15)
ax.set_ylabel('$y$', fontsize = 15)
ax.set_zlabel('$z$', fontsize = 15)
ax.set_title('3-D Self Avoiding Random Walk')
ax.plot3D(X2, Y2 ,Z2)

'''



#Linear Regression
df = pd.read_excel('3d_pivot_data.xlsx', sheetname=0) # can also index sheet by name or fetch all sheets
log10_nodes = df['log(length)'].tolist()
log10_Rg = df['log(Rg)'].tolist()


pt.figure()
m,b = np.polyfit(log10_nodes, log10_Rg, 1)
x_lin = np.linspace(0,4)
y_lin = m*x_lin + b
pt.title('Relationship Between $N$ and $R_g$ for 3D Linear Chains')
pt.xlabel('$\log_{10}(N)$')
pt.ylabel('$\log_{10}(R_g)$')
pt.plot(log10_nodes, log10_Rg, 'x')
pt.plot(x_lin, y_lin,'-')
pt.show()

print(m,b)


'''

#Plotting the end-to-end distance as a function of epsilon

df = pd.read_excel('3d_linear_chain_with_interactions_data.xlsx')

r = df['end to end distance'].tolist()
e = df['epsilon'].tolist()
logr = np.log10(r)
loge = np.log10(e)


pt.figure()
pt.title('The End-to-end Distance of a Linear Chain as a Function of $\epsilon$')
pt.xlabel('$\epsilon$', fontsize = 14)
pt.ylabel('$r^2$', fontsize = 14)
pt.plot(e, r, 'x')


pt.figure()
pt.title('The End-to-end Distance of a Linear Chain as a Function of $\epsilon$')
pt.xlabel('$\log_{10}(\epsilon)$', fontsize = 14)
pt.ylabel('$\log_{10}(r^2)$', fontsize = 14)
pt.plot(loge, logr, 'x')
m,b = np.polyfit(loge, logr, 1)
x_lin = np.linspace(-1,0.5)
y_lin1 = m*x_lin + b
pt.plot(x_lin, y_lin1)
'''