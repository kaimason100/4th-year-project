#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 10:37:15 2019

@author: kai
"""


import numpy as np
import matplotlib.pyplot as pt
import time
from mpl_toolkits import mplot3d
ax = pt.axes(projection='3d') 
import pandas as pd

#Define function to generate self avoiding random walk
def SAW(max_algs, X,Y,Z):
    global table, rejections, points,X2,Y2,Z2,overlap
    points = 6
    rejections = 0
    length = len(X[0])
    
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
    X2 = np.copy(X)
    Y2 = np.copy(Y)
    Z2 = np.copy(Z)
        
        
    #loop for number of times pivot applied
    algs = 0
    table = {}
    while algs < max_algs:
        table2 = dict(table)
        X3 = np.copy(X2)
        Y3 = np.copy(Y2)
        Z3 = np.copy(Z2)
        pivot = np.random.randint(1,length - 2)
        rand_matrix = np.random.randint(0,len(matrices)-1)
        trans_matrix = matrices[rand_matrix]
        
        #loop for applying pivot to end of walk
        j = pivot + 1
        selector = np.random.randint(0,6)
        while j < length:
            [X2[selector][j], Y2[selector][j], Z2[selector][j]] = trans_matrix.dot(([X2[selector][j] - X2[selector][pivot], Y2[selector][j] - Y2[selector][pivot], Z2[selector][j] - Z2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot], Z2[selector][pivot]] 
            j = j + 1
            
        #Check for self avoidance
        overlap = False
        table = {}
        for i in range(points):
            for j in range(length):
                table[X2[i][j],Y2[i][j],Z2[i][j]] = True    
        if len(table) < 6*length - 5:
            overlap = True
        
                
        #If not self avoiding then revert back to previous configuration
        if overlap:
            table = dict(table2)
            X2 = np.copy(X3)
            Y2 = np.copy(Y3)
            Z2 = np.copy(Z3)
            rejections = rejections + 1
        algs = algs + 1
        
#Define function to calculate radius of gyration
def Rg(x,y,z):
    global R_cm
    N = len(x)*len(x[0])
    R = np.zeros((N,3))
    x2,y2,z2 = [], [], []
    for i in range(len(x)):
        for j in range(len(x[0])):
            x2.append(x[i][j])
            y2.append(y[i][j])
            z2.append(z[i][j])
            
    R[:,0] = x2
    R[:,1] = y2
    R[:,2] = z2
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



#Set initial lengths and positions of arms        
temp_length = 30

t0 = time.time()
zeros = np.zeros(temp_length)
line = np.arange(temp_length)

X0 = [line, zeros, zeros, -line, zeros, zeros]
Y0 = [zeros, line, -line, zeros, zeros, zeros]
Z0 = [zeros, zeros, zeros, zeros, line, -line]


#plotting 3d star 
pt.figure()
ax.set_xlabel('$x$', fontsize = 15)
ax.set_ylabel('$y$', fontsize = 15)
ax.set_zlabel('$z$', fontsize = 15)
ax.set_title('A Six-armed Star Generated Using the Pivot Algorithm', fontsize = 13)
for i in range(points):
    ax.plot3D(X2[i], Y2[i] ,Z2[i])
t1=time.time()
print(t1-t0)




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

#Plotting equilibration graphs
pt.figure()
pt.title('Equilibrium of Six-Armed Star')
pt.xlabel('Number of pivots applied')   
pt.ylabel('$R_g$')     
pt.plot(x_axis,rsq)

#Acceptance rate
SAW(temp_length, 10000, X2, Y2, Z2)
print(rejections)






#Plot acceptance fraction
df = pd.read_excel('3d_star_data.xlsx', sheetname=0) # can also index sheet by name or fetch all sheets
acceptance_rate = df['acceptance_rate'].tolist()
arm_length = df['arm_length'].tolist()

pt.figure()
pt.title('')
pt.xlabel('')
pt.ylabel('')
pt.plot(np.log10(arm_length), np.log10(acceptance_rate), 'x')




