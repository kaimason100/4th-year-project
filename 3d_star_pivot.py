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
#ax = pt.axes(projection='3d') 
import pandas as pd


def SAW(max_algs, X,Y,Z):
    global table, rejections, points
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
    global X2,Y2,Z2
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
            
        
        overlap = False
        table = {}
        for i in range(points):
            for j in range(length):
                table[X2[i][j],Y2[i][j],Z2[i][j]] = True    
        if len(table) < 6*length - 5:
            overlap = True
        
                
     
        if overlap:
            table = dict(table2)
            X2 = np.copy(X3)
            Y2 = np.copy(Y3)
            Z2 = np.copy(Z3)
            rejections = rejections + 1
        algs = algs + 1
        
    
        
temp_length = 30

t0 = time.time()
zeros = np.zeros(temp_length)
line = np.arange(temp_length)

X0 = [line, zeros, zeros, -line, zeros, zeros]
Y0 = [zeros, line, -line, zeros, zeros, zeros]
Z0 = [zeros, zeros, zeros, zeros, line, -line]





#Equilibration
N = temp_length*6-5
x_axis=[]
rsq = []
for count in range(1,900):
    print(count)
    SAW(2*count, X0,Y0,Z0)
    x_axis.append(2*count)
    Rg_sq = 0
    R_cm = 0
    for (i,j,k) in table:
        R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2 + k**2))
    for (l,m,n) in table:
        Rg_sq = Rg_sq + np.sqrt((1/N)*(np.sqrt(l**2+m**2+n**2) - R_cm)**2)
    rsq.append(Rg_sq)

pt.figure()
pt.title('Equilibrium of Six-Armed Star')
pt.xlabel('Number of pivots applied')   
pt.ylabel('$R_g$')     
pt.plot(x_axis,rsq)

'''


#Acceptance rate
SAW(temp_length, 10000, X2, Y2, Z2)
print(rejections)



pt.figure()

ax.set_xlabel('$x$', fontsize = 15)
ax.set_ylabel('$y$', fontsize = 15)
ax.set_zlabel('$z$', fontsize = 15)
ax.set_title('A Six-armed Star Generated Using the Pivot Algorithm', fontsize = 13)
for i in range(points):
    ax.plot3D(X2[i], Y2[i] ,Z2[i])
t1=time.time()
print(t1-t0)




#Plot acceptance fraction
df = pd.read_excel('3d_star_data.xlsx', sheetname=0) # can also index sheet by name or fetch all sheets
acceptance_rate = df['acceptance_rate'].tolist()
arm_length = df['arm_length'].tolist()

pt.figure()
pt.title('')
pt.xlabel('')
pt.ylabel('')
pt.plot(np.log10(arm_length), np.log10(acceptance_rate), 'x')



'''

