#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 15:36:42 2018

@author: kai
"""
import pandas as pd
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as pt
ax = pt.axes(projection='3d')       



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

max_algs = 100 #number of times pivot algorithm applied 
length = 10 #length of initial line
x_axis = np.arange(0,length) #empty array to be the x axis
av_rsq2 = []
max_repeats = 5
repeats = 0
r = [] #Empty array which will store r values
X = np.arange(0,length) #initial line in x-axis
Y = np.zeros(length) #y-axis initially zeros
Z = np.zeros(length) #x-axis initially zeros
X2 = np.copy(X)
Y2 = np.copy(Y)
Z2 = np.copy(Z)
#loop for number of times pivot applied
algs = 0
while algs < max_algs:
    X3 = np.copy(X2)
    Y3 = np.copy(Y2)
    Z3 = np.copy(Z2)
    pivot = np.random.randint(1,length - 2)
    rand_matrix = np.random.randint(0,len(matrices)-1)
    trans_matrix = matrices[rand_matrix]
    j = pivot + 1
    #loop for applying pivot to end of walk
    while j < length -1:
        [X2[j], Y2[j], Z2[j]] = trans_matrix @ ([X3[j] - X3[pivot], Y3[j] - Y3[pivot], Z3[j] - Z3[pivot]])
        [X2[j], Y2[j], Z2[j]] = [X2[j] + X3[pivot], Y2[j] + Y3[pivot], Z2[j] + Z3[pivot]]
        j = j + 1
              
        #check for self avoidance
        k = 0 
        overlap = False
        while k < pivot:
            l = pivot + 1
            while l < length-1:
                if X2[k] == X2[l] and Y2[k] == Y2[l] and Z2[k] == Z2[l]:
                    overlap = True
                    break
                l = l + 1
            if overlap:
                break
            k = k + 1
            
        #if not self avoiding then revert back to config at beginning of loop
        if overlap:  
            X2 = np.copy(X3)
            Y2 = np.copy(Y3)
            Z2 = np.copy(Z3)
    algs = algs + 1
        
'''model = LinearRegression()

x_axis = x_axis.reshape(-1,1)
model.fit(x_axis, av_rsq3)'''

'''pt.figure()
#pt.xlim(0,1.3)
pt.plot(x_axis,av_rsq3, 'rx')
pt.plot(x_axis, model.predict(x_axis))
'''








    
#plot figures
'''pt.figure()
ax.set_xlabel('$x$', fontsize = 15)
ax.set_ylabel('$y$', fontsize = 15)
ax.set_zlabel('$z$', fontsize = 15)
ax.set_title('3D Self Avoiding Random Walk Using the Pivot Algorithm', fontsize = 15)
ax.plot3D(X2, Y2 ,Z2)

pt.figure()
pt.ylabel('$r$', fontsize = 15)
pt.xlabel('Number of Pivot Algorithm Repeats')
pt.plot(x_axis,r)'''


