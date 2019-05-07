#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 17:04:58 2018

@author: kai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 17:18:49 2018

@author: kai
"""

import numpy as np
import matplotlib.pyplot as pt
from scipy.spatial.distance import cdist



def SAW(length, max_algs, X,Y):
    global epsilon, neighbours
    neighbours = 0
    epsilon = 0.269
    #below are the 47 transformation matrices which can be applied
    matrices = np.array([[[0,-1],[1,0]],[[-1,0],[0,-1]],
                        [[0,1],[-1,0]], [[1,0],[0,-1]], [[-1,0], [0,1]], [[0,1], [1,0]], [[0,-1], [-1,0]]])
    global X2,Y2
    #X = np.arange(length) #initial line in x-axis
   # Y = np.zeros(length) #y-axis initially zeros
   # Z = np.zeros(length) #x-axis initially zeros
    X2 = np.copy(X)
    Y2 = np.copy(Y)
        
    #loop for number of times pivot applied
    algs = 0
    while algs < max_algs:
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
        #check for self avoidance
        k = 0 
        overlap = False
        while k < pivot:
            l = pivot 
            while l < length:
                if X2[k] == X2[l] and Y2[k] == Y2[l]:
                    overlap = True
                    break
                l = l + 1 
            if overlap:
                #print('overlap')
                break
            k = k + 1                          
            #if not self avoiding then revert back to config at beginning of loop
     
       
        if overlap == False:
            old_neighbours = neighbours      
            neighbours = 0 
            for i in range(0,length):
                for j in range(i+2,length):
                    if (X2[i] - X2[j])**2 + (Y2[i] - Y2[j])**2 == 1:
                        neighbours = neighbours + 1
            if neighbours < old_neighbours:
                acc = np.random.rand()
                if acc > np.exp((epsilon*(neighbours-old_neighbours))):
                    overlap = True
                    #print('rejection!!!!')
                    neighbours = old_neighbours
                    
        if overlap:
            X2 = np.copy(X3)
            Y2 = np.copy(Y3)
        algs = algs + 1
 

#Data collection

      
temp_length = 30
rsq = []
#x = np.arange(0,300)
X0 = np.arange(temp_length)
Y0 = np.zeros(temp_length)
'''for i in range(0,200):
    print(i)       
    SAW(temp_length, i, X0, Y0),k
    rsq.append(X2[temp_length - 1]**2 + Y2[temp_length - 1]**2)

pt.plot(x,rsq)
'''


SAW(temp_length, 300, X0, Y0)


for i in range(0, 10000):
    print(i)
    SAW(temp_length, 5, X2, Y2)
    rsq.append((X2[temp_length-1]**2 + Y2[temp_length-1]**2))
    
print(np.mean(rsq))
pt.plot(X2, Y2)


#plot of fraction of pivot attempts that are successful vs N in 2d and 3d, expecting power law?
#suggests log plot, Madras and Sokal, original pivot algorithm paper, READ!! journal of statistical physics

#expect as epsilon increases, walk become more collapsed
#is there a value of epsilon so flory exponent = 1/2???
#Papers by grassberger, READ!!!, see notes
