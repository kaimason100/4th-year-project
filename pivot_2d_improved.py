#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 18:40:53 2018

@author: kai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 11:43:58 2018

@author: kai
"""

import numpy as np
import matplotlib.pyplot as pt
import time


def SAW(length, max_algs, X,Y):
    global table
    
    #below are the 47 transformation matrices which can be applied
    matrices = np.array([[[0,-1],[1,0]],[[-1,0],[0,-1]],
                        [[0,1],[-1,0]], [[1,0],[0,-1]], [[-1,0], [0,1]], [[0,1], [1,0]], [[0,-1], [-1,0]]])
    global X2,Y2
    #X = np.arange(length) #initial line in x-axis
   # Y = np.zeros(length) #y-axis initially zeros
    X2 = np.copy(X)
    Y2 = np.copy(Y)
        
        
    #loop for number of times pivot applied
    algs = 0
    while algs < max_algs:
        table = {}
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
                break
            k = k + 1
            
            
            
        if overlap:
            X2 = np.copy(X3)
            Y2 = np.copy(Y3)
        algs = algs + 1
        
            
        
 


t0 = time.time()
SAW(4000, 100, np.arange(4000), np.zeros(4000))
t1=time.time()

print(t1-t0)






