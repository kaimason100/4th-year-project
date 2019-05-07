#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 20:59:42 2018

@author: kai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 13:30:39 2018

@author: kai
"""
#Import relevant packages
import numpy as np
import matplotlib.pyplot as pt


#Initialise variables and counter for number of walks
repeats = 0 
max_repeats = 5
num_walks = 5
av_rsq2 = []
max_length = 10
x_axis = np.arange(0,max_length)

#First loop, number of random walks. Initialising arrays and counter for number of repeats
while repeats < max_repeats:
    length_counter = 0
    rsq = []
    rsq2 = []
    av_rsq = []
    #Second loop, counter of length of walk
    while length_counter < max_length:
        N = 0
        #third loop, number of steps in walk, this loop generates one walk
        while N < num_walks:
            n = 0
            X = [0]
            Y = [0] 
            Z = [0]
            #fourth loop, grows the walk from nothing 
            while n < length_counter: 
                rand = np.random.randint(6)
                if rand == 0:
                    X.append(X[n] + 1)
                    Y.append(Y[n])
                    Z.append(Z[n])
                    n = n + 1
                
                if rand == 1:
                    X.append(X[n] -1)
                    Y.append(Y[n])
                    Z.append(Z[n])
                    n = n + 1
                
                if rand == 2:
                    X.append(X[n])
                    Y.append(Y[n] + 1)
                    Z.append(Z[n])
                    n = n + 1
                    
                if rand == 3:
                    X.append(X[n])
                    Y.append(Y[n] - 1)
                    Z.append(Z[n])
                    n = n + 1
                    
                if rand == 4:
                    X.append(X[n])
                    Y.append(Y[n])
                    Z.append(Z[n] + 1)
                    n = n + 1
                    
                if rand == 5:
                    X.append(X[n])
                    Y.append(Y[n])
                    Z.append(Z[n] - 1)
                    n = n + 1
                    
          
            #Append array storing r^2 values. 
            rsq.append(X[length_counter]**2 + Y[length_counter]**2)  
            N = N + 1
        rsq2.append(rsq[num_walks*length_counter:num_walks*length_counter+num_walks])        
        length_counter = length_counter + 1
    for i in range(0,len(rsq2)):
        av_rsq.append(np.mean(rsq2[i]))
    av_rsq2.append(av_rsq)
    repeats = repeats + 1


av_rsq3 = np.zeros((len(av_rsq2[0]), len(av_rsq2)))
for j in range(0,len(av_rsq2)):
    for i in range(0,len(av_rsq2[0])):
        av_rsq3[i][j] = av_rsq2[j][i]


    




