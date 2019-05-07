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
max_length = 100
x_axis = np.arange(0,max_length)


rsq_list = []


#First loop, number of random walks. Initialising arrays and counter for number of steps
while repeats < max_repeats:
    length_counter = 0
    rsq = []
    rsq2 = []
    av_rsq = []
    while length_counter < max_length:
        N = 0
        while N < num_walks:
            n=0
            X = [0]
            Y = [0] 
            
            #Second loop, number of steps per random walk
            #Append arrays with coordinates of next node in the random walk
            
            while n < length_counter:
                rand = np.random.randint(4)
                if rand == 0:
                    X.append(X[n] + 1)
                    Y.append(Y[n])
                    n = n + 1
                
                if rand == 1:
                    X.append(X[n] - 1)
                    Y.append(Y[n])
                    n = n + 1
                
                if rand == 2:
                    X.append(X[n])
                    Y.append(Y[n] + 1)
                    n = n + 1
                    
                elif rand == 3:
                    X.append(X[n])
                    Y.append(Y[n] - 1)
                    n = n + 1
                
            #Append array storing r^2 values and then calculate average of r^2, r = the distance from first to last node.
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

#print average of r^2
print(av_rsq3)












pt.figure()
pt.title('2D', fontsize = 30)
pt.xlabel('N', fontsize = 15)
pt.ylabel('r', fontsize = 15)
pt.plot(x_axis, np.sqrt(av_rsq))


pt.figure()
pt.plot(X,Y)
pt.xlabel('$x$',fontsize=15)
pt.ylabel('$y$', fontsize=15)
pt.title('2D Random Walk', fontsize=20)





