#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 11:17:35 2019

@author: kai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 12:25:12 2019

@author: kai
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 11:13:19 2019

@author: kai
"""

#Import packages 
import numpy as np
import matplotlib.pyplot as pt
import time
import pandas as pd
import xlrd, xlwt
from xlwt import Workbook
from xlutils.copy import copy as xl_copy
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import pool



num_cores = multiprocessing.cpu_count()
#Define function for SAW pivot algorithm

def SAW(max_algs, X,Y):
    #Define global variables to be called outside the function
    global neighbours, rejections, alt_rej, orig_rej, points, epsilon, table, table2, table3, type_of_step, overlap, y_selector, matrices, p
    #Initial nearest neighbours
    neighbours = 0
    #Strength of attractive interaction
    epsilon = 1.0
    #Prob. of original pivot step 
    p = 0.8
    
    #Define number of arms of star and length of each arm
    points = len(X)
    length = len(X[0])
    
    #below are the 7 transformation matrices which can be applied
    matrices = np.array([[[0,-1],[1,0]],[[-1,0],[0,-1]],
                        [[0,1],[-1,0]], [[1,0],[0,-1]], [[-1,0], [0,1]], [[0,1], [1,0]], [[0,-1], [-1,0]]])
    #Copy and store arrays 
    global X2,Y2, X3, Y3
    X2, Y2 = np.copy(X), np.copy(Y)
        
    algs = 0
    #Number of total, alternative and original pivot steps rejected initialised
    rejections = 0
    alt_rej = 0
    orig_rej = 0
    #Dictionary to store coordinates
    table = {}
    #Loop for pivot algorithm
    while algs < max_algs:
        type_of_step = 'none'
        #Variables used to count rejected steps
        orig, alt = False, False
        #Copy of table
        table2 = dict(table)
        #Copy of arrays before transformation
        X3, Y3 = np.copy(X2), np.copy(Y2)
        #Select random pivot point
        pivot = np.random.randint(1,length - 2)
        #Selectr random transformation matrix
        rand_matrix = np.random.randint(0,len(matrices)-1)
        trans_matrix = matrices[rand_matrix]
        #Random number used to determine type of pivot step
        trans_selector = np.random.random(1)
        #loop for applying pivot to end of walk
        j = pivot + 1
        #select arm of star
        arm_selector = np.random.randint(0,points)
        
        #Original pivot step 
        if trans_selector < p:
            while j < length:
                [X2[arm_selector][j], Y2[arm_selector][j]] = trans_matrix.dot(([X2[arm_selector][j] - X2[arm_selector][pivot], Y2[arm_selector][j] - Y2[arm_selector][pivot]])) + [X2[arm_selector][pivot], Y2[arm_selector][pivot]]      
                j = j + 1
                orig = True
                type_of_step = 'original'
        #Alternative pivot steps
        else:
            alt = True
            if np.random.random(1) <= 0.9:
                type_of_step = 'inversion'
                point1 = np.random.randint(1, int(length))
                point2 = np.random.randint(1, int(length))
                if point2 < point1:
                    temp = point1
                    point1 = point2
                    point2 = temp
                for i in range(point1, point2):
                    j = point1 + point2 - i 
                    [X2[arm_selector][i], Y2[arm_selector][i]] = [X3[arm_selector][point1] + X3[arm_selector][point2] - X3[arm_selector][j], Y3[arm_selector][point1] + Y3
                    [arm_selector][point2] - Y3[arm_selector][j]]
        
            else:           
                type_of_step = 'alt reflection'
                y_selector = np.random.randint(-temp_length,temp_length)
                for i in range(points):    
                    temp_pivots = []
                    for j in range(length):
                        if Y2[i][j] == y_selector:
                            temp_pivots.append(j)
                    if len(temp_pivots) != 0:
                        a = np.random.randint(0,len(temp_pivots))
                        pivot = temp_pivots[a]
                        k = pivot + 1
                        while k < length:
                            [X2[i][k], Y2[i][k]] = matrices[3].dot(([X2[i][k] - X2[i][pivot], Y2[i][j] - Y2[i][pivot]])) + [X2[i][pivot], Y2[i][pivot]]   
                            k = k + 1   
            
                        
                    
            
            
        
        #Check for self avoidance
        overlap = False
        table = {}
        #Storing coordinates in dictionary
        for i in range(points):
            for j in range(length):
                table[X2[i][j],Y2[i][j]] = True
        if len(table) < points*length:
            overlap = True
          
        #Including attractive interactions    
        if overlap == False:
            #Dictionary for attractively interacting nodes
            table3 = {}
            for i in range(points):
                for j in range(int(0.8*length), length):
                    table3[X2[i][j], Y2[i][j]] = True
            old_neighbours = neighbours      
            neighbours = 0 
            for (i, j) in table3:
                if (i+1, j) in table3:
                    neighbours = neighbours + 1
                if (i-1, j) in table3:
                    neighbours = neighbours + 1
                if (i, j+1) in table3:
                    neighbours = neighbours + 1
                if (i, j-1) in table3:
                    neighbours = neighbours + 1    
            neighbours = neighbours/2
            neighbours = neighbours -points*0.2*length + points
            
            #Monte carlo step to determine whether move accepted      
            if neighbours < old_neighbours:
                acc = np.random.rand()
                if acc > np.exp((epsilon*(neighbours-old_neighbours))):
                    overlap = True
                    neighbours = old_neighbours
        
        #If move rejected then revert back to previous config.
        if overlap == True:
            table = dict(table2)
            X2 = np.copy(X3)
            Y2 = np.copy(Y3)
            rejections = rejections + 1
            if alt:
                alt_rej = alt_rej + 1
            if orig:
                orig_rej = orig_rej +1
                
        algs = algs + 1

#Define funtion to caluclate radius of gyration
def Rg(x,y):
    global R_cm
    N = len(x)*len(x[0])
    R = np.zeros((N,2))
    x2,y2 = [], []
    for i in range(len(x)):
        for j in range(len(x[0])):
            x2.append(x[i][j])
            y2.append(y[i][j])
            
    R[:,0] = x2
    R[:,1] = y2
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





#Defining length of arm        
temp_length = 100

#Initialising arms
offset=1
box = 2
zeros = np.zeros(temp_length)
line = np.arange(box,temp_length+box)
X0 = [line,[offset]*temp_length, -line, line, -line, [-offset]*temp_length, [offset]*temp_length, [-offset]*temp_length]
Y0 = [[offset]*temp_length, line, [offset]*temp_length, [-offset]*temp_length, [-offset]*temp_length, line, -line, -line]


#Equilibration
x_axis=[]
rsq = []
points = 8
N = 8*temp_length

for count in range(1,1000):
    print(count)
    SAW(10*count, X0,Y0)
    x_axis.append(10*count)
    rsq.append(Rg(X2,Y2))
    
#Plotting equilibration 
pt.figure()
pt.title('Equilibrium of Eight-Armed Star')
pt.xlabel('Number of pivots applied')   
pt.ylabel('$R_g$')     
pt.plot(x_axis,rsq)



#Data collection of radius of gyration
N = 8*temp_length
Rg = 0
R_cm = 0


for count in range(5):
    rsq = []
    SAW(300000, X0,Y0)
    for i in range(100):
        SAW(10, X2 ,Y2)
        rsq.append(Rg(X2,Y2))
    print(np.mean(rsq))
print('p: ', p)

############################
#Data collection of arm lengths 

#create workbook
wb = Workbook()
sheet1 = wb.add_sheet('e=0.4') 

#Add sheets to existing workbook
wb = xlrd.open_workbook('8_armed_star_interactions_data.xls', formatting_info=True)
wb2 = xl_copy(wb)
sheet2 = wb2.add_sheet('e=1.0')

#Equilibrate
SAW(250000, X0, Y0)
#collect data
for i in range(1000):
    print(i)
    SAW(400, X2, Y2)
    for n in range(points):
        r_n = np.sqrt(X2[n][temp_length-1]**2 + Y2[n][temp_length-1]**2)
        #sheet1.write(points*i+n,3, r_n)
        sheet2.write(points*i+n,3, r_n)

#save to Excel spreadsheet
wb2.save('8_armed_star_interactions_data.xls') 


##############################
#Plotting eight armed star
pt.figure()
pt.title('An Eight-armed Star Generated Using the Pivot Algorithm')
pt.xlabel('$x$')
pt.ylabel('$y$')
for i in range(points):
    pt.plot(X2[i],Y2[i])
for i in range(points):
    pt.plot(X2[i],Y2[i])
pt.plot([-box,box,box,-box,-box], [box,box,-box,-box,box])

##################################     
#Acceptance Fraction data collection 
SAW(400000, X0, Y0)

for j in range(5):
    rej = []
    alternative_rejections = []
    original_rejections = []
    for i in range(10):
        SAW(1000, X2, Y2)
        rej.append(rejections)
        alternative_rejections.append(alt_rej)
        original_rejections.append(orig_rej)
    print('total rejections: ',np.mean(rej))
    print('alternative pivot step rejections: ', np.mean(alternative_rejections))
    print('original pivot step rejections: ', np.mean(original_rejections))


####################################