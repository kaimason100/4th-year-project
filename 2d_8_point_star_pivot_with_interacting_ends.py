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


#Defining length of arm        
temp_length = 100

#Initialising arms
offset=1
box = 2
zeros = np.zeros(temp_length)
line = np.arange(box,temp_length+box)
X0 = [line,[offset]*temp_length, -line, line, -line, [-offset]*temp_length, [offset]*temp_length, [-offset]*temp_length]
Y0 = [[offset]*temp_length, line, [offset]*temp_length, [-offset]*temp_length, [-offset]*temp_length, line, -line, -line]
'''

#Equilibration
t0 = time.time()
x_axis=[]
rsq = []
points = 8
N = 8*temp_length
table = {}
for i in range(points):
    for j in range(temp_length):
        table[X0[i][j],Y0[i][j]] = True
Rg = 0
R_cm = 0
for (i,j) in table:
    R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2))
for (l,m) in table:
    Rg = Rg + np.sqrt((1/N)*(np.sqrt(l**2+m**2) - R_cm)**2)



for count in range(1,5):
    print(count)
    SAW(100*count, X0,Y0)
    x_axis.append(100*count)
    if overlap == False:
        Rg = 0
        R_cm = 0
        for (i,j) in table:
            R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2))
        for (l,m) in table:
            Rg = Rg + np.sqrt((1/N)*(np.sqrt(l**2+m**2) - R_cm)**2)
        rsq.append(Rg)
    else: 
        rsq.append(Rg)
    
# =============================================================================
#     if Rg >=400 and count > 15:
#     #if type_of_step == 'alt reflection' and overlap == False:
#         print('Radius of gyration after: ', Rg)
#         print('Type of step: ', type_of_step)
#         print('Overlap?: ', overlap)
#         pt.figure()
#         pt.title('Before')
#         pt.xlabel('$x$')
#         pt.ylabel('$y$')
#         for i in range(points):
#             pt.plot(X3[i],Y3[i])
#         for i in range(points):
#             pt.plot(X3[i],Y3[i])
#         pt.plot([-box,box,box,-box,-box], [box,box,-box,-box,box])
#         pt.figure()
#         pt.title('After')
#         pt.xlabel('$x$')
#         pt.ylabel('$y$')
#         for i in range(points):
#             pt.plot(X2[i],Y2[i])
#         for i in range(points):
#             pt.plot(X2[i],Y2[i])
#         pt.plot([-box,box,box,-box,-box], [box,box,-box,-box,box])
# =============================================================================
    
t1 = time.time()
print(t1-t0)    
pt.figure()
pt.title('Equilibrium of Eight-Armed Star')
pt.xlabel('Number of pivots applied')   
pt.ylabel('$R_g$')     
pt.plot(x_axis,rsq)
print(rsq)
print('p', p)
print('epsilon: ', epsilon)



#Data collection of radius of gyration
SAW(300000, X0,Y0)
N = 8*temp_length
Rg = 0
R_cm = 0
for (i,j) in table:
    R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2))
for (l,m) in table:
    Rg = Rg + np.sqrt((1/N)*(np.sqrt(l**2+m**2) - R_cm)**2)

for count in range(5):
    rsq = []
    for i in range(100):
        SAW(10, X2 ,Y2)
        if overlap == False:
            Rg = 0
            R_cm = 0
            for (i,j) in table:
                R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2))
            for (l,m) in table:
                Rg = Rg + np.sqrt((1/N)*(np.sqrt(l**2+m**2) - R_cm)**2)
            rsq.append(Rg)
        else:
            rsq.append(Rg)
    print(np.mean(rsq))
print('p: ', p)
# =============================================================================
#     
#     if Rg >=400 and overlap == False and count != 0 and Rg != 0:
#     #if type_of_step == 'alt reflection' and overlap == False:
#         print('Radius of gyration after: ', Rg)
#         print('Type of step: ', type_of_step)
#         print('Overlap?: ', overlap)
#         pt.figure()
#         pt.title('Before')
#         for i in range(points):
#             pt.plot(X3[i],Y3[i])
#         for i in range(points):
#             pt.plot(X3[i],Y3[i])
#         pt.plot([-box,box,box,-box,-box], [box,box,-box,-box,box])
#         pt.figure()
#         pt.title('After')
#         for i in range(points):
#             pt.plot(X2[i],Y2[i])
#         for i in range(points):
#             pt.plot(X2[i],Y2[i])
#         pt.plot([-box,box,box,-box,-box], [box,box,-box,-box,box])
#     
# =============================================================================






#Data collection of arm lengths 

#create workbook
# =============================================================================
# wb = Workbook()
# sheet1 = wb.add_sheet('e=0.4') 
# =============================================================================

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

       
'''

#Plotting the distribution of arm lengths with interactions
    
df1 = pd.read_excel('8_armed_star_interactions_data.xls', sheetname=2)
arm_length1 = df1['arm_length'].tolist()

df2 = pd.read_excel('8_armed_star_interactions_data.xls', sheetname=3)
arm_length2 = df2['arm_length'].tolist()

print(np.mean(arm_length1), np.mean(arm_length2))



pt.figure()
pt.title('Distribution of Arm End-to-end Distances \n for an Eight-armed Star in 2-D')
pt.xlabel('Arm end-to-end distance')
pt.ylabel('Frequency')
pt.hist(arm_length1, normed=True, bins=30)

pt.figure()
pt.title('Distribution of Arm End-to-end Distances \n for an Eight-armed Star in 2-D')
pt.xlabel('Arm end-to-end distance')
pt.ylabel('Frequency')
pt.hist(arm_length2, normed=True, bins=30)




'''

#Plotting radius of gyration vs epsilon 

df3 = pd.read_excel('eight_armed_star_radius_of_gyration_data.xlsx', sheetname=0)
strength = df3['epsilon'].tolist()
radius_of_gyration = df3['radius of gyration'].tolist()


pt.scatter(strength, radius_of_gyration)




#Acceptance Fraction data collection and plotting

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

#Plotting acceptance fraction
#Import data
df4 = pd.read_excel('8_armed_star_interactions_data.xls', sheet_name = 0)
df5 = pd.read_excel('8_armed_star_interactions_data.xls', sheet_name = 1)


#Plotting for epsilon = 0.4
p_b_weak = df4['p before'].tolist()
acc_b_weak = df4['acceptance fraction before'].tolist()

p_a_weak = df4['p after'].tolist()
acc_a_weak = df4['acceptance fraction after'].tolist()


pt.figure()
pt.title('Acceptance Fraction for Eight-armed Star')
pt.xlabel('$p$')
pt.ylabel('Acceptance Fraction')
m,b = np.polyfit(p_b_weak, acc_b_weak, 1)
n,a = np.polyfit(p_a_weak, acc_a_weak, 1)
x_lin = np.linspace(0,1)
y_lin1 = m*x_lin + b
y_lin2 = n*x_lin + a
pt.plot(p_b_weak, acc_b_weak,'bx', label = 'Before equilibration')
pt.plot(p_a_weak, acc_a_weak,'rx', label = 'After equilibration')
pt.plot(x_lin, y_lin1)
pt.plot(x_lin, y_lin2)
pt.legend(loc = 'upper right')

#Plotting for epsilon = 1.0
p_b_strong = df5['p before'].tolist()
acc_b_strong = df5['acceptance fraction before'].tolist()

p_a_strong = df5['p after'].tolist()
acc_a_strong = df5['acceptance fraction after'].tolist()

pt.figure()
pt.title('Acceptance Fraction for Eight-armed Star')
pt.xlabel('$p$')
pt.ylabel('Acceptance Fraction')
m,b = np.polyfit(p_b_strong, acc_b_strong, 1)
n,a = np.polyfit(p_a_strong, acc_a_strong, 1)
x_lin = np.linspace(0,1)
y_lin1 = m*x_lin + b
y_lin2 = n*x_lin + a
pt.plot(p_b_weak, acc_b_strong,'bx', label = 'Before equilibration')
pt.plot(p_a_weak, acc_a_strong,'rx', label = 'After equilibration')
pt.plot(x_lin, y_lin1)
pt.plot(x_lin, y_lin2)
pt.legend(loc = 'upper right')


#Plotting percentage of original and alternative pivot moves rejected

alt_b_weak = df4['percentage of alt moves rejected before'].tolist()
orig_b_weak = df4['percentage of orig moves rejected before'].tolist()
alt_a_weak = df4['percentage of alt moves rejected after'].tolist()
orig_a_weak = df4['percentage of orig moves rejected after'].tolist()

alt_b_strong = df5['percentage of alt moves rejected before'].tolist()
orig_b_strong = df5['percentage of orig moves rejected before'].tolist()
alt_a_strong = df5['percentage of alt moves rejected after'].tolist()
orig_a_strong = df5['percentage of orig moves rejected after'].tolist()


pt.figure()
pt.title('A Graph Showing the Percentage of Unsuccessful \n Pivot Steps Before Equilibration')
pt.xlabel('$x$')
pt.ylabel('Percentage Unsuccessful')
pt.plot(p_b_weak, alt_b_weak, 'bx', label = 'Alternative pivot steps')
pt.plot(x_lin, y_lin1)
pt.plot(p_b_weak, orig_b_weak, 'rx', label = 'Original pivot steps')
pt.legend()


pt.figure()
pt.title('A Graph Showing the Percentage of Unsuccessful \n Pivot Steps After Equilibration')
pt.xlabel('$x$')
pt.ylabel('Percentage Unsuccessful')
pt.plot(p_b_weak, alt_a_weak, 'bx', label = 'Alternative pivot steps')
pt.plot(p_b_weak, orig_a_weak, 'rx', label = 'Original pivot steps')
pt.legend()


pt.figure()
pt.title('A Graph Showing the Percentage of Unsuccessful \n Pivot Steps Before Equilibration')
pt.xlabel('$x$')
pt.ylabel('Percentage Unsuccessful')
pt.plot(p_b_weak, alt_b_strong, 'bx', label = 'Alternative pivot steps')
pt.plot(p_b_weak, orig_b_strong, 'rx', label = 'Original pivot steps')
pt.legend()


pt.figure()
pt.title('A Graph Showing the Percentage of Unsuccessful \n Pivot Steps After Equilibration')
pt.xlabel('$x$')
pt.ylabel('Percentage Unsuccessful')
pt.plot(p_b_weak, alt_a_strong, 'bx', label = 'Alternative pivot steps')
pt.plot(p_b_weak, orig_a_strong, 'rx', label = 'Original pivot steps')
pt.legend()
'''

