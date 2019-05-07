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

import numpy as np
import matplotlib.pyplot as pt
import time
import pandas as pd
import xlrd, xlwt
from xlwt import Workbook
from xlutils.copy import copy as xl_copy

def SAW(max_algs, X,Y):
    global neighbours, rejections
    global table
    global points
    points = len(X)
    length = len(X[0])
    
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
    rejections = 0
    table = {}
    while algs < max_algs:
        table2 = dict(table)
        X3 = np.copy(X2)
        Y3 = np.copy(Y2)
        pivot = np.random.randint(1,length - 2)
        rand_matrix = np.random.randint(0,len(matrices)-1)
        trans_matrix = matrices[rand_matrix]
        #loop for applying pivot to end of walk
        j = pivot + 1
        selector = np.random.randint(0,points)
        while j < length:
            [X2[selector][j], Y2[selector][j]] = trans_matrix.dot(([X2[selector][j] - X2[selector][pivot], Y2[selector][j] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]      
            j = j + 1
            
        
        overlap = False
        table = {}
        for i in range(length):
            for j in range(points):
                table[X2[j][i],Y2[j][i]] = True    
        if len(table) < points*length:
            overlap = True
            
        
                
     
        if overlap:
            table = dict(table2)
            X2 = np.copy(X3)
            Y2 = np.copy(Y3)
            rejections = rejections + 1
        algs = algs + 1
        
temp_length = 150

t0 = time.time()
offset=1
box = 2
zeros = np.zeros(temp_length)
line = np.arange(box,temp_length+box)

X0 = [line,[offset]*temp_length, -line, line, -line, [-offset]*temp_length, [offset]*temp_length, [-offset]*temp_length]
Y0 = [[offset]*temp_length, line, [offset]*temp_length, [-offset]*temp_length, [-offset]*temp_length, line, -line, -line]







'''

N = temp_length*8

#Collecting data about radius of gyration 

rsq = []
#SAW(75000,X0,Y0)
for count in range(500):
    SAW(100, X2, Y2)
    Rg = 0
    R_cm = 0
    for (i,j) in table:
        R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2))
    for (l,m) in table:
        Rg = Rg + np.sqrt((1/N)*(np.sqrt(l**2+m**2) - R_cm)**2)
    print(Rg)









#Equilibration
x_axis=[]
rsq = []
for count in range(1,2000):
    print(count)
    x_axis.append(5*count)
    SAW(temp_length, 5*count, X0,Y0)
    Rg_sq = 0
    R_cm = 0
    for (i,j) in table:
        R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2))
    for (l,m) in table:
        Rg_sq = Rg_sq + np.sqrt((1/N)*(np.sqrt(l**2+m**2) - R_cm)**2)
    rsq.append(Rg_sq)
pt.figure()
pt.title('Equilibrium of Eight-Armed Star')
pt.xlabel('Number of pivots applied')   
pt.ylabel('$R_g$')     
pt.plot(x_axis,rsq)




#Data collection of arm lengths 


#create workbook
wb = Workbook()
sheet1 = wb.add_sheet('e=0.4') 


#Add sheets to existing workbook
wb = xlrd.open_workbook('star_data.xls', formatting_info=True)
wb3 = xl_copy(wb)
sheet3 = wb3.add_sheet('8_armed')

#Equilibrate
SAW(25000, X0, Y0)

#collect data
for i in range(1000):
    print(i)
    SAW(400, X2, Y2)
    for n in range(points):
        r_n = np.sqrt(X2[n][temp_length-1]**2 + Y2[n][temp_length-1]**2)
        #sheet1.write(points*i+n,3, r_n)
        sheet3.write(points*i+n,3, r_n)

#save to Excel spreadsheet
wb3.save('star_data.xls') 

#SAW(temp_length,1000, X0, Y0)
pt.figure()
pt.title('An Eight-armed Star Generated Using the Pivot Algorithm')
pt.xlabel('$x$')
pt.ylabel('$y$')
for i in range(points):
    pt.plot(X2[i],Y2[i])
pt.plot([-box,box,box,-box,-box], [box,box,-box,-box,box])
t1=time.time()
print(t1-t0)

'''

#Plotting the distribution of arm lengths without interactions


df = pd.read_excel('star_data.xls', sheetname=2) # can also index sheet by name or fetch all sheets
arm_lengths = df['arm_lengths'].tolist()
print(len(arm_lengths))



pt.figure()
pt.title('Distribution of Arm End-to-end Distances \n for an Eight-armed Star in 2-D')
pt.xlabel('Arm end-to-end distance')
pt.ylabel('Frequency')
pt.hist(arm_lengths, normed=True, bins=30)


'''

#acceptance fraction 
SAW(temp_length, 10000, X0, Y0)

for i in range(5):                                                              
    SAW(temp_length, 10000, X2, Y2)
    print(rejections)
    

#plotting acceptance fraction for eight armed star
df = pd.read_excel('2d_star_data.xlsx', sheetname=6) # can also index sheet by name or fetch all sheets
arm_length = df['log_10(length)'].tolist()
acc = df['log_10(acc)'].tolist()


pt.figure()
m,b = np.polyfit(arm_length, acc, 1)
x_lin = np.linspace(0,4)
y_lin = m*x_lin + b
pt.title('Acceptance Fraction for Eight-armed Star in 2D')
pt.xlabel('$\log_{10}(r)$')
pt.ylabel('$\log_{10}(f)$')
pt.plot(arm_length, acc, 'x')
pt.plot(x_lin, y_lin)
'''
















