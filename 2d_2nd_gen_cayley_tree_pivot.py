#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 13:49:47 2019

@author: kai
"""

import numpy as np
import matplotlib.pyplot as pt
import time
import pandas as pd


def SAW(max_algs, X,Y, X_branches1, Y_branches1, X_branches2, Y_branches2):
    global neighbours,rejections, table, epsilon, points, branches1, branches2
    length = temp_length
    branch_length1 = b1_length
    branch_length2 = b2_length
    neighbours = 0
    epsilon = 0.9
    points = 4
    branches1 = 8
    branches2 = 16
    
    #below are the 7 transformation matrices which can be applied
    matrices = np.array([[[0,-1],[1,0]],[[-1,0],[0,-1]],
                        [[0,1],[-1,0]], [[1,0],[0,-1]], [[-1,0], [0,1]], [[0,1], [1,0]], [[0,-1], [-1,0]]])
    global X2, Y2, X2_branches1, Y2_branches1, X2_branches2, Y2_branches2
    X2, X2_branches1, X2_branches2 = np.copy(X), np.copy(X_branches1), np.copy(X_branches2)
    Y2, Y2_branches1, Y2_branches2 = np.copy(Y), np.copy(Y_branches1), np.copy(Y_branches2)            
    
    #loop for number of times pivot applied
    rejections = 0
    table = {}
    for i in range(max_algs):      
        table2 = table
        X3, X3_branches1, X3_branches2 = np.copy(X2), np.copy(X2_branches1), np.copy(X2_branches2)
        Y3, Y3_branches1, Y3_branches2 = np.copy(Y2), np.copy(Y2_branches1), np.copy(Y2_branches2)
        rand_matrix = np.random.randint(0,len(matrices)-1)
        trans_matrix = matrices[rand_matrix]
        
        
        #loop for applying pivot to end of walk
        arm_or_branch = np.random.randint(0,3)
        if arm_or_branch == 0:
            selector = np.random.randint(0,points)
            pivot = np.random.randint(1,length - 2)
            j = pivot + 1
            k = 0
            l = 0
            while j < length:
                [X2[selector][j], Y2[selector][j]] = trans_matrix.dot(([X2[selector][j] - X2[selector][pivot], Y2[selector][j] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]   
                j = j + 1 
            while k < branch_length1:
                [X2_branches1[2*selector][k], Y2_branches1[2*selector][k]] = trans_matrix.dot(([X2_branches1[2*selector][k] - X2[selector][pivot], Y2_branches1[2*selector][k] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]
                [X2_branches1[2*selector+1][k], Y2_branches1[2*selector+1][k]] = trans_matrix.dot(([X2_branches1[2*selector+1][k] - X2[selector][pivot], Y2_branches1[2*selector+1][k] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]
                k = k + 1
            while l < branch_length2:
                [X2_branches2[4*selector][l], Y2_branches2[4*selector][l]] = trans_matrix.dot(([X2_branches2[4*selector][l] - X2[selector][pivot], Y2_branches2[4*selector][l] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]
                [X2_branches2[4*selector+1][l], Y2_branches2[4*selector+1][l]] = trans_matrix.dot(([X2_branches2[4*selector+1][l] - X2[selector][pivot], Y2_branches2[4*selector+1][l] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]
                [X2_branches2[4*selector+2][l], Y2_branches2[4*selector+2][l]] = trans_matrix.dot(([X2_branches2[4*selector+2][l] - X2[selector][pivot], Y2_branches2[4*selector+2][l] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]
                [X2_branches2[4*selector+3][l], Y2_branches2[4*selector+3][l]] = trans_matrix.dot(([X2_branches2[4*selector+3][l] - X2[selector][pivot], Y2_branches2[4*selector+3][l] - Y2[selector][pivot]])) + [X2[selector][pivot], Y2[selector][pivot]]
                l= l + 1
        
        if arm_or_branch == 1:
            selector = np.random.randint(0, branches1)
            pivot  = np.random.randint(1, branch_length1 - 2)
            j = pivot + 1
            k = 0
            while j < branch_length1:
                [X2_branches1[selector][j], Y2_branches1[selector][j]] = trans_matrix.dot(([X2_branches1[selector][j] - X2_branches1[selector][pivot], Y2_branches1[selector][j] - Y2_branches1[selector][pivot]])) + [X2_branches1[selector][pivot], Y2_branches1[selector][pivot]]
                j = j + 1
            while k < branch_length2:
                [X2_branches2[2*selector][k], Y2_branches2[2*selector][k]] = trans_matrix.dot(([X2_branches2[2*selector][k] - X2_branches1[selector][pivot], Y2_branches2[2*selector][k] - Y2_branches1[selector][pivot]])) + [X2_branches1[selector][pivot], Y2_branches1[selector][pivot]]
                [X2_branches2[2*selector+1][k], Y2_branches2[2*selector+1][k]] = trans_matrix.dot(([X2_branches2[2*selector+1][k] - X2_branches1[selector][pivot], Y2_branches2[2*selector+1][k] - Y2_branches1[selector][pivot]])) + [X2_branches1[selector][pivot], Y2_branches1[selector][pivot]]
                k = k + 1
                
        
        elif arm_or_branch == 2:
            selector = np.random.randint(0,branches2)
            pivot  = np.random.randint(1, branch_length2 - 2)
            j = pivot + 1
            while j < branch_length2:
                [X2_branches2[selector][j], Y2_branches2[selector][j]] = trans_matrix.dot(([X2_branches2[selector][j] - X2_branches2[selector][pivot], Y2_branches2[selector][j] - Y2_branches2[selector][pivot]])) + [X2_branches2[selector][pivot], Y2_branches2[selector][pivot]]
                j = j + 1
                
        
        
    
        table = {}
        overlap = False
        for i in range(points):
            for j in range(length):
                table[X2[i][j],Y2[i][j]] = True           
        for i in range(branches1):
            for j in range(branch_length1):
                table[X2_branches1[i][j], Y2_branches1[i][j]] = True 
        for i in range(branches2):
            for j in range(branch_length2):
                table[X2_branches2[i][j], Y2_branches2[i][j]] = True
                
        if len(table) < points*length + branches1*branch_length1 + branches2*branch_length2 - 27:
            overlap = True
            table = table2
 
        '''if overlap == False:
            old_neighbours = neighbours      
            neighbours = 0 
            for (i, j) in table2:
                if (i+1, j) in table2:
                    neighbours = neighbours + 1
                if (i-1, j) in table2:
                    neighbours = neighbours + 1
                if (i, j+1) in table2:
                    neighbours = neighbours + 1
                if (i, j-1) in table2:
                    neighbours = neighbours + 1
                    
            #neighbours = neighbours -2*points*end_length + 2*points
            #neighbours = neighbours/2
            
                    
            if neighbours < old_neighbours:
                acc = np.random.rand()
                if acc > np.exp((epsilon*(neighbours-old_neighbours))):
                    overlap = True
                    neighbours = old_neighbours'''
                    
        
                
        if overlap:
            X2, X2_branches1, X2_branches2 = np.copy(X3), np.copy(X3_branches1), np.copy(X3_branches2)
            Y2, Y2_branches1, Y2_branches2 = np.copy(Y3), np.copy(Y3_branches1), np.copy(Y3_branches2)
            rejections = rejections + 1        
        
        
        
temp_length = 100
b1_length = int(temp_length/4)
b2_length = int(b1_length/2)

t0 = time.time()
zeros = np.zeros(temp_length)
line = np.arange(temp_length)


branch1 = np.arange(b1_length)
b1_loc1 = [temp_length-1]*b1_length
b1_loc2 = [-temp_length+1]*b1_length


branch2_u = np.arange(b2_length) + (temp_length-1)
branch2_d = -np.arange(b2_length) + (temp_length-1)
b2_loc1 = [b1_length-1]*b2_length
b2_loc2 = [-b1_length+1]*b2_length

x_branches1 = [-branch1, branch1, b1_loc1, b1_loc1, branch1, -branch1, b1_loc2, b1_loc2]
y_branches1 = [b1_loc1, b1_loc1, branch1, -branch1, b1_loc2, b1_loc2, -branch1, branch1]



x_branches2 = [b2_loc2, b2_loc2, b2_loc1, b2_loc1, branch2_u, branch2_d, branch2_u, branch2_d,
               b2_loc1, b2_loc1, b2_loc2, b2_loc2, -branch2_u, -branch2_d, -branch2_u, -branch2_d]
y_branches2 = [branch2_u, branch2_d, branch2_u, branch2_d, b2_loc1, b2_loc1, b2_loc2, b2_loc2,
               -branch2_u, -branch2_d, -branch2_u, -branch2_d, b2_loc2, b2_loc2, b2_loc1, b2_loc1]


X0 = [zeros, line, zeros, -line] 
Y0 = [line, zeros, -line, zeros] 


     
#SAW(3000, X0, Y0, x_branches1, y_branches1, x_branches2, y_branches2)


#Plotting initial configuration
pt.figure()
pt.title('2nd Generation Cayley Tree', fontsize=15)
pt.xlabel('$x$',fontsize=11)
pt.ylabel('$y$', fontsize=11)
for i in range(points):
    pt.plot(X0[i],Y0[i])
    
for j in range(branches1):
    pt.plot(x_branches1[j], y_branches1[j])
    
for k in range(16):
    pt.plot(x_branches2[k], y_branches2[k])
    
        
#Plotting chain after pivot steps
pt.figure()
pt.title('2nd Generation Cayley Tree', fontsize=15)
pt.xlabel('$x$',fontsize=11)
pt.ylabel('$y$', fontsize=11)
for i in range(points):
    pt.plot(X2[i],Y2[i])
    
for j in range(branches1):
    pt.plot(X2_branches1[j], Y2_branches1[j])
    
for k in range(16):
    pt.plot(X2_branches2[k], Y2_branches2[k])
    
print(len(table))
  
#Equilibration
N=765
x_axis=[]
rsq = []
for count in range(1,2000):
    print(count)
    SAW(5*count, X0, Y0, x_branches1, y_branches1, x_branches2, y_branches2)
    Rg_sq = 0
    R_cm = 0
    x_axis.append(5*count)
    for (i,j) in table:
        R_cm = R_cm + (1/N)*np.sqrt((i**2 + j**2))
    for (l,m) in table:
        Rg_sq = Rg_sq + np.sqrt((1/N)*(np.sqrt(l**2+m**2) - R_cm)**2)
    rsq.append(Rg_sq)
pt.figure()
pt.title('Equilibrium of 2nd Generation Cayley Tree')
pt.xlabel('Number of pivots applied')   
pt.ylabel('$R_g$')     
pt.plot(x_axis,rsq)





