#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 12:11:10 2019

@author: kai
"""
import matplotlib.pyplot as pt
import numpy as np

sigma = 342*10**(-12)
epsilon = 128
x = np.linspace(0,10**-9,10000)
def V_a(r):
    return -(sigma/r)**6

def V_r(r):
    return (sigma/r)**12

def V(r):
    return (V_a(r) + V_r(r))

pt.figure()
pt.title('Lennard-Jones (12,6) Potential', fontsize = 15)
pt.ylim(-1,1)
pt.xlim(0,2.5)
pt.xlabel('Separation, $ r/ \sigma $', fontsize=11)
pt.ylabel('Potential, $V/4 \\alpha $', fontsize=11)
pt.plot(x/sigma, V_a(x), label = 'Attractive')
pt.plot(x/sigma, V_r(x), label = 'Repulsive')
pt.plot(x/sigma, V(x), label = 'Total')
pt.plot(x, 'k', linewidth = 1)
pt.legend(loc = 'upper right')


pt.figure()
pt.title('Attractive Potential Between Neighbouring Nodes', fontsize=13)
pt.xlabel('Separation (lattice spaces)', fontsize=11)
pt.ylabel('Potential (arbitrary units)', fontsize=11)
pt.yticks([])
pt.xlim(-0.1,5.1)
pt.ylim(-1.5,1)
pt.plot([1,2,3,4,5], [-1,0,0,0,0], 'ko')