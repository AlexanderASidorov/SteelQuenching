#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:09:17 2023

@author: aas
"""

import numpy as np
import pandas as pd
from scipy import integrate
from scipy.interpolate import splrep, splev, interp1d
from scipy.optimize import root
from alloy import Alloy
from sigmoidal import Sigmoidal
from factors import Factors
import matplotlib.pyplot as plt
from fraction import Fraction
from initial_data import Composition, t, T, delta_x, delta_T, delta_t



class Solve (Fraction):
    
    def __init__(self, Composition=Composition, t=t, T=T, delta_x=0.1, 
                 delta_T=delta_T, delta_t=delta_t):
        
        self.Factors=Factors()
        self.Sigmoidal=Sigmoidal()
        self.Fraction=Fraction(Composition, t, T)
        self.t_int=self.Fraction.t_int
        self.T_int=self.Fraction.T_int
        
        # for convenience reasons we return here all uncorrected phase fractions
        self.f_f=self.Fraction.f_f
        self.f_p=self.Fraction.f_p
        self.f_b=self.Fraction.f_b
        self.f_m=self.Fraction.f_m
        
        # fraction increments
        self.f_f_inc=self.get_fraction_increment()[0]
        self.f_p_inc=self.get_fraction_increment()[1]
        self.f_b_inc=self.get_fraction_increment()[2]
        self.f_m_inc=self.get_fraction_increment()[3]
        
        # First way of solving system of equations
        # self.f_01=self.solve_01()
        
        # Second way of solving system of equations
        self.f_02=self.solve_02()
        
        
        
    def get_fraction_increment (self):
        """
        Calculate phase fraction increment

        """        
        # increment of phase fraction
        f_f_inc = np.zeros(self.f_f.shape)
        f_p_inc = np.zeros(self.f_p.shape)
        f_b_inc = np.zeros(self.f_b.shape)
        f_m_inc = np.zeros(self.f_m.shape)
        
        f_f_inc[1:] = np.diff(self.f_f)
        f_p_inc[1:] = np.diff(self.f_p)
        f_b_inc[1:] = np.diff(self.f_b)
        f_m_inc[1:] = np.diff(self.f_m)
        
        return f_f_inc, f_p_inc, f_b_inc, f_m_inc
    
    def create_dataframe(self):
        """
        Pandas dataframe for collecting fraction data

        """      
        f = pd.DataFrame(columns=['t', 'T', 'ferrite', 'pearlite','bainite', 'martensite', 'austenite'])
        f['t'] = self.t_int
        f['T'] = self.T_int
        f.fillna(0, inplace=True)
        f.loc[0, 'ferrite'] = self.f_f[0]
        f.loc[0, 'pearlite'] = self.f_p[0]
        f.loc[0, 'bainite'] = self.f_b[0]
        f.loc[0, 'martensite'] = self.f_m[0]
        f.loc[0, 'austenite'] = 1. - self.f_f[0] - self.f_p[0] - self.f_b[0] - self.f_m[0]
        return f
    

    

    def solve_01 (self): # not required
        # Uncorrected phase fractions
        f_f=self.Fraction.f_f
        f_p=self.Fraction.f_p
        f_b=self.Fraction.f_b
        f_m=self.Fraction.f_m
        # time and temperature
        t=self.Fraction.t_int
        T=self.Fraction.T_int
        # phase increments
        f_f_inc = np.zeros(f_f.shape)
        f_p_inc = np.zeros(f_p.shape)
        f_b_inc = np.zeros(f_b.shape)
        f_m_inc = np.zeros(f_m.shape)
        
        f_f_inc[1:] = np.diff(f_f)
        f_p_inc[1:] = np.diff(f_p)
        f_b_inc[1:] = np.diff(f_b)
        f_m_inc[1:] = np.diff(f_m)

        # create dataframe for phase fractions
        f = pd.DataFrame(columns=['t', 'T', 'ferrite', 'pearlite',
                                  'bainite', 'martensite', 'austenite'])
        f['t'] = t
        f['T'] = T
        f.fillna(0, inplace=True)
        f.loc[0, 'ferrite'] = f_f[0]
        f.loc[0, 'pearlite'] = f_p[0]
        f.loc[0, 'bainite'] = f_b[0]
        f.loc[0, 'martensite'] = f_m[0]
        f.loc[0, 'austenite'] = 1. - f_f[0] - f_p[0] - f_b[0] - f_m[0]

        def f1(i, x, y, z, w):
            if f_f[i] < 1:
                return f.loc[i-1, 'ferrite'] + f_f_inc[i]*(1 - x - y - z - w)/(1 - f_f[i]) - x
            else:
                return f.loc[i-1, 'ferrite'] + f_f_inc[i]*(1 - y - z - w) - x

        def f2(i, x, y, z, w):
            if f_p[i] < 1:
                return f.loc[i-1, 'pearlite'] + f_p_inc[i]*(1 - x - y - z - w)/(1 - f_p[i]) - y
            else:
                return f.loc[i-1, 'pearlite'] + f_p_inc[i]*(1 - x - z - w) - y

        def f3(i, x, y, z, w): 
            return f.loc[i-1, 'bainite'] + f_b_inc[i]*(1 - x - y - w) - z

        def f4(i, x, y, z, w): 
            return f.loc[i-1, 'martensite'] + f_m_inc[i]*(1 - x - y - z) - w
        
        for i in range(len(f))[1:]:
            x0 = [f.loc[i-1, 'ferrite'], f.loc[i-1, 'pearlite'],
                  f.loc[i-1, 'bainite'], f.loc[i-1, 'martensite']]

            # Solves system of non-linear equations to get corrected phase fractions
            res = root(lambda x: [f1(i, *x), f2(i, *x), f3(i, *x), f4(i, *x)], x0=x0)

            if res.x[0]<0: f.loc[i, 'ferrite']=0 
            else: f.loc[i, 'ferrite'] = res.x[0]
            
            if res.x[1]<0: f.loc[i, 'pearlite']=0 
            else:f.loc[i, 'pearlite'] = res.x[1]
            
            if res.x[2]<0: f.loc[i, 'bainite']=0 
            else: f.loc[i, 'bainite'] = res.x[2]
            
            if res.x[3]<0: f.loc[i, 'martensite']=0 
            else: f.loc[i, 'martensite'] = res.x[3]
            
            if (1. - res.x.sum())<0: f.loc[i, 'austenite'] = 0
            else: f.loc[i, 'austenite'] = 1. - res.x.sum()
            
        return f.round(7)
    
    def solve_02 (self):
        # Uncorrected phase fractions
        f_f=self.Fraction.f_f
        f_p=self.Fraction.f_p
        f_b=self.Fraction.f_b
        f_m=self.Fraction.f_m
        
        # time and temperature
        t=self.Fraction.t_int
        T=self.Fraction.T_int
        
        # phase increments
        f_f_inc=self.f_f_inc
        f_p_inc=self.f_p_inc
        f_b_inc=self.f_b_inc
        f_m_inc=self.f_m_inc
        
        # variable for calculating corrected phase fractions:
        f=np.zeros([len(t),4])
        
        for i in range(len(t))[1:]:
        # coefficients of linear equations
            #######################
            # Bainite and Martensite coefficients
            c=f_b_inc[i]
            d=f_m_inc[i]
            #######################
        
        # Pearlite and Ferite coeficients depends on whether its uncorrected 
        # fraction >1 or not. There are 4 possible cases:
            #####################################################
            if f_f[i] < 1 and f_p[i]<1: # first possible case
                a=f_f_inc[i]/(1-f_f[i])
                b=f_p_inc[i]/(1-f_p[i])
                       
            # System of linear equations looks like A*X=B
            # so we define A:
            
                A=np.array([[a+1, a, a, a],
                            [b, b+1, b, b],
                            [c, c, 1, c],
                            [d, d, d, 1]])
            #####################################################
            if f_f[i] < 1 and 1<=f_p[i]: # second possible case
                a=f_f_inc[i]/(1-f_f[i])
                b=f_p_inc[i]
            
                A=np.array([[a+1, a, a, a],
                            [b, 1, b, b],
                            [c, c, 1, c],
                            [d, d, d, 1]])
          
            #####################################################
            if 1<=f_f[i] and f_p[i]<1: # third possible case
                a=f_f_inc[i]
                b=f_p_inc[i]/(1-f_p[i])
                       
                      
                A=np.array([[1, a, a, a],
                            [b, b+1, b, b],
                            [c, c, 1, c],
                            [d, d, d, 1]])
            
            ####################################################
            if 1<=f_f[i] and 1<=f_p[i]: # fourth possible case
                a=f_f_inc[i]
                b=f_p_inc[i]
                       
                      
                A=np.array([[1, a, a, a],
                            [b, 1, b, b],
                            [c, c, 1, c],
                            [d, d, d, 1]])
                  
        # and define B
            B = np.array([[f[i-1, 0]+a], [f[i-1, 1]+b], [f[i-1, 2]+c], [f[i-1, 3]+d]])
            # solve system of linear equations:
            Solution=np.linalg.solve(A,B)
            f[i] = np.transpose(Solution)
            for j in range (4):
                if f[i,j] <0: f[i,j]=0
        return f.round(7)

Steel=Solve()
# f_01=Steel.f_01
# f_02=Steel.f_02

# print(f_01.iloc[-1,:])
print('###################')
print(Steel.f_02[-1, :])
     


