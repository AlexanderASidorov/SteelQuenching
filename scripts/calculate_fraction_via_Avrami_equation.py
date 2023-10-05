#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 17:37:54 2023

@author: aas
"""

from plot import Plot
from solve_JMAK import JMAK
from read_data import Create_variables_from_Excel as create
from equilibria import Equilibria
from export import k_n_2Excel, k_2text, n_2text, eq_2text
from initial_data import t, T, delta_x, delta_T, delta_t

import os.path
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import splrep, splev, interp1d
import math
import matplotlib.pyplot as plt

class Calculate_fraction(Equilibria):
    
    def __init__(self, t=t, T=T, delta_x=delta_x, delta_T=delta_T, delta_t=delta_t):
        
        # Initial data (time, Temperature and steps)
        self.t=t
        self.T=T
        self.delta_t=delta_t
        self.delta_T=delta_T
        
        # Interpolation functions
        self.t2T=self.draw_thermal_cycle()[2]
        self.T2t=self.draw_thermal_cycle()[3]
        
        # Interpolated time and temperature
        self.t_int=self.draw_thermal_cycle()[0]
        self.T_int=self.draw_thermal_cycle()[1]
        
        # Phasa dataframes
        self.Ferrite = Equilibria().Ferrite
        self.Pearlite = Equilibria().Pearlite
        self.Bainite = Equilibria().Bainite
        
        # Temperatures of transformation start
        self.Ae3=self.get_points()[0]
        self.Ae1=self.get_points()[1]
        self.Bs=self.get_points()[2]
        self.Ms=self.get_points()[3]
        
    def draw_thermal_cycle(self):
        """
        Draw thermal cycle (cooling curve)

        """
        if len(self.t) > 3:
        # Fits T(t) by spline
            def t2T(t_): return splev(t_, splrep(self.t, self.T))
        else:
        # Uses linear interpolator
            t2T = interp1d(self.t, self.T)
            T2t = interp1d(self.T, self.t)
        
        # Defining time step
        TimeStep=delta_t
        
        # Defining time range
        t_int = np.arange(min(t), max(t)+TimeStep, TimeStep)
        
        # Defining Temperature range
        T_int = t2T(t_int)
        return t_int, T_int, t2T, T2t
    
    def get_points(self):
        Ae3=max(self.Ferrite.temp)
        Ae1=max(self.Pearlite.temp)
        Bs=max(self.Bainite.temp)
        Ms=min(self.Bainite.temp)
        return Ae3, Ae1, Bs, Ms
    
    #@staticmethod
    def coef2temp(self, dataframe):
        k2temp=interp1d(dataframe.temp, dataframe.k)
        n2temp=interp1d(dataframe.temp, dataframe.n)
        f_eq2temp=interp1d(dataframe.temp, dataframe.f_eq)
        return k2temp, n2temp, f_eq2temp
    
    @staticmethod    
    def S(x):
        """
        Sigmoidal function S(x)
        """
        return (x**(0.4*(1. - x))*(1. - x)**(0.4*x))
    
    @staticmethod
    def get_f_from_rate(tau_1, tau_2, constant):
        f=constant*Calculate_fraction.S
        return
        
    
        
#%%       
xdata=np.arange(0, 1.01, 0.01)     

ydata=Calculate_fraction.S(xdata)*100

fig01=plt.figure()

       
# curve
plt.plot(xdata, ydata)
             
# Plot settings
#plt.xscale("log")
plt.xlabel('Fraction')
plt.ylabel('S(x)')
plt.title('S(x) via fraction')
plt.grid()        
        
        




















