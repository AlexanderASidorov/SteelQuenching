#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 15:07:47 2023

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
from initial_data import Composition, t, T, delta_x, delta_T, delta_t




class Fraction (Factors):
    
    def __init__(self, Composition=Composition, t=t, T=T, delta_x=delta_x):
        
        self.Factors=Factors(Composition)
        self.Sigmoidal=Sigmoidal(delta_x)
        
        # for convenience reasons we return here all the transformation factors
        # and temperatures
        self.FC=self.Factors.Alloy.FC
        self.PC=self.Factors.Alloy.PC
        self.BC=self.Factors.Alloy.BC
        self.alpha_martensite=self.Factors.Alloy.alpha_martensite
        
        self.Ae1=self.Factors.Alloy.Ae1
        self.Ae3=self.Factors.Alloy.Ae3
        self.Bs=self.Factors.Alloy.Bs
        self.Ms=self.Factors.Alloy.Ms
        
        self.factor_f=self.Factors.factor_f
        self.factor_p=self.Factors.factor_p
        self.factor_b=self.Factors.factor_b
        
        # Initial data (time, Temperature and steps)
        self.t=t
        self.T=T
        self.delta_t=self.Factors.delta_t
        self.delta_T=self.Factors.delta_T
        
        # Interpolation functions
        self.t2T=self.draw_thermal_cycle()[2]
        self.T2t=self.draw_thermal_cycle()[3]
        
        # Interpolated time and temperature
        self.t_int=self.draw_thermal_cycle()[0]
        self.T_int=self.draw_thermal_cycle()[1]
        
        # Phase fractions
        
        #Ferrite
        self.filtr_f=self.get_ferrite_fraction ()[0]
        self.T__f=self.get_ferrite_fraction ()[1]
        self.nucleation_time_f=self.get_ferrite_fraction ()[2]
        self.f_f=self.get_ferrite_fraction ()[3]
        
        #Pearlite
        self.filtr_p=self.get_pearlite_fraction ()[0]
        self.T__p=self.get_pearlite_fraction ()[1]
        self.nucleation_time_p=self.get_pearlite_fraction ()[2]
        self.f_p=self.get_pearlite_fraction ()[3]
        
        #Bainite
        self.filtr_b=self.get_bainite_fraction ()[0]
        self.T__b=self.get_bainite_fraction ()[1]
        self.nucleation_time_b=self.get_bainite_fraction ()[2]
        self.f_b=self.get_bainite_fraction ()[3]
        
        # Martensite
        self.filtr_m=self.get_martensite_fraction ()[0]
        self.T__m=self.get_martensite_fraction ()[1]
        self.f_m=self.get_martensite_fraction ()[2]
        
        # Temparature gaps for each phase
        self.temp_f =self.Factors.temp_f
        self.temp_p =self.Factors.temp_p
        self.temp_b =self.Factors.temp_b
        
        # time gaps at each temperature for each phase
        self.tau_f=self.Factors.tau_f
        self.tau_p=self.Factors.tau_p
        self.tau_b=self.Factors.tau_b
        
        # Data frame with f_f, f_p, f_b and f_m results
        self.f_data=self.create_dataframe()
        
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
        TimeStep=self.Factors.delta_t
        
        # Defining time range
        t_int = np.arange(min(t), max(t)+TimeStep, TimeStep)
        
        # Defining Temperature range
        T_int = t2T(t_int)
        return t_int, T_int, t2T, T2t
    
    def get_ferrite_fraction (self):
        """
        Calculate Ferrite fraction based on thermal cycle and composition

        """
        # ititial data for calculation
        T=self.T_int
        t=self.t_int
        dt = self.delta_t
        
        # initialize necessary variables
        nucleation_time = np.full(t.shape, 0, dtype=float)
        f = np.full(T.shape, 0, dtype=float)
        
        
        filtr = (T < self.Ae3) & (T > self.Bs)
        
        if np.any(filtr):
            nucleation_time[filtr] = dt/self.Factors.get_ferrite_factor(T[filtr])[0]
            nucleation_time = nucleation_time.cumsum()
            if T[0] < self.Ae3:
                # This is the factor corresponding to the transformed fraction at t[0]
                 nucleation_time += min(t)/self.Factors.get_ferrite_factor(T[0])[0]
            filtr = (nucleation_time >= self.Sigmoidal.ymin) & (nucleation_time <= self.Sigmoidal.ymax)
        
            if np.any(filtr):
                f[filtr] = self.Sigmoidal.x_from_y(nucleation_time[filtr])
                f[nucleation_time < self.Sigmoidal.ymin] = 0
                f[nucleation_time > self.Sigmoidal.ymax] = 1
        
        #f=1
        return filtr, T, nucleation_time, f
    
    def get_pearlite_fraction (self):
        """
        Calculate Pearlite fraction based on thermal cycle and composition

        """
        # ititial data for calculation
        T=self.T_int
        t=self.t_int
        dt = self.delta_t
        
        # initialize necessary variables
        nucleation_time = np.full(t.shape, 0, dtype=float)
        f = np.full(T.shape, 0, dtype=float)
        
        
        filtr = (T < self.Ae1) & (T > self.Bs)
        
        if np.any(filtr):
            nucleation_time[filtr] = dt/self.Factors.get_pearlite_factor(T[filtr])[0]
            nucleation_time = nucleation_time.cumsum()
            if T[0] < self.Ae3:
                # This is the factor corresponding to the transformed fraction at t[0]
                  nucleation_time += min(t)/self.Factors.get_pearlite_factor(T[0])[0]
            filtr = (nucleation_time >= self.Sigmoidal.ymin) & (nucleation_time <= self.Sigmoidal.ymax)
        
            if np.any(filtr):
                f[filtr] = self.Sigmoidal.x_from_y(nucleation_time[filtr])
                f[nucleation_time < self.Sigmoidal.ymin] = 0
                f[nucleation_time > self.Sigmoidal.ymax] = 1
        
        #f=1
        return filtr, T, nucleation_time, f
    
    def get_bainite_fraction (self):
        """
        Calculate Bainite fraction based on thermal cycle and composition

        """
        # ititial data for calculation
        T=self.T_int
        t=self.t_int
        dt = self.delta_t
        
        # initialize necessary variables
        nucleation_time = np.full(t.shape, 0, dtype=float)
        f = np.full(T.shape, 0, dtype=float)
        
        
        filtr = (T < self.Bs) & (T > self.Ms)
        
        if np.any(filtr):
            nucleation_time[filtr] = dt/self.Factors.get_bainite_factor(T[filtr])[0]
            nucleation_time = nucleation_time.cumsum()
            if T[0] < self.Bs:
                # This is the factor corresponding to the transformed fraction at t[0]
                  nucleation_time += min(t)/self.Factors.get_bainite_factor(T[0])[0]
            filtr = (nucleation_time >= self.Sigmoidal.ymin) & (nucleation_time <= self.Sigmoidal.ymax)
        
            if np.any(filtr):
                f[filtr] = self.Sigmoidal.x_from_y(nucleation_time[filtr])
                f[nucleation_time < self.Sigmoidal.ymin] = 0
                f[nucleation_time > self.Sigmoidal.ymax] = 1
        
        #f=1
        return filtr, T, nucleation_time, f
        
    def get_martensite_fraction (self):
       """
       Calculate Martensite fraction based on thermal cycle and composition

       """
       # ititial data for calculation
       T=self.T_int
       t=self.t_int
       dt = self.delta_t
       
       # initialize necessary variables
       nucleation_time = np.full(t.shape, 0, dtype=float)
       f = np.full(T.shape, 0, dtype=float)
       
       
       filtr = T < self.Ms
       
       if np.any(filtr):
           f[filtr] = 1 - np.exp(-self.alpha_martensite*(self.Ms - T[filtr]))       
             
       #f=1
       return filtr, T, f
        
    def create_dataframe(self):
        """
        Pandas dataframe for collecting fraction data

        """      
        f = pd.DataFrame(columns=['t', 'T', 'f_f', 'f_p','f_b', 'f_m'])
        f['t'] = self.t_int
        f['T'] = self.T_int
        f.fillna(0, inplace=True)
        f['f_f'] = self.f_f
        f['f_p'] = self.f_p
        f['f_b'] = self.f_b
        f['f_m'] = self.f_m
                
        return f        
        
        
        
        
        