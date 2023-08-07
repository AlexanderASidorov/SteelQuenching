#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 12:51:46 2023

@author: aas
"""

import numpy as np
from scipy import integrate
from scipy.interpolate import splrep, splev, interp1d
from scipy.optimize import root
from alloy import Alloy
from sigmoidal import Sigmoidal
from factors import Factors
import matplotlib.pyplot as plt
from initial_data import Composition, t, T, delta_x, delta_T, delta_t
from fraction import Fraction


#Composition=Composition

class Plot (Fraction):
    
    def __init__(self):
        
        self.Factors=Factors()
        self.Fraction=Fraction()
        
        self.T_int=self.Fraction.T_int
        self.t_int=self.Fraction.t_int
               
        #Temperature gaps
        self.temp_f=self.Factors.temp_f
        self.temp_p=self.Factors.temp_p
        self.temp_b=self.Factors.temp_b
        
        
        # Time
        # Start Curves
        self.start_f=self.Factors.tau_f[:, 0]
        self.start_p=self.Factors.tau_p[:, 0]
        self.start_b=self.Factors.tau_b[:, 0]
         
        # End Curves
        self.end_f=self.Factors.tau_f[:, -1]
        self.end_p=self.Factors.tau_p[:, -1]
        self.end_b=self.Factors.tau_b[:, -1]
        
        # Plot TTT diagram
        self.TTT=self.PlottingTTT()


    def PlottingTTT(self):
        fig01=plt.figure()
       
        # Ferrite curves
        plt.plot(self.start_f, self.temp_f,'--', label="Start Ferrite")
        plt.plot(self.end_f, self.temp_f, '--' ,label="End Ferrite")
        # Pearlite curves
        plt.plot(self.start_p, self.temp_p, label="Start Perlite")
        plt.plot(self.end_p, self.temp_p, label="End Perlite")
        # Bainite curves
        plt.plot(self.start_b, self.temp_b,'-.', label="Start Bainite")
        plt.plot(self.end_b, self.temp_b,'-.', label="End Bainite")
        
        # Thermal cycle
        plt.plot(self.t_int, self.T_int, label="Thermal curve", linewidth=3.0)
        
        
        # Plot settings
        plt.xscale("log")
        plt.xlabel('Time, sec')
        plt.ylabel(u'Temperature, \u00B0C')
        plt.title('Time Temperature Transformation')
        plt.legend()
        plt.grid()
        return fig01
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        