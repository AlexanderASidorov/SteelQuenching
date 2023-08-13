#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 12:51:46 2023

@author: aas
"""

# import numpy as np
# from scipy import integrate
# from scipy.interpolate import splrep, splev, interp1d
# from scipy.optimize import root
# from alloy import Alloy
# from sigmoidal import Sigmoidal
from factors import Factors
import matplotlib.pyplot as plt
# from initial_data import Composition, t, T, delta_x, delta_T, delta_t
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
        
        
        #########################################
        # for TTTT
        #########################################
        # Start Curves
        self.start_f=self.Factors.tau_f[:, 0]
        self.start_p=self.Factors.tau_p[:, 0]
        self.start_b=self.Factors.tau_b[:, 0]
         
        # End Curves
        self.end_f=self.Factors.tau_f[:, -1]
        self.end_p=self.Factors.tau_p[:, -1]
        self.end_b=self.Factors.tau_b[:, -1]
        ##########################################
        # Plot TTT diagram
        self.TTT=self.plottingTTT()
        
        
        ###########################################
        # for f_f, f_p, f_b and f_m
        ###########################################
        self.f_f=self.Fraction.f_f
        self.f_p=self.Fraction.f_p
        self.f_b=self.Fraction.f_b
        self.f_m=self.Fraction.f_m
        ###########################################
        # Plot uncorrected fracture
        #self.f_uncorrected=self.plotting_f_uncorrected()
    
    def plottingTTT(self):
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
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        return fig01
    
    def plotting_f_uncorrected(self):
        fig02=plt.figure()
       
        # Ferrite curves
        plt.plot(self.t_int, self.f_f, '--', label="f_f")
        plt.plot(self.t_int, self.f_p, label="f_p")
        plt.plot(self.t_int, self.f_b, label="f_b")
        plt.plot(self.t_int, self.f_m, '-.', label="f_f")
      
        # Plot settings
        #plt.xscale("log")
        plt.xlabel('Time, sec')
        plt.ylabel('Fraction')
        plt.title('Uncorrected fractures')
        plt.legend()
        plt.grid()
        return fig02
    

        
    def plotting_phase_change(data):
        # Plotting data with X time and Temperature  
        
        # X data
        time = data.loc[:,'time']
        Temperature = data.loc[:, 'Temperature']
        
        # Y data
        Ferrite = data.loc[:,'Ferrite']
        Pearlite = data.loc[:,'Pearlite']
        Bainite = data.loc[:,'Bainite']
        Martensite = data.loc[:,'Martensite']
        Austenite = data.loc[:, 'Austenite']
        
        # The main plot is Phase Fraction - time, so we build all the curves
        # in this coordinates
        fig03, ax1 = plt.subplots()
                      
        # Phase change curves
        ax1.plot(time, Ferrite, '--', label="Ferrite")
        ax1.plot(time, Pearlite, label="Pearlite")
        ax1.plot(time, Bainite, label="Bainite")
        ax1.plot(time, Martensite, label="Martensite")
        ax1.plot(time, Austenite, label="Austenite")
      
        # X settings for time
        ax1.set_xlabel('Time, sec')
        ax1.set_ylabel('Phase Fraction')
        ax1.grid()
        ax1.set(xlim=(0, max(time)))
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        # add an axis to show the correspondent Temperature
        ax2 = ax1.twiny()
        # Add some extra space for the second axis at the bottom
        fig03.subplots_adjust(bottom=0.1)

        # Move twinned axis ticks and label from top to bottom
        ax2.xaxis.set_ticks_position("bottom")
        ax2.xaxis.set_label_position("bottom")
        # Offset the twin axis below the host
        ax2.spines["bottom"].set_position(("axes", -0.175))
        
        
        # X settings for Temperature
        ax2.set_xlabel(u'Temperature, \u00B0C')
        ax2.set(xlim=(max(Temperature), min(Temperature)))
      
        
        
        return fig03     
    
    def plotting_comparison (data_01, data_02, phase):
        # Plotting data with X time and Temperature  
        
        # Y data
        fraction_caclulated = data_01.loc[:, phase]
        fraction_FV = data_02.loc[:, phase]
        
        # X data
        time_caclulated = data_01.loc[:, 'time']
        time_FV = data_02.loc[:, 'time']
        
        fig04, ax1 = plt.subplots()
                      
        # Phase change curves
        ax1.plot(time_caclulated, fraction_caclulated, label="Calculated")
        ax1.plot(time_FV, fraction_FV, label="FlowVison")
       
        # X settings for time
        ax1.set_xlabel('Time, sec')
        ax1.set_ylabel('Phase Fraction')
        ax1.grid()
        ax1.set(xlim=(0, max(time_FV)))
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        plt.title(phase)
             
        
        
        return fig04    
        
        
   
        
        
        
        
        
        
        