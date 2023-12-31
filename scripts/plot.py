#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 12:51:46 2023

@author: aas
"""

import numpy as np
# from scipy import integrate
# from scipy.interpolate import splrep, splev, interp1d
# from scipy.optimize import root
# from alloy import Alloy
# from sigmoidal import Sigmoidal
from factors import Factors
import matplotlib.pyplot as plt
from initial_data import t, T
from fraction import Fraction
from solve_JMAK import JMAK


#Composition=Composition

class Plot (Fraction):
    
    def __init__(self):
        
        
        
        self.Factors=Factors()
        self.Fraction=Fraction()
        
        
        ###########################################
        # for martensite rate
        ###########################################
        self.rate_martensite=self.Fraction.rate_m
        
        
        
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
        # Plot Martensite growth rate
        #self.rate_m_plot=self.plotting_rate_martensite()
        
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
        plt.xlabel('Time, sec', fontsize=16)
        plt.ylabel(u'Temperature, \u00B0C', fontsize=16)
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
        ax1.set_xlabel('Time, sec', fontsize=16)
        ax1.set_ylabel('Phase Fraction', fontsize=16)
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
        ax2.set_xlabel(u'Temperature, \u00B0C', fontsize=16)
        ax2.set(xlim=(max(Temperature), min(Temperature)))
      
        # calculate and plot cooling rate
        Cooling_rate=(T[0]-T[1])/t[1]
        plt.figtext(0.15, 0.45, 
                    'Cooling rate = ' + str(Cooling_rate) + "$^\circ$C/sec ",
                    fontsize=12)
        
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
        
    def plotting_rate_martensite(self):
        fig05=plt.figure()
        
        rate_m=self.rate_martensite
        time=self.t_int
       
        # Ferrite curves
        plt.plot(time, rate_m, label="Martensite growth rate")
             
        # Plot settings
        #plt.xscale("log")
        plt.xlabel('Time, sec')
        plt.ylabel('Rate')
        plt.title('Martensite growth rate')
        plt.legend()
        plt.grid()
        return fig05    
    
    
    @staticmethod
    def plot_f(Phase, Temp):
        '''
        Plot phase change at constant temperature based on calculated k and n
        
        '''
        # Search for index whete Temperature == Temp
        i=Phase[Phase['temp']==Temp].index.values[0]
        t_s=Phase.time_s[i]
        t_e=Phase.time_e[i]
        #time=[Phase.time_s[i], Phase.time_01[i], Phase.time_05[i], Phase.time_09[i], Phase.time_e[i]]
        time = np.linspace(t_s, t_e, 25)
        tau=Phase.tau[i]
        n=Phase.n[i]
        #f=JMAK.equation(time, k, n)
        f=JMAK.equation_incubation_tau(time, tau, n, t_s)
        
        
        
        fig06=plt.figure()

       
        # Ferrite curves
        plt.plot(time, f)
             
        # Plot settings
        #plt.xscale("log")
        plt.xlabel('Time, sec')
        plt.ylabel('Phase fraction')
        plt.title('Phase change at Temperature = ' + str(Temp) + '$^\circ$C')
        #plt.legend()
        plt.grid()
        return fig06
    
    @staticmethod
    def plot_TTT_imp(Ferrite, Pearlite, Bainite):
        
        '''
        Plot phase change at constant temperature based on calculated k and n
        '''
        fig07=plt.figure()
       
        # Ferrite curves
        plt.plot(Ferrite.time_s, Ferrite.temp,'--', label="Start Ferrite")
        plt.plot(Ferrite.time_e, Ferrite.temp, '--' ,label="End Ferrite")
        # Pearlite curves
        plt.plot(Pearlite.time_s, Pearlite.temp, label="Start Perlite")
        plt.plot(Pearlite.time_e, Pearlite.temp, label="End Perlite")
        # Bainite curves
        plt.plot(Bainite.time_s, Bainite.temp,'-.', label="Start Bainite")
        plt.plot(Bainite.time_e, Bainite.temp, '-.', label="End Bainite")
        
        # Plot settings
        plt.xscale("log")
        plt.xlabel('Time, sec')
        plt.ylabel(u'Temperature, \u00B0C')
        plt.title('Time Temperature Transformation')
        plt.legend()
        plt.grid()
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        
        
        return fig07
    
    
    
    @staticmethod
    def plot_isothermal_transformation(time, fraction, Temp):
        fig08=plt.figure()

       
        # Ferrite curves
        plt.plot(time, fraction)
             
        # Plot settings
        #plt.xscale("log")
        plt.xlabel('Time, sec')
        plt.ylabel('Phase fraction')
        plt.title('Phase change at Temperature = ' + str(Temp) + '$^\circ$C')
        #plt.legend()
        plt.grid()
        return fig08        
        
        
        
       
        
        
        

        
   
        
        
        
        
        
        
        