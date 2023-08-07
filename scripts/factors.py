#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:56:51 2023

@author: aas
"""

import numpy as np
from scipy import integrate
from scipy.interpolate import splrep, splev, interp1d
from scipy.optimize import root
from alloy import Alloy
from sigmoidal import Sigmoidal
from initial_data import Composition, t, T, delta_x, delta_T, delta_t




class Factors(Sigmoidal, Alloy):
        
    Q=27500*4.184
    # n1=0.32
    # n2=3
    K=273.15
    R=8.314459
  
        
    def __init__(self, Composition=Composition, delta_T=delta_T, delta_t=delta_t):
        
        self.Alloy=Alloy(Composition)
        self.Sigmoidal=Sigmoidal()
        self.delta_x=self.Sigmoidal.delta_x
        self.delta_T= delta_T
        self.delta_t=delta_t
                       
        # Factor for each phase
        self.factor_f=self.get_ferrite_factor()[0]
        self.factor_p=self.get_pearlite_factor()[0]
        self.factor_b=self.get_bainite_factor()[0]
        
        # Temperature gap for each phase
        self.temp_f=self.get_ferrite_factor()[1]
        self.temp_p=self.get_pearlite_factor()[1]
        self.temp_b=self.get_bainite_factor()[1]
        
        
        # Denumerator
        #self.Denumerator=self.get_ferrite_factor()[2]
        
        # Time necessary for transform into X fraction
        self.tau_f=self.get_tau_f()
        self.tau_p=self.get_tau_p()
        self.tau_b=self.get_tau_b()
    
    
  
    
    def get_ferrite_factor(self, T=None):
       # define Temperature gap
       Ts=round(self.Alloy.Ae3, 0)
       Te=round(self.Alloy.Bs, 0)
       
       if T is None:
           T=np.arange(Ts-1, Te-1, -self.delta_T)
       
       Numerator=self.Alloy.FC
       FirstTerm=2**(self.Alloy.n1_F*self.Alloy.gs)
       SecondTerm=(Ts - T)**self.Alloy.n2_F
       ThirdTerm=np.exp(-self.Q/(self.R*(T + self.K)))
       Denumerator=FirstTerm*SecondTerm*ThirdTerm
       return  Numerator/Denumerator, T #, Denumerator 

        
    def get_pearlite_factor(self, T=None):
       # define Temperature gap
       Ts=round(self.Alloy.Ae1, 0)
       Te=round(self.Alloy.Bs, 0)
       
       if T is None:
           T=np.arange(Ts-1, Te-1, -self.delta_T)
       
       Numerator=self.Alloy.PC
       FirstTerm=2**(self.Alloy.n1_P*self.Alloy.gs)
       SecondTerm=(Ts - T)**self.Alloy.n2_P
       ThirdTerm=np.exp(-self.Q/(self.R*(T + self.K)))
       Denumerator=FirstTerm*SecondTerm*ThirdTerm
       return  Numerator/Denumerator, T
    
    def get_bainite_factor(self, T=None):
       # define Temperature gap
       Ts=round(self.Alloy.Bs, 0)
       Te=round(self.Alloy.Ms, 0)
       
       if T is None:
           T=np.arange(Ts-1, Te-1, -self.delta_T)

       Numerator=self.Alloy.BC
       FirstTerm=2**(self.Alloy.n1_B*self.Alloy.gs)
       SecondTerm=(Ts - T)**self.Alloy.n2_B
       ThirdTerm=np.exp(-self.Q/(self.R*(T + self.K)))
       Denumerator=FirstTerm*SecondTerm*ThirdTerm
       return  Numerator/Denumerator, T
   
    def get_factor(self, Ts, Te, Coef_comp, n1, n2, gs):
        """
        Function for calculation factor for users's defined temperatures and
        phases
        """
        T=np.arange(Ts-1, Te-1, -self.delta_T)
        Numerator=Coef_comp
        FirstTerm=2**(n1*gs)
        SecondTerm=(Ts - T)**n2
        ThirdTerm=np.exp(-self.Q/(self.R*(T + self.K)))
        Denumerator=FirstTerm*SecondTerm*ThirdTerm
        return  Numerator/Denumerator
        
    
    
    def get_tau_f(self):
        # define Temperature gap
        tau_f=np.outer(self.factor_f, self.Sigmoidal.y)
        return tau_f
    
    
    def get_tau_p(self):
        # define Temperature gap
        tau_p=np.outer(self.factor_p, self.Sigmoidal.y)
        return tau_p
    
    def get_tau_b(self):
        # define Temperature gap
        tau_b=np.outer(self.factor_b, self.Sigmoidal.y)
        return tau_b
    

    
    
    
    
    
    
    
    
    
    
    
    
      