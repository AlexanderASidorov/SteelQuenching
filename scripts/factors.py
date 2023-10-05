#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:56:51 2023

@author: aas
"""

import numpy as np
import pandas as pd
from scipy import integrate
from scipy.interpolate import splrep, splev, interp1d
from scipy.optimize import root
from alloy import Alloy
from sigmoidal import Sigmoidal
from initial_data import Composition, t, T, delta_x, delta_T, delta_t




class Factors(Alloy):
        
    Q1=27500*4.184
    Q2=27500*4.184
    Q3=37000*4.184
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
        
                       
        # Factor for each phase using S function
        self.factor_fS=self.get_ferrite_factor_S()[0]
        self.factor_pS=self.get_pearlite_factor_S()[0]
        self.factor_bS=self.get_bainite_factor_S()[0]
        
        # Factor for each phase using I function
        self.factor_fI=self.get_ferrite_factor_I()[0]
        self.factor_pI=self.get_pearlite_factor_I()[0]
        self.factor_bI=self.get_bainite_factor_I()[0]
        
        if self.Sigmoidal.f==self.Sigmoidal.S:
            # if we choose to use S(x) function:
            self.factor_f=self.get_ferrite_factor_S()[0]
            self.factor_p=self.get_pearlite_factor_S()[0]
            self.factor_b=self.get_bainite_factor_S()[0]
            
            self.get_ferrite_factor=self.get_ferrite_factor_S
            self.get_pearlite_factor=self.get_pearlite_factor_S
            self.get_bainite_factor=self.get_bainite_factor_S
        
        if self.Sigmoidal.f==self.Sigmoidal.I:
            # if we choose to use I(x) function:
            self.factor_f=self.get_ferrite_factor_I()[0]
            self.factor_p=self.get_pearlite_factor_I()[0]
            self.factor_b=self.get_bainite_factor_I()[0]        
        
            self.get_ferrite_factor=self.get_ferrite_factor_I
            self.get_pearlite_factor=self.get_pearlite_factor_I
            self.get_bainite_factor=self.get_bainite_factor_I        
        
        # Temperature gap for each phase
        self.temp_f=self.get_ferrite_factor_S()[1]
        self.temp_p=self.get_pearlite_factor_S()[1]
        self.temp_b=self.get_bainite_factor_S()[1]
        
        
        
        # Time necessary for transform into X fraction
        self.tau_f=self.get_tau_f()
        self.tau_p=self.get_tau_p()
        self.tau_b=self.get_tau_b()
    
#############################################################################
##########Functions for the case we use S(x)#################################
#############################################################################
    def get_ferrite_factor_S(self, T=None):
       # define Temperature gap
       Ts=round(self.Alloy.Ae3, 0)
       Te=round(self.Alloy.Bs, 0)
       
       if T is None:
           T=np.arange(Ts-1, Te-1, -self.delta_T)
       
       Numerator=self.Alloy.FC_S()
       FirstTerm=2**(self.Alloy.n1_F*self.Alloy.gs)
       SecondTerm=(Ts - T)**self.Alloy.n2_F
       ThirdTerm=np.exp(-self.Q1/(self.R*(T + self.K)))
       Denumerator=FirstTerm*SecondTerm*ThirdTerm
       return  Numerator/Denumerator, T #, Denumerator 

        
    def get_pearlite_factor_S(self, T=None):
       # define Temperature gap
       Ts=round(self.Alloy.Ae1, 0)
       Te=round(self.Alloy.Bs, 0)
       
       if T is None:
           T=np.arange(Ts-1, Te-1, -self.delta_T)
       else: T=T
       
       Numerator=self.Alloy.PC_S()
       FirstTerm=2**(self.Alloy.n1_P*self.Alloy.gs)
       SecondTerm=(Ts - T)**self.Alloy.n2_P
       ThirdTerm=np.exp(-self.Q1/(self.R*(T + self.K)))
       Denumerator=FirstTerm*SecondTerm*ThirdTerm
       return  Numerator/Denumerator, T
 
    
    def get_bainite_factor_S(self, T=None):
       # define Temperature gap
       Ts=round(self.Alloy.Bs, 0)
       Te=round(self.Alloy.Ms, 0)
       
       if T is None:
           T=np.arange(Ts-1, Te-1, -self.delta_T)

       Numerator=self.Alloy.BC_S()
       FirstTerm=2**(self.Alloy.n1_B*self.Alloy.gs)
       SecondTerm=(Ts - T)**self.Alloy.n2_B
       ThirdTerm=np.exp(-self.Q1/(self.R*(T + self.K)))
       Denumerator=FirstTerm*SecondTerm*ThirdTerm
       return  Numerator/Denumerator, T

#############################################################################
##########Functions for the case we use I(x)#################################
#############################################################################    
    
    def Df(self, T): # see formula 8 in Saunders publiction 
        return np.exp(-self.Q2/(self.R*(T + self.K)))
    
    def Dp(self, T): # see formula 9 in Saunders publiction
        Dp1=1./(np.exp(-self.Q1/(self.R*(T + self.K))))
        Dp2=(0.5*self.Alloy.Mo)/(np.exp(-self.Q3/(self.R*(T + self.K))))
        return (Dp1+Dp2)**(-1)
    
    def Db(self, T): # see formula 8 in Saunders publiction 
        return np.exp(-self.Q1/(self.R*(T + self.K)))
    

    def get_ferrite_factor_I(self, T=None):
       # define Temperature gap
       Ts=round(self.Alloy.Ae3, 0)
       Te=round(self.Alloy.Bs, 0)
       
       if T is None:
           T=np.arange(Ts-1, Te-1, -self.delta_T)
       
       Numerator=self.Alloy.FC_I()
       FirstTerm=6*2**((self.Alloy.gs)/8)
       SecondTerm=(Ts - T)**self.Alloy.n2_F
       ThirdTerm=self.Df(T)
       Denumerator=FirstTerm*SecondTerm*ThirdTerm
       return  Numerator/Denumerator, T #, Denumerator     
    
    def get_pearlite_factor_I(self, T=None):
       # define Temperature gap
       Ts=round(self.Alloy.Ae1, 0)
       Te=round(self.Alloy.Bs, 0)
       
       if T is None:
           T=np.arange(Ts-1, Te-1, -self.delta_T)
       else: T=T
       
       Numerator=self.Alloy.PC_I()
       FirstTerm=6*2**((self.Alloy.gs)/8)
       SecondTerm=(Ts - T)**self.Alloy.n2_P
       ThirdTerm=self.Dp(T)
       Denumerator=FirstTerm*SecondTerm*ThirdTerm
       return  Numerator/Denumerator, T    

    def get_bainite_factor_I(self, T=None):
       # define Temperature gap
       Ts=round(self.Alloy.Bs, 0)
       Te=round(self.Alloy.Ms, 0)
       
       if T is None:
           T=np.arange(Ts-1, Te-1, -self.delta_T)

       Numerator=self.Alloy.BC_I()
       FirstTerm=6*2**((self.Alloy.gs)/8)
       SecondTerm=(Ts - T)**self.Alloy.n2_B
       ThirdTerm=self.Db(T)
       Denumerator=FirstTerm*SecondTerm*ThirdTerm
       return  Numerator/Denumerator, T    
    
   
    # def get_factor(self, Ts, Te, Coef_comp, n1, n2, gs):
    #     """
    #     Function for calculation factor for users's defined temperatures and
    #     phases
    #     """
    #     T=np.arange(Ts-1, Te-1, -self.delta_T)
    #     Numerator=Coef_comp
    #     FirstTerm=2**(n1*gs)
    #     SecondTerm=(Ts - T)**n2
    #     ThirdTerm=np.exp(-self.Q/(self.R*(T + self.K)))
    #     Denumerator=FirstTerm*SecondTerm*ThirdTerm
    #     return  Numerator/Denumerator
        
    
    
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
    

    
# #%%
# Factors_=Factors()
# i=50
# f=Sigmoidal.S
# tau001_01=Factors_.tau_p[i, 0]
# tau099_01=Factors_.tau_p[i, -1]
# ratio_01=tau099_01/tau001_01



 
# Factor_Pearlite=Factors_.factor_pS[i]
# Current_Temperature=Factors_.temp_p[i]
# tau001_02=Factor_Pearlite*Sigmoidal.Integral(f, 0, 0.01)
# tau099_02=Factor_Pearlite*Sigmoidal.Integral(f, 0, 0.99)
# ratio_02=tau099_02/tau001_02


# Sigmoidal_=Sigmoidal()
# tau099_=Factor_Pearlite*Sigmoidal_.y[-1]
# tau099_outer=(np.outer(Factor_Pearlite, Sigmoidal_.y))[-1]  
# tau099_outer=tau099_outer[-1]
    
    
    
    
    
    
    
      