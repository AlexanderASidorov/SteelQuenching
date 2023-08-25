#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:25:54 2023

@author: aas
"""
import numpy as np
from scipy import integrate
from scipy.interpolate import splrep, splev, interp1d
from initial_data import Composition, t, T, delta_x, delta_T, delta_t, FunctionType

class Sigmoidal():
    """
    Class for transformtion of sigmoidal function S(x)
    Requires delta_x that is descritisation step for Phase Fracture. Maximum 
    delta_x = 0.1.
    """
    def __init__(self, delta_x=delta_x, FunctionType=FunctionType):
        # integration step
        if delta_x <= 0.1:
            self.delta_x=delta_x
        else: self.delta_x=0.1
        self.x = self.DevidingX ()
        self.xmin=min(self.x)
        self.xmax=max(self.x)
        
        # Integration of function S(X) or I(x)
        # (depends on S or I in brackets)
        if FunctionType == 'I(x)':
            self.f=self.I
        else:
            self.f=self.S
        self.y=self.Integration(self.f)[0]
        self.ymin=self.Integration(self.f)[1]
        self.ymax=self.Integration(self.f)[2]
        
        # Inversing Integration function
        self.x_from_y=self.Integration_inv()[0]
        self.y_from_x=self.Integration_inv()[1]
        
        # Inversing S(x)
        self.S_inv=self.inv(self.S)
        # Inversing I(x)
        self.I_inv=self.inv(self.I)    
    
    def DevidingX (self):
        """
        Creating linspace for phase fraction X 
        """
        return np.linspace(0.001, 0.999, int((0.999-0.001)/self.delta_x))
    
    @staticmethod
    def S(x):
        """
        Sigmoidal function S(x)
        """
        return 1./(x**(0.4*(1. - x))*(1. - x)**(0.4*x))
    
    @staticmethod
    def I(x):
        """
        Sigmoidal function I(x)
        """
        return 1./(x**(2.*(1. - x)/3.)*(1. - x)**(2.*x/3.)) 
            
    def inv(self, f):
        """
        Invariant of Sigmoidal function S(x) or I(x)
        """
        y=f(self.x)
        x=self.x
        inv=interp1d(y, x)
        return inv
 
    def Integration (self, f):
        """
        Integration of S(x) or I(x)
        """
        y=np.zeros(len(self.x))
        for i in range (0, (len(y))):
            y[i] = integrate.quad(f, self.x[0]*0.1, self.x[i])[0]
        return y, min(y), max(y)       
    
    def Integration_inv(self):
        """
        Invariant of Integration of S(x) or I(x)
        """
        x_from_y=interp1d(self.y, self.x, kind= 'cubic')
        y_from_x=interp1d(self.x, self.y, kind='cubic')
        return x_from_y, y_from_x
    
    
 
    
 
    # def S_inv(self):
    #     """
    #     Invariant of Sigmoidal function S(x)
    #     """
    #     y=self.S(self.x)
    #     x=self.x
    #     S_inv=interp1d(y, x)
    #     return S_inv
      
    # def I_inv(self):
    #     """
    #     Invariant of Sigmoidal function S(x)
    #     """
    #     y=self.I(self.x)
    #     x=self.x
    #     I_inv=interp1d(y, x)
    #     return I_inv
    
    
    # def Integration (self):
    #     """
    #     Integration of S(x)
    #     """
    #     y=np.zeros(len(self.x))
    #     for i in range (0, (len(y))):
    #         y[i] = integrate.quad(self.S, self.x[0]*0.1, self.x[i])[0]
    #     return y, min(y), max(y)
    
    # def Integration_inv(self):
    #     """
    #     Invariant of Integration of S(x)
    #     """
    #     x_from_y=interp1d(self.y, self.x, kind= 'cubic')
    #     y_from_x=interp1d(self.x, self.y, kind='cubic')
    #     return x_from_y, y_from_x
        
    
    
                           