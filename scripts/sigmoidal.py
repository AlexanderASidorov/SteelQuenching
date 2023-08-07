#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:25:54 2023

@author: aas
"""
import numpy as np
from scipy import integrate
from scipy.interpolate import splrep, splev, interp1d
from initial_data import Composition, t, T, delta_x

class Sigmoidal():
    """
    Class for transformtion of sigmoidal function S(x)
    Requires delta_x that is descritisation step for Phase Fracture. Maximum 
    delta_x = 0.1.
    """
    def __init__(self, delta_x=delta_x):
        # integration step
        if delta_x <= 0.1:
            self.delta_x=delta_x
        else: self.delta_x=0.1
        self.x = self.DevidingX ()
        self.xmin=min(self.x)
        self.xmax=max(self.x)
        
        # Integration of S(X)
        self.y=self.Integration()[0]
        self.ymin=self.Integration()[1]
        self.ymax=self.Integration()[2]
        
        # Inversing Integration function
        self.x_from_y=self.Integration_inv()[0]
        self.y_from_x=self.Integration_inv()[1]
        
        # Inversing S(x)
        self.S_inv=self.S_inv()
    
    
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
    
    
    def S_inv(self):
        """
        Invariant of Sigmoidal function S(x)
        """
        y=self.S(self.x)
        x=self.x
        S_inv=interp1d(y, x)
        return S_inv
        
    
    def Integration (self):
        """
        Integration of S(x)
        """
        y=np.zeros(len(self.x))
        for i in range (0, (len(y))):
            y[i] = integrate.quad(self.S, self.x[0]*0.1, self.x[i])[0]
        return y, min(y), max(y)
    
    def Integration_inv(self):
        """
        Invariant of Integration of S(x)
        """
        x_from_y=interp1d(self.y, self.x, kind= 'cubic')
        y_from_x=interp1d(self.x, self.y, kind='cubic')
        return x_from_y, y_from_x
        
    
    
                           