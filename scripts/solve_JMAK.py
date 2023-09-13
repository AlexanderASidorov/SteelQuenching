#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:51:43 2023

@author: aas
"""

import numpy as np
from scipy.optimize import curve_fit
import math
#from read_data import Create_variables_from_Excel as create

# Pearlite=create()._Pearlite
# Ferrite=create()._Ferrite
# Bainite=create()._Bainite


class JMAK ():
        
    #def __init__(self):
       #self._dataframe = dataframe
   
        
    @staticmethod
    def equation(t, k, n):
        """
        base formula for Avrami equation
        """   
        fraction=1-np.exp(-1/k*(t**n))
        return fraction
    
    @staticmethod
    def set_solving_method(case):
        ### method TRF
        if case == 0:
            bounds=([0.1, 1.1], [np.inf, 4])
            guess=[10e+20, 4.0]
            method='trf'
        if case == 1:
            bounds=([0.1, 1.1], [np.inf, 4])
            guess=[1.1, 4.0]
            method='trf'
        if case == 2:
            bounds=([0.1, 1.1], [np.inf, 4])
            guess=[1.1, 1.1]
            method='trf' 
        if case == 3:
            bounds=([0.1, 1.1], [np.inf, 4])
            guess=[10e+20, 4.0]
            method='trf'         
        
        ### method LM
        if case == 4:
            bounds=([0.1, 1.1], [np.inf, 4])
            guess=[10e+20, 4.0]
            method='lm'
        if case == 5:
            bounds=([0.1, 1.1], [np.inf, 4])
            guess=[1.1, 4.0]
            method='lm'
        if case == 6:
            bounds=([0.1, 1.1], [np.inf, 4])
            guess=[1.1, 1.1]
            method='lm'
        if case == 7:
            bounds=([0.1, 1.1], [np.inf, 4])
            guess=[10e+20, 4.0]
            method='lm'
        else: 
            if case>7:
                print('!!!!!!!!!!case number can be integer 0 -7 only')
            
        return bounds, guess, method
    
    @staticmethod
    def Fit_k_n(t001, t099, initial_guess=None):
        """
        function that calculates k and n caeficients of JMAK_system for a single
        temperature point based on t start and t end.
        """ 
        
        #######
        # As we don't know the midle points between t001 and t099 we create
        # them artificially
        #########
        t050=0.5*(t099-t001)
        t045=0.45*(t099-t001)
        t055=0.55*(t099-t001)
        t010=0.25*(t099-t001)
        t090=0.75*(t099-t001)
        # Create x and y data for curve_fit        
        xdata=np.array([t001, t010, t045, t050, t055, t090, t099])
        ydata=np.array([0.01, 0.1, 0.45, 0.5, 0.55, 0.9, 0.99])
        # Create bounds and guesses for curve_fit
        #p0=np.array(guess)
        #b0=bounds
        # Error=np.ones(4)
        # k=np.zeros(4)
        # n=np.zeros(4)
        Error=1000
        i=-1
        while Error<=0 or Error > 0.5 and i <=6:
            i=i+1
            if i==0 or i==1 or i==2 or i==3:
                bound, guess, meth = JMAK.set_solving_method(i)
                if initial_guess is not None:
                    guess=initial_guess
                popt, pcov = curve_fit(JMAK.equation, xdata, ydata, 
                                       p0=guess, bounds=bound, method=meth, 
                                       maxfev=5000)
                k=popt[0]
                n=popt[1]
                errors=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
                Error=errors[3]
                # print('#########')
                # print(pcov)
                # print('#########')
            
                         
              
            if i==4 or i==5 or i==6 or i==7:
                bound, guess, meth = JMAK.set_solving_method(i)
                if initial_guess is not None:
                    guess=initial_guess
                popt, pcov = curve_fit(JMAK.equation, xdata, ydata, 
                                       p0=guess, method=meth, maxfev=5000)
                k=popt[0]
                n=popt[1]
                errors=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
                Error=errors[3]
                # print('#########')
                # print(pcov)
                # print('#########')
            
        return k, n, Error, guess

    @staticmethod
    def solve_k_n(t001, t099):
        ###################################################
        # here we solve system of two logorythmic equations
        # 0.01=1-exp(-1/k*(t001**n)) and
        # 0.99=1-exp(-1/k*(t099**n)) and
        # for n we solve the equation a**x=b
        ####################################
        b=np.log(1-0.99)/np.log(1-0.01)
        a=t099/t001
        n=math.log(b,a)
        k=1/(-(np.log(1-0.01))/(t001**n))
        return k, n

        
            
    
    
    