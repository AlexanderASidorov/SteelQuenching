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
    def equation_classic(t, b, n):
        """
        base formula for Avrami equation
        """   
        fraction=1-np.exp(-b*t**n)
        return fraction

    @staticmethod
    def equation_classic_tau(t, tau, n):
        """
        base formula for Avrami equation
        """   
        fraction=1-np.exp(-((t/tau)**n))
        return fraction



        
    @staticmethod
    def equation_incubation_tau(t, tau, n, t0):
        """
        base formula for Avrami equation
        """   
        fraction=1-np.exp(-(((t-t0)/tau)**n))
        return fraction
    
       
    @staticmethod
    def equation_incubation_b(t, b, n, t0):
        """
        base formula for Avrami equation
        """   
        fraction=1-np.exp(-b*(t-t0)**n)
        return fraction
    
    
    
    
    @staticmethod
    def equation_derivative(f, K, n, f_eq):
        """
        derivative of Avrami equation with equilibrium fraction
        """      
        fraction_rate=((f_eq-f)/K)*n*(-np.log(1-f/f_eq))**(1-1/n)
        return fraction_rate
    
    @staticmethod
    def set_solving_method(case, t001, t099):
        ### method TRF
        if case == 0:
            bounds=([0., 1.1], [np.inf, 4])
            guess=JMAK.b_n_tau(t001, t099)[0:2]
            method='trf'
        if case == 1:
            bounds=([0., 1.1], [np.inf, 4])
            guess=[10e+20, 4.0]
            method='trf'
        if case == 2:
            bounds=([0., 1.1], [np.inf, 4])
            guess=[1.1, 4.0]
            method='trf'
        if case == 3:
            bounds=([0., 1.1], [np.inf, 4])
            guess=[1.1, 1.1]
            method='trf' 
        if case == 4:
            bounds=([0., 1.1], [np.inf, 4])
            guess=[10e+20, 4.0]
            method='trf'         
        
        ### method LM
        if case == 5:
            bounds=([0., 1.1], [np.inf, 4])
            guess=JMAK.b_n_tau(t001, t099)[0:2]
            method='lm'        
        if case == 6:
            bounds=([0., 1.1], [np.inf, 4])
            guess=[10e+20, 4.0]
            method='lm'
        if case == 7:
            bounds=([0., 1.1], [np.inf, 4])
            guess=[1.1, 4.0]
            method='lm'
        if case == 8:
            bounds=([0., 1.1], [np.inf, 4])
            guess=[1.1, 1.1]
            method='lm'
        if case == 9:
            bounds=([0., 1.1], [np.inf, 4])
            guess=[10e+20, 4.0]
            method='lm'
        else: 
            if case>9:
                print('!!!!!!!!!!case number can be integer 0-9 only')
            
        return bounds, guess, method
    
    @staticmethod
    def x_t_interpalation(t001, t099):
        """
        Creating points inbetween t001 and t099
        """
             
        t050=0.5*(t099-t001)
        t045=0.45*(t099-t001)
        t055=0.55*(t099-t001)
        t010=0.25*(t099-t001)
        t090=0.75*(t099-t001)
        # Create x and y data for curve_fit        
        xdata=np.array([t001, t010, t045, t050, t055, t090, t099])
        ydata=np.array([0.01, 0.1, 0.45, 0.5, 0.55, 0.9, 0.99])
        return xdata, ydata         
    
    
    @staticmethod
    def Fit_b_n(t001, t099, initial_guess=None):
        """
        function that calculates k and n caeficients of JMAK_system for a single
        temperature point based on t start and t end.
        """ 
        
        #######
        # As we don't know the midle points between t001 and t099 we create
        # them artificially
        #########
        xdata=JMAK.x_t_interpalation(t001, t099)[0]
        ydata=JMAK.x_t_interpalation(t001, t099)[1]
        
        def equation_incu_b(t, b, n, t0=t001):
            """
            base formula for Avrami equation
            """   
            fraction=1-np.exp(-b*(t-t0)**n)
            return fraction
        
        Error=1000
        i=-1
        while Error<=0 or Error > 0.5 and i <=8:
            i=i+1
            if i==0 or i==1 or i==2 or i==3 or i==4:
                bound, guess, meth = JMAK.set_solving_method(i, t001, t099)
                if initial_guess is not None:
                    guess=initial_guess
                popt, pcov = curve_fit(equation_incu_b, xdata, ydata, 
                                       p0=guess, bounds=bound, method=meth, 
                                       maxfev=5000)
                b=popt[0]
                n=popt[1]
                errors=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
                Error=errors[3]
                # print('#########')
                # print(pcov)
                # print('#########')
            
                         
              
            if i==5 or i==6 or i==7 or i==8 or i==9:
                bound, guess, meth = JMAK.set_solving_method(i, t001, t099)
                if initial_guess is not None:
                    guess=initial_guess
                popt, pcov = curve_fit(equation_incu_b, xdata, ydata, 
                                       p0=guess, method=meth, maxfev=5000)
                b=popt[0]
                n=popt[1]
                errors=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
                Error=errors[3]
                # print('#########')
                # print(pcov)
                # print('#########')
                
        tau=np.exp(-(np.log(b)/n))  
        return tau, n, Error, guess, b, i

        
    @staticmethod
    def b_n_tau(t001, t099, x1=0.01, x2=0.99):
        """
        ########################################################
        # here we solve system of two logorythmic equations
        # x1=1-exp(-((1/tau)*t001)**n)) and
        # x2=1-exp(-((1/tau)*t099)**n)) and
        #
        # by default x1=0.01 and x2=0.99 whis corresponds to the
        # start and end curves in TTT
        #
        # we can write this equations in the following form
        # x1=1-exp(-(1/tau)**n*t001**n)
        # x2=1-exp(-(1/tau)**n*t099**n)
        # and denominate (1/tau)**n as b
        # so tau = exp(-ln(b)/n)
        # 
        # after it we can write the system in the following way
        # which is a classical look of Avrami equation
        # x1=1-exp(-b*t001**n)
        # x2=1-exp(-b*t099**n)
        #
        # for details of the following handling with this equations
        # see the video of inimitable Taylor Sparks
        # https://www.youtube.com/watch?v=N3KGw5rF6oM
        #
        ########################################################
        """
        nominator = np.log(np.log(1/(1-x1)))-np.log(np.log(1/(1-x2)))
        denominator = np.log(t001)-np.log(t099)
        n=nominator/denominator
        #n=2
        b=np.exp(np.log(np.log(1/(1-x1)))-n*np.log(t001))
        tau=np.exp(-(np.log(b)/n))
        
        return b, n, tau

    @staticmethod
    def b_n_tau_incub(t001, t099, x001=0.01, x099=0.99):
        """
        ########################################################
        # here we solve system of two logorythmic equations
        # x099=1-exp(-((1/tau)*(t099-t001))**n)) and
        # x050=1-exp(-((1/tau)*(t050-t001))**n)) and
        #
        # by default x001=0.01 and x099=0.99 whis corresponds to the
        # start and end curves in TTT and x050 is the middle point with
        # phase fraction = 50%
        #
        # we can write this equations in the following form
        # x099=1-exp(-(1/tau)**n*(t099-t001)**n)
        # x050=1-exp(-(1/tau)**n*(t050-t001)**n)
        # and denominate (1/tau)**n as b
        # so tau = exp(-ln(b)/n)
        # 
        # after it we can write the system in the following way
        # which is close to a classical look of Avrami equation
        # x099=1-exp(-b*(t099-t001)**n)
        # x050=1-exp(-b*(t050-t001)**n)
        #
        ########################################################
        """
        x050=0.5
        t050=(t099-t001)/2
        
        nominator = np.log(np.log(1/(1-x099)))-np.log(np.log(1/(1-x050)))
        denominator = np.log(t099-t001)-np.log(t050-t001)
        n=nominator/denominator
        #n=2
        b=np.exp(np.log(np.log(1/(1-x099)))-n*np.log(t099-t001))
        tau=np.exp(-(np.log(b)/n))
        
        return b, n, tau




# #%%
# t001=10000000
# t099=t001*18.66
# t050=(t099-t001)/2

 
# b, n, tau = JMAK.b_n_tau_incub(t001, t099, x001=0.01, x099=0.99)

            

# # tau, n, Error, guess, b, i= JMAK.Fit_b_n(t001, t099)
# # x001=JMAK.equation_incubation_b(t001, b_, n_, t001)
# # x099=JMAK.equation_incubation_b(t099, b_, n_, t001)
# # x050=JMAK.equation_incubation_b(t050, b_, n_, t001)


# x001_=JMAK.equation_incubation_tau(t001, tau, n, t001)
# x099_=JMAK.equation_incubation_tau(t099, tau, n, t001)
# x050_=JMAK.equation_incubation_tau(t050, tau, n, t001)






