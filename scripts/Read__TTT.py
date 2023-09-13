#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 12:24:16 2023

@author: aas
"""

import numpy as np
import os.path
import pandas as pd
from scipy.optimize import least_squares, curve_fit
import matplotlib.pyplot as plt
from solve_JMAK import JMAK


class Read_TTT ():
    
    def __init__(self, filename='TTT.xlsx', filepath='../data'):
       self.__filename = filename
       self.__filepath = filepath
       ######################################
       ######################################
       self.__dataTTT=self.get_TTT()[0]
       self.__curve_names=self.get_TTT()[1]
       ######################################
       ######################################
       self.__Ferrite=self.get_variables()[0]
       self.__Pearlite=self.get_variables()[1]
       self.__Bainite=self.get_variables()[2]
       ######################################
       ######################################
       self.Ferrite_plus_Pearlite_region=self.extr_comp_temp()
       self.Equilibria=self.get_equilibria()[4]


    
    
    @staticmethod
    def JMAK(t, k, n):
        """
        base formula for Avrami equation
        """   
        fraction=1-np.exp(-1/k*(t**n))
        return fraction
 
    @staticmethod
    def set_method(case):
        if case == 1:
            bounds=([1.1, 1.1], [np.inf, 4])
            guess=[1.1, 4.0]
            method='trf'
        if case == 2:
            bounds=([1.1, 1.1], [np.inf, 4])
            guess=[10e+10, 1.1]
            method='trf'
        if case == 3:
            bounds=([1.1, 1.1], [np.inf, 4])
            guess=[10e+10, 1.1]
            method='lm'
        if case == 4:
            bounds=([1.1, 1.1], [np.inf, 4])
            guess=[1.1, 4]
            method='lm'
        else: print('case number can be 1 ,2, 3 or 4 only')
            
        return bounds, guess, method
 
    @staticmethod
    def Fit_k_n(t001, t099):
        """
        function that calculates k and n caeficients of JMAK_system for a single
        temperature point based on t start and t end.
        """ 
        
        #######
        # As we don't know the midle points between t001 and t099 we create
        # them artificially
        #########
        t050=0.5*(t099-t001)
        t010=0.25*(t099-t001)
        t090=0.75*(t099-t001)
        # Create x and y data for curve_fit        
        xdata=np.array([t001, t010, t050, t090, t099])
        ydata=np.array([0.01, 0.1, 0.5, 0.9, 0.99])
        # Create bounds and guesses for curve_fit
        #p0=np.array(guess)
        #b0=bounds
        i=0
        for i in range (1, 4):
            bounds, guess, method = set_method(i)
        
            popt, pcov = curve_fit(Read_TTT.JMAK, xdata, ydata, 
                               p0=guess, bounds=bounds, method=method, maxfev=5000)
       
        k=popt[0]
        n=popt[1]
        errors=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
        Error=errors[3]
        case=1
        if Error > 0 and Error < 0.01: 
            pass
        # if the error is too big we try another methode
        else:
            bounds, guess, method = Read_TTT.set_methode(2)
            popt, pcov = curve_fit(Read_TTT.JMAK, xdata, ydata, 
                                   p0=guess, bounds=bounds, method=method, maxfev=50000)
            k=popt[0]
            n=popt[1]
            errors=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
            Error=errors[3]
            case=2
            if Error > 0 and Error < 0.01: 
                pass
            # if the error is still too big we try the third one
            else:
                popt, pcov = curve_fit(Read_TTT.JMAK, xdata, ydata, 
                                       p0=[1.1, 1.1], method='lm', maxfev=50000)
                k=popt[0]
                n=popt[1]
                errors=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
                Error=errors[3]
                case=3
                
        if Error > 0.01:
            print('ATTENTION!!!!! Error is too big!!!!!!!')
            print(t001)
            print(t099)
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

        return k, n, Error, case      
 

    
    def get_TTT(self):
        """
        create dataframe from Excel file and a list of curve names
        """
        os.chdir(self.__filepath) # moving to the directory data
        data_TTT=pd.read_excel(self.__filename, index_col=0) #create dataframe with data
        curve_names = data_TTT.columns.values.tolist()
        return data_TTT, curve_names
    
    
    def get_variables(self):
        """
        create a separate dataframe for each phase from the common one
        """
        data_TTT=self.__dataTTT
        curve_names=self.__curve_names
        data={}
        data[curve_names[0]]=data_TTT[curve_names[0]].to_numpy(dtype=float)
        for i in range (1, len(curve_names)):
            name=curve_names[i]
            dataframe=pd.DataFrame (data = [data[curve_names[0]], 
                                            data_TTT[name]]).transpose()
            dataframe=dataframe.dropna(axis=0)
            temp=dataframe.to_numpy(dtype=float)[:, 0]
            time=dataframe.to_numpy(dtype=float)[:, 1]
            data[name + ' temp']=temp
            data[name + ' time']=time
            
        # split the dictionary to separate variables
        # by default we assume that we have Ferrite, Pearlite, Bainite
        # the "resе" variable for the case we have something in addition
        a, b, c, d, e, f, g, h, i, j, k, l, m, *rest = data.values()
        Temperature = a
        temp_f_s, time_f_s, temp_f_e, time_f_e = b, c, d, e
        temp_p_s, time_p_s, temp_p_e, time_p_e = f, g, h, i
        temp_b_s, time_b_s, temp_b_e, time_b_e = j, k, l, m
        
        # add points for better sollution when we define k and n in 
        # JMAK equation
        # 10% point
        coef1=0.25
        
        time_f_01=coef1*(time_f_e-time_f_s)
        time_p_01=coef1*(time_p_e-time_p_s)
        time_b_01=coef1*(time_b_e-time_b_s)
        
        # 50% point
        time_f_05=0.5*(time_f_e-time_f_s)
        time_p_05=0.5*(time_p_e-time_p_s)
        time_b_05=0.5*(time_b_e-time_b_s)
        
        # 90% point
        coef2=0.75
        time_f_09=coef2*(time_f_e-time_f_s)
        time_p_09=coef2*(time_p_e-time_p_s)
        time_b_09=coef2*(time_b_e-time_b_s)    
        
     
         
        ######################################################################
        # Create dataframe with values for Ferrite 
        Ferrite=pd.DataFrame(data=[temp_f_s, time_f_s, time_f_01, time_f_05,
                                   time_f_09, time_f_e]).transpose()
        Ferrite.columns = ['temp', 'time_s', 'time_01', 'time_05', 'time_09','time_e']

        ######################################################################
        # Create dataframe with values for Pearlite
        Pearlite=pd.DataFrame(data=[temp_p_s, time_p_s, time_p_01, time_p_05, 
                                    time_p_09,time_p_e]).transpose()
        Pearlite.columns = ['temp', 'time_s', 'time_01', 'time_05', 'time_09','time_e']
        
        ######################################################################
        # Create dataframe with values for Bainite
        Bainite=pd.DataFrame(data=[temp_b_s, time_b_s, time_b_01, time_b_05, 
                                   time_b_09, time_b_e]).transpose()
        Bainite.columns = ['temp', 'time_s', 'time_01', 'time_05', 'time_09','time_e']    
            
        return Ferrite, Pearlite, Bainite 
    
    def extr_comp_temp(self):
        dataframe=self.__dataTTT
        new_dataframe=pd.DataFrame(columns=dataframe.columns)
        for i in range (len(dataframe.index)):
            raw=dataframe.loc[i, :].values # copy dataframe raw into a single numpy raw
            raw = raw[~np.isnan(raw)] # delete all NaN values
            # chech it there are more then 3 values in a raw (if there are more then
            # three then we have a competition between phase changes)        
            if len (raw)>3:
                new_dataframe.loc[i]=dataframe.loc[i]
        new_dataframe=new_dataframe.dropna(axis='columns')
        return new_dataframe
    
    def get_equilibria(self):
        new_dataframe=self.Ferrite_plus_Pearlite_region
        Temperature=new_dataframe.iloc[:,0].to_numpy()
        #elta=0.1
        t_eq=np.zeros(len(new_dataframe.index))
        f=np.zeros(len(new_dataframe.index))
        f_equilibria=np.ones(len(new_dataframe.index))
        p_equilibria=np.ones(len(new_dataframe.index))
        for i in range (len(new_dataframe.index)):
            tau001_1, tau099_1, tau001_2, tau099_2 = new_dataframe.iloc[i, 1:5]
            k1, n1 = Read_TTT.Fit_k_n (tau001_1, tau099_1)[0:2]
            k2, n2 = Read_TTT.Fit_k_n (tau001_2, tau099_2)[0:2]
            t_s=min(tau001_1, tau001_2)
            t_e=min(tau099_1, tau099_2)
            t_av=t_e-(t_e-t_s)*0.5
            function = lambda t: Read_TTT.JMAK(t, k1, n1) + Read_TTT.JMAK(t, k2, n2)
            f[i]=function(t_av)
            j=0
            while f[i] <= 0.99 or f[i] > 1:
                if f[i]<= 0.99:
                    j=j+1
                    #print(j)
                    t_s=t_av
                    t_e=t_e
                    t_av=t_e-(t_e-t_s)*0.5
                    f[i]=function(t_av)
                    t_eq[i]=t_av
                if f[i]>=1:
                    j=j+1
                    #print(j)
                    t_s=t_s
                    t_e=t_av
                    t_av=t_e-(t_e-t_s)*0.5
                    f[i]=function(t_av)
                    t_eq[i]=t_av
            f_equilibria[i]=Read_TTT.JMAK(t_eq[i], k1, n1)
            p_equilibria[i]=Read_TTT.JMAK(t_eq[i], k2, n2)
        Equilibria=pd.DataFrame(data=[Temperature, f_equilibria, 
                                      p_equilibria]).transpose()
        Equilibria.columns=['Temperature', 'Ferrite', 'Pearlite']
        
        return f_equilibria, p_equilibria, t_eq, f, Equilibria
 
    
 
    
 
    
 
    
 
    
 
    
 