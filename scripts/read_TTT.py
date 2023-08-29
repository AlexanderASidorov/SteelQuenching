#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 19:02:57 2023

@author: aas
"""
import numpy as np
import os.path
import pandas as pd
from scipy.optimize import least_squares, curve_fit
import matplotlib.pyplot as plt
#from read_TTT import *




def JMAK(t, k, n):
    """
    base formula for Avrami equation
    """    
        
    fraction=1-np.exp(-1/k*(t**n))
    return fraction


def Fit_k_n_single (t001, t099):
    """
    function that calculates k and n caeficients of JMAK_system for a single
    temperature point based on t start and t end.
    """ 
    
    bounds=([1.1, 1.1], [np.inf, 4])
    guess=np.array([1.1, 4], dtype=float)
    t050=0.5*(t099-t001)
    t010=0.25*(t099-t001)
    t090=0.75*(t099-t001)
    xdata=np.array([t001, t010, t050, t090, t099])
    ydata=np.array([0.01, 0.1, 0.5, 0.9, 0.99])
    popt, pcov = curve_fit(JMAK, xdata, ydata, p0=guess, bounds=bounds, 
                           method='trf')
    k=popt[0]
    n=popt[1]
    errors=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
    if errors[3] > 0 and errors[3] < 0.01: 
        pass
    else:
        popt, pcov = curve_fit(JMAK, xdata, ydata, p0=[10e+10, 4], method='lm')
        k=popt[0]
        n=popt[1]
        errors=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
    return k, n

#%%

def readTTT(filename):
    """
    create dataframe from Excel file and a list of curve names
    """
    os.chdir('../data') # moving to the directory data
    data_TTT=pd.read_excel(filename, index_col=0) #create dataframe with data
    curve_names = data_TTT.columns.values.tolist()
    return data_TTT, curve_names

# data_TTT=readTTT('TTT.xlsx')[0]
# curve_names=readTTT('TTT.xlsx')[1]

def create_variables(filename):
    """
    create a separate dataframe for each phase from the common one
    """
    
    
    data_TTT=readTTT(filename)[0]
    curve_names=readTTT(filename)[1]
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
        
    return data_TTT, Ferrite, Pearlite, Bainite 
    
data_TTT, Ferrite, Pearlite, Bainite = create_variables('TTT.xlsx')



#%%
# Temp=739
# index01=Ferrite.temp[Ferrite.temp==Temp].index.tolist()[0]
# index02=Pearlite.temp[Pearlite.temp==Temp].index.tolist()[0]
# k1, n1 =  Fit_k_n_single (Ferrite.time_s.iloc[index01], Ferrite.time_e.iloc[index01])
# k2, n2 =  Fit_k_n_single (Pearlite.time_s.iloc[index02], Pearlite.time_e.iloc[index02])
# print(k1, n1)
# print(k2, n2)








def extracting_competitor_temp(dataframe):
    new_dataframe=pd.DataFrame(columns=dataframe.columns)
    for i in range (len(dataframe.index)):
        raw=dataframe.loc[i, :].values # copy dataframe raw into a sinfle numpy raw
        raw = raw[~np.isnan(raw)] # delete all NaN values
        # chech it there are more then 3 values in a raw (if there are more then
        # three then we have a competition between phase changes)        
        if len (raw)>3:
            new_dataframe.loc[i]=dataframe.loc[i]
    new_dataframe=new_dataframe.dropna(axis='columns')
    return new_dataframe
#new_dataframe=extracting_competitor_temp(data_TTT)

def calculate_fraction_equilibria(dataframe):
    new_dataframe=extracting_competitor_temp(dataframe)
    Temperature=new_dataframe.iloc[:,0].to_numpy()
    #elta=0.1
    t_eq=np.zeros(len(new_dataframe.index))
    f=np.zeros(len(new_dataframe.index))
    f_equilibria=np.ones(len(new_dataframe.index))
    p_equilibria=np.ones(len(new_dataframe.index))
    for i in range (len(new_dataframe.index)):
        tau001_1, tau099_1, tau001_2, tau099_2 = new_dataframe.iloc[i, 1:5]
        k1, n1 = Fit_k_n_single (tau001_1, tau099_1)
        k2, n2 = Fit_k_n_single (tau001_2, tau099_2)
        t_s=min(tau001_1, tau001_2)
        t_e=min(tau099_1, tau099_2)
        t_av=t_e-(t_e-t_s)*0.5
        function = lambda t: JMAK(t, k1, n1) + JMAK(t, k2, n2)
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
        f_equilibria[i]=JMAK(t_eq[i], k1, n1)
        p_equilibria[i]=JMAK(t_eq[i], k2, n2)
    Equilibria=pd.DataFrame(data=[Temperature, f_equilibria, 
                                  p_equilibria]).transpose()
    Equilibria.columns=['Temperature', 'Ferrite', 'Pearlite']
    
    return f_equilibria,p_equilibria, t_eq, f, Equilibria

#Equilibria = calculate_fraction_equilibria(data_TTT)

# def fraction_equilibria_at_temp(dataframe, temp):
#     new_dataframe=extracting_competitor_temp(dataframe).reset_index(drop=True)
#     for i in range (len(new_dataframe.Temperature)):
#         if new_dataframe.Temperature[i]==temp:
#             Index=new_dataframe.index[i]
    
#     tau001_1, tau099_1, tau001_2, tau099_2 = new_dataframe.iloc[Index, 1:5]
#     k1, n1 = Fit_k_n_single (tau001_1, tau099_1)
#     k2, n2 = Fit_k_n_single (tau001_2, tau099_2)
#     print(k1, n1)
#     print(k2, n2)
#     t_s=min(tau001_1, tau001_2)
#     t_e=min(tau099_1, tau099_2)
#     t_av=t_e-(t_e-t_s)*0.5
#     function = lambda t: JMAK(t, k1, n1) + JMAK(t, k2, n2)
#     f=function(t_av)
#     j=0
#     while f <= 0.99 or f > 1:
#         if f<= 0.99:
#             j=j+1
#             #print(j)
#             t_s=t_av
#             t_e=t_e
#             t_av=t_e-(t_e-t_s)*0.5
#             f=function(t_av)
#             t_eq=t_av
#         if f>=1:
#             j=j+1
#             #print(j)
#             t_s=t_s
#             t_e=t_av
#             t_av=t_e-(t_e-t_s)*0.5
#             f=function(t_av)
#             t_eq=t_av
#         f_equilibria=JMAK(t_eq, k1, n1)
#         p_equilibria=JMAK(t_eq, k2, n2)
#         Avrami_coef=[k1, k2, n1, n2]
      
#     return new_dataframe, f_equilibria, p_equilibria, t_eq, Avrami_coef

#new_dataframe, f_equilibria, p_equilibria, t_eq, Avrami_coef = fraction_equilibria_at_temp(data_TTT, 739)

#%%





def Fit_k_n(Dataframe):
    """
    function that iterates the whole dataframe and solves JMAK_system for each
    temperature step with fit_curve function.
    input:
    Dataframe = dataframe vaiable that we created via  create_variables
    name - string (name of the phase that we create a variable for)
    """ 
    bounds=([1.1, 1.1], [np.inf, 4])
    time_s=Dataframe['time_s'].to_numpy(dtype=float)
    time_e=Dataframe['time_e'].to_numpy(dtype=float)
    time_05=Dataframe['time_05'].to_numpy(dtype=float)
    time_01=Dataframe['time_01'].to_numpy(dtype=float)
    time_09=Dataframe['time_09'].to_numpy(dtype=float)
    guess=np.full([len(time_s), 2], [1.1, 4], dtype=float)
    k = np.zeros(len(time_s), dtype=float)
    #k=np.float128(k)
    n = np.zeros(len(time_s))
    i=len(time_s)
    errors=np.zeros([len(time_s), 4])
    while i >=1:
        i=i-1
        t001=time_s[i]
        t010=time_01[i]
        t050=time_05[i]
        t090=time_09[i]
        t099=time_e[i]
        xdata=np.array([t001, t010, t050, t090, t099])
        ydata=np.array([0.01, 0.1, 0.5, 0.9, 0.99])
        popt, pcov = curve_fit(JMAK, xdata, ydata, p0=guess[i], bounds=bounds, 
                               method='trf')
        k[i]=popt[0]
        n[i]=popt[1]
        errors[i]=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
        
        # let the initial guess for the next step be the result of the current one
        if i>0: guess[i-1]=[k[i], n[i]]
        else: pass
        
        # the trf mothod can't find the solution for temperature with a long
        # incubation period, so for this case we swithc to 'lm' mothode 
        if errors[i, 3] > 0 and errors[i, 3] < 0.01: 
            pass
        else:
            popt, pcov = curve_fit(JMAK, xdata, ydata, p0=[10e+10, 1.1], method='lm')
            k[i]=popt[0]
            n[i]=popt[1]
            errors[i]=[pcov[0,0], pcov[0,1], pcov[1,0], pcov[1,1]]
        #print(pcov)
    Dataframe['k']=k
    Dataframe['n']=n
    Dataframe['Error']=errors[:, 3]
    Dataframe['Equilibria']=np.ones(len(time_s))
    ######################
    ######################
    ######################
    return Dataframe  #errors #k, n, guess, time_s, time_e

Pearlite=Fit_k_n(Pearlite)
Ferrite=Fit_k_n(Ferrite)
Bainite=Fit_k_n(Bainite) 

def add_equilibria(dataframe, name):
    Equilibria = calculate_fraction_equilibria(data_TTT)[-1]
    if name == 'Pearlite':
        for i in range (len(Equilibria.Temperature)):
            for j in range (len(dataframe.temp)):
                if Equilibria.Temperature[i]==dataframe.temp[j]:
                    dataframe.Equilibria[j]=Equilibria.Pearlite[i]
                else: pass
    if name == 'Ferrite':
         for i in range (len(Equilibria.Temperature)):
             for j in range (len(dataframe.temp)):
                 if Equilibria.Temperature[i]==dataframe.temp[j]:
                     dataframe.Equilibria[j]=Equilibria.Ferrite[i]
                 else: pass
    else: pass
    return dataframe
                    
Pearlite = add_equilibria(Pearlite, 'Pearlite')
Ferrite = add_equilibria(Ferrite,'Ferrite')





def plot_fit_results(Phase, index):
    """
    function that plots the fit Fit_k_n results for the defined Phase and defined 
    index
    """ 
    i=index
    Dataframe=Phase
    time=np.linspace(Dataframe['time_s'][i], Dataframe['time_e'][i], 25)
    X=JMAK(time, Dataframe['k'][i], Dataframe['n'][i])
    X2=np.array([0.01, 0.1, 0.5, 0.9, 0.99])
    t2=np.array([Dataframe['time_s'][i], Dataframe['time_01'][i], 
                 Dataframe['time_05'][i], Dataframe['time_09'][i], 
                 Dataframe['time_e'][i] ])
    fig01=plt.figure()
           
    # Ferrite curves
    plt.plot(time, X, label="X recived")
    #plt.plot(time, X1, label="X arbitrary")
    plt.plot(t2, X2,'*' , label="X from dataframe")
    # Plot settings
    #plt.xscale("log")
    plt.xlabel('Time, sec')
    plt.ylabel('Fraction')
    plt.title('Isothermal curve')
    plt.legend()
    plt.grid()
    return fig01

#Figure01= plot_fit_results(Pearlite, 0)

def plot_two_curves(Phase1, Phase2, Temp):
    """
    function that plots two isothermal curves for the defined temperature. The aim
    is to visialise the case when Ferrite and Pearlit competites each other
    """
    
    # defining indexes for the defind temperature
    index01=Phase1.temp[Phase1.temp==Temp].index.tolist()[0]
    index02=Phase2.temp[Phase2.temp==Temp].index.tolist()[0]
    
    # defining time stretch for the plot
    time1=np.linspace(Phase1['time_s'][index01], Phase1['time_e'][index01], 25)
    time2=np.linspace(Phase2['time_s'][index02], Phase2['time_e'][index02], 25)
    # extracting coefficients
    k1=Phase1.k[index01]
    k2=Phase2.k[index02]
    n1=Phase1.n[index01]
    n2=Phase2.n[index02]
    
    # defining phase fraction via Avrami equation
    X1=JMAK(time1, k1, n1)
    X2=JMAK(time2, k2, n2)
    
    fig01=plt.figure()
           
    # Ferrite curves
    plt.plot(time1, X1, label="Phase 01")
    plt.plot(time2, X2, label="Phase 02")
    # Plot settings
    #plt.xscale("log")
    plt.xlabel('Time, sec')
    plt.ylabel('Fraction')
    plt.title('Isothermal curves')
    plt.legend()
    plt.grid()
    return fig01

#Figure01= plot_two_curves(Pearlite, Ferrite, 739)








    


# #%%

# def JMAK_system(z, t001, t010, t050, t090, t099):
#     """
#     function for solvin system of two Avrami equations
#     the first one for 0.01 phase fraction, the second one for 0.99
#     """    
#     k=z[0]
#     n=z[1]
#     solution = np.empty(5)
#     solution[0]=JMAK(t001, k, n)-0.01
#     solution[1]=JMAK(t010, k, n)-0.1
#     solution[2]=JMAK(t050, k, n)-0.5
#     solution[3]=JMAK(t090, k, n)-0.9
#     solution[4]=JMAK(t099, k, n)-0.99
#     return solution


# def Solve_k_n(Dataframe):
#     """
#     function that iterates the whole dataframe and solves JMAK_system for each
#     temperature step with least_squares. Works not fine at points with a long
#     incubation period (near Ae1, Ae3 and Bs)
#     """      
    
#     bounds=([1.1, 1.1], [np.inf, 4])
#     time_s=Dataframe['time_s'].to_numpy(dtype=float)
#     time_e=Dataframe['time_e'].to_numpy(dtype=float)
#     time_05=Dataframe['time_05'].to_numpy(dtype=float)
#     time_01=Dataframe['time_01'].to_numpy(dtype=float)
#     time_09=Dataframe['time_09'].to_numpy(dtype=float)
#     guess=np.full([len(time_s), 2], [5, 1.1], dtype=float)
#     k = np.zeros(len(time_s), dtype=float)
#     k=np.float128(k)
#     n = np.zeros(len(time_s))
#     i=len(time_s)
#     while i >=1:
#         i=i-1
#         t001=time_s[i]
#         t010=time_01[i]
#         t050=time_05[i]
#         t090=time_09[i]
#         t099=time_e[i]
#         zguess=guess[i]
#         z = least_squares(JMAK_system, zguess, bounds=bounds,  
#                           args=(t001, t010, t050, t090, t099)).x
#         k[i]=z[0]
#         n[i]=z[1]
#         if i>0: guess[i-1, 0]=k[i]
#         else: pass
#     Dataframe['k']=k
#     Dataframe['n']=n
#     return Dataframe #k, n, guess, time_s, time_e
 
# Pearlite=Solve_k_n(Pearlite)
# Ferrite=Solve_k_n(Ferrite)
# Bainite=Solve_k_n(Bainite) 






# i=30
# tau001_1, tau099_1, tau001_2, tau099_2 = new_dataframe.iloc[i, 1:5]
# k1, n1 = Fit_k_n_single (tau001_1, tau099_1)
# k2, n2 = Fit_k_n_single (tau001_2, tau099_2)

# t=tau_eq[i]
# Equilibrium = JMAK(t, k1, n1) + JMAK(t, k2, n2)
# Phase_Ferrite=JMAK(t, k1, n1)
# PHASE_Perlite=JMAK(t, k2, n2)







