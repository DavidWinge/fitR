#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 10:07:59 2020

This is the diode class. 
Updated and cleaned up to package in the fitWdark project 28 Oct 2021.

@author: dwinge
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import lambertw

class Diode:
    """
    Class holding and fitting data to different types of diode equations.
    1) Fit a log-slope to the ideal diode equation
    2) Fit data to a modified diode equation including a series resistance
    
    """
    
    k_boltzmann = 8.617e-5 # eV/K
    kTroom = k_boltzmann*300 # eV
    
    def __init__(self, J0=1e-10, n=1.0, kT=kTroom, area=1.0):
        """
        Constructor with optional parameters. Please make sure that you use 
        consistent units for say, current and area to match your physical 
        situation.

        Parameters
        ----------
        J0 : float, optional
            Saturation current. The default is 1e-10.
        n : float, optional
            Ideality factor. The default is 1.0.
        kT : float, optional
            Thermal voltage. The default is kTroom.
        area : float, optional
            Area convert current to current density. The default is 1.0.

        Returns
        -------
        None.

        """
        self.J0 = J0
        self.n = n
        self.kT = kT
        self.area = area
        
    def add_data(self,V,J) :
        """
        Add data to the Diode class.

        Parameters
        ----------
        V : array
            Voltage values.
        J : array
            Current values.

        Returns
        -------
        None.

        """
        # sort the data as a security measure
        sort_idxs = np.argsort(V)
        V = V[sort_idxs]
        J = J[sort_idxs]
        self.v_space = V
        self.j_space = J
        
    def add_data_I(self,V,I,area) :
        """
        Add data and scale by area to current density

        Parameters
        ----------
        V : array
            Voltage values.
        I : array
            Current values.
        area : float
            Device area.

        Returns
        -------
        None.

        """
        self.add_data(V,I/area)
        
    def diode_eq(self,V,J0,n) :
        """
        Ideal (Schockley) diode equation

        Parameters
        ----------
        V : float, array
            Voltage valeues, independent variable.
        J0 : float
            Saturation current.
        n : float
            Ideality factor.

        Returns
        -------
        float, array
            Diode current.

        """
        return J0 * (np.exp(V/n/self.kT)-1.)
    
    def log_diode_eq(self,V,J0,n) :
        """
        Logarithmic version of diode_eq

        """
        return np.log(abs(J0)) + np.log(abs(np.exp(V/n/self.kT)-1.))
                         
    def diode_Rseries_eq(self,V,J0,n,R,V0) :
        """
        Modified diode equation to include a series resistance R as well as a
        constant bias drop V0.
        
        Parameters
        ----------
        V : float, array
            Voltage valeues, independent variable.
        J0 : float
            Saturation current.
        n : float
            Ideality factor.
        R : float
            Series resistance.
        V0 : float
            Constant bias drop.

        Returns
        -------
        float, array
            Diode current.

        """
        J0=abs(J0)
        alpha = J0*R/n/self.kT
        beta = (V-V0)/n/self.kT
        z = alpha*np.exp(alpha+beta)
        # z should be strictly positive which allows us to use the principal 
        # branch of the Lambert function (no optional arguments to lambertw)
        y = (1./alpha)*lambertw(z)-1.
        return J0*y.real # if complex, use real part
    
    def diode_Rseries_log_eq(self,V,J0,n,R,V0) :
        """
        Logarithmic version of diode_Rseries_eq
        
        """
        J0=abs(J0)
        alpha = J0*R/n/self.kT
        beta = (V-V0)/n/self.kT
        z = alpha*np.exp(alpha+beta)
        # z is strictly positive which allows us to use the principal branch
        # of the Lambert function
        y = (1./alpha)*lambertw(z)-1.
        return np.log(J0*y.real)
           
    def fit_diode_Rseries(self,J0_guess=1e-16,n_guess=2.0,R_guess=1e5,V0_guess=0.4,
                   Vmin=None,Vmax=None, fit_log=False):
        """
        Fits data to the function diode_Rseries_eq. 

        Parameters
        ----------
        J0_guess : float, optional
            Starting guess. The default is 1e-16.
        n_guess : float, optional
            Starting guess. The default is 2.0.
        R_guess : float, optional
            Starting guess. The default is 1e5.
        V0_guess : float, optional
            Starting guess. The default is 0.4.
        Vmin : float, optional
            Lowest volrage point to include in fit. The default is None.
        Vmax : float, optional
            Highest voltage point to include in fit. The default is None.
        fit_log : bool, optional
            Fit the logarithmic data. The default is False.

        Returns
        -------
        popt : array
            Values for best fit.
        V_fit : array
            Voltage values used in the fit.

        """
        # create mask for fitting the data
        if (Vmin and Vmax) :
            v_fit_mask = ((self.v_space < Vmax) * (self.v_space > Vmin))
        elif (Vmin) :
            v_fit_mask = (self.v_space > Vmin)
        elif (Vmax) :
            v_fit_mask = (self.v_space < Vmax)
        else :
            v_fit_mask = (self.v_space != None)
                
        # initial guess for the parameters
        p0 = np.array([J0_guess,n_guess,R_guess,V0_guess])
        
        if fit_log :
            fit_func = self.diode_Rseries_log_eq
            ydata = np.log(self.j_space[v_fit_mask])
        else :  
            fit_func = self.diode_Rseries_eq()
            ydata = self.j_space[v_fit_mask]
                
        popt, pcov = curve_fit(fit_func, 
                               self.v_space[v_fit_mask], 
                               ydata, 
                               p0=p0,
                               maxfev=5000)                     
        
        return popt, self.v_space[v_fit_mask]
        
    def fit_diode(self,J0_guess=1e-16,n_guess=2.0,Vmin=None,Vmax=None):
        """
        Fits data to the ideal diode function diode_eq.

        Parameters
        ----------
        J0_guess : float, optional
            Starting guess. The default is 1e-16.
        n_guess : float, optional
            Starting guess. The default is 2.0.
        Vmin : float, optional
            Lowest volrage point to include in fit. The default is None.
        Vmax : float, optional
            Highest voltage point to include in fit. The default is None.

        Returns
        -------
        popt : array
            Values for best fit.
        V_fit : array
            Voltage values used in the fit.

        """
        # create mask for fitting the data
        if (Vmin and Vmax) :
            v_fit_mask = ((self.v_space < Vmax) * (self.v_space > Vmin))
        elif (Vmin) :
            v_fit_mask = (self.v_space > Vmin)
        elif (Vmax) :
            v_fit_mask = (self.v_space < Vmax)
        else :
            v_fit_mask = (self.v_space != None)
                
            
        p0 = np.array([J0_guess,n_guess])
        popt, pcov = curve_fit(self.log_diode_eq,
                               self.v_space[v_fit_mask],
                               np.log(abs(self.j_space[v_fit_mask])),
                               p0)
        
        return popt, self.v_space[v_fit_mask]
                    
    def get_ideal_slope(self,Vmin,Vmax) :
        """
        Use a line to calculate the slope between Vmin and Vmax

        Parameters
        ----------
        Vmin : float
            Lowest voltage to include.
        Vmax : float
            Highest value to include.

        Returns
        -------
        n, float
            Resulting ideality factor.
        line_seg : array
            Two points: (Vmin,Jmin) and (Vmax, Jmax).

        """
        # Find the point corresponding to Vmin
        for idx, val in enumerate(self.v_space) :
            if (val > Vmin) :
                # call by data explicitly here
                Vmin_idx = idx
                break
        # Find the point corresponding to Vmax
        for idx, val in enumerate(self.v_space) :    
            if (val > Vmax) :
                Vmax_idx = idx-1
                break
        # Now get the slope of the two points
        slope = (np.log10(self.j_space[Vmax_idx])-np.log10(self.j_space[Vmin_idx]))
        slope = slope / (self.v_space[Vmax_idx]-self.v_space[Vmin_idx])
        
        line_seg = np.array([[self.v_space[Vmin_idx],self.v_space[Vmax_idx]],
                             [self.j_space[Vmin_idx],self.j_space[Vmax_idx]]])
        
        return 1.0/slope/59.6e-3, line_seg
    
