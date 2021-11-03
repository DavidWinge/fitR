#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 12:21:12 2021

Post-processing analysis of the fitted data

HOW TO RUN:
    python fitR_analyze.py fitR_result.csv
    
OUTPUT:
    Augmented .xlsx files (copy files is standard option)
    Table of results in following format:
        Filename, dir, J0, n, Rs, V0
        first.xlsx dir1, 1.0, 1.0, 1.0, 1.0

TODO: Consider using pythons argparser for a nice user interface and 
--help option. Can also be used to supply settings, like --plot and --verbos

@author: dwinge
"""

import pandas as pd
import numpy as np
import sys

e_const = 1.602e-19 

def calculate_R(L,D,mu,N) :
    area = np.pi*(D/2)**2
    return L/area/N/e_const/1e-4/mu
    
def doping_from_R(L,D,mu,R):
    area = np.pi*(D/2)**2
    return L/area/R/e_const/1e-4/mu

def analyze(filename):
    df = pd.read_csv(filename)
    
    # Read a few parameters from fitWanalyze.config
    Lp,Li,Ln,D,mup,mui,mun,Np,Nn = np.loadtxt('fitR_analyze.config',unpack=True)
    
    # Estimate the resistance values of the n+ and p+ regions
    Rp = calculate_R(Lp,D,mup,Np)  
    Rn = calculate_R(Ln,D,mun,Nn)
    Ri = df['Rs']-Rp-Rn
    
    # Estimate the doping based on the resitance on the i/p-/n- region
    Ni = doping_from_R(Li,D,mui,Ri)
    
    df['Ni'] = Ni
    
    return df

def plot_Ni(df) :
    # Find the different samples
    unique_samples = np.unique(df['Sample'])
    
    # Organize the sample measurment in dictionaries
    samples = {}
    # Save the statistics of each sample in dictionaries
    stats = {}
    for sample in unique_samples :
        samples[sample] = df['Ni'][df['Sample']==sample]
        # Get the mean and standard deviation
        stats[sample] = samples[sample].describe()
        
    # Show mean and standard deviation
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()  
    
    # Tried a few options to visualize the error bars
    ax.errorbar(range(len(samples)),
                [stats[key]['mean'] for key in stats],
                yerr=[stats[key]['std'] for key in stats],
                linestyle='',
                marker='s',
                ms=10,
                elinewidth=2,
                #ecolor='black',
                capsize=6,
                capthick=3)

    # Name the data from the sample names
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(stats.keys())
    
    ax.set_xlabel('Sample label')
    ax.set_ylabel(r'Estimated doping density (1/cm$^3$)')
    
    plt.tight_layout()
    plt.savefig('estimated_mid_segment_doping.png')
    plt.show()      

if __name__ == '__main__' :
    # if file is run as a script
    if len(sys.argv) > 1 :
        df = analyze(sys.argv[1])
        #df = analyze('/home/dwinge/Projects/Solarcells_darkIV/fitWresult.csv')
        plot_Ni(df)
        df.to_csv('fitR_analyzedresult.csv',float_format='%.4g',index=False)
        
    else :
        print('No input given, please provide fitR_result.csv')