#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 10:47:22 2021

Top-level script to control the analysis script fitR_file.py
Takes a list of directories as input, searches those directories for .xlsx files,
fits an IV curve which is appended to the existing data.

A table of the resulting parameters is written in the top level where the 
script is run.

HOW TO RUN:
    python fitR_directories.py dir1 dir2 ... dir3
    
OUTPUT:
    Augmented .xlsx files (copy files is standard option)
    Table of results in following format:
        Filename, dir, J0, n, Rs, V0
        first.xlsx dir1, 1.0, 1.0, 1.0, 1.0

TODO: Consider using pythons argparser for a nice user interface and 
--help option. Can also be used to supply settings, like --plot and --verbos

@author: dwinge
"""

import sys
import os
import fitR_file
import pandas as pd
import numpy as np

def search_dirs(basedir, datadirs) :
    file_list = []
    path_list = []
    folder_names = []
    for d in datadirs :
        #print(os.listdir(d))
        for file in os.listdir(d) :
            print(file)
            if file.endswith('.xlsx'):
                # This works for both Windows and Linux systems
                path_list.append(os.path.join(basedir,d,file))
                file_list.append(file)
                folder_names.append(d)
            
    return file_list, path_list, folder_names

if __name__ == '__main__':
    # if this file is executed as a script
    # Where are we?
    basedir = os.getcwd()
  
    if len(sys.argv) > 1 :
        # save both names and paths
        file_list, path_list, folder_names = search_dirs(basedir, sys.argv[1:])

        # Fit the data by looping over files and save best fit parameters
        res = np.zeros((len(path_list),4))
    
        for k, file in enumerate(path_list) :
            res[k] = fitR_file.fitR(file)

        df = pd.DataFrame(data=res,index=file_list, columns=['J0','n','Rs','V0'])
        # append also a the folder name as a Sample label
        df.insert(0,'Sample',folder_names)
    
        # Write to .csv file
        df.to_csv('fitR_result.csv', sep=',', float_format='%.4g',index_label='Measurement')
    
    else :
        print('No input arguments read, please provide the folder names')
