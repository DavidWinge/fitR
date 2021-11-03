#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 15:03:42 2021

Fits the data in an .xlsx file to a modified diode equation. 
Generates output as a .xlsx file, either a new copy or augmented input file.

RUN AS SCRIPT:
    python fitR_file.py filename.xlsx

TODO: Consider using pythons argparser for a nice user interface and 
--help option. Can also be used to supply settings, like --plot and --verbose

Dependiencies: numpy, scipy, pandas, xlsxwriter

@author: dwinge
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import diode_local
import sys

def fitR(filename):
    # Read the config file
    Vheader, Iheader, J0, n0, R0, V0, Vstart, Vflip, make_copy, make_plot = np.genfromtxt('fitR_file.config',dtype='str')
    
    # Convert certain input from string format to floats
    parameter_guess = [J0,n0,R0,V0]
    parameter_guess = np.array([float(x) for x in parameter_guess])
    Vstart = float(Vstart)
    
    # Read the xlsx file using pandas, use this for more than one sheet
    # excel_reader = pd.ExcelFile(filename)  
    # Parse the excel into a DataFrame, simple way for a single sheet 
    df = pd.read_excel(filename)
    
    # Structure on how to do this for all sheets
    # =============================================================================
    # for sheet in excel_reader.sheet_names:
    #     sheet_df = excel_reader.parse(sheet)
    #     append_df = to_update.get(sheet)
    # 
    #     if append_df is not None:
    #         sheet_df = pd.concat([sheet_df, append_df], axis=1)
    # 
    #     sheet_df.to_excel(excel_writer, sheet, index=False)
    # =============================================================================
    
    # Preprocessing, find out if a sign change is necessary for the voltage column
    def sum_derivative(a) :
        sum = 0
        for k in range(1,len(a)) :
            sum += a[k]-a[k-1]    
        return sum
    
    # Create a copy of the original data frame:
    df_org = df.copy(deep=True)
    # This is True if voltage should be flipped
    negative_voltage = sum_derivative(df_org[Vheader].to_numpy())<0
    # Flip current and voltage
    if negative_voltage :
        df[Vheader] = df[Vheader].apply(lambda x : -x)
        df[Iheader] = df[Iheader].apply(lambda x : -x)
        
    # Use the diode class to analyze the data
    my_diode = diode_local.Diode()
    my_diode.add_data(df[Vheader].to_numpy(), df[Iheader].to_numpy())
    # Fit to the log of the data, this is much more accurate for J0, n and V0
    log_opt, V_fit = my_diode.fit_diode_Rseries(*parameter_guess,Vmin=Vstart, fit_log=True)
    
    J = log_opt[0]
    n = log_opt[1]
    R = log_opt[2]
    V = log_opt[3]
    
    print('*********************************\n' + 
          f'File: {filename}\n' +
          'Result from the fit:\n' +
          f'Voltage and current flipped: {Vflip}\n'+
          'Fitted values:\n' +
          f'J0={J:.4g} A \n' +
          f'n={n:.2f} \n' + 
          f'Rseries={R:.4g} Ohm \n'+
          f'V0={V:.2f} V')
          
    # Plot the result
    if make_plot.lower()=='true' :
        figurename = filename[:-5]+'_fitR.png'
        fig, ax = plt.subplots()
        ax.plot(df[Vheader],df[Iheader].apply(abs),label='data')
        ax.plot(df[Vheader],my_diode.diode_Rseries_eq(df[Vheader].to_numpy(), J, n, R, V), label='best fit')
        ax.set_ylabel('Current (A)')
        ax.set_xlabel('Voltage (V)')
        ax.set_yscale('log')
        plt.tight_layout()
        plt.savefig(figurename)

    # Construct new DataFrame to append
    # Fill the unfitted part with zeros
    data_append = np.zeros_like(df[Vheader].to_numpy())
    data_append[-len(V_fit):] = my_diode.diode_Rseries_eq(V_fit, J, n, R, V)
    # Name the header
    Fheader = f'J0={J:.4g} A, n={n:.2f}, Rs={R:.4g} Ohm, V0={V:.2f} V'
    # Create new DataFrame
    df_append = pd.DataFrame(data_append,columns=[Fheader])
    
    # Append to existning Excel sheet
    if Vflip.lower()=='true' :
        sheet_df = pd.concat([df, df_append], axis=1)
    else :  
        sheet_df = pd.concat([df_org, df_append], axis=1)
            
    # use xlsxwrite in order to not destroy formatting
    if make_copy.lower()=='true' :
        filename_write = filename[:-5] + '_fitR.xlsx'
    else:
        filename_write = filename
            
    excel_writer = pd.ExcelWriter(filename_write, engine='xlsxwriter')
    sheet_df.to_excel(excel_writer, startrow=1, header=False, index=False)
    
    # Get the xlsxwriter workbook and worksheet objects.
    workbook  = excel_writer.book
    worksheet = excel_writer.sheets['Sheet1']
    
    # Add a header format.
    header_format = workbook.add_format({
        'bold': False,
        'text_wrap': False,
        'border': 0})
    
    # Write the column headers with the defined format.
    for col_num, value in enumerate(sheet_df.columns.values):
        worksheet.write(0, col_num + 0, value, header_format)
        
    excel_writer.save()
    
    return J,n,R,V
    
if __name__ == "__main__" :
    # if file is run as script
    filename = sys.argv[-1]
    #filename = '/home/dwinge/Projects/Solarcells_darkIV/13557/AX13557_NW1_dark.xlsx'
    res = fitR(filename)

    
    
