fitR package 
**************************************
Version: 0.1, 2021-11-03
Author: David Winge, david.winge@sljus.lu.se
**************************************
Description:

A collection of python scripts that fits a modified form 
of the diode equation to pin-diode IV data. It takes
as input .xlsx files, adds an additional column containing
the fitted function values.
If used with many measurment files at the same time, 
it outputs a fitR_results.csv file that contains the fitted 
parameters for each input file. This can be used 
for subsequent analysis by the fitR_analyze.py script,
as an example.

The package contains example data in a folder structure
that is representative for measurements on different
samples. 

The equation that the script tries to fit is given as
J(V)=J0*(1/alpha * W(alpha*exp(q(V-V0+J0*Rs)/(nkT)))-1)
where 
alpha=qJ0*Rs/(nkT)
and J0,n,Rs,V0 are the fit parameters saturation current,
ideality factor, series resistance and voltage offset, 
respectively. W is the Lambert function evaluated on its 
principle branch. 

The J(V) function is the solution of the Schockley diode 
equation in the presence of a series resistance, where J(V)
has been isolated on the left hand side. 

In the fitR_analyze.py script, the doping density is linked 
to the resistance in each segment using
Rs=1/(N*e*mu) L/A
where N is the doping density, mu the mobility and L and A
the lenght and area of the segment, respectively.

Dependiencies: 

numpy, scipy, pandas, xlsxwriter

How to run:

Scripts are run either from terminal (see examples below) or from a IDE of 
some sort (very minor changes to the code are necessary in that case).

Example usage of the provided scripts on the example 
data folders in the repository:

Single file:
    python fitR_file.py 13557/AX13557_NW1_dark.xlsx
    # outputs 13557/AX13557_NW1_dark_fitR.xlsx as default
    # outputs 13557/AX13557_NW1_dark_fitR.png as an option
    
Several directories:
    python fitR_directories.py 13557 13558 13559
    # outputs fitR_result.csv

Analyze result:
    python fitR_analyze.py fitR_result.csv
    # outputs fitR_analyzedresult.csv
    # outputs estimated_mid_segment_doping.png
    
Config files needed:

fitR_file.py and fitR_directories.py need fitR_file.config in the run-directory.
fitR_analyze.py needs fitR_analyze.config in the run-directory. 
For sample-specific instructions, put a copy of fitR_analyze.config in the sample 
directory with the specific parameters. The script will look for it.

Generated output:

fitR_file.py augments an .xlsx file with the fitted function values.
It has an option to produce a plot of the fitted data.
fitR_directories.py outputs a fitR_result.csv file containing the fit
parameters for all files. 
fitR_analyze.py produces an output file fitR_analyzedresult.csv with an
estimate of the doping in the i-region of the pin-diode. It produces a
plot of the result.

Features to be added:

diode_local: should be able to catch exceptions if they occur.

help text: each script should have a --help option to give a 
short documentation.
