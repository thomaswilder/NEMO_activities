#!/usr/bin/env python3
'''
Author: T Wilder
Date: 03/11/2023
Description: This script (calc_cdfmoy.py) runs the python function cdfmoy from pyCDFTOOLS to compute mean NEMO data over a set period.
Method: 
'''

import re
import os
import glob
from pyCDFTOOLS.cdfmoy import cdfmoy

# directory
user_path = "/gws/nopw/j04/terrafirma/twilder/"
data_directory = "u-cz312/data/"

# model and grid specifics
exp = "cz312o"
grid = "U"

# list of years to compute
years = [str(year) for year in range(1977,1983)]
print(years)

# loop through each year
for year in years:
    print(f'starting year {year}')
    # search for files with a specified year in e.g. 1985 and assign them to list
    # year = "1980"
    search_pattern = r"nemo_" + exp + "_1m_" + year + "*" + grid + ".nc"
    
    matching_files_path = glob.glob(os.path.join(user_path + data_directory, search_pattern))
    
    matching_files = []
    for file_path in matching_files_path:
        filename = os.path.basename(file_path)
        matching_files.append(filename)
        print(filename)
    
    print(matching_files)
        
    # create filename cdfmoy_exp_year_grid-*
    cdfmoy_filename = "cdfmoy_" + exp + "_" + year + "_grid-" + grid 
    
    print(cdfmoy_filename)
        
    # compute mean using file list and new filename
    var_name = ["uo", "thkcello"]
    kwargs = {"lprint": False, "grid": "U"}
    output = cdfmoy(user_path + data_directory, matching_files, var_name, cdfmoy_filename, **kwargs)
    





