#!/usr/bin/env python3
'''
Author: T Wilder
Date: 03/11/2023
Description: This script (calc_transport.py) 
Method:
'''

import os
import glob
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from cftime import date2num, num2date
from pyCDFTOOLS.cdftransport import cdftransport

def write_to_netcdf(filename, time, data, variable_name, units):
    """
    Write data to a NetCDF file.

    Inputs: 
        - filename
        - time (in units of "hours since 1800-01-01")
        - data
        - variable name and its units

    Outputs:
        - netcdf file of transport for chosen experiment
        
    """
    try:
        with nc.Dataset(filename, 'w', format='NETCDF4') as dataset:
            # Create a dimension for time
            dataset.createDimension('time', len(time))
            
            # Create a variable for time
            time_var = dataset.createVariable('time', 'f8', ('time',))
            time_var.units = "hours since 1800-01-01"
            time_var[:] = time
            
            # Create a variable for the time series data
            data_var = dataset.createVariable(variable_name, 'f4', ('time',))
            data_var.units = units
            data_var[:] = data

    except Exception as e:
        return f'Error: {str(e)}'

    

# search directory for cdfmoy filenames
# directory
user_path = "/gws/nopw/j04/terrafirma/twilder/"
data_directory = "u-da643/data/"

# model and grid specifics
# year = "1977"
exp = "da643o"
grid = "_grid-U"

# search all years
search_pattern = r"cdfmoy_" + exp + "_*" + grid + ".nc"
matching_files_path = glob.glob(os.path.join(user_path + data_directory, search_pattern))

matching_files = []    
for file_path in matching_files_path:
    filename = os.path.basename(file_path)
    matching_files.append(filename)
    print(filename)

# list of years to compute
# years = [str(year) for year in range(1977,1979)]

print(matching_files)

# compute transport
var_name = "uo"
kwargs = {"kt": 0, "lprint": False, "drake": True}

transport = np.zeros(len(matching_files))
# transport = np.zeros(3)
# for i in range(3):
for i in range(len(matching_files)):
    transport[i] = cdftransport(user_path + data_directory, matching_files[i], var_name, **kwargs)
print(transport)

# save to netcdf file

# set up date coordinates
dates = [datetime(1977,6,1) + n*relativedelta(years=+1) for n in range(len(matching_files))]
# dates = [datetime(1977,6,1) + n*relativedelta(years=+1) for n in range(3)]
time = date2num(dates,"hours since 1800-01-01")
print(dates)
print(date2num(dates,"hours since 1800-01-01"))     


write_to_netcdf("transport_" + exp + ".nc", time, transport, "acc_transport", "Sv")











