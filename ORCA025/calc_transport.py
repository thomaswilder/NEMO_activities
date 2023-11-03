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
from pyCDFTOOLS.cdftransport import cdftransport

# def write_to_netcdf(filename, time, data, variable_name, units):
#     """
#     Write data to a NetCDF file

#     Inputs:
        
#     """
#     with nc.Dataset(filename, 'w', format='NETCDF4') as dataset:
#         # Create a dimension for time
#         dataset.createDimension('time', len(time))
        
#         # Create a variable for time
#         time_var = dataset.createVariable('time', 'f8', ('time',))
#         time_var.units = 'seconds since 1970-01-01 00:00:00'
#         time_var[:] = time
        
#         # Create a variable for the time series data
#         data_var = dataset.createVariable(variable_name, 'f4', ('time',))
#         data_var.units = units
#         data_var[:] = data
    

# search directory for cdfmoy filenames
# directory
user_path = "/gws/nopw/j04/terrafirma/twilder/"
data_directory = "u-cz312/data/"

# model and grid specifics
# year = "1977"
exp = "cz312o_"
grid = "_grid-U"

# search all years
search_pattern = r"cdfmoy_" + exp + "*" + grid + ".nc"
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
for i in range(len(matching_files)):
    transport[i] = cdftransport(user_path + data_directory, matching_files[i], var_name, **kwargs)
print(transport)

# # save to netcdf file
# dates = [datetime(1976,2,1) + n*relativedelta(months=+1) for n in range(len(files))]
# times[:] = date2num(dates,"hours since 1800-01-01")





