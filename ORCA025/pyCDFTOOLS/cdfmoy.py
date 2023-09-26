#!/usr/bin/env python3
# 
# An adapted version of cdfmoy that employs xarray dask to compute the mean
# and mean of squared velocities
#
# PROGRAM cdfmoy
#   !!======================================================================
#   !!                     ***  PROGRAM  cdfmoy  ***
#   !!=====================================================================
#   !!  ** Purpose : Compute mean values for all the variables in a bunch
#   !!               of cdf files given as argument
#   !!               Store the results on a 'similar' cdf file.
#   !!
#   !!  ** Method  : Also store the mean squared values for the nn_sqdvar
#   !!               variables belonging to cn_sqdvar(:), than can be changed 
#   !!               in the nam_cdf_names namelist if wished.
#   !!               Optionally order 3 moments for some variables can be
#   !!               computed.
#   !!
#   !! History : 2.0  : 11/2004  : J.M. Molines : Original code
#   !!         : 2.1  : 06/2007  : P. Mathiot   : Modif for forcing fields
#   !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
#   !!                  04/2015  : S. Leroux    : add nomissincl option
#   !!         : 4.0  : 03/2017  : J.M. Molines  


import xarray as xr
import os
import re

def cdfmoy(data_directory, filename, var_name, new_filename, **kwargs):
    """
    Computes:
    1) mean of a supplied variable
    2) mean of a variable squared

    Inputs:
      data_dir   = string for data directory
      filename   = string of a single or list of filenames on the same grid
      var_name   = string for variable name
      new_filename = a filename that could be of form cdfmoy*

    Optional arguments:
      lprint   = True   print out xarray dataset
      chunk    = chunk size for xarray-dask (not working yet)
      surf     = compute only surface values (not implemented yet)
      grid     = which grid, either U, V. (T not set up yet)

    Returns:
      two netcdf files cdfmoy*grid-U.nc and cdfmoy*grid-U2.nc


    History:
      Sep 2023     First write up in Python.
    """
    
    # some defaults for optional keyword arguments
    opt_dic = {"lprint": False, "chunk": 6, "surf": False, "grid": "U"}

    # overwrite the options by cycling through the input dictionary
    for key in kwargs:
        opt_dic[key] = kwargs[key]
    
    # read in data
    if type(filename) == list: # filename is a list of files with single/multiple time indices
        
        if opt_dic["grid"] == "U":
            
            common_prefix = os.path.commonprefix(filename)
            files_to_read = os.path.join(data_directory, common_prefix + '*grid-U.nc') # we assume file format
            
        elif opt_dic["grid"] == "V":
            
            common_prefix = os.path.commonprefix(filename)
            files_to_read = os.path.join(data_directory, common_prefix + '*grid-V.nc') # we assume file format
            
        ds = xr.open_mfdataset(files_to_read, parallel=True, chunks={"time_counter": 1})
        
        if opt_dic["lprint"]:
            print(ds)
        
    elif len(filename) < 2: # filename is one file with multiple time indices
        
        ds = xr.open_dataset(data_directory + filename, chunks={"time_counter": "chunk"})
        if opt_dic["lprint"]:
            print(ds)
    

    print(var_name)
    # mean velocity data
    var = ds[var_name].mean(dim="time_counter")
    print('mean computed')
    # print(var)

    # square velocity data and mean
    var_sqd = (ds[var_name] ** 2).mean(dim="time_counter")
    var_sqd = var_sqd.rename(var_name + "2")
    
    # write data to netcdf
    var.to_netcdf(data_directory + new_filename + ".nc")
    var_sqd.to_netcdf(data_directory + new_filename + "2.nc")
    
    return print('variables saved to netcdf')