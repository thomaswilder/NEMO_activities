#!/usr/bin/env python3
#
# Function adapted from cdftools/cdfzonalmean.f90 to compute the zonal mean
#
# !!======================================================================
# !!                     ***  PROGRAM  cdfzonalmean  ***
# !!=====================================================================
# !!  ** Purpose : Compute the zonal mean of a file
# !!
# !!  ** Method  : In this program the 'zonal' mean is in fact a mean 
# !!               along the I coordinate. 
# !!
# !! History : 2.1  : 11/2005  : J.M. Molines : Original code
# !!                : 06/2007  : P. Mathiot   : adaptation for 2D files
# !!           3.0  : 05/2011  : J.M. Molines : Doctor norm + Lic.
# !!         : 4.0  : 03/2017  : J.M. Molines  
# !!----------------------------------------------------------------------

import netCDF4 as nc
import numpy as np
import math


def cdfzonalmean(data_dir, file, mesh_mask, var, **kwargs):
    '''
    Computes:
        1) the zonal mean along I coordinates for a chosen variable

    Needs associated mesh_mask.nc file as input.

    Inputs:
        data_dir = string for data directory
        file     = string for file e.g. U| V| T| W| F
        var      = string for variable name

    Optional arguments:
        lprint   = print out variable names in netcdf file
        kt       = index of timestep to evaluate
        C-point  = choose which point variable lies on c-grid e.g. U| V| T| W| F
        basin    = {"glo", "so"} - glo = global, so = southern ocean
        o        = output filename, default set to zonal_mean.nc

    Returns:
        a netcdf file with a zonal mean at each latitude coordinate and depth level

    History:
        Nov 23'      Initial inception, configured for southern ocean and T point

    '''


    # some default optional keyword arguments
    opt_dic = {'kt'      : 0,
               'lprint'  : False,
               'C-point' : 'T',
               'basin'   : 'glo',
               'o'       : 'zonal_mean.nc'} 

    # overwrite the options by cycling through the input dictionary
    for key in kwargs:
        opt_dic[key] = kwargs[key]
    
    # pull data from mask  
    with nc.Dataset(mesh_mask, 'r') as cn_mask:
        if opt_dic["lprint"]:
            print(cn_mask)
        if opt_dic["basin"]=="so":
            if opt_dic['C-point'] == 'T':
                e1t   = cn_mask.variables["e1t"][0, 159:549, :]
                e2t   = cn_mask.variables["e2t"][0, 159:549, :]
                zmask = cn_mask.variables["tmask"][0, :, 159:549, :]
            nav_lat = cn_mask.variables['nav_lat'][159:549,400]
        elif opt_dic["basin"]=="glo":
            if opt_dic['C-point'] == 'T':
                e1t   = cn_mask.variables["e1t"][0, 159:1000, :]
                e2t   = cn_mask.variables["e2t"][0, 159:1000, :]
                zmask = cn_mask.variables["tmask"][0, :, 159:1000, :]
            nav_lat = cn_mask.variables['nav_lat'][159:1000,400]
        else:
            raise Exception("no basin chosen for zonal averaging")
    
    # get data from input file
    with nc.Dataset(data_dir + file, 'r') as cf_tfil:
        if opt_dic["lprint"]:
            print(cf_tfil)
        npi  = len(cf_tfil.dimensions["x"])
        if opt_dic["basin"]=="so":
            npj  = len(range(159,549))
        elif opt_dic["basin"]=="glo":
            npj  = len(range(159,1000))
            # print(npj)
        else:
            raise Exception("no basin chosen for zonal averaging")
        npk     = len(cf_tfil.dimensions["deptht"])
        npt     = len(cf_tfil.dimensions["time_counter"])
        if opt_dic["basin"]=="so":
            zv      = cf_tfil.variables[var][opt_dic["kt"],:,159:549,:]
        elif opt_dic["basin"]=="glo":
            zv      = cf_tfil.variables[var][opt_dic["kt"],:,159:1000,:]
        # replace masked values with zero
        if np.ma.is_masked(zv):
            zv = np.ma.filled(zv,0)
            # print(zv)
        if opt_dic['C-point'] == 'T':
            depth = cf_tfil.variables['deptht'][:]

    # area of cell
    if opt_dic['C-point'] == 'T':
        dl_surf = e1t * e2t
    
    # initilise variables
    dzomean = np.zeros((npj,npk))
    darea = np.zeros((npj,npk))
    
    # main computing loop
    for jk in range(npk):
        for jj in range(npj):
            darea[jj, jk] = np.sum(dl_surf[jj, :] * zmask[jk, jj, :])
            dzomean[jj, jk] = np.sum(dl_surf[jj, :] * (zmask[jk, jj, :] * zv[jk, jj, :]))
            # assign zero values as nan for plotting purposes
            if dzomean[jj, jk] == 0:
                dzomean[jj, jk] = math.nan
            # do computations if area is not zero
            if darea[jj, jk] != 0:
                dzomean[jj, jk] = dzomean[jj, jk] / darea[jj, jk]

    
    # assign to new variable for writing to netcdf file
    zo = np.zeros((1,npk,npj,1))
    for jk in range(npk):
        for jj in range(npj):
            zo[0, jk, jj, 0] = dzomean[jj, jk]
    
    # write data to netcdf
    with nc.Dataset(opt_dic['o'], 'w', format='NETCDF4') as dataset:
        # Create a dimension
        dataset.createDimension("time_counter", 1)
        dataset.createDimension("depth", npk)
        dataset.createDimension("y", npj)
        dataset.createDimension("x", 1)
        
        # Create dimension variables
        depth_var = dataset.createVariable('depth', 'f4', ('depth',))
        depth_var.units = "m"
        depth_var[:] = depth
        nav_lat_variable = dataset.createVariable("nav_lat", "f4", ("y", "x"))
        nav_lat_variable[:] = nav_lat
        
        # Create a variable for the data
        data_var = dataset.createVariable(f"zo_{var}_{opt_dic['basin']}", 'f4', ("time_counter", "depth", "y", "x"))
        data_var.units = "n/a"
        data_var[:] = zo
    
        
    return print('data written to netcdf')