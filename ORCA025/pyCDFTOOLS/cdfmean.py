#!/usr/bin/env python3
#
# Function adapted from cdfmean.F90 to compute the horizontal average of fields.
# n.b. this function does not maintain the same level of functionality as cdfmean.F90.

# PROGRAM cdfmean
#  !!======================================================================
#  !!                     ***  PROGRAM  cdfmean  ***
#  !!=====================================================================
#  !!  ** Purpose : Compute the Mean Value over the ocean or part of the
#  !!               ocean (spatial mean).
#  !!
#  !!  ** Method  : mean= sum( V * e1 *e2 * e3 *mask )/ sum( e1 * e2 * e3 *mask ))
#  !!               Partial cell version
#  !!
#  !! History : 2.1  : 10/2005  : J.M. Molines : Original code
#  !!         : 2.1  : 07/2009  : R. Dussin    : Netcdf output
#  !!           3.0  : 01/2011  : J.M. Molines : Doctor norm + Lic.
#  !!         : 4.0  : 03/2017  : J.M. Molines  
#  !!----------------------------------------------------------------------

from netCDF4 import Dataset
import numpy as np
import netCDF4 as nc
import os

def cdfmean(data_dir, file, var, mesh_mask, **kwargs):

    """
    Computes:
        1) the horizontal average of an input field for all vertical levels,
    
    n.b. the input field has 3 dimensions: depth, x, y. 

    Inputs:
        data_dir  = string for data directory
        file      = string for file e.g. U| V| T| W| F
        var       = string for variable name
        mesh_mask = string for mesh mask filename

    Optional arguments:
        lprint   = print out variable names in netcdf file
        kt       = index of timestep to evaluate
        C-point  = choose which point variable lies on c-grid e.g. U| V| T| W| F
        basin    = {"global", "southern_ocean"}

    Returns:
        horizontal average on each z-level -> dhw[t,u,v,w], depth[t,u,v,w]

    History:
        May 2023     Original inception T. Wilder
        Nov 23'      Added time counter to dlt and option of sub-basins

    """

    # some default optional keyword arguments
    opt_dic = {'kt'      : 0,
               'lprint'  : True,
               'C-point' : 'T',
               'basin'   : 'global'} 

    # overwrite the options by cycling through the input dictionary
    for key in kwargs:
        opt_dic[key] = kwargs[key]

    # open files and pull data out
    if opt_dic['C-point'] == 'T':
        cf_tfil = Dataset(data_dir + file)
        if opt_dic["lprint"]:
            print(cf_tfil)
        # npiglo  = len(cf_tfil.dimensions["x"])
        # npjglo  = len(cf_tfil.dimensions["y"])
        npk     = len(cf_tfil.dimensions["deptht"])
        if opt_dic["basin"]=="global":
            dlt     = cf_tfil.variables[var][opt_dic["kt"], :, :, :]
        elif opt_dic["basin"]=="southern_ocean":
            dlt     = cf_tfil.variables[var][opt_dic["kt"], :, 159:549, :]
            print(dlt.shape)
        else:
            raise Exception("neither global or southern_ocean basin has been chosen")
        deptht  = cf_tfil.variables["deptht"][:]
        cf_tfil.close()
        # is the array a masked array, if so, get data.
        if isinstance(dlt, np.ma.MaskedArray):
            dlt = np.ma.getdata(dlt)
    else:
        raise Exception("You must choose a valid C-point for your input file... Or U|V|F need adding.")

    cn_mask = Dataset(mesh_mask)
    if opt_dic['C-point'] == 'T':
        if opt_dic["basin"]=="global":
            e1t   = cn_mask.variables["e1t"][0, :, :]
            e2t   = cn_mask.variables["e2t"][0, :, :]
            tmask = cn_mask.variables["tmask"][0, :, :, :]
        elif opt_dic["basin"]=="southern_ocean":
            e1t   = cn_mask.variables["e1t"][0, 159:549, :]
            e2t   = cn_mask.variables["e2t"][0, 159:549, :]
            tmask = cn_mask.variables["tmask"][0, :, 159:549, :]
            print(e1t.shape)
            print(e2t.shape)
            print(tmask.shape)
    cn_mask.close()

    # Allocate variables
    if opt_dic['C-point'] == 'T':
        e1e2 = np.zeros((npk)) # domain area at each jk level
        dsum = np.zeros((npk)) # integral of all values in domain at each jk level
        dhwt = np.zeros((npk)) # horizontal average at each jk level of a T-point variable
        print(dhwt.shape)

    # Begin computations
    if opt_dic['C-point'] == 'T':
        for jk in range(npk-1):
            e1e2[jk] = np.sum( e1t[:, :] * e2t[:, :]  * tmask[jk, :, :] )
            dsum[jk] = np.sum( dlt[jk, :, :] * e1t[:, :] * e2t[:, :] * tmask[jk, :, :] )
            dhwt[jk] = dsum[jk]/e1e2[jk]

    # if opt_dic['C-point'] == 'T':
    #     return ( dhwt, deptht )
    # else:
    #     print("nothing to output!")
    
    # write data to netcdf file
    base_name, extension = os.path.splitext(file)
    with nc.Dataset(data_dir + base_name + "_horizontal_average_" + var + ".nc", 'w', format='NETCDF4') as dataset:
        # Create a dimension for time
        dataset.createDimension('deptht', len(deptht))
        
        # Create a variable for time
        depth_var = dataset.createVariable('deptht', 'f8', ('deptht',))
        depth_var.units = "m"
        depth_var[:] = deptht
        
        # Create a variable for the data
        data_var = dataset.createVariable(var + "_hor_av", 'f4', ('deptht',))
        data_var.units = "n/a"
        data_var[:] = dhwt
        
    return print("data written to netcdf file")
