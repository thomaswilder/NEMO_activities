#!/usr/bin/env python3
#
# Function adapted from CDFTOOLS/cdftransport.f90 to compute the zonal transport.
# In addition, cdftransport.py also computes the vertical heat transport.
# n.b. this function may need the addition of bathymetric level `mbkt`.
#
#  !!======================================================================
#  !!                     ***  PROGRAM  cdftransport  ***
#  !!=====================================================================
#  !!  ** Purpose : Compute Transports across a section. 
#  !!               By default, mass (Sv) and  heat(PW)/salt(kT/s) transports
#  !!               are computed unless -noheat option is used (mass 
#  !!               transport only).
#  !!
#  !!  ** Method  : The begining and end point of the section are given in 
#  !!               term of F-points index. A broken line joining successive
#  !!               F-points is defined between the begining and end point
#  !!               of the section. Therefore each segment between F-points
#  !!               is either a zonal or meridional segment corresponding to
#  !!               V or U velocity component. Doing so, the volume conservation
#  !!               is ensured as velocities are not interpolated, and stay 
#  !!               on the native model grid. 
#  !!                 The section name and the begin/end point of a section are
#  !!               read from standard input, till 'EOF' is given as section
#  !!               name. This make possible to give a bunch of sections in 
#  !!               an ASCII files and use the < redirection.
#  !!            SIGN CONVENTION : The transport is positive when the flow cross
#  !!               the section to the right, negative otherwise. This depends
#  !!               on the sense the section is described.  With this convention
#  !!               The algebric sum of transports accross sections forming a 
#  !!               closed area is 0. 
#  !!            OPTIONS :
#  !!               -full   : full step case
#  !!               -noheat : only mass transport is computed.
#  !!               -time   : specify the time frame to be used
#  !!               -zlimit : transports can be computed in different depth layers
#  !!                         defined by their depth limit
#  !!            REQUIREMENT :
#  !!               mesh-mask file are required in the current directory.
#  !!            
#  !!
#  !! History : 2.1  : 01/2005  : J.M. Molines : Original code
#  !!           2.1  : 07/2009  : R. Dussin : add cdf output
#  !!           2.1  : 01/2010  : M.A. Balmaseda : Change integration signs 
#  !!                             so that the transport across a segment is 
#  !!                             independent of the chosen trajectory.
#  !!           3.0  : 04/2011  : J.M. Molines : Doctor norm + Lic.
#  !!         :  4.0  : 03/2017  : J.M. Molines
#  !!----------------------------------------------------------------------
#  !!----------------------------------------------------------------------

from netCDF4 import Dataset
import logging
# import numba
from numpy import zeros

import os
os.environ['HDF5_USE_FILE_LOCKING']='FALSE'

logging.basicConfig(filename='output.log', level=logging.ERROR)

def cdftransport(data_dir, file, var, **kwargs):

  """
  Computes:
  1) the zonal transport of an abitrary section by supplying a U-file and u-var,
  2) the vertical heat transport by also supplying a W-file and w-var.

  Needs associated mesh_mask.nc file in the same data folder

  Inputs:
    data_dir = string for data directory
    file   = string for file with u in or a list of files [U, T, W]
    var    = string for u variable name or a list of variables [u, t, w]

  Optional arguments:
    lprint   = True   print out variable names in netcdf file
    kt       = index of timestep to evaluate
    lheat    = Compute heat transport
    lvert    = Supplying vertical velocity for vertical heat transport

  Returns:
    total zonal transport, vertical heat transport, depth on T-point => (dvoltrpsum, vht, deptht)


  History:
    -        First write up in Python.
    May      Added condition to check if input data is masked array, if so, get data to do numpy computations on.
  """

  w_file = t_file = w_var = t_var = None

  if type(file) == list:

    if len(file) < 2:

      raise ValueError('file lists must have 2 or more elements.')

    # horizontal heat transport if only u_ and t_ files supplied
    elif len(file) == 2 and 'vovecrtz' not in var:
      u_file, t_file = file
      u_var,  t_var  = var

    # horizontal and vertical heat transport if all files supplied
    elif len(file) == 3:
      u_file, t_file, w_file = file
      u_var,  t_var,  w_var  = var

  # this is the default state
  else:
    u_file = file
    u_var  = var


  # define some global variables
  rau0 = 1000.0 # density of pure water (kg/m^3)
  rcp = 4000.0  # heat capacity (J/kg/K)

  # some defaults for optional keyword arguments
  opt_dic = {"kt"     : 0,
             "lprint" : False,
             "lheat"  : False,
             "lvert"  : False}

  # overwrite the options by cycling through the input dictionary
  for key in kwargs:
    opt_dic[key] = kwargs[key]

  # open some files and pull variables out
  cf_ufil = Dataset(data_dir + u_file)
  if opt_dic["lprint"]:
    print(cf_ufil)
  npiglo  = len(cf_ufil.dimensions["x"])
  npjglo  = len(cf_ufil.dimensions["y"])
  npk     = len(cf_ufil.dimensions["depthu"])
  dlu     = cf_ufil.variables[u_var][opt_dic["kt"], :, :, :]
  # is the array a masked array, if so, get data.
  if isinstance(dlu, np.ma.MaskedArray):
    dlu = np.ma.getdata(dlu)
  if opt_dic["lheat"]:
#    t_file = kwargs.get('t_file')
#    t_var = kwargs.get('t_var')
    cf_tfil = Dataset(data_dir + t_file)
    print(cf_tfil)
    dlt     = cf_tfil.variables[t_var][opt_dic["kt"], :, :, :]
    deptht  = cf_tfil.variables["deptht"][:]
    cf_tfil.close()
    if isinstance(dlt, np.ma.MaskedArray):
      dlt = np.ma.getdata(dlt)
  if opt_dic["lvert"]:
#    w_file = kwargs.get('w_file')
#    w_var = kwargs.get('w_var')
    cf_wfil = Dataset(data_dir + w_file)
    dlw     = cf_wfil.variables[w_var][opt_dic["kt"], :, :, :]
    cf_wfil.close()
    if isinstance(dlw, np.ma.MaskedArray):
      dlw = np.ma.getdata(dlw)
  cf_ufil.close()

#  val = dlt[0,100,300]
#  logging.debug(f"a temperature value is {val}")
#  val = dlw[0,100,300]
#  logging.debug(f"a vertical velocity value is {val}")

  mask = 'mesh_mask.nc'
  cn_mask = Dataset(data_dir + mask)
  if opt_dic["lprint"]:
    print(cn_mask)
  e2u     = cn_mask.variables["e2u"][0, :, :]
  e3u_0   = cn_mask.variables["e3u_0"][0, :, :, :]
  gdepw   = cn_mask.variables["gdepw_1d"][0, :]
  if opt_dic["lvert"]:
    e1t     = cn_mask.variables["e1t"][0, :, :]
    e2t     = cn_mask.variables["e2t"][0, :, :]
    e3t_0   = cn_mask.variables["e3t_0"][0, :, :, :]
  cn_mask.close()


  

  # Begin calculations
  # Allocate variables
  dwku = zeros((npk, npjglo, npiglo))
  dtrpu = zeros((npjglo, npiglo))
  if opt_dic["lheat"]:
    dlwt = zeros((npk, npjglo, npiglo))
    dwkwt = zeros((npk, npjglo, npiglo))
    dtrpwt = zeros((npk))

    # compute temperature flux at every grid point
    for jk in range(npk):
      for jj in range(npjglo):
        for ji in range(npiglo):
          if jk == 0:
            dlwt[jk, jj, ji] = 0.0
          elif jk < npk-1:
            dlwt[jk, jj, ji] = dlw[jk, jj, ji] * 0.5*( dlt[jk, jj, ji] + dlt[jk+1, jj, ji] ) * rau0 * rcp
          else:
            dlwt[jk, jj, ji] = 0.0

    val = dlwt[5,100,300]
    logging.debug(f"a temperature flux value is {val}")
        
  # find volume transport at each cell boundary (e.g. dwku = dlu*e2u*e3u)
  # and integrate vertically (dtrpu = dtrpu + dwku)
  for jk in range(npk):
    for jj in range(npjglo):
      for ji in range(npiglo):
        dwku[jk, jj, ji] = dlu[jk, jj, ji] * e2u[jj, ji] * e3u_0[jk, jj, ji]
        dtrpu[jj, ji] = dtrpu[jj, ji] + dwku[jk, jj, ji]
        if opt_dic["lheat"]:
          dwkwt[jk, jj, ji] = dlwt[jk, jj, ji] * e1t[jj, ji] * e2t[jj, ji]
          
  # integrate over domain for vertical heat transport
  if opt_dic["lheat"]:
    for jk in range(npk):
      for jj in range(npjglo):
        for ji in range(npiglo):
          if jk == 0 or jk == npk-1:
            dtrpwt[jk] = 0.0
          else:
            dtrpwt[jk] = dtrpwt[jk] + ( dwkwt[jk, jj, ji] * e3t_0[jk, jj, ji] )



  # choose a section to take the transport at (e.g. ji = 100)
  # dvoltrpsum = dvoltrpsum + dtrpu
  # direction of integration is positive for south to north (line 797 in cdftransport.f90)
  dvoltrpsum = 0.0
  for jj in range(npjglo):
    dvoltrpsum = dvoltrpsum + dtrpu[jj, 100]/1e+6

  if opt_dic["lheat"]:
    return (dvoltrpsum, dtrpwt, deptht)
  else:
    return (dvoltrpsum)
