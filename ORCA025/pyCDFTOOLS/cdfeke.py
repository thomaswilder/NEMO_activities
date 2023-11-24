#!/usr/bin/env python3

# A adapted cdfeke ...

# PROGRAM cdfeke
#   !!======================================================================
#   !!                     ***  PROGRAM cdfeke   ***
#   !!=====================================================================
#   !!  ** Purpose : Compute Eddy Kinetic Energy
#   !!
#   !!  ** Method  : Use gridU gridU2, gridV gridV2 files produced by
#   !!               cdfmoy. Velocities are interpolated both on T points
#   !!               and the variance is computed. If -mke option is used
#   !!               the program also outputs MKE field
#   !!
#   !! History : pre  : 11/2004  : J.M. Molines : Original code
#   !!           2.1  : 04/2005  : J.M. Molines : use modules
#   !!           3.0  : 12/2010  : J.M. Molines : Doctor norm + Lic.
#   !!         : 4.0  : 03/2017  : J.M. Molines
#   !!----------------------------------------------------------------------

from netCDF4 import Dataset
import numpy as np


def cdfeke(data_dir, filenames, var_names, mesh_mask, **kwargs):

    '''
    Computes:
        1) eddy and mean kinetic energy

    Needs associated mesh_mask.nc file in working directory.

    Inputs:
        data_dir   = string for data directory
        filenames  = string for file e.g. U| V| T| W| F
        var_names  = string for variable name
        mesh_mask  = string for mesh mask filename

    Optional arguments:
        lprint   = print out variable names in netcdf file
        kt       = index of timestep to evaluate
        eke/mke  = choose whether to calculate eke or mke
        cdfmoy   = if True, each variable has their own file
        o        = output filename, default set to eke.nc

    Returns:
        a netcdf file with eke/mke

    History:
        2023      Initial inception
        Nov 23'   Added default output filename
    '''

    
    # some defaults for optional keyword arguments
    opt_dic = {"kt": 0, 
               "lprint": False, 
               "eke": True, 
               "mke": False,
               "cdfmoy": False,
               "o"     : "eke.nc"}

    # overwrite the options by cycling through the input dictionary
    for key in kwargs:
        opt_dic[key] = kwargs[key]

    if type(filenames) == list:
        
        if opt_dic["cdfmoy"]:

            if len(filenames) < 2:
                raise ValueError("Error: file list must have at least 2 elements!")
    
            elif len(filenames) == 2:  # mke
                u_file, v_file = filenames
                u_var, v_var = var_names
    
            elif len(filenames) == 3:
                raise ValueError("file list cannot have 3 elements!")
    
            elif len(filenames) == 4:  # eke
                u_file, v_file, u2_file, v2_file = filenames
                u_var, v_var, u2_var, v2_var = var_names
                
        else:
            if len(filenames) < 2:
                raise ValueError("Error: file list must have at least 2 elements!")
            
            elif len(filenames) == 2 and len(var_names) == 4: # eke
                u_file, v_file = filenames
                u_var, v_var, u2_var, v2_var = var_names
            
            elif len(filenames) == 2 and len(var_names) == 2: # mke
                u_file, v_file = filenames
                u_var, v_var = var_names
            
    else:
        raise ValueError("Error: input files should be type cdfmoy or nemo output.")

    # import masks from mesh_mask file
    cn_mask = Dataset(mesh_mask)
    umask = cn_mask.variables["umask"][0, :, :, :]
    vmask = cn_mask.variables["vmask"][0, :, :, :]
    cn_mask.close()

    if opt_dic["cdfmoy"]:
        # open some files and pull variables out
        cf_ufil = Dataset(data_dir + u_file)
        if opt_dic["lprint"]:
            print(cf_ufil)
        # npiglo = len(cf_ufil.dimensions["x"])
        # npjglo = len(cf_ufil.dimensions["y"])
        npk = len(cf_ufil.dimensions["depthu"])
        nav_lon = cf_ufil.variables["nav_lon"][:, :]
        nav_lat = cf_ufil.variables["nav_lat"][:, :]
        uc = cf_ufil.variables[u_var][opt_dic["kt"],:, :, :]
        depthu = cf_ufil.variables["depthu"][:]
        # replace masked values with zero
        if np.ma.is_masked(uc):
            uc = np.ma.filled(uc,0)
        npjglo = np.size(uc, 1)
        npiglo = np.size(uc, 2)
        cf_ufil.close()
        
        cf_u2fil = Dataset(data_dir + u2_file)
        if opt_dic["lprint"]:
            print(cf_u2fil)
        u2 = cf_u2fil.variables[u2_var][opt_dic["kt"],:, :, :]
        # replace masked values with zero
        if np.ma.is_masked(u2):
            u2 = np.ma.filled(u2,0)
        cf_u2fil.close()
        
        cf_vfil = Dataset(data_dir + v_file)
        if opt_dic["lprint"]:
            print(cf_vfil)
        vc = cf_vfil.variables[v_var][opt_dic["kt"],:, :, :]
        # replace masked values with zero
        if np.ma.is_masked(vc):
            vc = np.ma.filled(vc,0)
        cf_vfil.close()
        
        cf_v2fil = Dataset(data_dir + v2_file)
        if opt_dic["lprint"]:
            print(cf_v2fil)
        v2 = cf_v2fil.variables[v2_var][opt_dic["kt"],:, :, :]
        # replace masked values with zero
        if np.ma.is_masked(v2):
            v2 = np.ma.filled(v2,0)
        cf_v2fil.close()
        
    else:
        # open some files and pull variables out
        cf_ufil = Dataset(data_dir + u_file)
        if opt_dic["lprint"]:
            print(cf_ufil)
        # npiglo = len(cf_ufil.dimensions["x"])
        # npjglo = len(cf_ufil.dimensions["y"])
        npk = len(cf_ufil.dimensions["depthu"])
        nav_lon = cf_ufil.variables["nav_lon"][:, :]
        nav_lat = cf_ufil.variables["nav_lat"][:, :]
        uc = cf_ufil.variables[u_var][opt_dic["kt"],:, :, :]
        # replace masked values with zero        
        if np.ma.is_masked(uc):
            uc = np.ma.filled(uc,0)
        if opt_dic["eke"]:
            u2 =cf_ufil.variables[u2_var][opt_dic["kt"],:,:,:]
            # replace masked values with zero
            if np.ma.is_masked(u2):
                u2 = np.ma.filled(u2,0)
        depthu = cf_ufil.variables["depthu"][:]
        npjglo = np.size(uc, 1)
        npiglo = np.size(uc, 2)
        cf_ufil.close()
        
        cf_vfil = Dataset(data_dir + v_file)
        if opt_dic["lprint"]:
            print(cf_vfil)
        vc = cf_vfil.variables[v_var][opt_dic["kt"],:, :, :]
        # replace masked values with zero
        if np.ma.is_masked(vc):
            vc = np.ma.filled(vc,0)
        if opt_dic["eke"]:
            v2 =cf_vfil.variables[v2_var][opt_dic["kt"],:,:,:]
            # replace masked values with zero
            if np.ma.is_masked(v2):
                v2 = np.ma.filled(v2,0)
        cf_vfil.close()
        
    # # mask the data
    # uc = uc * umask
    # vc = vc * vmask
    # u2 = u2 * umask
    # v2 = v2 * vmask
        

    # begin calculation
    if opt_dic["eke"]:
        eke = np.zeros((1, npk, npjglo, npiglo))
        for jk in range(npk):
            for jj in range(1, npjglo):
                for ji in range(1, npiglo):
                    value = 0.5 * (0.5 * (
                        (u2[jk, jj, ji] - uc[jk, jj, ji] * uc[jk, jj, ji])
                        + (u2[jk, jj, ji - 1] - uc[jk, jj, ji - 1] * uc[jk, jj, ji - 1])
                    ) + 0.5 * (
                        (v2[jk, jj, ji] - vc[jk, jj, ji] * vc[jk, jj, ji])
                        + (v2[jk, jj - 1, ji] - vc[jk, jj - 1, ji] * vc[jk, jj - 1, ji])
                    ))
                    eke[0,jk,jj,ji] = value


    if opt_dic["mke"]:
        rmke = np.zeros((1, npk, npjglo, npiglo))
        for jk in range(npk):
            for jj in range(1, npjglo):
                for ji in range(1, npiglo):
                    value = 0.5 * (
                        0.5
                        * (
                            uc[jk, jj, ji] * uc[jk, jj, ji]
                            + uc[jk, jj, ji - 1] * uc[jk, jj, ji - 1]
                        )
                        + 0.5
                        * (
                            vc[jk, jj, ji] * vc[jk, jj, ji]
                            + vc[jk, jj - 1, ji] * vc[jk, jj - 1, ji]
                        )
                    )
                    rmke[0,jk,jj,ji] = value
                    
    print(data_dir + opt_dic["o"])              

    # create netcdf file
    netcdf_file = Dataset(data_dir + opt_dic["o"], "w")

    # Create dimensions for the netCDF file.
    netcdf_file.createDimension("time_counter", 1)
    netcdf_file.createDimension("deptht", npk)
    netcdf_file.createDimension("y", npjglo)
    netcdf_file.createDimension("x", npiglo)

    # Create variables for the netCDF file.
    if opt_dic["eke"]:
        eke_variable = netcdf_file.createVariable("eke", "f4", ("time_counter", "deptht", "y", "x"))
    if opt_dic["mke"]:
        rmke_variable = netcdf_file.createVariable("rmke", "f4", ("time_counter", "deptht", "y", "x"))
    nav_lat_variable = netcdf_file.createVariable("nav_lat", "f4", ("y", "x"))
    nav_lon_variable = netcdf_file.createVariable("nav_lon", "f4", ("y", "x"))
    deptht = netcdf_file.createVariable("deptht", "f4", ("deptht"))

    # Write the data to the netCDF file.
    if opt_dic["eke"]:
        eke_variable[:] = eke
    if opt_dic["mke"]:
        rmke_variable[:] = rmke
    nav_lat_variable[:] = nav_lat
    nav_lon_variable[:] = nav_lon
    deptht[:] = depthu # U and T depth levels are the same.

    # Close the netCDF file.
    netcdf_file.close()
    
    
    return print('data written to netcdf')
