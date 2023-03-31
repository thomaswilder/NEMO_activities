# A modified version of `gen_nemo_unagi_fields.py`. March 2023.
# Changes made to temperature profile.
# -------------------------------------------------------
# 30 Jan 2019
# generates the NEMO horizontal grid and bathymetry file

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from copy import deepcopy
import xarray as xr

#%%
# create the output data class

def gen_bathy_meter(jpiglo, jpjglo, e1, e2, bathy, config, filename):

    # open a new netCDF file for writing.
    ncfile = Dataset(filename, "w", format = "NETCDF4") 
    ncfile.title = "homebrew bathymetry file for %s" % config

    # create the dimensions.
    ncfile.createDimension("x", jpiglo)
    ncfile.createDimension("y", jpjglo)
    ncfile.createDimension("t", None)

    # first argument is name of variable, 
    # second is datatype,
    # third is a tuple with the names of dimensions.

    lon_netcdf = ncfile.createVariable("nav_lon", np.dtype("float32").char, ("x"), fill_value = False)
    lon_netcdf[:] = e1
    lon_netcdf.units = "m"
    lon_netcdf.long_name = "x"

    lat_netcdf = ncfile.createVariable("nav_lat", np.dtype("float32").char, ("y"), fill_value = False)
    lat_netcdf[:] = e2
    lat_netcdf.units = "m"
    lat_netcdf.long_name = "y"

    time_netcdf = ncfile.createVariable("time_counter", np.dtype("float32").char, ("t"), fill_value = False)
    time_netcdf[:] = 0.0
    time_netcdf.units = "s"

    bathy_netcdf = ncfile.createVariable("Bathymetry", np.dtype("float64").char, ("t", "y", "x"), fill_value = False)
    bathy_netcdf[:] = bathy
    bathy_netcdf.units = "m"

    # close the file.
    ncfile.close()

    print("*** SUCCESS writing example file %s" % filename)
    
def gen_forcing(jpiglo, jpjglo, lonT, latT, utau, vtau, qtot, qsr, emp, filename):

    # open a new netCDF file for writing.
    ncfile = Dataset(filename , "w", format = "NETCDF4") 
    ncfile.title = "homebrew forcing file for modEEL"

    # create the dimensions.
    ncfile.createDimension("x", jpiglo)
    ncfile.createDimension("y", jpjglo)
    ncfile.createDimension("t", None)

    # first argument is name of variable, 
    # second is datatype,
    # third is a tuple with the names of dimensions.

    lon_netcdf = ncfile.createVariable("nav_lon", np.dtype("float32").char, ("y", "x"), fill_value = False)
    lon_netcdf[:] = lonT
    lon_netcdf.units = "m"
    lon_netcdf.long_name = "x"

    lat_netcdf = ncfile.createVariable("nav_lat", np.dtype("float32").char, ("y", "x"), fill_value = False)
    lat_netcdf[:] = latT
    lat_netcdf.units = "m"
    lat_netcdf.long_name = "y"

    utau_netcdf = ncfile.createVariable("utau", np.dtype("float64").char, ("t", "y", "x"), fill_value = False)
    utau_netcdf[:] = utau
    utau_netcdf.units = "N m-2"

    vtau_netcdf = ncfile.createVariable("vtau", np.dtype("float64").char, ("t", "y", "x"), fill_value = False)
    vtau_netcdf[:] = vtau
    vtau_netcdf.units = "N m-2"

    qtot_netcdf = ncfile.createVariable("qtot", np.dtype("float64").char, ("t", "y", "x"), fill_value = False)
    qtot_netcdf[:] = qtot
    qtot_netcdf.units = "W m-2"

    qsr_netcdf = ncfile.createVariable("qsr", np.dtype("float64").char, ("t", "y", "x"), fill_value = False)
    qsr_netcdf[:] = qsr
    qsr_netcdf.units = "W m-2"

    emp_netcdf = ncfile.createVariable("emp", np.dtype("float64").char, ("t", "y", "x"), fill_value = False)
    emp_netcdf[:] = emp
    emp_netcdf.units = "m s-1"

    # close the file.
    ncfile.close()

    print("*** SUCCESS writing example file %s!" % filename)
    
def gen_istate(x, y, z, toce, soce, sst, sss, filename):

    # open a new netCDF file for writing.
    ncfile = Dataset(filename, "w", format = "NETCDF4") 
    ncfile.title = "homebrew initial state file for modEEL"

    # create the dimensions.
    ncfile.createDimension("x", jpiglo)
    ncfile.createDimension("y", jpjglo)
    ncfile.createDimension("z", jpkglo)
    ncfile.createDimension("t", None)

    # first argument is name of variable, 
    # second is datatype,
    # third is a tuple with the names of dimensions.

    lon_netcdf = ncfile.createVariable("nav_lon", np.dtype("float32").char, ("y", "x"), fill_value = False)
    lon_netcdf[:] = x
    lon_netcdf.units = "m"
    lon_netcdf.long_name = "x"

    lat_netcdf = ncfile.createVariable("nav_lat", np.dtype("float32").char, ("y", "x"), fill_value = False)
    lat_netcdf[:] = y
    lat_netcdf.units = "m"
    lat_netcdf.long_name = "y"

    lat_netcdf = ncfile.createVariable("nav_lev", np.dtype("float32").char, ("z"), fill_value = False)
    lat_netcdf[:] = z
    lat_netcdf.units = "m"
    lat_netcdf.long_name = "z"

    toce_netcdf = ncfile.createVariable("toce", np.dtype("float64").char, ("t", "z", "y", "x"), fill_value = False)
    toce_netcdf[:] = toce
    toce_netcdf.units = "C"

    soce_netcdf = ncfile.createVariable("soce", np.dtype("float64").char, ("t", "z", "y", "x"), fill_value = False)
    soce_netcdf[:] = soce
    soce_netcdf.units = "g kg-1"

    sst_netcdf = ncfile.createVariable("sst", np.dtype("float64").char, ("t", "y", "x"), fill_value = False)
    sst_netcdf[:] = sst
    sst_netcdf.units = "C"

    sss_netcdf = ncfile.createVariable("sss", np.dtype("float64").char, ("t", "y", "x"), fill_value = False)
    sss_netcdf[:] = sss
    sss_netcdf.units = "g kg-1"

    # close the file.
    ncfile.close()

    print("*** SUCCESS writing example file %s!" % filename)

#%%
# EEL configuration, flat bottom + mound

jpiglo = 83
jpjglo = 242

reso_m = 2000.0
e1 = reso_m * np.arange(0, jpiglo)
e2 = reso_m * np.arange(0, jpjglo)
Lx, Ly = e1[-1], e2[-1]

# EEL has 4000m depth, also put two walls in
bathy = 4000.0 * np.ones((1, jpjglo, jpiglo))

xx, yy = np.meshgrid(e1, e2)
mound = 1500.0 * np.exp(   -(  ( (yy - Ly / 2.0) ** 2 + (xx - Lx / 2.0) ** 2 ) / 35e3 ** 2  ))
plt.subplot(1, 2, 1)
plt.contourf(xx, yy, mound)
plt.colorbar()
plt.axis("equal")
plt.tight_layout()

bathy[0, :, :] -= mound
plt.subplot(1, 2, 2)
plt.contourf(xx, yy, bathy[0, :, :])
plt.colorbar()
plt.axis("equal")
plt.tight_layout()

bathy[0,  0, :] = 0.0
bathy[0, -1, :] = 0.0

gen_bathy_meter(jpiglo, jpjglo, e1, e2, bathy, "EEL", "bathy_meter_EEL.nc")


#%%
# UNAGI configuration, based on SO channel of Dave Munday
# 9000km long, 2400km wide, 3000m deep

# relevant numbers
# res  [nx short]  nx    ny      nz
# 100   40         90    24 + 2  30 + 1
#  50   80        180    48 + 2  30 + 1
#  25  160        360    96 + 2  30 + 1
#  15             600   160 + 2  30 + 1
#  10  400        900   240 + 2  30 + 1

# bathy_filename = "bathy_meter_UNAGI_R050.nc"
bathy_filename = "bathy_meter.nc"

jpiglo = 90
jpjglo = 24 + 2 # add to grid points on
jpkglo = 30 + 1 # add the bottom level in
reso_m = 100.0e3
e1 = reso_m * np.arange(0, jpiglo)
e2 = reso_m * np.arange(0, jpjglo)
Lx, Ly = e1[-1], e2[-1]

# UNAGI has 3000m depth, also put two walls in
oce_depth = 3000.0
bathy = oce_depth * np.ones((1, jpjglo, jpiglo))

l_ridge = True
ridge_H = 1500.0
ridge_L = 500.0e3 # half width of ridge

if l_ridge:
    ridge_x1d = np.zeros(jpiglo)
    for ji in range(jpiglo):
        if (e1[ji] > -ridge_L - 100.0e3 + Lx / 2.0) & (e1[ji] < -100.0e3 + Lx / 2.0):
            ridge_x1d[ji] = 0.5 * ridge_H * (1.0 + np.cos(np.pi * (e1[ji] + Lx / 2.0 + 100.0e3) / ridge_L))
        elif (e1[ji] > 100.0e3 + Lx / 2.0) & (e1[ji] < ridge_L + 100.0e3 + Lx / 2.0):
            ridge_x1d[ji] = 0.5 * ridge_H * (1.0 + np.cos(np.pi * (e1[ji] + Lx / 2.0 - 100.0e3) / ridge_L))
        elif (e1[ji] >= -100.0e3 + Lx / 2.0) & (e1[ji] <= 100.0e3 + Lx / 2.0):
            ridge_x1d[ji] = ridge_H
else:
    ridge_x1d = oce_depth * np.zeros(jpiglo)

bathy[0, :, :] -= ridge_x1d[np.newaxis, :]
bathy[0,  0, :] = 0.0
bathy[0, -1, :] = 0.0
    
plt.contourf(e1 / 1e3, e2 / 1e3, bathy[0, :, :])
plt.colorbar()

gen_bathy_meter(jpiglo, jpjglo, e1, e2, bathy, "UNAGI", bathy_filename)

#%% Flat bottom single layer channel model for intial testing stage

# bathy_filename = "bathy_meter_UNAGI_R050.nc"
bathy_filename = "IDEAL_bump/bathy_meter.nc"

jpiglo = 600 # 1000 km in long
jpjglo = 200 + 2 # add to grid points on # 2000 km latitude
jpkglo = 30 + 1 # add the bottom level in
reso_m = 10.0e3 # 10 km resolution
e1 = reso_m * np.arange(0, jpiglo)
e2 = reso_m * np.arange(0, jpjglo)
Lx, Ly = e1[-1], e2[-1]

# UNAGI has 3000m depth, also put two walls in
oce_depth = 3000.0
bathy = oce_depth * np.ones((1, jpjglo, jpiglo))

l_ridge = True
ridge_H = 1500.0 # no ridge for barotropic setup
ridge_L = 1000.0e3 # half width of ridge

if l_ridge:
    ridge_x1d = np.zeros(jpiglo)
    for ji in range(jpiglo):
        if (e1[ji] > -ridge_L - 100.0e3 + Lx / 2.0) & (e1[ji] < -100.0e3 + Lx / 2.0):
            ridge_x1d[ji] = 0.5 * ridge_H * (1.0 + np.cos(np.pi * (e1[ji] + Lx / 2.0 + 100.0e3) / ridge_L))
        elif (e1[ji] > 100.0e3 + Lx / 2.0) & (e1[ji] < ridge_L + 100.0e3 + Lx / 2.0):
            ridge_x1d[ji] = 0.5 * ridge_H * (1.0 + np.cos(np.pi * (e1[ji] + Lx / 2.0 - 100.0e3) / ridge_L))
        elif (e1[ji] >= -100.0e3 + Lx / 2.0) & (e1[ji] <= 100.0e3 + Lx / 2.0):
            ridge_x1d[ji] = ridge_H
else:
    ridge_x1d = oce_depth * np.zeros(jpiglo)

bathy[0, :, :] -= ridge_x1d[np.newaxis, :]
bathy[0,  0, :] = 0.0
bathy[0, -1, :] = 0.0

plt.contourf(e1 / 1e3, e2 / 1e3, bathy[0, :, :],16)
plt.colorbar()

gen_bathy_meter(jpiglo, jpjglo, e1, e2, bathy, "IDEAL", bathy_filename)


#%%
# use the generated bathy_meter.nc in DOMAINcfg to generate a domaincfg file
# the use that domaincfg file to generate the initial state and forcing files

# no meridional wind here so just put everything on T points
filename = "IDEAL_bump/domain_cfg.nc"
data = Dataset(filename)
jpiglo = data.variables["jpiglo"][:]
jpjglo = data.variables["jpjglo"][:]
jpkglo = data.variables["jpkglo"][:]
lonV   = data.variables["glamv"][0, :, :]
latV   = data.variables["gphiv"][0, :, :]
lonT   = data.variables["glamt"][0, :, :]
latT   = data.variables["gphit"][0, :, :]
e1t    = data.variables["e1t"][0, 0, 0]
z      = data.variables["nav_lev"][:]
data.close()

# have a sinusoidally varying wind but zero everything else
Ly     = 2000.0e3
Ly_mid = (Ly - e1t) / 2.0 # take into account the slight offset of the T/Vgrid
Lz     = 3000

tau0 = 0.2 * 1.0
utau = np.zeros((1, jpjglo, jpiglo))
utau[0, :, :] = tau0 * 0.5 * ( 1.0 + np.cos(2.0 * np.pi * (latT * 1.0e3 - Ly_mid) / Ly) )
utau[0, 0, :] = 0.0
utau[0, -1, :] = 0.0

vtau = np.zeros((1, jpjglo, jpiglo))
qtot = np.zeros((1, jpjglo, jpiglo))
qsr  = np.zeros((1, jpjglo, jpiglo))
emp  = np.zeros((1, jpjglo, jpiglo))

plt.plot(utau[0, :, 0], latT, 'rx-')
plt.grid()

gen_forcing(jpiglo, jpjglo, lonT, latT, utau, vtau, qtot, qsr, emp, filename.replace("IDEAL_bump/domain_cfg", "IDEAL_bump/forcing"))

#%%
# Choose the amplitude of the temperature variations
dtheta = 15.0

# Choose the e-folding scale for the stratification.
z0 = 1000.0

# Generate the idealised temperature stratification, used as both an initial
# condition and for the restoring temperature in the sponge regions.

toce = np.zeros((1, jpkglo, jpjglo, jpiglo))
soce = 35.0 * np.ones((1, jpkglo, jpjglo, jpiglo))

sss = deepcopy(soce[:, 0, :, :])

yy, zz = np.meshgrid(latT[:, 0], -z)
yy *= 1.0e3

# Linear gradient at surface with exponential decay at depth and 0oC at
# the southern boundary (similar to Abernathey et al., 2011)
for ji in range(jpiglo):
    toce[0, :, :, ji] =   ( ( 1.0 + ( dtheta * ( yy / Ly )  ) ) * (
                                    ( np.exp( zz / z0) - np.exp(-Lz / z0) ) 
                                    / ( 1.0 - np.exp( -Lz / z0 ) )
                                    ) 
                        )
    
# for ji in range(jpiglo):
#     toce[0, :, :, ji] =   (  (
#                                ( np.exp( zz / z0) - np.exp(-Lz / z0) ) 
#                                / ( 1.0 - np.exp( -Lz / z0 ) )
#                                ) 
#                     )

# because of the grid there is some offset, add it back on
# toce -= toce[0, 0, 0, 0]

# add some noise
# toce += np.random.normal(0, 0.01, (1, jpkglo, jpjglo, jpiglo))

# % Make sure there is no water colder than 0.5oC.
toce[toce < 0.25] = 0.25
    
# Pick out the surface temperature for the restoring condition.
sst = deepcopy(toce[:, 0, :, :])

# note the first and last y point is going to be set to masked out
plt.contourf(yy / 1.0e3, zz, toce[0, :, :, 1], np.arange(0, 18, 1), cmap = "RdBu_r")
lines = plt.contour(yy / 1.0e3, zz, toce[0, :, :, 1], np.arange(2, 18, 2), colors = "w")
plt.clabel(lines, fmt = r"$%i\ {}^\circ \mathrm{C}$", colors = 'w')

gen_istate(lonT, latT, z, toce, soce, sst, sss, filename.replace("IDEAL_bump/domain_cfg", "IDEAL_bump/state"))


#%%
# MITgcm
rho0 = 1035
alpT = 2.0e-4
Tref = 2.5
print("MITgcm value = %.6f" % (rho0 * (1 - alpT * (15 - Tref))))
print("MITgcm value = %.6f" % (rho0 * (1 - alpT * (0  - Tref))))

# NEMO
rho0 = 1026
alpT = 2.0e-4
Tref = 10
print("NEMO  value = %.6f" % (rho0 * (1 - alpT * (15 - Tref))))
print("NEMO  value = %.6f" % (rho0 * (1 - alpT * (0  - Tref))))

#%%
# testing code for generating diffkr

filename = "UNAGI_R100_domcfg.nc"
data = Dataset(filename)
jpiglo = data.variables["jpiglo"][:]
jpjglo = data.variables["jpjglo"][:]
jpkglo = data.variables["jpkglo"][:]
data.close()
reso_m = 100.0e3 # 100 km resolution
e2 = reso_m * np.arange(0, jpjglo)

avt0 = 1.0e-5
avtf = 5.0e-3
L_sponge = 300.0e3
Ly   = 2400.0e3 # 2400km is where the wall SHOULD be at

diffkr = avt0 + (0.5 * avtf * ( 1 + np.cos( np.pi * (e2 - e2[-2]) / L_sponge ) )
               - 0.5 * avt0 * ( 1 + np.cos( np.pi * (e2 - e2[-2]) / L_sponge ) )
                )
diffkr = np.where(e2 > 2000.0e3, diffkr, avt0)

amp_factor = 1 + (0.5 * 500 * ( 1 + np.cos( np.pi * (e2 - e2[-2]) / L_sponge ) )
               -  0.5       * ( 1 + np.cos( np.pi * (e2 - e2[-2]) / L_sponge ) )
                )
amp_factor = np.where(e2 > 2000.e3, amp_factor, 1.0)

plt.plot(diffkr, e2 / 1e3, 'bx-')
plt.plot(amp_factor * avt0, e2 / 1e3, 'ro-')
plt.plot([0, avtf], [2000, 2000], 'k--')
plt.grid()


# test = (1. + 0.5 * 500.0 * (    1. + np.cos(  np.pi * ( yT - 2300. ) / 300.  )   )
#            - 0.5  * ( 1. + np.cos(  np.pi * ( yT - 2300. ) / 300.  )   )
#        )
spacing = 15
y = np.arange(-spacing, 2400 + spacing, spacing)

test = (0. + 0.5 * 500.0 * ( 1. + np.cos(  np.pi * ( y - y[-2] ) / 300.  )   )
           - 0.5         * ( 1. + np.cos(  np.pi * ( y - y[-2] ) / 300.  )   )
       )
test = np.where(y >= 2100, test, 0.0)
test *= 1.0e-5

test1 = (0. + 0.5 * 250.0 * ( 1. + np.cos(  np.pi * ( y - y[-2] ) / 300.  )   )
            - 0.5         * ( 1. + np.cos(  np.pi * ( y - y[-2] ) / 300.  )   )
       )
test1 = np.where(y >= 2100, test1, 0.0)
test1 *= 1.0e-5

test2 = (0. + 0.5 * 500.0 * ( 1. + np.cos(  np.pi * ( y - 2000 ) / 150.  )   )
            - 0.5         * ( 1. + np.cos(  np.pi * ( y - 2000 ) / 150.  )   )
       )
test2 = np.where(y > 1850, test2, 0.0)
test2 = np.where(y <= 2000, test2, 0.0)
test2 *= 1.0e-5

plt.plot(test, y, 'ro-', test1, y, 'g^-', test2, y, 'bx-')
plt.plot([0, 0.005], [2400, 2400], 'k--')

int1 = np.trapz(test[(y < 2400)], y[(y < 2400)])
int2 = np.trapz(test1[(y < 2400)], y[(y < 2400)])
int3 = np.trapz(test2[(y <= 2000)], y[(y <= 2000)])

print("int1 = %.4f, int2 = %.4f, int3 = %.4f" % (int1, int2, int3) )


#%%
# Dave's configuration, based on SO channel of Dave Munday
# 4000km long, 2000km wide, 3000m deep

jpiglo = 40
jpjglo = 20 + 2 # add to grid points on
jpkglo = 30 + 1 # add the bottom level in
reso_m = 100.0e3 # 100 km resolution
e1 = reso_m * np.arange(0, jpiglo)
e2 = reso_m * np.arange(0, jpjglo)
Lx, Ly = e1[-1], e2[-1]

# UNAGI has 3000m depth, also put two walls in
oce_depth = 3000.0
bathy = oce_depth * np.ones((1, jpjglo, jpiglo))

l_ridge = True
ridge_H = 1500.0
ridge_L = 400.0e3 # half width of ridge

if l_ridge:
    ridge_x1d = np.zeros(jpiglo)
    for ji in range(jpiglo):
        if (e1[ji] > -ridge_L - 100.0e3 + Lx / 2.0) & (e1[ji] < -100.0e3 + Lx / 2.0):
            ridge_x1d[ji] = 0.5 * ridge_H * (1.0 + np.cos(np.pi * (e1[ji] + Lx / 2.0 + 100.0e3) / ridge_L))
        elif (e1[ji] > 100.0e3 + Lx / 2.0) & (e1[ji] < ridge_L + 100.0e3 + Lx / 2.0):
            ridge_x1d[ji] = 0.5 * ridge_H * (1.0 + np.cos(np.pi * (e1[ji] + Lx / 2.0 - 100.0e3) / ridge_L))
        elif (e1[ji] >= -100.0e3 + Lx / 2.0) & (e1[ji] <= 100.0e3 + Lx / 2.0):
            ridge_x1d[ji] = ridge_H
else:
    ridge_x1d = oce_depth * np.zeros(jpiglo)

bathy[0, :, :] -= ridge_x1d[np.newaxis, :]
bathy[0,  0, :] = 0.0
bathy[0, -1, :] = 0.0
    
plt.contourf(e1 / 1e3, e2 / 1e3, bathy[0, :, :])
plt.colorbar()

gen_bathy_meter(jpiglo, jpjglo, e1, e2, bathy, "DAVE", "bathy_meter_DAVE_R100.nc")


#%%
# use the generated bathy_meter.nc in DOMAINcfg to generate a domaincfg file
# the use that domaincfg file to generate the initial state and forcing files

# no meridional wind here so just put everything on T points
filename = "domcfg_R100_DAVE.nc"
data = Dataset(filename)
jpiglo = data.variables["jpiglo"][:]
jpjglo = data.variables["jpjglo"][:]
jpkglo = data.variables["jpkglo"][:]
lonV   = data.variables["glamv"][0, :, :]
latV   = data.variables["gphiv"][0, :, :]
lonT   = data.variables["glamt"][0, :, :]
latT   = data.variables["gphit"][0, :, :]
e1t    = data.variables["e1t"][0, 0, 0]
z      = data.variables["nav_lev"][:]
data.close()

# have a sinusoidally varying wind but zero everything else
Ly     = 2000.0e3
Ly_mid = (Ly - e1t) / 2.0 # take into account the slight offset of the T/Vgrid
Lz     = 3000

tau0 = 0.2
utau = np.zeros((1, jpjglo, jpiglo))
utau[0, :, :] = tau0 * 0.5 * ( 1.0 + np.cos(2.0 * np.pi * (latT * 1.0e3 - Ly_mid) / Ly) )
utau[0, 0, :] = 0.0
utau[0, -1, :] = 0.0

vtau = np.zeros((1, jpjglo, jpiglo))
qtot = np.zeros((1, jpjglo, jpiglo))
qsr  = np.zeros((1, jpjglo, jpiglo))
emp  = np.zeros((1, jpjglo, jpiglo))

plt.plot(utau[0, :, 0], latT, 'rx-')
plt.grid()

gen_forcing(jpiglo, jpjglo, lonT, latT, utau, vtau, qtot, qsr, emp, "forcing_R100_DAVE.nc")


#%%
# Choose the amplitude of the temperature variations
dtheta = 15.0

# Choose the e-folding scale for the stratification.
z0 = 1000.0

# Generate the idealised temperature stratification, used as both an initial
# condition and for the restoring temperature in the sponge regions.

toce = np.zeros((1, jpkglo, jpjglo, jpiglo))
soce = 35.0 * np.ones((1, jpkglo, jpjglo, jpiglo))

sss = deepcopy(soce[:, 0, :, :])

yy, zz = np.meshgrid(latT[:, 0], -z)
yy *= 1.0e3

# Linear gradient at surface with exponential decay at depth and 0oC at
# the southern boundary (similar to Abernathey et al., 2011)
for ji in range(jpiglo):
    toce[0, :, :, ji] = (dtheta * (  ( yy / Ly )
                                   * ( np.exp( zz / z0) - np.exp(-Lz / z0) ) 
                                   / ( 1.0 - np.exp( -Lz / z0 ) )
                                   ) 
                        )

# because of the grid there is some offset, add it back on
toce -= toce[0, 0, 0, 0]

# add some noise
toce += np.random.normal(0, 0.05, (1, jpkglo, jpjglo, jpiglo))

# % Make sure there is no water colder than 0.5oC.
toce[toce < 0.25] = 0.25
    
# Pick out the surface temperature for the restoring condition.
sst = deepcopy(toce[:, 0, :, :])

# note the first and last y point is going to be set to masked out
plt.contourf(yy / 1.0e3, zz, toce[0, :, :, 1], np.arange(0, 16, 1), cmap = "RdBu_r")
lines = plt.contour(yy / 1.0e3, zz, toce[0, :, :, 1], np.arange(2, 16, 2), colors = "w")
plt.clabel(lines, fmt = r"$%i\ {}^\circ \mathrm{C}$", colors = 'w')

gen_istate(lonT, latT, z, toce, soce, sst, sss, "state_R100_DAVE.nc")




