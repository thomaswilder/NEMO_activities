#!/usr/bin/env python3
'''
Author
    T Wilder

Description 
    plot_moc.py makes use of nemo_toolkit/data_view.py to plot a map of Meridional Overturning Circulation data computed using cdfmoc (depth levels).

History:
    Nov 23'   Initial inception
'''

import cmocean
import numpy as np
import matplotlib.pyplot as plt
from nemo_toolkit import data_view

viscosity_used = ["Biharm", "2D Leith", "QG Leith"]
exp = ["cy516", "da643", "da506"]

variable_name = "zomsfatl"

for experiment, viscosity in zip(exp, viscosity_used):

    # Create an instance of MapsVisualiser
    visualiser = data_view.MapsVisualiser()
    
    filepath = f"/home/users/twilder/Python/ORCA025/data/moc_{experiment}o_2007.nc"
    # for the depth
    filepath2 = f"/gws/nopw/j04/terrafirma/twilder/u-cy516/data/nemo_cy516o_thetao_con_2007_grid-T.nc"

    # load data by netcdf
    data = visualiser.load_netcdf_data(filepath, variable_name)
    depth = visualiser.load_netcdf_data(filepath2, "deptht")
    lat = visualiser.load_netcdf_data(filepath, "nav_lat")

    levels = np.arange(-8,22,2)
    visualiser.plot_transect_data(
        lat[:,0], depth[:], data[0,:,:,0], f"{variable_name}, u-{experiment}_2007", 
        levels, cmocean.cm.balance, "both", 0, 0
    )


    # save individual figures
    plt.savefig(f"{variable_name}_{experiment}_2007.png",dpi=100, bbox_inches='tight')
