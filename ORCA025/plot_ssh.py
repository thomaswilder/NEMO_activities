#!/usr/bin/env python3
'''
Author
    T Wilder

Description 
    plot_ssh.py makes use of nemo_toolkit/data_view.py to plot a map of sea surface height data.

History:
    Nov 23'   Initial inception
'''

import cmocean
import numpy as np
from nemo_toolkit import data_view

viscosity_used = ["Biharm", "2D Leith", "QG Leith"]
exp = ["cy516", "da643", "da506"]


for experiment, viscosity in zip(exp, viscosity_used):

    # Create an instance of MapsVisualiser
    visualiser = data_view.MapsVisualiser()
    
    filepath = f"/gws/nopw/j04/terrafirma/twilder/u-{experiment}/data/nemo_{experiment}o_ssh_2007_grid-T.nc"
    
    print(visualiser.preview_iris_data(filepath))
    
    variable_name = "sea_surface_height_above_geoid"
    
    # load data
    data = visualiser.load_iris_data(filepath, variable_name)
    print(data)
    
    # transform data
    new_variable = visualiser.transform_iris_data_to_proj(data)
    print(new_variable)
    
    # plot data
    levels = np.arange(-2,2.2,0.2)
    visualiser.plot_map_iris_data(new_variable, f"{variable_name}, {viscosity}_2007", levels, cmocean.cm.balance, "both", f"{variable_name}_u-{experiment}_2007")
