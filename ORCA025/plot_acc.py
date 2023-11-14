#!/usr/bin/env python3
'''
Author
    T Wilder

Description 
    plot_acc.py makes use of nemo_toolkit/data_view.py to plot a time series of transport data in a netcdf format.

History:
    Nov 23'   Initial inception
'''

import matplotlib.pyplot as plt
from cftime import num2date
import datetime
import numpy as np

from nemo_toolkit import data_view


# Create an instance of TimeSeriesVisualiser
visualiser = data_view.TimeSeriesVisualiser()
        
# read in data and prep for plotting
data_directory = "/home/users/twilder/Python/ORCA025/"

time = visualiser.load_netcdf_data(data_directory + "transport_da506o.nc", "time")
acc_qg = visualiser.load_netcdf_data(data_directory + "transport_da506o.nc", "acc_transport")
acc_2d = visualiser.load_netcdf_data(data_directory + "transport_da643o.nc", "acc_transport")
acc_bh = visualiser.load_netcdf_data(data_directory + "transport_cy516o.nc", "acc_transport")

data_series = acc_qg, acc_2d, acc_bh
labels = "Leith_QG", "Leith_2D", "Biharm"
title = "Transport through Drake passage"
y_label = "Sv"

# convert time into datetime
time_cftime = num2date(time,"hours since 1800-01-01")
time = []
for i in range(len(time_cftime)):
    time.append(datetime.datetime.strptime(str(time_cftime[i]),'%Y-%m-%d %H:%M:%S'))

# making plot
visualiser.plot_time_series(time, data_series, labels, title, y_label)
plt.savefig('acc_transport_test_test.png',dpi=100, bbox_inches='tight')