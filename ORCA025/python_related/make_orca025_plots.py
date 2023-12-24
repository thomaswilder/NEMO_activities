#!/usr/bin/env python3
'''
Author: T Wilder
Date: 06/11/2023
Description: Plots various data fields from ORCA025 model runs and analysis. User has options to choose what they want to plot.
Method:
'''

# eventually make this into a package that the user can call from a simple script

import matplotlib.pyplot as plt
from netCDF4 import Dataset
from cftime import num2date
import datetime
import numpy as np
import logging


logging.basicConfig(filename='output.log', level=logging.DEBUG)

class DataVisualiser:
    def __init__(self):
        self.fig, self.ax = plt.subplots()

    def load_netcdf_data(self, file_path, variable_name):
        # Open the NetCDF file
        with Dataset(file_path, 'r') as ds:
            # Access the variable by name
            data = ds.variables[variable_name][:]
        return data

    def plot_time_series(self, time, data_series, labels, title, y_label):
        for data in data_series:
            self.ax.plot(time, data)
        self.ax.set_title(title)
        self.ax.set_xlabel("Time")
        self.ax.set_ylabel(y_label)
        self.ax.grid(True)
        self.ax.legend(labels) 
        plt.show()


# Create an instance of DataVisualizer
visualiser = DataVisualiser()
        
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
plt.savefig('acc_transport_orca025.png',dpi=100, bbox_inches='tight')





