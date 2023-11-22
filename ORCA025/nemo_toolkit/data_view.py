#!/usr/bin/env python3
'''
Author
    T Wilder

Description
    data_view.py defines a class and functions for a user to easily load and plot data from their nemo orca experiments.

History
    Nov 23'   Initial inception


'''

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import iris
import iris.analysis.cartography
import iris.plot as iplt

class TimeSeriesVisualiser:
    '''
    A class DataVisualiser for visualising nemo model output.

    Methods:
        __init__(): Initializes an instance with a matplotlib figure and axis.

        load_netcdf_data(file_path, variable_name): Loads data from a NetCDF file.

        plot_time_series(time, data_series, labels, title, y_label): Plots time series data.
    '''
    def __init__(self):
        '''
        Initializes an instance with a matplotlib figure and axis.
        '''
        self.fig, self.ax = plt.subplots()

    def load_netcdf_data(self, file_path, variable_name):
        '''
        Opens a NetCDF file and retrieves data based on the specified variable name.

        Parameters:
            file_path (str): The path to the NetCDF file.
            variable_name (str): The name of the variable to retrieve.

        Returns:
            numpy.ndarray: The loaded data.
        '''
        # Open the NetCDF file
        with Dataset(file_path, 'r') as ds:
            # Access the variable by name
            data = ds.variables[variable_name][:]
        return data

    def plot_time_series(self, time, data_series, labels, title, y_label):
        '''
        Plots time series data.

        Parameters:
            time (numpy.ndarray/datetime): Time values for the x-axis.
            data_series (list of numpy.ndarray): List of data arrays to plot.
            labels (list of str): List of labels for the legend.
            title (str): Title for the plot.
            y_label (str): Label for the y-axis.
        '''
        for data in data_series:
            self.ax.plot(time, data)
        self.ax.set_title(title)
        self.ax.set_xlabel("Time")
        self.ax.set_ylabel(y_label)
        self.ax.grid(True)
        self.ax.legend(labels) 
        plt.show()

class MapsVisualiser:
    '''
    A class for visualizing maps using Iris, Cartopy, and NetCDF4.

    Methods:
        __init__(num_rows=1, num_cols=1, width=12, height=6): Initializes an instance with a matplotlib figure and axis.
        
        preview_iris_data(filepath): Prints information about the Iris data loaded from a file.

        load_iris_data(filepath, variable, **kwargs): Loads Iris data from a file, with optional keyword arguments.

        load_netcdf_data(file_path, variable_name): Loads data from a NetCDF file.

        transform_iris_data_to_proj(data): Transforms Iris data to a target projection.

        plot_map_iris_data(data, title, levels, cmap_name, extend, row, col): Plots Iris data on a map.

        plot_transect_data(xpoints, depth, data, title, levels, cmap_name, extend, row, col): Plots transect data.
    '''
    
    def __init__(self, num_rows=1, num_cols=1, width = 12, height = 6):
        '''
        Initializes an instance with a matplotlib figure and axis.

        n.b. currently only configured for maps and transects, not yet working for odd and even row and cols e.g. 2 by 1.

        Parameters:
            num_rows (int): Number of rows in the subplot grid.
            num_cols (int): Number of columns in the subplot grid.
            width    (int): Width of figure in inches (default 12)
            height   (int): Height of figure in inches (default 6)
        '''
        self.num_rows = num_rows
        self.num_cols = num_cols
        self.width = width
        self.height = height
        self.fig_maps, self.ax_maps = plt.subplots(num_rows, num_cols, figsize=(width, height),
                                        subplot_kw={'projection': ccrs.PlateCarree()})
        # 
        
        # If there is only one subplot, convert it to a 2D array for consistency
        if num_rows == 1 and num_cols == 1:
            self.ax_maps = np.array([[self.ax_maps]])

        # Create a separate subplot for transects
        self.fig_transect, self.ax_transect = plt.subplots(num_rows, num_cols, figsize=(width, height))

        # If there is only one subplot, convert it to a 2D array for consistency
        if num_rows == 1 and num_cols == 1:
            self.ax_transect = np.array([[self.ax_transect]])


    def preview_iris_data(self, filepath):
        '''
        Prints information about the Iris data loaded from a file.

        Parameters:
            filepath (str): The path to the Iris data file.
        '''
        print(iris.load(filepath))

    def load_iris_data(self, filepath, variable, **kwargs):
        '''
        Loads Iris data from a file, with optional keyword arguments.

        Parameters:
            filepath (str): The path to the Iris data file.
            variable (str): The name of the variable to load.
            **kwargs: Choose surface, 150 m, or 1000 m, or None.

        Returns:
            iris.cube.Cube: The loaded Iris data.
        '''
        # default keyword arguments
        opt_dic = {"level": None}
        # overwrite the options by cycling through the input dictionary
        for key in kwargs:
            opt_dic[key] = kwargs[key]
        # logging.info("overwritten default keyword arguments")
        # load full dataset
        data = iris.load(filepath, variable)[0]
        # subset
        if opt_dic["level"]=="surface":
            subset_data = data[0,0]
        elif opt_dic["level"]=="150 m":
            subset_data = data[0,28]
        elif opt_dic["level"]=="1000 m":
            subset_data = data[0,47]
        else:
            # probably 3D field like ssh
            subset_data = data[0]

        return subset_data

    def create_iris_cube(self, data_in, depth, lat, lon, data_name, unit):
        """
        Creates an Iris cube from input data with specified dimensions.
    
        Parameters:
        - data_in (numpy.ndarray): Input data for the cube.
        - depth (numpy.ndarray): Vertical levels for the cube (e.g., depth in meters).
        - lat (numpy.ndarray): Latitude values for the cube in degrees North.
        - lon (numpy.ndarray): Longitude values for the cube in degrees East.
        - data_name (str): Long name for the data, used as the cube's long_name attribute.
        - unit (str): Units for the data, used as the cube's units attribute.
    
        Returns:
        iris.cube.Cube: An Iris cube.
    
        Example:
        >>> data = np.random.rand(10, 5, 5)  # Replace with actual data
        >>> depth_levels = np.linspace(0, 100, 10) or None if surface variable.
        >>> latitudes = np.linspace(-90, 90, 5)
        >>> longitudes = np.linspace(-180, 180, 5)
        >>> cube = create_iris_cube(data, depth_levels, latitudes, longitudes, 'Temperature', 'K')
        """


        # Create a cube
        cube = iris.cube.Cube(data_in, long_name=data_name, units=unit)

        if depth is not None:
            # Add dimension coordinates
            vertical_levels_coord = iris.coords.DimCoord(depth, standard_name='depth', units='m')
            cube.add_dim_coord(vertical_levels_coord, 0)
        
        # Add auxiliary coordinates
        latitude_coord = iris.coords.AuxCoord(lat, standard_name='latitude', units='degrees_north')
        longitude_coord = iris.coords.AuxCoord(lon, standard_name='longitude', units='degrees_east')
        cube.add_aux_coord(latitude_coord, (1, 2))
        cube.add_aux_coord(longitude_coord, (1, 2))

        return cube

    def load_netcdf_data(self, file_path, variable_name):
        '''
        Loads data from a NetCDF file.

        Parameters:
            file_path (str): The path to the NetCDF file.
            variable_name (str): The name of the variable to load.

        Returns:
            numpy.ndarray: The loaded NetCDF data.
        '''
        # Open the NetCDF file
        with Dataset(file_path, 'r') as ds:
            # Access the variable by name
            data = ds.variables[variable_name][:]
        return data

    def transform_iris_data_to_proj(self, data):
        '''
        Transforms Iris data to a target projection.

        Parameters:
            data (iris cube): Iris data to transform.

        Returns:
            iris cube: Transformed Iris data.
        '''
        new_data, extent = iris.analysis.cartography.project(
            data, ccrs.PlateCarree(), nx=1207, ny=1442
        )
        # logging.info("data transformed into projection")
        return new_data

    def plot_map_iris_data(self, data, units, title, levels, cmap_name, extend, row, col):
        '''
        Plots Iris data on a map, either as a single figure, or in a subplot.

        Parameters:
            data (iris cube): Iris data to plot.
            units (str): units of data to be plotted.
            title (str): Title for the plot.
            levels (list or array): Contour levels for the plot.
            cmap_name (str): Name of the colormap.
            extend (str): Colormap extension (e.g., "max", "both").
            row, col (int): position of figure in subplot.
        '''
        # Define the xticks for longitude
        self.ax_maps[row, col].set_xticks(np.arange(-180, 181, 60), crs=ccrs.PlateCarree())
        lon_formatter = cticker.LongitudeFormatter()
        self.ax_maps[row, col].xaxis.set_major_formatter(lon_formatter)
    
        # Define the yticks for latitude
        self.ax_maps[row, col].set_yticks(np.arange(-90, 91, 30), crs=ccrs.PlateCarree())
        lat_formatter = cticker.LatitudeFormatter()
        self.ax_maps[row, col].yaxis.set_major_formatter(lat_formatter)

        # Select the subplot
        self.ax_maps[row, col].coastlines()
        self.ax_maps[row, col].add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
                   
        # Set limits
        self.ax_maps[row, col].set_global()
    
        # set title
        self.ax_maps[row, col].set_title(title)
        
        # plot with Iris plot
        plt.sca(self.ax_maps[row,col])
        cf = iplt.contourf(data, levels=levels, cmap=cmap_name, extend=extend)
        # colorbar
        cb = plt.colorbar(cf, ax=self.ax_maps[row, col], shrink=0.75)
        cb.ax.set_title(f"{units}")

        plt.show()
    
        return self.fig_maps

    def plot_transect_data(self, xpoints, depth, data, title, levels, 
                               cmap_name, extend, row, col):
        '''
        Plots transect data.

        Parameters:
            xpoints (array-like): X-axis values for the transect.
            depth (array-like): Depth values for the transect.
            data (array-like): Data values for the transect.
            title (str): Title for the plot.
            levels (list or array): Contour levels for the plot.
            cmap_name (str): Name of the colormap.
            extend (str): Colormap extension (e.g., "max", "both").
            row, col (int): position of figure in subplot.
        '''
        # plot with contourf
        plt.sca(self.ax_transect[row,col])
        cf = self.ax_transect[row,col].contourf(xpoints, depth, data, 
                                       levels=levels, cmap=cmap_name, extend=extend)
        # colorbar
        cb = plt.colorbar(cf, ax=self.ax_transect[row, col], shrink=0.75)

        # Flip the y-axis
        self.ax_transect[row, col].invert_yaxis()

        # set title
        self.ax_transect[row, col].set_title(title)

        # set x and y labels
        self.ax_transect[row, col].set_xlabel("Latitude (deg)")
        self.ax_transect[row, col].set_ylabel("Depth (m)")

                               
        plt.show()

        return self.fig_transect





