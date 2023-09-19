#!/usr/bin/env python

# See Bard chat for details.
# When plotting data, convert time data to datetime using num2date, see below.

import os
import re
import numpy as np
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from cftime import date2num, num2date
from netCDF4 import Dataset
from pyCDFTOOLS.cdftransport import cdftransport


def get_date_from_filename(filename):
	"""
	Some description 
	"""
	date_str = filename.split("-")[0] # splits filename at "-" and takes the first string
	return date_str[-8:] # returns the final 8 characters corresponding to date

def get_U_files(directory_path):
	"""
	Some description 
	"""
	files = os.listdir(directory_path)

	filtered_files = []
	for file in files:
		match = re.search(r"grid-U.nc$", file) # search filename for string
		if match: # if string in filename then match and append to list
			filtered_files.append(file)

	return filtered_files

def get_acc_transport(directory_path, data_variable, kwargs):
	"""
	Some description 
	"""
	files = get_U_files(directory_path)
	files.sort(key=get_date_from_filename) # sorts filenames using get_date function

	acc_transport = np.zeros((len(files)))
	for i in range(len(files)):
		acc_transport[i] = cdftransport(directory_path, files[i], data_variable, **kwargs)

	return acc_transport



def main():
	user_path = "/gws/nopw/j04/terrafirma/twilder/"

	directories = ["u-cz312/data/", "u-cy517/data/", "u-cy516/data/"]

	# find length of file list to intialise numpy array
	files = get_U_files(user_path + directories[0])

	acc_transport = np.zeros((len(files),3))
	for j, directory in enumerate(directories):
		print(f"starting cdftransport for directory {directory}")
		acc_transport[:,j] = get_acc_transport(user_path + directory, "uo", {"kt": 0, "lprint": False, "drake": True})

	# make into a function called write_to_netcdf()
	print("writing data to netcdf")
	# create netcdf file that stores acc transport data
	dataset = Dataset("acc_transport.nc", "w")

	# create dimensions
	level = dataset.createDimension("level", 1)
	time = dataset.createDimension("time", len(files)) # make this into a date string?

	# create variables
	levels = dataset.createVariable("level", "i4", ("level",))
	times = dataset.createVariable("time", "f8", ("time",))
	acc_qg = dataset.createVariable("acc_qg", "f4", ("time","level",))
	acc_2d = dataset.createVariable("acc_2d", "f4", ("time","level",))
	acc_bh = dataset.createVariable("acc_bh", "f4", ("time","level",))

	# set up date coordinates
	dates = [datetime(1976,2,1) + n*relativedelta(months=+1) for n in range(len(files))]

	times[:] = date2num(dates,"hours since 1800-01-01")
	print(times)
	print(num2date(times[:],"hours since 1800-01-01"))

	# set up attributes
	acc_qg.standard_name = "acc_transport_qg"
	acc_qg.long_name = "ACC transport in QG Leith"
	acc_qg.units = "Sv"

	acc_2d.standard_name = "acc_transport_2d"
	acc_2d.long_name = "ACC transport in 2D Leith"
	acc_2d.units = "Sv"

	acc_bh.standard_name = "acc_transport_bh"
	acc_bh.long_name = "ACC transport in Biharmonic"
	acc_bh.units = "Sv"

	# populate netcdf file
	acc_qg[:] = acc_transport[:,0]
	acc_2d[:] = acc_transport[:,1]
	acc_bh[:] = acc_transport[:,2]


	# close file
	dataset.close()

		

if __name__ == "__main__":
	main()


