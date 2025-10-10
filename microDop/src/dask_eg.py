import xarray as xr
# import dask # loads automatically with xarray

# Define the path to your NetCDF files.
# You can use wildcards (e.g., *.nc) to select multiple files.
# file_path = '~/data/ATOMIC/microDop/reprocessed20210224/*.nc'
file_path = '../data/*.nc'

# Open the multiple NetCDF files as a single Xarray Dataset using Dask.
# The 'chunks' argument specifies how the data should be divided into Dask chunks.
# For example, {'time': 10} chunks the data along the 'time' dimension into blocks of 10.
# Adjust chunk sizes based on your data and available memory for optimal performance.

# ds = xr.open_mfdataset(file_path, chunks={'time': 10})
#  What's a block of 10? Does dask default chunk size perform better?
#  How do we determine a good chunk size?

#  http://docs.xarray.dev/en/latest/user-guide/dask.html
# ds = xr.open_mfdataset(file_path, parallel=True, engine="h5netcdf")
#  h5netcdf should be faster, but ncdump -k shows it's classic not h5 netcdf
#ds = xr.open_mfdataset(file_path, 
#        parallel=True,
#        combine='nested', 
#        concat_dim='__new_dim__', 
#        decode_times=False)
# this doesn't finish for minutes.
# Too many files (881), bad chunking in files, bad lazy loading?

ds_combined = xr.open_mfdataset(file_path,
                                combine='nested',        # handle nonmonotonic time variable
                                concat_dim='time_index', #
                                decode_times=False,      #
                                join="outer",
                                parallel=True,
                                chunks={'time': 8000, 'range': 250})
# rechunks in 2x10^6-element chunks

# At this point, 'ds' is an Xarray Dataset where the data variables are backed by Dask arrays.
# Computations on 'ds' will be performed in parallel and out-of-core using Dask.
# For example, to compute the mean of a variable 'temperature':
# mean_temp = ds['temperature'].mean(dim='time').compute()
# The '.compute()' method triggers the Dask computation and loads the result into memory.

# You can inspect the dataset to see its structure and Dask array information
print(ds)
