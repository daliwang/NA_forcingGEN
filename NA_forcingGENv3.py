import os,sys
import glob
import math
import netCDF4 as nc
import numpy as np
from time import process_time
from datetime import datetime, timedelta

def is_leap_year(year):
    """ Check if a given year is a leap year. """
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)

version_date = "03312025"

def forcing_1dNA(input_path, file, var_name, period, time, output_path):
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    source_file = input_path + '/'+ file
    src = nc.Dataset(source_file, 'r', format='NETCDF3_64BIT_DATA')

    # Parse the filename into separate substrings
    parts = file.split('.')
    parts.insert(3, '1d')
    # Join the parts back together to form the new filename
    new_filename = '.'.join(parts)

    total_rows = src.dimensions['x'].size
    total_cols = src.dimensions['y'].size
    total_time = src.dimensions['time'].size

    # remove the leap day 
    if total_time == 232:
            total_time = 224

    #print('total timesteps is :' + str(total_timesteps))
    if time == -1:
        time = total_time

#    data = src[var_name][0:time, :, :]
    x_dim = src['x'][...]
    y_dim = src['y'][...]
    latxy = src['lat'][:]
    lonxy = src['lon'][:]

    # Sample code for time variable transformation:
    time_variable = src['time'][:]  # Read full time variable

    # remove the leap day
    if (len(time_variable)) == 232:
        time_variable = time_variable[0:224]
    tunit = src.variables['time'].getncattr('units')  # Expected format: "days since 1950-01-01 00:00:00"
    print(tunit[10:])

    t0 = datetime.strptime(tunit[11:], '%Y-%m-%d %H:%M:%S')  # Parse the reference time

        # Initialize an output array for the days of the current month
    days_in_month = np.zeros_like(time_variable, dtype=float)
    months = np.zeros_like(time_variable, dtype=int)  # Store month indices

    # Reference handling for each year
    for i, days_since in enumerate(time_variable):
        # Calculate the full year based on total days since t0
        current_date = t0 + timedelta(days=float(days_since))  # Current date from base date

        # Get the year, month, and the day of the month
        year = current_date.year
        month = current_date.month
        day = current_date.day

        # Store month and day information
        months[i] = month
        #days_in_month[i] = day
        # Instead of storing days, calculate times for 8 intervals per day
        days_in_month[i] = (day-1)  + (1.5 + (i%8) *3 )/24

    # Prepare the output time units in the required format
    new_time_unit = f"days since {int(year)}-{str(month).zfill(2)}-01 00:00:00"

    # Example for how to update the data_time output
    data_time = days_in_month  # This now holds the day of the current month

    # Get mask and create gridID, latxy, lonxy first

    #create land mask
    mask = src[var_name][0:1, :, :]   # data is in (time, Y, X) format
    mask = np.where(~np.isnan(mask), 1, np.nan)

    #create gridIDs
    total_gridcells = total_rows * total_cols
    grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)
    grid_ids = grid_ids.reshape(total_cols,total_rows)

    # create an flattened list of land gridID and reduce the size of gridIDs array
    grid_ids = np.multiply(mask,grid_ids)
    grid_ids = grid_ids[~np.isnan(grid_ids)]

    latxy = np.multiply(mask,latxy)
    latxy = latxy[~np.isnan(latxy)]
    lonxy = np.multiply(mask,lonxy)
    lonxy = lonxy[~np.isnan(lonxy)]

    # convert local grid_id_lists into an array
    grid_id_arr = np.array(grid_ids)

    #data_arr = np.array(FSDS_list[i])
    #data_arr = np.array(data)
    lonxy_arr= np.array(lonxy)
    latxy_arr= np.array(latxy)

    dst_name = output_path + '/clmforc.Daymet4.1km.1d.' + var_name + '.' + period +'.nc'

    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    dst = nc.Dataset(dst_name, 'w', format='NETCDF3_64BIT_DATA')
    dst.title = var_name + '('+period+') creted from '+ input_path +' on ' +version_date

    # create the gridIDs, lon, and lat variable
    x = dst.createDimension('time', time)
    x = dst.createDimension('ni', grid_id_arr.size)
    x = dst.createDimension('nj', 1)
    x = dst.createDimension('x_dim', len(x_dim))
    x = dst.createDimension('y_dim', len(y_dim))

    w_nc_var = dst.createVariable('gridID', np.int32, ('nj','ni'))

    #  set variable attributes
    w_nc_var.long_name = "gridId in the NA domain" ;
    w_nc_var.decription = "Covers all land and ocean gridcells, with #0 at the upper left corner of the domain" ;
    dst.variables['gridID'][...] = grid_id_arr.reshape(grid_id_arr.size,1)

    # Create variables for x_dim and y_dim in the output NetCDF file
    x_dim_var = dst.createVariable('x_dim', np.float32, ('x_dim',), zlib=True, complevel=5)
    y_dim_var = dst.createVariable('y_dim', np.float32, ('y_dim',), zlib=True, complevel=5)

    # Assign values to the x and y variables
    dst.variables['x_dim'][...] = x_dim
    dst.variables['y_dim'][...] = y_dim

    # Set attributes for x and y variables (optional)
    x_dim_var.long_name = "x coordinate"
    x_dim_var.units = "m: east to west"  # Adjust units if necessary
    y_dim_var.long_name = "y coordinate"
    y_dim_var.units = "m: north to south"  # Adjust units if necessary

    # Copy the variables from the source to the target
    for name, variable in src.variables.items():

        start = process_time()
        print("Working on varibale: "+ name + " dimensions: " + str(variable.dimensions))
       
        if variable.datatype == np.int32:
            fill_value = -9999  # or any other value that you want to use to represent missing data
        if variable.datatype == np.float32:
            fill_value = np.float32(-9999)  # or any other value that you want to use to represent missing data

        # Check if the last two dimensions are y and x and save the data 
 
        if len(variable.dimensions) == 3 and variable.dimensions[-2:] == ('y', 'x'):
            chunk_size = 30  # Number of time steps to process at a time
            landcells = len(grid_ids)

            # Create the variable in the output NetCDF file
            w_nc_var = dst.createVariable(name, np.float32, ('time', 'nj', 'ni'), zlib=True, complevel=5)

            for start in range(0, time, chunk_size):
                end = min(start + chunk_size, time)  # Ensure we don't exceed the total time steps
                print(f"Processing time steps {start} to {end}")

                # Read a chunk of data from the source file
                data_chunk = src[name][start:end, :, :]

                # Extract the data over land grid cells
                data_chunk = np.multiply(mask, data_chunk)
                data_chunk = data_chunk[~np.isnan(data_chunk)]
                data_chunk = np.reshape(data_chunk, (end - start, landcells))

                # Write the chunk to the output file
                dst.variables[name][start:end, 0, :] = data_chunk

            # Copy attributes for the variable
            for attr_name in variable.ncattrs():
                if attr_name != '_FillValue':  # Skip the _FillValue attribute
                    dst[name].setncattr(attr_name, variable.getncattr(attr_name))
 
        if (name == 'time'):
            dvname = 'time'
            w_nc_var = dst.createVariable(dvname, np.float32, ('time'), fill_value=fill_value, zlib=True, complevel=5)
            # Ensure data_time matches the 'time' dimension
            data_time = data_time[:time]  # Slice data_time to match the 'time' dimension

            dst.variables[dvname][...] = data_time
            for attr_name in variable.ncattrs():
                if 'units' in attr_name:
                    dst[dvname].units = new_time_unit
                elif 'calendar' in attr_name:
                    dst[dvname].calendar = "no_leap"                     
                else:
                    dst[dvname].setncattr(attr_name, variable.getncattr(attr_name))
 
        if (name == 'lat'):
            dvname = 'LATIXY'
            w_nc_var = dst.createVariable(dvname, np.float64, ('nj','ni'), zlib=True, complevel=5)
            dst.variables[dvname][...] = latxy_arr
            for attr_name in variable.ncattrs():
                dst[dvname].setncattr(attr_name, variable.getncattr(attr_name))

        if (name == 'lon'):
            dvname = 'LONGXY'
            w_nc_var = dst.createVariable(dvname, np.float64, ('nj','ni'))
            dst.variables[dvname][...] = lonxy_arr
            for attr_name in variable.ncattrs():
                dst[dvname].setncattr(attr_name, variable.getncattr(attr_name))

    src.close()  # close the source file 
    dst.close()  # close the new file        
    

def main():
    args = sys.argv[1:]

    if len(sys.argv) != 4  or sys.argv[1] == '--help':  # sys.argv includes the script name as the first argument
        print("Example use: python NA_forcingGEN.py <input_path> <output_path> <time steps>")
        print(" <input_path>: path to the 1D source data directory")
        print(" <output_path>:  path for the 1D AOI forcing data directory")
        print(" <time steps>: timesteps to be processed or -1 (all time series)")
        print(" The code converts 2D NA forcing inot  1D NA forcing")              
        exit(0)

    input_path = args[0]
    output_path = args[1]
    time = int(args[2])

    # Iterate over all subdirectories in the input directory
    for root, dirs, files in os.walk(input_path):
        for file in files:
            # Check if the file ends with '.nc'
            if file.endswith('.nc'):
                # Parse the filename into separate substrings
                parts = file.split('.')
                var_name = parts[3]
                period = parts[4]         
                print('processing '+ var_name + '(' + period + ') in the file ' + file )
                # Create the corresponding subfolder in the output directory
                new_dir = os.path.join(output_path, os.path.relpath(root, input_path))
                os.makedirs(new_dir, exist_ok=True)
                start = process_time()
                # Copy the file to the new location
                print(root, new_dir)
                forcing_1dNA(root, file, var_name, period, time, new_dir)
                end = process_time()
                print("Generating 1D forcing data takes {}".format(end-start))


if __name__ == '__main__':
    main()
    
