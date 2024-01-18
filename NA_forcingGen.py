
# data_partition module for batch processing
# based on array_split and function definition

import os,sys 
import netCDF4 as nc
import numpy as np
from time import process_time
from datetime import datetime

# Get current date
current_date = datetime.now()
# Format date to mmddyyyy
formatted_date = current_date.strftime('%m%d%Y')

def forcing_save_1dNA(input_path, file, var_name, period, time, output_path):
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    source_file = input_path + '/'+ file
    src = nc.Dataset(source_file, 'r', format='NETCDF4')

    total_rows = src.dimensions['x'].size
    total_cols = src.dimensions['y'].size
    total_time = src.dimensions['time'].size

    #print('total timesteps is :' + str(total_timesteps))
    if time == -1:
        time = total_time

    data = src[var_name][0:time, :, :] # read (time, y, x) format

    #create land mask
    mask = data[0]    # data is in (time, Y, X) format
    mask = np.where(~np.isnan(mask), 1, np.nan)

    #create gridIDs
    total_gridcells = total_rows * total_cols
    grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)
    grid_ids = grid_ids.reshape(total_cols,total_rows)

    # create an flattened list of land gridID and reduce the size of gridIDs array
    grid_ids = np.multiply(mask,grid_ids)
    grid_ids = grid_ids[~np.isnan(grid_ids)]

    # extract the data over land gridcells
    landcells = len(grid_ids)    
    data = np.multiply(mask, data)
    data = data[~np.isnan(data)]
    data = np.reshape(data,(time,landcells))

    # convert local grid_id_lists into an array
    grid_id_arr = np.array(grid_ids)

    #data_arr = np.array(FSDS_list[i])
    data_arr = np.array(data)
    dst_name = output_path + '/clmforc.Daymet4.1km.1d.' + var_name + '.' + period +'.nc'

    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    dst = nc.Dataset(dst_name, 'w', format='NETCDF4')
    dst.title = var_name + '('+period+') creted from '+ input_path +' on ' +formatted_date

    # create the gridIDs, lon, and lat variable
    x  = dst.createDimension('gridcell', grid_id_arr.size)
    x = dst.createDimension('time', time)

    w_nc_var = dst.createVariable('gridID', np.int32, ('gridcell',))
    #  set variable attributes
    w_nc_var.long_name = "gridId in the NA domain" ;
    w_nc_var.decription = "Covers all land and ocean gridcells, with #0 at the upper left corner of the domain" ;
    dst.variables['gridID'][:] = grid_id_arr.reshape(grid_id_arr.size)

    # Copy the variables from the source to the target
    for name, variable in src.variables.items():
        if (name == var_name):
            w_nc_var = dst.createVariable(var_name, variable.datatype, ('time', 'gridcell'), zlib=True, complevel=5)
            dst.variables[var_name][:] =data_arr.reshape(time,grid_id_arr.size)
            for attr_name in variable.ncattrs():
                dst[name].setncattr(attr_name, variable.getncattr(attr_name))
        
    src.close()  # close the source file 
    dst.close()  # close the new file        
    
def get_files(input_path):
    print(input_path)
    files = os.listdir(input_path) 

    files.sort() 

    file_no =0

    files_nc = [f for f in files if (f[-2:] == 'nc')] 
    print("total " + str(len(files_nc)) + " files need to be processed")
    return files_nc

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
    
    files_nc = get_files(input_path)

    for f in files_nc: 
        var_name = f[20:-11]
        period = f[-10:-3]
        print('processing '+ var_name + '(' + period + ') in the file ' + f )
        start = process_time() 
        forcing_save_1dNA(input_path, f, var_name, period, time, output_path)
        end = process_time()
        print("Generating 1D forcing data takes {}".format(end-start))

if __name__ == '__main__':
    main()
    
