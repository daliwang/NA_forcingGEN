
# data_partition module for batch processing
# based on array_split and function definition

import os,sys
import netCDF4 as nc
import numpy as np
import pandas as pd
from time import process_time
from datetime import datetime

# Get current date
current_date = datetime.now()
# Format date to mmddyyyy
formatted_date = current_date.strftime('%m%d%Y')

def AOI_forcing_save_1d(input_path, file, AOI, AOI_points, var_name, period, time, output_path):
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    source_file = input_path + '/'+ file
    src = nc.Dataset(source_file, 'r', format='NETCDF4')

    #total_gridcell = src.dimensions['gridcell'].size
    total_time = src.dimensions['time'].size

    #print('total timesteps is :' + str(total_timesteps))
    if time == -1:
        time = total_time

    #read gridIDs
    grid_ids = src['gridID']    # gridID for NA

    # read forcing data
    data = src[var_name][0:time, :] # read (time, gridcell) format
    print(data.shape)
    # create the mask for data subsetting using AOI points
    # TODO
    #AOI_mask = np.where(~np.isnan(mask), 1, np.nan)
    #AOI_landcells = len(AOI_points)

    # Create a boolean mask
    AOI_mask = np.isin(grid_ids, AOI_points)

    # You can apply the same mask to this array
    masked_data = data[:, AOI_mask]

    # extract the data over land gridcells into a new data_arr for netcdf  
    data_arr = np.copy(masked_data)
    print(masked_data.shape)
    # convert local grid_id_lists into an array
    #grid_id_arr = np.array(grid_ids)

    dst_name = output_path + '/'+ AOI + '_clmforc.Daymet4.1km.1d.' + var_name + '.' + period +'.nc'

    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    dst = nc.Dataset(dst_name, 'w', format='NETCDF4')
    dst.title = var_name + '('+period+') creted from '+ input_path +' on ' +formatted_date

    # create the gridIDs, lon, and lat variable
    x = dst.createDimension('gridcell', len(AOI_points))
    x = dst.createDimension('time', time)

    w_nc_var = dst.createVariable('gridID', np.int32, ('gridcell',))
    #  set variable attributes
    w_nc_var.long_name = "gridId in the "+AOI+ " domain" ;
    w_nc_var.decription = "Land gridcells within the "+ AOI +" domain" ;
    dst.variables['gridID'][:] = AOI_points

    # Copy the variables from the source to the target
    for name, variable in src.variables.items():
        if (name == var_name):
            print(variable.datatype)
            w_nc_var = dst.createVariable(var_name, np.float32, ('time', 'gridcell'))
            dst.variables[var_name][:] =data_arr
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
        print("Example use: python AOI_forcingGEN.py <points_file_name> <AOI_domain_name>")
        print(" <input_path>: path to the 1D source data directory")
        print(" <output_path>:  path for the 1D AOI forcing data directory")
        print(" <points_file_name>:  <AOI>_gridID.csv")
        exit(0)

    input_path = args[0]
    output_path = args[1]
    AOI_gridID_file = args[2]
    AOI='AKSP'
    time = 1

    files_nc = get_files(input_path)

    #AOI_gridcell_file = AOI+'_gridID.csv'  # user provided gridcell IDs
    df = pd.read_csv(AOI_gridID_file, sep=",", skiprows=1, names = ['gridID'])
    #read gridIds
    AOI_points = np.array(df['gridID'])

    for f in files_nc: 
        var_name = f[23:-11]
        period = f[-10:-3]
        print(var_name, period)
        print('processing '+ var_name + '(' + period + ') in the file ' + f )
        start = process_time() 
        #forcing_save_1dNA(input_path, f, var_name, period, time, output_path)
        AOI_forcing_save_1d(input_path, f, AOI, AOI_points, var_name, period, time, output_path)
        end = process_time()
        print("Generating 1D forcing data for "+AOI+ " domain takes {}".format(end-start))

if __name__ == '__main__':
    main()
    