
import os,sys 
import netCDF4 as nc
import numpy as np
from time import process_time
from datetime import datetime

# Get current date
current_date = datetime.now()
# Format date to mmddyyyy
formatted_date = current_date.strftime('%m%d%Y')

def forcing_save_1dTES(input_path, file, var_name, period, time, output_path):
    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    source_file = input_path + '/'+ file
    src = nc.Dataset(source_file, 'r', format='NETCDF4')
    
    # Parse the filename into separate substrings
    parts = file.split('.')
    parts.insert(3, '1d')
    # Join the parts back together to form the new filename
    new_filename = '.'.join(parts)
    
    total_rows = src.dimensions['lon'].size
    total_cols = src.dimensions['lat'].size
    total_time = src.dimensions['time'].size

    #print('total timesteps is :' + str(total_timesteps))
    if time == -1:
        time = total_time

    lat= src['lat'][:]
    lon= src['lon'][:]

    lonxy, latxy = np.meshgrid(lon, lat)

    print(lonxy.shape)

    # time in 'days since ' current yyyy-mm-01-01 00:00:00
    data_time = src['time'][0:time] # read (time, y, x) format
    tunit = src.variables['time'].getncattr('units')  # Kao's data is with datetime of leap-year 
    t0=str(tunit.lower()).strip('days since')
    t0=datetime.strptime(t0,'%Y-%m-%d %X')
    iyr = t0.year + np.floor(data_time/365.0)
    data_time0 = datetime.strptime(str(int(iyr[0]))+'-01-01','%Y-%m-%d')
    data_time0 = (data_time0-t0).total_seconds()/86400.0  # days from t0 at the beginning of year
    iday = data_time - data_time0   # now in DOY of leapyear format, will be re-filled below
    imm = np.zeros_like(iday)
    mdoy = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
    for m in range(1,12):
        tpts = np.where((iday>mdoy[m-1]) & (iday<=mdoy[m]))
        if (len(tpts[0])>0): 
            imm[tpts] = m             # in MM, may be 1 day off for leap-year
            iday[tpts]= iday[tpts]-mdoy[m-1]     # day of current month 
 
    data_time=iday          # in days of current year-month
    tunit = tunit.replace(str(t0.year).zfill(4)+'-', str(int(iyr[0])).zfill(4)+'-')
    tunit = tunit.replace('-'+str(t0.month).zfill(2)+'-', '-'+str(int(imm[0])).zfill(2)+'-')
    if(tunit.endswith(' 00') and not tunit.endswith(' 00:00:00')):
        tunit=tunit+':00:00'
    
    # Get mask and create gridID, latxy, lonxy first
    
    for name, variable in src.variables.items():
        
        # Check if the last two dimensions are lat and lon
        # Get the mask and create gridIDs, latxy, longxy
        
        if (variable.dimensions[-2:] == ('lat', 'lon')):
            
            data = src[name][0:1, :, :]
            #create a land mask
            mask = data[0]    # data is in (time, Y, X) format
            mask = np.where(~np.isnan(mask), 1, np.nan)

            #create gridIDs
            total_gridcells = total_rows * total_cols
            grid_ids = np.linspace(0, total_gridcells-1, total_gridcells, dtype=int)
            grid_ids = grid_ids.reshape(total_cols,total_rows)

             # create a flattened list of land gridID
            grid_ids = np.multiply(mask,grid_ids)
            grid_ids = grid_ids[~np.isnan(grid_ids)]

            print(lonxy.shape, latxy.shape, mask.shape)
            latxy = np.multiply(mask,latxy)
            latxy1 = latxy[~np.isnan(latxy)]
            lonxy = np.multiply(mask,lonxy)
            lonxy1 = lonxy[~np.isnan(lonxy)]

            # convert local grid_id_lists into an array
            grid_id_arr = np.array(grid_ids)

            lonxy_arr= np.array(lonxy1)
            latxy_arr= np.array(latxy1)
    
    
    dst_name = output_path + '/'+ new_filename

    # Open a new NetCDF file to write the data to. For format, you can choose from
    # 'NETCDF3_CLASSIC', 'NETCDF3_64BIT', 'NETCDF4_CLASSIC', and 'NETCDF4'
    dst = nc.Dataset(dst_name, 'w', format='NETCDF4')
    dst.title = var_name + '('+period+') creted from '+ input_path +' on ' +formatted_date

    # create the gridIDs, lon, and lat variable
    x = dst.createDimension('time', time)
    x = dst.createDimension('ni', grid_id_arr.size)
    x = dst.createDimension('nj', 1)

    w_nc_var = dst.createVariable('gridID', np.int32, ('nj','ni'))
    #  set variable attributes
    w_nc_var.long_name = "gridId in the NA domain" ;
    w_nc_var.decription = "Covers all land and ocean gridcells, with #0 at the upper left corner of the domain" ;
    dst.variables['gridID'][...] = grid_id_arr.reshape(grid_id_arr.size,1)

    # Copy the variables from the source to the target

    
    for name, variable in src.variables.items():
        #if (name == var_name):
        start = process_time()
        print("Working on varibale: "+ name + " dimensions: " + str(variable.dimensions))
        
        # Check if the last two dimensions are lat and lon and save the data 
        
        if (variable.dimensions[-2:] == ('lat', 'lon') ) :
            
            data = src[name][0:time, :, :]
            
            # extract the data over land gridcells
            landcells = len(grid_ids)    
            data = np.multiply(mask, data)
            data = data[~np.isnan(data)]
            data = np.reshape(data,(time,landcells))

            data_arr = np.array(data)
            
            w_nc_var = dst.createVariable(name, np.float32, ('time', 'nj', 'ni'), zlib=True, complevel=5)
            dst.variables[name][:] =data_arr.reshape(time,grid_id_arr.size)
            for attr_name in variable.ncattrs():
                dst[name].setncattr(attr_name, variable.getncattr(attr_name))
        
        if (name == 'time'):
            dvname = 'time'
            w_nc_var = dst.createVariable(dvname, np.float32, ('time'))
            dst.variables[dvname][...] = data_time
            for attr_name in variable.ncattrs():
                if 'units' in attr_name:
                    dst[dvname].units = tunit
                else:
                    dst[dvname].setncattr(attr_name, variable.getncattr(attr_name))

        if (name == 'lat'):
            dvname = 'LATIXY'
            w_nc_var = dst.createVariable(dvname, np.float64, ('nj','ni'))
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
    

'''
#input_path = "/gpfs/wolf2/cades/cli185/proj-shared/Daymet_GSWP3_4KM_TESSFA/netcdf/2014/TPHWL3Hrly/"
#input_path = "/gpfs/wolf2/cades/cli185/proj-shared/Daymet_GSWP3_4KM_TESSFA/"
#input_path = "/gpfs/wolf2/cades/cli185/proj-shared/Daymet_GSWP3_4KM_TESSFA/netcdf/2014/Solar3Hrly/"
#output_path = "./temp/"
#time = 1

files_nc = get_files(input_path)

for f in files_nc: 
    if (not f.startswith('clmforc')): continue
    parts = f.split('.')
    common_string = parts[0]
    project=parts[1]
    res = parts[2]
    var_name = parts[3]
    period = parts[4]
    print('processing '+ var_name + '(' + period + ') in the file ' + f )
    start = process_time() 
    forcing_save_1dTES(input_path, f, var_name, period, time, output_path)
    end = process_time()
    print("Generating 1D forcing data takes {}".format(end-start))
'''



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
                forcing_save_1dTES(root, file, var_name, period, time, new_dir)
                end = process_time()
                print("Generating 1D forcing data takes {}".format(end-start))


if __name__ == '__main__':
    main()
    
