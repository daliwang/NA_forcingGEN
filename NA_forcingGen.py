
# data_partition module for batch processing
# based on array_split and function definition

import os,sys
import glob
import math
import netCDF4 as nc
import numpy as np
from time import process_time
from datetime import datetime
from pyproj import Proj
from pyproj import Transformer
from pyproj import CRS

try:
    from mpi4py import MPI
    HAS_MPI4PY=True
except ImportError:
    HAS_MPI4PY=False

# Get current date
current_date = datetime.now()
# Format date to mmddyyyy
formatted_date = current_date.strftime('%m-%d-%Y')

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

    data = src[var_name][0:time, :, :]
    #latxy= src['lat'][:,:]   # lat/lon in original nc seems NOT matching with its projected y/x
    #lonxy= src['lon'][:,:]   # So we do trust more in its y/x data
    x_dim = src['x'][...]
    y_dim = src['y'][...]
    #Proj4: +proj=lcc +lon_0=-100 +lat_1=25 +lat_2=60 +k=1 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m  +no_defs
    geoxy_proj_str = "+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"
    geoxyProj = CRS.from_proj4(geoxy_proj_str)
    # EPSG: 4326
    # Proj4: +proj=longlat +datum=WGS84 +no_defs
    lonlatProj = CRS.from_epsg(4326) # in lon/lat coordinates
    Txy2lonlat = Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)
    grid_x, grid_y = np.meshgrid(x_dim,y_dim)
    lonxy,latxy = Txy2lonlat.transform(grid_x,grid_y)


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
    for m in range(1,13):
        tpts = np.where((iday>mdoy[m-1]) & (iday<=mdoy[m]))
        if (len(tpts[0])>0): 
            imm[tpts] = m             # in MM, may be 1 day off for leap-year
            iday[tpts]= iday[tpts]-mdoy[m-1]     # day of current month 
 
    data_time=iday          # in days of current year-month
    tunit = tunit.replace(str(t0.year).zfill(4)+'-', str(int(iyr[0])).zfill(4)+'-')
    tunit = tunit.replace('-'+str(t0.month).zfill(2)+'-', '-'+str(int(imm[0])).zfill(2)+'-')
    if(tunit.endswith(' 00') and not tunit.endswith(' 00:00:00')):
        tunit=tunit+':00:00'
 
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

    latxy = np.multiply(mask,latxy)
    latxy = latxy[~np.isnan(latxy)]
    lonxy = np.multiply(mask,lonxy)
    lonxy = lonxy[~np.isnan(lonxy)]

    # convert local grid_id_lists into an array
    grid_id_arr = np.array(grid_ids)

    #data_arr = np.array(FSDS_list[i])
    data_arr = np.array(data)
    lonxy_arr= np.array(lonxy)
    latxy_arr= np.array(latxy)
    
    
    dst_name = output_path + '/clmforc.Daymet4.1km.1d.' + var_name + '.' + period +'.nc'

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
        if (name == var_name):
            w_nc_var = dst.createVariable(var_name, np.float32, ('time', 'nj', 'ni'), zlib=True, complevel=5)
            dst.variables[var_name][:] =data_arr.reshape(time,grid_id_arr.size)
            for attr_name in variable.ncattrs():
                dst[name].setncattr(attr_name, variable.getncattr(attr_name))
        
        if (name == 'time'):
            dvname = 'time'
            w_nc_var = dst.createVariable(dvname, np.float64, ('time'), zlib=True, complevel=5)
            dst.variables[dvname][...] = data_time
            for attr_name in variable.ncattrs():
                if 'units' in attr_name:
                    dst[dvname].units = tunit
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
    
def get_files(input_path, ncheader='clmforc'):
    print(input_path+ncheader)
    #files = os.listdir(input_path)
    files = glob.glob("%s*.%s" % (inputpath+ncheader,'nc'))
    
    files.sort() 

    #files_nc = [f for f in files if (f[-2:] == 'nc')] 
    print("total " + str(len(files)) + " files need to be processed")
    return files

def main():
    args = sys.argv[1:]

    if len(sys.argv) != 4  or sys.argv[1] == '--help':  # sys.argv includes the script name as the first argument
        print("Example use: python NA_forcingGEN.py <input_path> <output_path> <time steps>")
        print(" <input_path>: path to the 2D source data directory")
        print(" <output_path>:  path for the 1D forcing data directory")
        print(" <time steps>: timesteps to be processed or -1 (all time series)")
        print(" The code converts 2D NA forcing inot  1D NA forcing")              
        exit(0)

    input_path = args[0]
    if not input_path.endswith('/'): input_path=input_path+'/'
    output_path = args[1]
    if not output_path.endswith('/'): output_path=output_path+'/'
    time = int(args[2])
    
    files_nc = get_files(input_path)
    n_files = len(files_nc)
    #
#------------------------------------------------------------------------------------------
# mpi implementation - simply round-robin 'n_files' over cpu_cores
    if HAS_MPI4PY:
        mycomm = MPI.COMM_WORLD
        myrank = mycomm.Get_rank()
        mysize = mycomm.Get_size()
    else:
        mycomm = 0
        myrank = 0
        mysize = 1

    ng = math.floor(n_files/mysize)
    ng_rank = np.full([mysize], int(1))
    ng_rank = np.cumsum(ng_rank)*ng
    xg = int(math.fmod(n_files, mysize))
    xg_rank = np.full([mysize], int(0))
    if xg>0: xg_rank[:xg]=1
    ng_rank = ng_rank + np.cumsum(xg_rank) - 1        # ending file index, starting 0, for each rank
    ng0_rank = np.hstack((0, ng_rank[0:mysize-1]+1))  # starting file index, starting 0, for each rank

    for i in range(ng0_rank[myrank], ng_rank[myrank]+1):
        f = files_nc[i]
        if (not f.split('/')[-1].startswith('clmforc')): continue
        var_name = f[20:-11]
        period = f[-10:-3]
        if (myrank==0): print('processing '+ var_name + '(' + period + ') in the file ' + f )
        start = process_time() 
        forcing_save_1dNA(input_path, f.split('/')[-1], var_name, period, time, output_path)
        end = process_time()
        if (myrank==0): print("Generating 1D forcing data takes {}".format(end-start))

if __name__ == '__main__':
    main()
    
