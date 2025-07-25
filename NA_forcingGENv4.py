import os
from datetime import datetime
from time import process_time

import netCDF4 as nc
import numpy as np
from pyproj import Transformer, CRS

def get_transformer():
    """Return a reusable LCC-to-lonlat Transformer."""
    geoxy_proj_str = (
        "+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 "
        "+x_0=0 +y_0=0 +R=6378137 +f=298.257223563 +units=m +no_defs"
    )
    geoxyProj = CRS.from_proj4(geoxy_proj_str)
    lonlatProj = CRS.from_epsg(4326)
    return Transformer.from_proj(geoxyProj, lonlatProj, always_xy=True)

def calculate_offsets(x_dim, y_dim):
    """Calculate half-grid offsets for x and y."""
    x_offset = abs(x_dim[1] - x_dim[0]) / 2
    y_offset = abs(y_dim[1] - y_dim[0]) / 2
    return x_offset, y_offset

def calculate_vertices(Txy2lonlat, XC, YC, x_offset, y_offset):
    """Calculate the four vertices of grid cells."""
    x_offsets = [-x_offset, x_offset, x_offset, -x_offset]
    y_offsets = [y_offset, y_offset, -y_offset, -y_offset]
    xv, yv = [], []
    for dx, dy in zip(x_offsets, y_offsets):
        x, y = Txy2lonlat.transform(XC + dx, YC + dy)
        xv.append(x)
        yv.append(y)
    return np.array(xv), np.array(yv)

def calculate_area_arcrad(xv, yv):
    """Calculate area of grid cells in arc radians squared."""
    xv_rad = np.radians(xv)
    yv_rad = np.radians(yv)
    lon0, lon1, lon2, lon3 = xv_rad
    lat0, lat1, lat2, lat3 = yv_rad
    E1 = (
        np.sin(lat0) * np.sin(lat1) * np.cos(lon1 - lon0)
        + np.sin(lat1) * np.sin(lat2) * np.cos(lon2 - lon1)
        + np.sin(lat2) * np.sin(lat0) * np.cos(lon0 - lon2)
    )
    E1 = np.arccos(E1)
    E2 = (
        np.sin(lat0) * np.sin(lat2) * np.cos(lon2 - lon0)
        + np.sin(lat2) * np.sin(lat3) * np.cos(lon3 - lon2)
        + np.sin(lat3) * np.sin(lat0) * np.cos(lon0 - lon3)
    )
    E2 = np.arccos(E2)
    spherical_excess = E1 + E2 - np.pi
    return abs(spherical_excess)

def create_projection_variable(nc_file, name, **kwargs):
    """Create a projection variable in a NetCDF file."""
    var = nc_file.createVariable(name, np.short)
    for key, value in kwargs.items():
        setattr(var, key, value)

def create_variable(nc_file, name, dtype, dims, data, **kwargs):
    """Create a NetCDF variable."""
    var = nc_file.createVariable(name, dtype, dims, zlib=True, complevel=5)
    for key, value in kwargs.items():
        setattr(var, key, value)
    var[...] = data

def get_masked_arrays(data, lat, lon):
    """Return masked arrays for landmask, landfrac, grid_id, lat, lon, XC, YC."""
    mask = ~np.isnan(data)
    masked = np.where(mask)
    landmask = mask.astype(np.int32)
    landfrac = mask.astype(np.float32)
    return mask, masked, landmask, landfrac

def domain_1dNA(output_path, grid_data, Txy2lonlat, x_offset, y_offset):
    """Save 1D domain data for the Daymet NA region."""
    formatted_date = datetime.now().strftime('%y%m%d')
    data = grid_data["data"][0]
    lat, lon = grid_data["lat"], grid_data["lon"]
    x_dim, y_dim = grid_data["x_dim"], grid_data["y_dim"]

    mask, masked, landmask, landfrac = get_masked_arrays(data, lat, lon)
    total_rows, total_cols = data.shape
    grid_ids = np.arange(total_rows * total_cols).reshape(data.shape)
    grid_id_arr = grid_ids[masked]
    landmask_arr = landmask[masked]
    landfrac_arr = landfrac[masked]
    lat_arr, lon_arr = lat[masked], lon[masked]

    XC, YC = np.meshgrid(x_dim, y_dim)
    XC_arr, YC_arr = XC[masked], YC[masked]
    area_arr = np.full(grid_id_arr.shape, 1.0)

    xv, yv = calculate_vertices(Txy2lonlat, XC_arr, YC_arr, x_offset, y_offset)
    area_arcrad2_arr = calculate_area_arcrad(xv, yv)

    file_name = os.path.join(output_path, f'domain.lnd.Daymet_NA.1km.1d.c{formatted_date}.nc')
    print(f"Saving 1D domain file: {file_name}")

    with nc.Dataset(file_name, 'w', format='NETCDF4') as dst:
        dst.title = '1D domain file for the Daymet NA region'
        create_projection_variable(
            dst, 'lambert_conformal_conic',
            grid_mapping_name="lambert_conformal_conic",
            longitude_of_central_meridian=-100.0,
            latitude_of_projection_origin=42.5,
            false_easting=0.0,
            false_northing=0.0,
            standard_parallel=[25.0, 60.0],
            semi_major_axis=6378137.0,
            inverse_flattening=298.257223563
        )
        dst.createDimension('ni', grid_id_arr.size)
        dst.createDimension('nj', 1)
        dst.createDimension('nv', 4)
        create_variable(dst, 'gridID', np.int32, ('nj', 'ni'), grid_id_arr, long_name='Grid ID in the NA domain')
        create_variable(dst, 'xc_lcc', np.float32, ('nj', 'ni'), XC_arr, long_name='x_coordinate (LCC) of grid cell center', units='m')
        create_variable(dst, 'yc_lcc', np.float32, ('nj', 'ni'), YC_arr, long_name='y_coordinate (LCC) of grid cell center', units='m')
        create_variable(dst, 'xc', np.float32, ('nj', 'ni'), lon_arr, long_name='Longitude of grid cell center', units='degrees_east')
        create_variable(dst, 'yc', np.float32, ('nj', 'ni'), lat_arr, long_name='Latitude of grid cell center', units='degrees_north')
        create_variable(dst, 'xv', np.float32, ('nv', 'nj', 'ni'), xv, long_name='Longitude of grid cell vertices', units='degrees_east')
        create_variable(dst, 'yv', np.float32, ('nv', 'nj', 'ni'), yv, long_name='Latitude of grid cell vertices', units='degrees_north')
        create_variable(dst, 'mask', np.int32, ('nj', 'ni'), landmask_arr, long_name='Land mask', units='unitless')
        create_variable(dst, 'frac', np.float32, ('nj', 'ni'), landfrac_arr, long_name='Land fraction', units='unitless')
        create_variable(dst, 'area', np.float32, ('nj', 'ni'), area_arr, long_name='Area of grid cells (LCC)', units='km^2')
        create_variable(dst, 'area_arcrad', np.float32, ('nj', 'ni'), area_arcrad2_arr, long_name='Area of grid cells (radian)', units='radian^2')

def domain_2dNA(output_path, grid_data, Txy2lonlat, x_offset, y_offset):
    """Save 2D domain data for the Daymet NA region."""
    formatted_date = datetime.now().strftime('%y%m%d')
    data = grid_data["data"][0]
    lat, lon = grid_data["lat"], grid_data["lon"]
    x_dim, y_dim = grid_data["x_dim"], grid_data["y_dim"]
    total_rows, total_cols = data.shape

    XC, YC = np.meshgrid(x_dim, y_dim)
    xv, yv = calculate_vertices(Txy2lonlat, XC, YC, x_offset, y_offset)
    area = np.full(data.shape, 1.0)
    area_arcrad2 = calculate_area_arcrad(xv, yv)
    landmask = (~np.isnan(data)).astype(int)
    landfrac = landmask.astype(float)

    file_name = os.path.join(output_path, f'domain.lnd.Daymet_NA.1km.2d.c{formatted_date}.nc')
    print(f"Saving 2D domain file: {file_name}")

    with nc.Dataset(file_name, 'w', format='NETCDF4') as dst:
        dst.title = '2D domain file for the Daymet NA region'
        create_projection_variable(
            dst, 'lambert_conformal_conic',
            grid_mapping_name="lambert_conformal_conic",
            longitude_of_central_meridian=-100.0,
            latitude_of_projection_origin=42.5,
            false_easting=0.0,
            false_northing=0.0,
            standard_parallel=[25.0, 60.0],
            semi_major_axis=6378137.0,
            inverse_flattening=298.257223563
        )
        dst.createDimension('ni', total_cols)
        dst.createDimension('nj', total_rows)
        dst.createDimension('nv', 4)
        create_variable(dst, 'x_lcc', np.float32, ('ni',), x_dim, long_name='x_coordinate of grid cell center', units='m')
        create_variable(dst, 'y_lcc', np.float32, ('nj',), y_dim, long_name='y_coordinate of grid cell center', units='m')
        create_variable(dst, 'xc', np.float32, ('nj', 'ni'), lon, long_name='Longitude of grid cell center', units='degrees_east')
        create_variable(dst, 'yc', np.float32, ('nj', 'ni'), lat, long_name='Latitude of grid cell center', units='degrees_north')
        create_variable(dst, 'xv', np.float32, ('nv', 'nj', 'ni'), xv, long_name='Longitude of grid cell vertices', units='degrees_east')
        create_variable(dst, 'yv', np.float32, ('nv', 'nj', 'ni'), yv, long_name='Latitude of grid cell vertices', units='degrees_north')
        create_variable(dst, 'mask', np.int32, ('nj', 'ni'), landmask, long_name='Land mask (1 means land)', units='unitless')
        create_variable(dst, 'frac', np.float32, ('nj', 'ni'), landfrac, long_name='Land fraction', units='unitless')
        create_variable(dst, 'area', np.float32, ('nj', 'ni'), area, long_name='Area of grid cells (LCC)', units='km^2')
        create_variable(dst, 'area_arcrad', np.float32, ('nj', 'ni'), area_arcrad2, long_name='Area of grid cells (radian)', units='radian^2')

def main():
    input_path = './'
    file_name = 'example_data.nc'
    output_path = input_path

    with nc.Dataset(file_name, 'r', format='NETCDF4') as src:
        x_dim, y_dim = src['x'][:], src['y'][:]
        data = src['TBOT'][0:1, :, :]
        lon, lat = src['lon'][:, :], src['lat'][:, :]
        domain_data = {
            "lon": lon,
            "lat": lat,
            "x_dim": x_dim,
            "y_dim": y_dim,
            "data": data,
        }

    Txy2lonlat = get_transformer()
    x_offset, y_offset = calculate_offsets(domain_data["x_dim"], domain_data["y_dim"])

    start = process_time()
    domain_2dNA(output_path, domain_data, Txy2lonlat, x_offset, y_offset)
    print(f"2D domain data saved in {process_time() - start:.2f} seconds")

    start = process_time()
    domain_1dNA(output_path, domain_data, Txy2lonlat, x_offset, y_offset)
    print(f"1D domain data saved in {process_time() - start:.2f} seconds")

if __name__ == '__main__':
    main()