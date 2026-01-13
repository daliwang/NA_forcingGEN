import os
import sys
from datetime import datetime
from collections import defaultdict
from multiprocessing import Pool, get_context
import gc

# Globals for worker processes
G_VALID_MASK_FLAT = None
G_VAR_NAME = None
G_LANDCELLS = None

def _worker_init(valid_mask_flat_local, var_name_local, landcells_local):
    global G_VALID_MASK_FLAT, G_VAR_NAME, G_LANDCELLS
    G_VALID_MASK_FLAT = valid_mask_flat_local
    G_VAR_NAME = var_name_local
    G_LANDCELLS = int(landcells_local)

def _read_one(args):
    idx, fpath, steps_this_file = args
    with nc.Dataset(fpath, "r") as src:
        vname = G_VAR_NAME if G_VAR_NAME in src.variables else find_data_variable(src)
        if vname is None:
            return idx, None
        arr = src[vname][0:steps_this_file, :, :]  # (t,y,x)
        # Vectorized masking over time dimension to avoid Python loops
        arr = np.array(arr, dtype=np.float32, copy=False)  # ensure float32 without extra copy
        t, ny, nx = arr.shape
        flat_mask = G_VALID_MASK_FLAT
        # reshape to (t, ny*nx) then select masked columns; returns contiguous array
        out = arr.reshape(t, ny * nx)[:, flat_mask]
        return idx, out
import netCDF4 as nc
import numpy as np


def parse_daily_filename(filename: str):
    """
    Parse Daymet_ERA5 daily filename like:
      clmforc.Daymet4.1km.<VAR>.<YYYYMMDD>.nc
    Returns (variable_name, year, month, day)
    """
    parts = filename.split(".")
    if len(parts) < 5:
        return None
    try:
        variable_name = parts[3]
        yyyymmdd = parts[4]
        year = int(yyyymmdd[0:4])
        month = int(yyyymmdd[4:6])
        day = int(yyyymmdd[6:8])
        return variable_name, year, month, day
    except Exception:
        return None


def is_leap_day(year: int, month: int, day: int) -> bool:
    return month == 2 and day == 29


def find_data_variable(src: nc.Dataset):
    """
    Identify the main 3D data variable with dims ('time','y','x').
    Returns the variable name or None if not found.
    """
    for name, var in src.variables.items():
        if hasattr(var, "dimensions"):
            dims = var.dimensions
            if len(dims) == 3 and dims[0] == "time" and dims[1] == "y" and dims[2] == "x":
                return name
    return None


def build_land_mask_and_coords(first_file_path: str, data_var_name: str):
    """
    Build the 2D land mask (1 for valid, NaN for invalid), grid indices, and 1D lat/lon arrays.
    Uses the first timestep of the first daily file for mask determination (consistent with v3).
    """
    with nc.Dataset(first_file_path, "r") as src:
        x_dim = src["x"][...]
        y_dim = src["y"][...]
        lat2d = src["lat"][:, :]
        lon2d = src["lon"][:, :]

        var_obj = src[data_var_name]
        # Detect fill/missing value if present
        fill_val = None
        if "_FillValue" in var_obj.ncattrs():
            fill_val = var_obj.getncattr("_FillValue")
        elif "missing_value" in var_obj.ncattrs():
            fill_val = var_obj.getncattr("missing_value")
        # Use first timestep for mask
        sample = var_obj[0, :, :]  # shape (y, x)
        if np.ma.isMaskedArray(sample):
            valid = ~sample.mask
            sample = sample.filled(np.nan)
        else:
            sample = np.array(sample)
            valid = ~np.isnan(sample)
        if fill_val is not None:
            valid = np.logical_and(valid, sample != fill_val)
        mask2d = np.where(valid, 1.0, np.nan).reshape(lat2d.shape)

        total_rows = lat2d.shape[0]  # y
        total_cols = lat2d.shape[1]  # x
        total_gridcells = total_rows * total_cols
        grid_ids = np.arange(total_gridcells, dtype=np.int32).reshape(total_rows, total_cols)

        # masked arrays flattened
        masked_grid_ids = (mask2d * grid_ids).ravel()
        masked_grid_ids = masked_grid_ids[~np.isnan(masked_grid_ids)].astype(np.int32)

        masked_lat = (mask2d * lat2d).ravel()
        masked_lat = masked_lat[~np.isnan(masked_lat)].astype(np.float64)

        masked_lon = (mask2d * lon2d).ravel()
        masked_lon = masked_lon[~np.isnan(masked_lon)].astype(np.float64)

        return (
            mask2d,
            masked_grid_ids,
            masked_lat,
            masked_lon,
            x_dim,
            y_dim,
            total_rows,
            total_cols,
        )


def build_time_vector(year: int, month: int, days_in_month: int, time_limit: int | None):
    """
    Build time coordinate as fractional days since YYYY-MM-01 00:00:00.
    Each day contributes 24 hourly timesteps centered on the hour:
      (day_index) + (0.5 + hour) / 24
    """
    times = []
    for day in range(1, days_in_month + 1):
        for hour in range(24):
            times.append((day - 1) + (0.5 + hour) / 24.0)
    if time_limit is not None and time_limit >= 0:
        times = times[:time_limit]
    new_time_unit = f"days since {year:04d}-{month:02d}-01 00:00:00"
    return np.array(times, dtype=np.float32), new_time_unit


def write_monthly_file(
    output_dir: str,
    variable_name: str,
    year: int,
    month: int,
    daily_files: list[str],
    mask2d: np.ndarray,
    grid_ids_1d: np.ndarray,
    lat_1d: np.ndarray,
    lon_1d: np.ndarray,
    x_dim: np.ndarray,
    y_dim: np.ndarray,
    time_limit: int | None,
    max_daily_files: int | None,
    intra_workers: int | None = None,
):
    """
    Create monthly 1D file for the specified variable and (year, month),
    reading and stacking daily files (each with 24 time steps).
    """
    os.makedirs(output_dir, exist_ok=True)
    dst_name = os.path.join(
        output_dir, f"clmforc.Daymet4.1km.1d.{variable_name}.{year:04d}-{month:02d}.nc"
    )

    # Determine number of days to process (skip leap day)
    filtered_files = []
    for fpath in sorted(daily_files):
        parsed = parse_daily_filename(os.path.basename(fpath))
        if parsed is None:
            continue
        _, y, m, d = parsed
        if is_leap_day(y, m, d):
            continue
        filtered_files.append(fpath)
    if max_daily_files is not None and max_daily_files >= 0:
        filtered_files = filtered_files[:max_daily_files]

    num_days = len(filtered_files)
    if num_days == 0:
        return

    # Time vector (24 hours per day)
    total_time_steps = num_days * 24
    if time_limit is not None and time_limit >= 0:
        total_time_steps = min(total_time_steps, time_limit)

    netcdf_format = "NETCDF3_64BIT_DATA"
    dst = nc.Dataset(dst_name, "w", format=netcdf_format)
    try:
        print(f"[{variable_name} {year:04d}-{month:02d}] preparing output: days={num_days}, total_steps={total_time_steps}")
        try:
            sys.stdout.flush()
        except Exception:
            pass
        dst.title = (
            f"{variable_name} ({year:04d}-{month:02d}) created from Daymet_ERA5 on {datetime.now():%Y-%m-%d}"
        )

        # Dimensions
        dst.createDimension("time", total_time_steps)
        dst.createDimension("ni", grid_ids_1d.size)
        dst.createDimension("nj", 1)
        dst.createDimension("x_dim", len(x_dim))
        dst.createDimension("y_dim", len(y_dim))

        # Static vars
        grid_var = dst.createVariable("gridID", np.int32, ("nj", "ni"))
        grid_var.long_name = "gridId in the NA domain"
        grid_var.decription = (
            "Covers gridcells with valid forcing; #0 at the upper left corner of the domain"
        )
        dst.variables["gridID"][...] = grid_ids_1d.reshape(1, grid_ids_1d.size)

        # compression settings: disabled for NETCDF3, enabled for NETCDF4
        comp_args = {} if netcdf_format.startswith("NETCDF3") else dict(zlib=True, complevel=5)

        x_dim_var = dst.createVariable("x_dim", np.float32, ("x_dim",), **comp_args)
        y_dim_var = dst.createVariable("y_dim", np.float32, ("y_dim",), **comp_args)
        x_dim_var.long_name = "x coordinate"
        x_dim_var.units = "m: east to west"
        y_dim_var.long_name = "y coordinate"
        y_dim_var.units = "m: north to south"
        dst.variables["x_dim"][...] = x_dim
        dst.variables["y_dim"][...] = y_dim

        lat_var = dst.createVariable("LATIXY", np.float64, ("nj", "ni"), **comp_args)
        lon_var = dst.createVariable("LONGXY", np.float64, ("nj", "ni"), **comp_args)
        dst.variables["LATIXY"][...] = lat_1d.reshape(1, lat_1d.size)
        dst.variables["LONGXY"][...] = lon_1d.reshape(1, lon_1d.size)

        # Time variable
        # Determine days in month excluding leap day
        month_days = len(filtered_files)
        time_values, new_time_unit = build_time_vector(year, month, month_days, time_limit)
        time_var = dst.createVariable("time", np.float32, ("time",), **comp_args)
        time_var.units = new_time_unit
        time_var.calendar = "no_leap"
        dst.variables["time"][...] = time_values

        # Data variable (set fill_value explicitly for downstream tools)
        data_var = dst.createVariable(
            variable_name,
            np.float32,
            ("time", "nj", "ni"),
            fill_value=np.float32(-9999.0),
            **comp_args,
        )

        # Parallelized read/transform per day; sequential write to NetCDF
        dst_time_index = 0
        landcells = grid_ids_1d.size
        valid_bool_flat = (~np.isnan(mask2d)).ravel()

        # Determine data variable name once from first file
        with nc.Dataset(filtered_files[0], "r") as src_probe:
            probe_var = variable_name if variable_name in src_probe.variables else find_data_variable(src_probe)
            if probe_var is None:
                return

        # Optional append mode: close file after define phase and reopen per-batch for writes
        # Default OFF for performance consistency with NETCDF3
        append_mode = os.environ.get("APPEND_MODE", "0") == "1"
        if append_mode:
            try:
                dst.sync()
            except Exception:
                pass
            dst.close()
            dst = None

        tasks = []
        remaining_total = total_time_steps
        for i, fpath in enumerate(filtered_files):
            if remaining_total <= 0:
                break
            steps = min(24, remaining_total)
            tasks.append((i, fpath, steps))
            remaining_total -= steps

        # Batch size (days) is tunable via env; default 5 days
        batch_days = int(os.environ.get("WRITE_BATCH_DAYS", "5"))
        batch_steps_target = max(1, batch_days) * 24
        # Optional periodic sync interval in steps (0 = never)
        sync_interval = int(os.environ.get("SYNC_INTERVAL_STEPS", "0"))

        if intra_workers is None or intra_workers <= 1:
            # Serial path with batching
            global G_VALID_MASK_FLAT, G_VAR_NAME, G_LANDCELLS
            G_VALID_MASK_FLAT = valid_bool_flat
            G_VAR_NAME = probe_var
            G_LANDCELLS = landcells

            t_idx = 0
            while t_idx < len(tasks):
                sub = []
                steps_in_batch = 0
                while t_idx < len(tasks) and (steps_in_batch + tasks[t_idx][2]) <= batch_steps_target:
                    sub.append(tasks[t_idx])
                    steps_in_batch += tasks[t_idx][2]
                    t_idx += 1
                if not sub:
                    sub.append(tasks[t_idx])
                    steps_in_batch += tasks[t_idx][2]
                    t_idx += 1
                # print the start time of the batch to the log file
                start_time = datetime.now()
                start_time_str = start_time.strftime("%Y-%m-%d %H:%M:%S")
                
                first_idx = sub[0][0]
                last_idx = sub[-1][0]
                first_name = os.path.basename(filtered_files[first_idx])
                last_name = os.path.basename(filtered_files[last_idx])
                print(f"[{variable_name} {year:04d}-{month:02d}] batch {first_idx+1}-{last_idx+1} reading: {first_name} ... {last_name} (steps={steps_in_batch}) started at {start_time_str}")
                try:
                    sys.stdout.flush()
                except Exception:
                    pass

                batch_chunks = []
                for idx, fpath, steps in sub:
                    _, chunk = _read_one((idx, fpath, steps))
                    if chunk is None:
                        continue
                    batch_chunks.append(chunk)
                if batch_chunks:
                    batch_arr = np.vstack(batch_chunks)
                    if append_mode:
                        with nc.Dataset(dst_name, "r+") as dsta:
                            dsta.variables[variable_name][dst_time_index : dst_time_index + batch_arr.shape[0], 0, :] = batch_arr
                            # Optional per-batch sync; default off
                            if int(os.environ.get("APPEND_SYNC", "0")) == 1:
                                dsta.sync()
                    else:
                        dst.variables[variable_name][dst_time_index : dst_time_index + batch_arr.shape[0], 0, :] = batch_arr
                        # defer sync to end or next periodic boundary
                    dst_time_index += batch_arr.shape[0]
                    # free memory eagerly
                    del batch_arr
                    batch_chunks.clear()
                    del batch_chunks
                    try:
                        gc.collect()
                    except Exception:
                        pass
                # end-of-batch timestamp and duration
                end_time = datetime.now()
                dur_s = (end_time - start_time).total_seconds()
                end_time_str = end_time.strftime("%Y-%m-%d %H:%M:%S")
                print(f"[{variable_name} {year:04d}-{month:02d}] batch {first_idx+1}-{last_idx+1} finished at {end_time_str} (elapsed {dur_s:.1f}s)")
                try:
                    sys.stdout.flush()
                except Exception:
                    pass

                if (not append_mode) and (sync_interval > 0) and (dst_time_index % sync_interval == 0):
                    try:
                        dst.sync()
                    except Exception:
                        pass
        else:
            # Parallel path with batching (spawn)
            ctx = get_context("spawn")
            with ctx.Pool(
                processes=intra_workers,
                initializer=_worker_init,
                initargs=(valid_bool_flat, probe_var, landcells),
                maxtasksperchild=1,
            ) as pool:
                t_idx = 0
                while t_idx < len(tasks):
                    sub = []
                    steps_in_batch = 0
                    while t_idx < len(tasks) and (steps_in_batch + tasks[t_idx][2]) <= batch_steps_target:
                        sub.append(tasks[t_idx])
                        steps_in_batch += tasks[t_idx][2]
                        t_idx += 1
                    if not sub:
                        sub.append(tasks[t_idx])
                        steps_in_batch += tasks[t_idx][2]
                        t_idx += 1

                    # batch start timestamp
                    start_time = datetime.now()
                    start_time_str = start_time.strftime("%Y-%m-%d %H:%M:%S")

                    first_idx = sub[0][0]
                    last_idx = sub[-1][0]
                    first_name = os.path.basename(filtered_files[first_idx])
                    last_name = os.path.basename(filtered_files[last_idx])
                    print(f"[{variable_name} {year:04d}-{month:02d}] batch {first_idx+1}-{last_idx+1} reading (parallel): {first_name} ... {last_name} (steps={steps_in_batch}) started at {start_time_str}")
                    try:
                        sys.stdout.flush()
                    except Exception:
                        pass

                    results = pool.map(_read_one, sub, chunksize=1)
                    batch_chunks = [chunk for _, chunk in results if chunk is not None]
                    if batch_chunks:
                        batch_arr = np.vstack(batch_chunks)
                        if append_mode:
                            with nc.Dataset(dst_name, "r+") as dsta:
                                dsta.variables[variable_name][dst_time_index : dst_time_index + batch_arr.shape[0], 0, :] = batch_arr
                                if int(os.environ.get("APPEND_SYNC", "0")) == 1:
                                    dsta.sync()
                        else:
                            dst.variables[variable_name][dst_time_index : dst_time_index + batch_arr.shape[0], 0, :] = batch_arr
                        dst_time_index += batch_arr.shape[0]
                        del batch_arr
                        batch_chunks.clear()
                        del batch_chunks
                        try:
                            gc.collect()
                        except Exception:
                            pass
                    # end-of-batch timestamp and duration
                    end_time = datetime.now()
                    dur_s = (end_time - start_time).total_seconds()
                    end_time_str = end_time.strftime("%Y-%m-%d %H:%M:%S")
                    print(f"[{variable_name} {year:04d}-{month:02d}] batch {first_idx+1}-{last_idx+1} finished at {end_time_str} (elapsed {dur_s:.1f}s)")
                    try:
                        sys.stdout.flush()
                    except Exception:
                        pass

                    if (not append_mode) and (sync_interval > 0) and (dst_time_index % sync_interval == 0):
                        try:
                            dst.sync()
                        except Exception:
                            pass
        # Final flush before close (only when file remained open)
        if not append_mode and dst is not None:
            try:
                dst.sync()
            except Exception:
                pass
    finally:
        if not append_mode and dst is not None:
            try:
                dst.close()
            except Exception:
                pass
        # Mark completion for this monthly file to enable upstream rsync+prune
        try:
            with open(dst_name + ".done", "w") as _f:
                _f.write("ok\n")
        except Exception:
            # Best-effort; lack of marker only affects asynchronous pruning
            pass


def collect_daily_files_by_var_month(input_root: str, var_filter: str | None = None, include_months: set[int] | None = None):
    """
    Walk input_root recursively and collect daily files by (variable, year, month).
    Returns dict: keys -> (variable, year, month), values -> list of file paths
    """
    groups = defaultdict(list)
    for root, _, files in os.walk(input_root):
        for fname in files:
            if not fname.endswith(".nc"):
                continue
            parsed = parse_daily_filename(fname)
            if parsed is None:
                continue
            variable_name, year, month, _ = parsed
            if var_filter is not None and variable_name != var_filter:
                continue
            if include_months is not None and month not in include_months:
                continue
            groups[(variable_name, year, month)].append(os.path.join(root, fname))
    return groups


def main():
    args = sys.argv[1:]
    if (len(args) < 2) or (args and args[0] == "--help"):
        print("Usage: python NA_ERAforcingGEN.py <input_root> <output_root> [time_limit] [max_daily_files_per_month] [variable_name] [months_csv] [intra_workers]")
        print("  <input_root>: Daymet_ERA5/netcdf_NA/<year> or Daymet_ERA5/netcdf_NA")
        print("  <output_root>: where monthly 1D files will be saved")
        print("  [time_limit]: optional int, limit total timesteps per monthly file (e.g., 1 for testing)")
        print("  [max_daily_files_per_month]: optional int, limit days processed per month (e.g., 1 for testing)")
        print("  [variable_name]: optional string, process only this variable (e.g., TBOT)")
        print("  [months_csv]: optional comma-separated months to include (e.g., 01,02,03)")
        sys.exit(0)

    input_root = args[0]
    output_root = args[1]
    time_limit = int(args[2]) if len(args) >= 3 else None
    max_daily_files_per_month = int(args[3]) if len(args) >= 4 else None
    var_filter = args[4] if len(args) >= 5 else None
    months_csv = args[5] if len(args) >= 6 else None
    intra_workers = int(args[6]) if len(args) >= 7 else None
    include_months = None
    if months_csv:
        try:
            include_months = set(int(m) for m in months_csv.split(",") if m.strip())
        except ValueError:
            include_months = None

    groups = collect_daily_files_by_var_month(input_root, var_filter=var_filter, include_months=include_months)
    if not groups:
        print(f"No valid daily files found under {input_root}")
        return

    # Process each variable-month group
    for (variable_name, year, month), files in sorted(groups.items()):
        # Determine mask and coordinates from the first available file
        first_file = sorted(files)[0]
        with nc.Dataset(first_file, "r") as src0:
            # Prefer explicit variable; fallback to detection
            data_var_name = variable_name if variable_name in src0.variables else find_data_variable(src0)
            if data_var_name is None:
                print(f"Skip {first_file}: cannot locate data variable")
                continue
        (
            mask2d,
            grid_ids_1d,
            lat_1d,
            lon_1d,
            x_dim,
            y_dim,
            _rows,
            _cols,
        ) = build_land_mask_and_coords(first_file, data_var_name)

        # Output directory per variable/year to keep files organized
        out_dir = os.path.join(output_root, variable_name, f"{year:04d}")
        write_monthly_file(
            out_dir,
            variable_name,
            year,
            month,
            files,
            mask2d,
            grid_ids_1d,
            lat_1d,
            lon_1d,
            x_dim,
            y_dim,
            time_limit,
            max_daily_files_per_month,
            intra_workers,
        )

    # Clean termination for batch wrapper
    try:
        sys.stdout.flush()
        sys.stderr.flush()
    except Exception:
        pass
    os._exit(0)

if __name__ == "__main__":
    main()


