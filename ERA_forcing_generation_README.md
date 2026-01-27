## NA forcing generation (single year)

This directory contains the single-year Daymet_ERA5 forcing generation workflow.
It converts daily Daymet_ERA5 NetCDFs into monthly 1D forcing files for a chosen year.

### Inputs

- Daily NetCDFs under `Daymet_ERA5/netcdf_NA/<YEAR>/`
  - Expected pattern: `clmforc.Daymet4.1km.<VAR>.<YYYYMMDD>.nc`
- Script entrypoint: `NA_ERAforcingGEN.py`

### Outputs

- Monthly forcing files written to:
  - `/gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/NA_cases_data/NADaymet_ERA5_TEST_<YEAR>/forcing/<VAR>/<YEAR>/`

### Direct python usage

```
python NA_ERAforcingGEN.py <input_root> <output_root> [time_limit] [max_days] [variable_name] [months_csv] [intra_workers]
```

Example:

```
/gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/python_test_env/conda_envs/testvenv/bin/python \
  NA_ERAforcingGEN.py ./Daymet_ERA5/netcdf_NA/2014 \
  /gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/NA_cases_data/NADaymet_ERA5_TEST_2014_monthly \
  1 1 TBOT
```

### Slurm single-year multi-variable workflow

Use `NA_ERAforcingGEN_sinlgeyear_multivariables_local_nnode.sub` to process multiple variables for one year.
Each variable is assigned to a node; within each node, 4 tasks split months.

Example submit:

```
sbatch -N 3 -t 04:00:00 \
  --export=ALL,YEAR=2024,VARS_STR="TBOT FLDS FSDS QBOT PRECTmms WIND PSRF",TIME_LIMIT=-1,PROCS=4,INTRA_WORKERS=2 \
  NA_ERAforcingGEN_sinlgeyear_multivariables_local_nnode.sub
```

Key environment variables (with defaults in the job script):

- `YEAR` (default `2024`)
- `VARS_STR` (default `TBOT FLDS FSDS QBOT PRECTmms WIND PSRF`)
- `PROCS` (default `4`, allowed `4` or `6`)
- `INTRA_WORKERS` (default `2`)
- `TIME_LIMIT` (default `-1` for full hourly series)
- `MAX_DAYS` (default `-1` for all non-leap days)
- Performance knobs: `WRITE_BATCH_DAYS`, `APPEND_MODE`, `SYNC_INTERVAL_STEPS`, `APPEND_SYNC`

The job writes to node-local scratch first and rsyncs to the final GPFS output.

### Convenience submission helpers

- `command.txt` contains example `sbatch` commands for single-year runs.
- `submit_years_2004_2022.sh` submits a year-by-year chain (afterok dependency) if you need a sequence.

