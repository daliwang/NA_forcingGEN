# NA_forcingGEN

## Overview

`NA_forcingGEN` is a Python-based tool designed to convert 2D North American (NA) forcing data into 1D format. This tool processes NetCDF files, specifically handling Daymet4 data, to generate forcing data for climate and land surface models like CLM (Community Land Model).

## Purpose

The primary purpose of this tool is to transform gridded 2D forcing data over North America into a 1D format, focusing on land grid cells. This conversion facilitates the use of the data in specific modeling frameworks by reducing dimensionality while preserving essential meteorological variables.

## Repository Structure

- `AOI_forcingGEN.py`: Script for generating forcing data for a specific Area of Interest (AOI).
- `NA_forcingGEN.py`: Base script for North American forcing data conversion.
- `NA_forcingGENv3.py`: Version 3 of the NA forcing data generator with enhancements.
- `NA_forcingGENv4.py`: Version 4 of the NA forcing data generator with further improvements.

## Usage

To use `NA_forcingGEN`, you need to run the script with appropriate command-line arguments to specify input and output paths along with the time steps to process.

### Prerequisites

- Python 3.x
- Required libraries: `os`, `sys`, `glob`, `math`, `netCDF4`, `numpy`

### Running the Script

```bash
python NA_forcingGENv3.py <input_path> <output_path> <time_steps>
```

- `<input_path>`: Path to the directory containing the 2D source NetCDF data files.
- `<output_path>`: Path where the converted 1D forcing data will be saved.
- `<time_steps>`: Number of time steps to process, or `-1` to process all time series data.

### Example

```bash
python NA_forcingGENv3.py /path/to/input/data /path/to/output/data -1
```

This command processes all time steps from the input directory and saves the 1D forcing data to the specified output directory.

## Data Processing

The script handles the following tasks:
- Reads 2D NetCDF files containing meteorological variables.
- Removes leap days if present (e.g., converting 232 days to 224 days).
- Converts time units to a format relative to the start of the current month.
- Applies a land mask to focus on land grid cells, creating 1D arrays for grid IDs, latitude, and longitude.
- Processes data in chunks to manage memory usage for large datasets.
- Outputs the processed data into new NetCDF files with a 1D structure.

## Version Information

- **Version Date**: March 31, 2025 (as per `NA_forcingGENv3.py`)

## Notes

- Ensure that the input NetCDF files are in a compatible format (`NETCDF3_64BIT_DATA`).
- The script processes files in subdirectories, maintaining the directory structure in the output path.
- Performance metrics are printed for each file conversion to monitor processing time.

## Contact

For questions or contributions, please contact the repository maintainer or open an issue in the repository. 