A Python package for computing opacity cross-sections, using line lists from VALD/Kurucz, ExoMol or HITRAN/HITEMP. 

# Usage
Give the input python-file as the first positional argument when calling `main.py` (see examples of the input files in `input/`).

## Computing cross-sections
Give the temporary output-file as an argument when calculating the cross-sections. This file will be placed in the `tmp_output_dir`-directory specified in the input-file. 
```
python main.py input/XXX -cs --tmp_output_file cross.hdf5
```
The pressure-temperature points can also be defined on the command line (i.e. different from the input-file), which helps to parallelise computations. 
```
python main.py input/XXX -cs --P 1.0 --tmp_output_file cross.hdf5
```
will save a file (`cross_P1.0.hdf5`) at a pressure of 1 bar for the temperatures specified in the input-file. The `--T` argument can be used similarly for a temperature point.

### (Computing multiple ExoMol `.trans-files`)
If the line list is separated into multiple `.trans`-files, run between `i_min=5` and `i_max=10` with:
```
python main.py input/XXX -cs -i_min 5 -i_max 10 --tmp_output_file 'cross{}.hdf5'
```
This will save the temporary outputs in 'cross5.hdf5', 'cross6.hdf5', etc. 


## Combining cross-sections
To combine output-files in the `tmp_output_dir`-directory into a `final_output_file`:
```
python main.py input/XXX -s
```
This will try to combine different PT-grids, but will only work if a rectangular grid can be made. You could expand the PT-grid of an existing `final_output_file` by placing it in the `tmp_output_dir`, alongside more PT-points, and running the command above. 

## Converting to pRT2-format
```
python main.py input/XXX --convert_to_pRT2
```
will create the `short_stream/` opacity-directory for petitRADTRANS version 2. 

## Plotting
```
python main.py input/XXX --plot
```
or (to plot the pRT-converted opacities):
```
python figures.py
```

## Downloading data
```
python main.py input/XXX -d
```
- **VALD**: downloads the tab-separated NIST energy levels. Transitions need to be requested from the VALD database directly (https://www.astro.uu.se/valdwiki).
- **Kurucz**: downloads the tab-separated NIST energy levels and transitions (`gf*.pos`) from Robert Kurucz' website.
- **ExoMol**: downloads the `*.json` definition file, all transitions (`*.trans.bz2`), states (`*.states.bz2`) and partition function (`*.pf`) (see `inputs/example_ExoMol_NaH.py`).
- **HITRAN/HITEMP**: downloads files in the `urls` list, specified in the input-file (see `inputs/example_HITEMP.py`). 