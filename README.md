# pyROX: Rapid Opacity X-sections

A Python package for computing opacity cross-sections using line lists from Kurucz, ExoMol or HITRAN/HITEMP. In addition, the package supports calculations of collision-induced absorption (CIA) coefficients from HITRAN or Borysow. 

# Usage
Give the python-script with configuration parameters as the first positional argument when calling `main.py` (see examples of these input-files in [`examples/`](examples/)).

Important steps are run with the following arguments:
- `--download` (or `-d`): Download data from the specified database.
- `--calculate` (or `-c`): Calculate and save the temporary outputs (i.e. cross-sections, CIA-coefficients).
- `--save` (or `-s`): Combine the temporary outputs into a single grid and save to a final output-file. 

Given as:
```
python main.py <config_file> -d -c -s
```
