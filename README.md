[![Documentation Status](https://readthedocs.org/projects/py-rox/badge/?version=latest)](https://py-rox.readthedocs.io/en/latest/)
[![License](https://img.shields.io/github/license/samderegt/pyROX)](https://github.com/samderegt/pyROX/blob/main/LICENSE)
[![Python Tests](https://github.com/samderegt/pyROX/actions/workflows/python-tests.yml/badge.svg)](https://github.com/samderegt/pyROX/actions/workflows/python-tests.yml)

# pyROX: Rapid Opacity X-sections
A Python package for computing opacity cross-sections and collision-induced absorption coefficients to be used in applications of (exo)-planetary and (sub)-stellar atmospheres. pyROX currently supports calculating opacities from the Kurucz, ExoMol, HITRAN/HITEMP and Borysow databases.

## Documentation
The full documentation for pyROX, including an installation guide and tutorial, can be found on [https://py-rox.readthedocs.io/en/latest/](https://py-rox.readthedocs.io/en/latest/). 

## Installation
To install pyROX, clone this repository and install via
```bash
git clone https://github.com/samderegt/pyROX.git
cd pyROX
pip install .
```

## Usage
After installation, pyROX can be called directly from the terminal using the `pyROX` command. The first argument should be a should be a python-script with configuration parameters (see [`examples/`](examples/)). 

Important steps can be run with the following arguments:
- `--download` (or `-d`): Download data from the specified database.
- `--calculate` (or `-c`): Calculate and save the temporary outputs (i.e. cross-sections or CIA-coefficients).
- `--save` (or `-s`): Combine the temporary outputs into a single grid and save to a final output-file. 

Example usage:
```bash
pyROX <config_file> --download --calculate --save
```

## License and Acknowledgements
pyROX is available under the [MIT license](LICENSE).