[![Documentation Status](https://readthedocs.org/projects/py-rox/badge/?version=latest)](https://py-rox.readthedocs.io/en/latest/)
[![License](https://img.shields.io/github/license/samderegt/pyROX)](https://github.com/samderegt/pyROX/blob/main/LICENSE)
[![Python Tests](https://github.com/samderegt/pyROX/actions/workflows/python-tests.yml/badge.svg)](https://github.com/samderegt/pyROX/actions/workflows/python-tests.yml)
[![JOSS Status](https://joss.theoj.org/papers/8f2ed067c20eaa7090b1dc1bfbbe0f54/status.svg)](https://joss.theoj.org/papers/8f2ed067c20eaa7090b1dc1bfbbe0f54)

# pyROX: Rapid Opacity X-sections
A Python package for computing opacity cross-sections and collision-induced absorption coefficients for applications in models of (exo)-planetary and (sub)-stellar atmospheres. pyROX currently supports line opacity calculations from the ExoMol, HITRAN, HITEMP, and Kurucz databases. Collision-induced absorption coefficients can also be computed using the HITRAN and Borysow databases. 

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

## Contributing
If you have suggestions for new features or found a potential bug, please [open an issue](https://github.com/samderegt/pyROX/issues). If you would like to implement a new feature, please read the [contribution guidelines](https://py-rox.readthedocs.io/en/latest/contributing.html).

## License and Acknowledgements
pyROX is available under the [MIT license](LICENSE).