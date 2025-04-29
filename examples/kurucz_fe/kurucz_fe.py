# Basic information on database and species
database = 'kurucz' # Can be ['exomol', 'hitran', 'hitemp', 'kurucz']
species  = 'fe'      # Species name
# mass = 55.845     # Will be read from atoms_info.csv

# Input/output-directories
input_data_dir  = './examples/kurucz_fe/input_data/'
output_data_dir = './examples/kurucz_fe/'

# Downloading can be done automaticall from the Kurucz database
# urls = []

# Input-data files
files = dict(
    transitions = f'{input_data_dir}/gf2600.all',
    states      = f'{input_data_dir}/NIST_states.tsv',
)

import numpy as np
# Pressure and temperature grids
P_grid = np.logspace(-3, 2, 6) # [bar]
T_grid = np.array([1000,2000]) # [K]

# Wavenumber grid
wave_min = 0.3; wave_max = 28.0 # [um]
delta_nu = 0.01 # [cm^-1]

# Pressure-broadening information
# perturber_info = {} # Fe is an atom, defaults to Sharp & Burrows (2007) description

# Function with arguments gamma_V [cm^-1], and P [bar]
wing_cutoff = lambda *_: 25

# Metadata to be stored in pRT3's .h5 file
pRT3_metadata = dict(
    DOI = 'http://kurucz.harvard.edu/atoms/', # DOI of the data
    mol_name = 'Fe',                           # Using the right capitalisation
    isotopologue_id = {'Fe':26},               # Atomic number of each element
)