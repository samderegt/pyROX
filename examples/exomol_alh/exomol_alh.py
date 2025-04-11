# Basic information on database and species
database = 'exomol' # Can be ['exomol', 'hitran', 'hitemp', 'kurucz']
species  = 'alh'    # Species name
mass = 27.98948     # Can be found in *.def.json file

# Input/output-directories
input_data_dir  = './examples/exomol_alh/input_data/'
output_data_dir = './examples/exomol_alh/'

# Instructions to download from ExoMol database
urls = [
    'https://www.exomol.com/db/AlH/27Al-1H/AloHa/27Al-1H__AloHa.def.json', 
    'https://www.exomol.com/db/AlH/27Al-1H/27Al-1H__H2.broad', 
    'https://www.exomol.com/db/AlH/27Al-1H/27Al-1H__He.broad', 
]

# Input-data files
files = dict(
    transitions = f'{input_data_dir}/27Al-1H__AloHa.trans.bz2',
    states      = f'{input_data_dir}/27Al-1H__AloHa.states.bz2',
    partition_function = f'{input_data_dir}/27Al-1H__AloHa.pf',
)

import numpy as np
# Pressure and temperature grids
P_grid = np.logspace(-3, 2, 6) # [bar]
T_grid = np.array([1000,2000]) # [K]

# Wavenumber grid
wave_min = 0.3; wave_max = 28.0 # [um]
delta_nu = 0.01 # [cm^-1]
adaptive_nu_grid = True # Switch to sparser wavenumber grid for high broadening?

# Pressure-broadening information
perturber_info = dict(
    H2 = dict(VMR=0.85, file=f'{input_data_dir}/27Al-1H__H2.broad'), # Read from file
    He = dict(VMR=0.15, file=f'{input_data_dir}/27Al-1H__He.broad'),
)

# Line-strength cutoffs
global_cutoff = 1e-45 # [cm^-1 / (molecule cm^-2)]
local_cutoff  = 0.25

# Function with arguments gamma_V [cm^-1], and P [bar]
wing_cutoff = lambda gamma_V, P: 25 if P<=200 else 100 # Gharib-Nezhad et al. (2024)

# Metadata to be stored in pRT3's .h5 file
pRT3_metadata = dict(
    DOI = '10.1093/mnras/stad3802',      # DOI of the data
    mol_name = 'AlH',                    # Using the right capitalisation
    linelist = 'AloHa',                  # Line-list name, used in .h5 filename
    isotopologue_id = {'Al':27, 'H':1},  # Atomic number of each element
)