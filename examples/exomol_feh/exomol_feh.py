# Basic information on database and species
database = 'exomol' # Can be ['exomol', 'hitran', 'hitemp', 'kurucz]
species  = 'feh'    # Species name
mass = 56.942762    # Can be found in *.def.json file

# Input/output-directories
input_data_dir  = './examples/exomol_feh/input_data/'
output_data_dir = './examples/exomol_feh/'

# Instructions to download from ExoMol database
urls = [
    'https://www.exomol.com/db/FeH/56Fe-1H/MoLLIST/56Fe-1H__MoLLIST.def.json', 
]

# Input-data files
files = dict(
    transitions = f'{input_data_dir}/56Fe-1H__MoLLIST.trans.bz2',
    states      = f'{input_data_dir}/56Fe-1H__MoLLIST.states.bz2',
    partition_function = f'{input_data_dir}/56Fe-1H__MoLLIST.pf',
)

import numpy as np
# Pressure and temperature grids
P_grid = 10**np.array([-3.,0.,1.]) # [bar]
T_grid = np.array([1000,2000])   # [K]

# Wavenumber grid
wave_min = 0.3; wave_max = 28.0 # [um]
delta_nu = 0.01 # [cm^-1]
adaptive_nu_grid = True # Switch to sparser wavenumber grid for high broadening?

# Pressure-broadening information
perturber_info = dict(
    H2 = dict(VMR=0.85, gamma=0.07, n=0.5), # gamma = [cm^-1]
    He = dict(VMR=0.15, gamma=0.07, n=0.5),
)

# Line-strength cutoffs
global_cutoff = 1e-45 # [cm^-1 / (molecule cm^-2)]
local_cutoff  = 0.25

# Function with arguments gamma_V [cm^-1], and P [bar]
wing_cutoff = lambda gamma_V, P: 25 if P<=200 else 100 # Gharib-Nezhad et al. (2024)