# Basic information on database and species
database = 'hitemp' # Can be ['exomol', 'hitran', 'hitemp', 'kurucz]
species  = 'co'     # Species name
mass = 28.01

# [HITRAN-specific]: .par-file includes all isotopologues
isotope_idx = 1
isotope_abundance = 9.86544e-1

# Input/output-directories
input_data_dir  = './input_data/'
output_data_dir = './'

# Instructions to download from HITEMP database
urls = [
    # Transitions (see https://hitran.org/hitemp/)
    'https://hitran.org/hitemp/data/bzip2format/05_HITEMP2019.par.bz2', 

    # Partition function (see https://hitran.org/docs/iso-meta/)
    'https://hitran.org/data/Q/q26.txt', 

    # Broadening parameters (from ExoMol)
    'https://www.exomol.com/db/CO/12C-16O/12C-16O__H2.broad', 
    'https://www.exomol.com/db/CO/12C-16O/12C-16O__He.broad', 
]

# Input-data files
files = dict(
    transitions        = f'{input_data_dir}/05_HITEMP2019.par.bz2', 
    partition_function = f'{input_data_dir}/q26.txt', 
)

import numpy as np
# Pressure and temperature grids
P_grid = np.logspace(-5,2,8) # [bar]
T_grid = np.array([500,1000,2000,3000]) # [K]

# Wavenumber grid
wave_min = 0.3; wave_max = 28.0 # [um]
resolution = 1e6
adaptive_nu_grid = False # Doesn't work for fixed resolution

# Pressure-broadening information
perturber_info = dict(
    H2 = dict(VMR=0.85, file=f'{input_data_dir}/12C-16O__H2.broad'), # Read from a file
    He = dict(VMR=0.15, file=f'{input_data_dir}/12C-16O__He.broad'), 
)

# Line-strength cutoffs
global_cutoff = 1e-45 # [cm^-1 / (molecule cm^-2)]
local_cutoff  = 0.25

# Function with arguments gamma_V [cm^-1], and P [bar]
wing_cutoff = lambda gamma_V, P: 25 if P<=200 else 100 # Gharib-Nezhad et al. (2024)

# Metadata to be stored in pRT3's .h5 file
pRT3_metadata = dict(
    DOI = '10.1088/0067-0049/216/1/15', # DOI of the data
    mol_name = 'CO',                    # Using the right capitalisation
    isotopologue_id = {'C':12, 'O':16}, # Atomic number of each element
)