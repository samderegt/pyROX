# Basic information on database and species
database = 'exomol' # Can be ['exomol', 'hitran', 'hitemp', 'kurucz']
species  = 'h2s'    # Species name
mass = 33.987721    # Can be found in *.def.json file

# Input/output-directories
input_data_dir  = './input_data/'
output_data_dir = './'

# Instructions to download from ExoMol database
urls = [
    'https://www.exomol.com/db/H2S/1H2-32S/AYT2/1H2-32S__AYT2.def.json', 
]

# Input-data files
files = dict(
    #transitions = f'{input_data_dir}/1H2-32S__AYT2.trans.bz2',
    transitions = [
        f'{input_data_dir}/1H2-32S__AYT2__{nu_min:05d}-{nu_min+1000:05d}.trans.bz2'
        for nu_min in range(0,35000,1000)
    ], 
    states      = f'{input_data_dir}/1H2-32S__AYT2.states.bz2',
    partition_function = f'{input_data_dir}/1H2-32S__AYT2.pf',
)

import numpy as np
# Pressure and temperature grids
P_grid = 10**np.array([-1.,0.]) # [bar]
T_grid = np.array([1000,2000])   # [K]

# Wavenumber grid
wave_min = 0.3; wave_max = 28.0 # [um]
delta_nu = 0.01 # [cm^-1]
adaptive_nu_grid = True

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

# Metadata to be stored in pRT3's .h5 file
pRT3_metadata = dict(
    DOI = ['10.1093/mnras/stw1133','10.1016/j.jqsrt.2018.07.012'], # DOI of the data
    mol_name = 'H2S',                    # Using the right capitalisation
    linelist = 'AYT2',                   # Line-list name, used in .h5 filename
    isotopologue_id = {'H2':1, 'S':32},  # Atomic number of each element
)