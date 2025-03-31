import numpy as np

database = 'exomol'
species  = 'h2o'
mass = 18.010565 # (in .json file)
isotopologue_id = {'H2':1, 'O':16}

# Instructions to download from ExoMol database
urls = [
    'https://www.exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL.def.json', 
    'https://www.exomol.com/db/H2O/1H2-16O/1H2-16O__H2.broad', 
    'https://www.exomol.com/db/H2O/1H2-16O/1H2-16O__He.broad',
    'https://www.exomol.com/db/H2O/1H2-16O/POKAZATEL/1H2-16O__POKAZATEL__04300-04400.trans.bz2'
]

# Input/output-directories
input_data_dir  = './examples/exomol_h2o/input_data/'
output_data_dir = './examples/exomol_h2o/'

files = dict(
    transitions = f'{input_data_dir}/1H2-16O__POKAZATEL__04300-04400.trans.bz2', 
    states      = f'{input_data_dir}/1H2-16O__POKAZATEL.states.bz2',

    partition_function = f'{input_data_dir}/1H2-16O__POKAZATEL.pf',
)

# Pressure-broadening information
perturber_info = dict(
    H2 = dict(
        VMR=0.85, file=f'{input_data_dir}/1H2-16O__H2.broad', # read from file
        ),
    He = dict(
        VMR=0.15, file=f'{input_data_dir}/1H2-16O__He.broad', # read from file
        ),
)

#P_grid = np.logspace(-5,2,8) # [bar]   # can be given in cmd, one point at a time
#T_grid = np.array([1000,2000])   # [K]     # can be given in cmd, one point at a time

wave_min = 2.299; wave_max = 2.304 # [um]
delta_nu = 0.01 # [cm^-1]

# Switch to sparser wavenumber grid for high broadening?
adaptive_nu_grid = True

# Line-strength cutoffs
local_cutoff  = 0.35
#global_cutoff = 1e-45

# gamma_V [cm^-1], P [bar]
wing_cutoff = lambda _, P: 25 if P<=200 else 100 # Gharib-Nezhad et al. (2024) DEFAULT