import numpy as np

database = 'kurucz'
species  = 'k'
# isotopologue_id = {'H2':1, 'O':16}

# Input/output-directories
input_data_dir  = './examples/kurucz_k/input_data/'
output_data_dir = './examples/kurucz_k/'

files = dict(
    transitions = f'{input_data_dir}/gf1900.pos', 
    states = f'{input_data_dir}/NIST_states.tsv',
)

# Pressure-broadening information
P_grid = np.logspace(-5,2,8) # [bar]
T_grid = np.array([1000,2000])   # [K]

# Wavenumber/wavelength grid
wave_min = 1.0/3.0; wave_max = 50.0 # [um]
delta_nu = 0.01 # [cm^-1]

wing_cutoff = lambda *_: 1500