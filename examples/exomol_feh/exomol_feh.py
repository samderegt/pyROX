import numpy as np

database = 'exomol'
species  = 'feh'
mass = 56.942762 # (in .json file)

# Instructions to download from ExoMol database
urls = [
    'https://www.exomol.com/db/FeH/56Fe-1H/MoLLIST/56Fe-1H__MoLLIST.def.json', 
]

# Input/output-directories
input_data_dir  = './examples/exomol_feh/input_data/'
output_data_dir = './examples/exomol_feh/'

files = dict(
    transitions = f'{input_data_dir}/56Fe-1H__MoLLIST.trans.bz2',
    states      = f'{input_data_dir}/56Fe-1H__MoLLIST.states.bz2',

    partition_function = f'{input_data_dir}/56Fe-1H__MoLLIST.pf',
)

# Pressure-broadening information
from pyROX.utils import Broaden_Gharib_Nezhad_ea_2021 as B
perturber_info = dict(
    H2 = dict(
        #VMR=0.85, gamma=B('FeH').gamma_H2, n=0.5
        VMR=0.85, gamma=0.07, n=0.5
        ),
    He = dict(
        #VMR=0.15, gamma=B('FeH').gamma_He, n=0.5
        VMR=0.15, gamma=0.07, n=0.5
        ),
)

#P_grid = np.logspace(-5,2,8) # [bar]
P_grid = 10**np.array([-5.,-1.]) # [bar]
T_grid = np.array([2000,3000])   # [K]

#wave_min = 1.0/3.0; wave_max = 50.0 # [um]
wave_min = 0.3; wave_max = 28.0 # [um]
#delta_nu = 0.01 # [cm^-1]
#delta_wave = 1e-6 # [um]
resolution = 1e6

# Switch to sparser wavenumber grid for high broadening?
#adaptive_nu_grid = True
adaptive_nu_grid = False

# Line-strength cutoffs
local_cutoff  = 0.25
global_cutoff = 1e-45

# gamma_V [cm^-1], P [bar]
#wing_cutoff = lambda _, P: 25 if P<=200 else 100 # Gharib-Nezhad et al. (2024) DEFAULT

pRT3_metadata = dict(
    DOI = '10.1016/j.jqsrt.2019.106687', # DOI of the data
    mol_name = 'FeH',
    linelist = 'MoLLIST',
    isotopologue_id = {'Fe':56, 'H':1}, 
)