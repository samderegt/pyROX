import numpy as np

database = 'cia_hitran'
species  = 'h2-h2'

# Instructions to download
urls = [
    'https://hitran.org/data/CIA/H2-H2_2011.cia', 
    'https://hitran.org/data/CIA/supplementary/H2-H2_eq_2018.cia', 
    #'https://hitran.org/data/CIA/supplementary/H2-H2_norm_2018.cia', 
]
# Input/output-directories
input_data_dir  = f'./h2-h2_HITRAN/input_data/'
output_data_dir = f'./h2-h2_HITRAN/'

files = dict(
    cia = [
        # Whichever file is last in the list will be prioritised at shared grid points
        (f'{input_data_dir}/H2-H2_2011.cia',    lambda T: (T>=200)&(T<=3000), lambda nu: (nu>=20)&(nu<=10000)),
        (f'{input_data_dir}/H2-H2_eq_2018.cia', lambda T: (T>=40)&(T<=400),   lambda nu: (nu>=0)&(nu<=2400)), 
    ], 
)

tmp_output_file = 'cia_coeffs_{}.hdf5'

# Wavenumber/wavelength grid
wave_min = 1.0/3.0; wave_max = 250.0 # [um]
delta_nu = 10 # [cm^-1]

pRT3_metadata = dict(
    DOI = ['10.1021/jp109441f', '10.3847/1538-4365/aaa07a'], 
    mol_mass = [2.0,2.0],
    mol_name = ['H2','H2'],
)