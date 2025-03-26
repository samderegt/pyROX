import numpy as np

database = 'cia_borysow'
species  = 'h2-h2'

# Input/output-directories
input_data_dir  = f'./h2-h2_Borysow/input_data/'
output_data_dir = f'./h2-h2_Borysow/'

files = dict(
    cia = [
        # File, temperature range, wavenumber range
        (f'{input_data_dir}/final_CIA_LT.dat', lambda T: (T>=60)&(T<=350),    lambda nu: (nu>=0)&(nu<=15000)),
        (f'{input_data_dir}/final_CIA_HT.dat', lambda T: (T>=400)&(T<=1000),  lambda nu: (nu>=0)&(nu<=17000)),
        (f'{input_data_dir}/CIA.H2H2.Yi',      lambda T: (T>=1000)&(T<=7000), lambda nu: (nu>=20)&(nu<=20000)),
    ], 
)

tmp_output_file = 'cia_coeffs_{}.hdf5'

# Wavenumber/wavelength grid
wave_min = 1.0/3.0; wave_max = 250.0 # [um]
delta_nu = 10 # [cm^-1]

pRT3_metadata = dict(
    DOI = ['10.1051/0004-6361:20020555', '10.1016/S0022-4073(00)00023-6'], 
    mol_mass = [2.0,2.0],
    mol_name = ['H2','H2'],
)