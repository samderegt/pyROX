# Basic information on database and species
database = 'cia_hitran' # Can be ['cia_hitran', 'cia_borysow']
species  = 'h2-h2'      # Species name

# Input/output-directories
input_data_dir  = './examples/cia_hitran_h2_h2/input_data/'
output_data_dir = './examples/cia_hitran_h2_h2/'

# Instructions to download from HITRAN database
urls = [
    'https://hitran.org/data/CIA/H2-H2_2011.cia', 
    'https://hitran.org/data/CIA/supplementary/H2-H2_eq_2018.cia', # Extends to lower temperatures
]

# Wavenumber grid
wave_min = 0.3; wave_max = 250.0 # [um]
delta_nu = 10. # [cm^-1] # Can be sparse

# Extents of the CIA data
T_mask_1 = lambda T: (T>=400) & (T<=3000)
nu_mask_1 = lambda nu: (nu>=20) & (nu<=10000) # [cm^-1]

T_mask_2 = lambda T: (T>=40) & (T<=400)
nu_mask_2 = lambda nu: (nu>=0) & (nu<=2400) # [cm^-1]

# Input-data files
files = dict(
    # Last file in the list is prioritised at shared grid points
    cia = [
        (f'{input_data_dir}/H2-H2_2011.cia', T_mask_1, nu_mask_1),
        (f'{input_data_dir}/H2-H2_eq_2018.cia', T_mask_2, nu_mask_2),
    ]
)

# Metadata to be stored in pRT3's .h5 file
pRT3_metadata = dict(
    DOI = ['10.1021/jp109441f', '10.3847/1538-4365/aaa07a'], 
    mol_mass = [2.0,2.0],
    mol_name = ['H2','H2'],
)