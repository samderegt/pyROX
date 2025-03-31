import numpy as np

database = 'hitemp'

species  = 'co'
mass = 28.01
isotope_idx = 1 # !! HITEMP specific !! .par-file includes all isotopologues
isotope_abundance = 9.86544e-1
#isotopologue_id = {'C':12, 'O':16}

# Instructions to download from HITRAN/HITEMP database
urls = [
    # Transitions (see https://hitran.org/hitemp/)
    'https://hitran.org/hitemp/data/bzip2format/05_HITEMP2019.par.bz2', 

    # Partition function (see https://hitran.org/docs/iso-meta/)
    'https://hitran.org/data/Q/q26.txt', 

    # Broadening parameters (from ExoMol)
    'https://www.exomol.com/db/CO/12C-16O/12C-16O__H2.broad', 
    'https://www.exomol.com/db/CO/12C-16O/12C-16O__He.broad', 
]

# Input/output-directories
input_data_dir  = './examples/hitemp_co/input_data/'
output_data_dir = './examples/hitemp_co/'

files = dict(
    transitions = f'{input_data_dir}/05_HITEMP2019.par.bz2', 
    partition_function = f'{input_data_dir}/q26.txt', 
)

perturber_info = dict(
    H2 = dict(
        VMR=0.85, file=f'{input_data_dir}/12C-16O__H2.broad', # read from file
        ), 
    He = dict(
        VMR=0.15, file=f'{input_data_dir}/12C-16O__He.broad', # read from file
        ), 
)

P_grid = np.logspace(-5,2,8) # [bar]
#T_grid = np.array([100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500,5000])   # [K]     # can be given in cmd, one point at a time
T_grid = np.array([1000,2000]) # [K]

wave_min = 1.0; wave_max = 5.0 # [um]
delta_nu = 0.01 # [cm^-1]

# Switch to sparser wavenumber grid for high broadening?
adaptive_nu_grid = True

# Line-strength cutoffs
local_cutoff  = 0.25
global_cutoff = 1e-45

# gamma_V [cm^-1], P [bar]
#wing_cutoff = lambda _, P: 25 if P<=200 else 100 # Gharib-Nezhad et al. (2024) DEFAULT
