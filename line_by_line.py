import numpy as np
from pandas import read_fwf, read_csv

import pathlib

from cross_sections import CrossSections
from utils import sc
import utils

class LineByLine(CrossSections):
    """
    Base class for line-by-line cross-sections.
    """

    def __init__(self, config, **kwargs):

        # Initialise the CrossSections parent class
        super().__init__(config, **kwargs)

    def _read_from_config(self, config):
        
        # Read the common parameters
        super()._read_from_config(config)
        
        # Read the parameters for line-by-line calculations
        self._read_mass()
        self._read_pressure_broadening_info()
        self._read_equation_of_state()
        #self._read_partition_function()

        self.N_lines_in_chunk = getattr(config, 'N_lines_in_chunk', 10_000_000)

        '''
        self.T_0 = getattr(conf, 'T_0', 296.)
        self.q_0 = np.interp(self.T_0, self.Q[:,0], self.Q[:,1])

        # Create wavenumber grid
        self._set_nu_grid(conf)
        self.adaptive_nu_grid = getattr(conf, 'adaptive_nu_grid', False)
        # Finer resolution to determine which lines to keep
        self.delta_nu_local_cutoff = getattr(conf, 'delta_nu_local_cutoff', 1e-3) # [cm^-1]
        
        # (P,T)-grid to compute cross-sections on
        self.P_grid = np.asarray(conf.P_grid) * 1e5 # [bar] -> [Pa]
        self.T_grid = np.asarray(conf.T_grid)
        self.N_PT = len(self.P_grid) * len(self.T_grid)

        # Cross-sections array
        self.sigma = np.zeros(
            (len(self.nu_grid), len(self.P_grid), len(self.T_grid))
            )

        # Line-wing cutoff
        self.wing_cutoff = getattr(
            conf, 'wing_cutoff', 
            lambda _, P: 25 if P<=200 else 100 # Gharib-Nezhad et al. (2024)
            )
        #self.cutoff_max = getattr(conf, 'cutoff_max', 25)
        self.cutoff_max = getattr(conf, 'cutoff_max', np.inf)

        # Line-strength cutoffs
        self.local_cutoff  = getattr(conf, 'local_cutoff', None)
        self.global_cutoff = getattr(conf, 'global_cutoff', None)
        '''

    def _read_equation_of_state(self):
        """
        Read the equation of state.
        """
        # Mean molecular weight of the atmosphere
        self.mean_mass = np.sum(
            [info['VMR']*info['mass'] for info in self.pressure_broadening_info.values()]
        )
        print(f'  Mean molecular weight using perturber info: {self.mean_mass/sc.amu:.2f} amu')

        EOS_table = getattr(self.config, 'EOS_table', None)
        if EOS_table is None:
            # Assume ideal gas
            self.get_number_density = lambda P, T: P/(sc.k*T)
            return

        # Load equation-of-state table and convert into SI
        from pandas import read_csv
        EOS_table = read_csv(EOS_table, delim_whitespace=True, skiprows=1)
        T = np.unique(EOS_table['T[K]'])
        P = np.unique(EOS_table['P[dyne_cm-2]'])
        P *= 1e-1 # [dyne cm^-2] -> [Pa]            

        mass_density_grid = EOS_table.pivot(
            index='P[dyne_cm-2]', columns='T[K]', values='rho[g_cm-3]'
            ).values
        number_density_grid = mass_density_grid*1e3 / self.mean_mass # [g cm^-3] -> [m^-3]
        
        # Create an interpolation function for the density grid
        from scipy.interpolate import RegularGridInterpolator
        interp_func = RegularGridInterpolator(
            points=(P, T), values=number_density_grid, method='cubic', 
            bounds_error=False, fill_value=None
            )

        # Define the number density function
        self.get_number_density = lambda P, T: interp_func([P, T])

    def _read_mass(self):
        """
        Read the mass of the species.
        """
        self.mass = getattr(self.config, 'mass', None)
        if (self.mass is None) and (self.database in ['hitran', 'exomol']):
            raise ValueError('Mass of the species must be provided in the configuration file.')
        elif (self.mass is None) and (self.database in ['vald', 'kurucz']):
            # Atomic species, read from table
            self.mass = self.atoms_info.loc[self.species, 'mass']
        self.mass *= sc.amu # [kg]

    def _read_pressure_broadening_info(self):
        """
        Read the pressure broadening parameters from the configuration file.
        """

        perturber_info = getattr(self.config, 'perturber_info', None)
        if perturber_info is None:
            if self.database in ['hitran','exomol']:
                raise ValueError('Perturber information must be provided in the configuration file.')

            # For atomic lines, use polarisabilities (Eq. 23, Sharp & Burrows 2007)
            perturber_info = dict(
                H2 = dict(VMR=0.85, alpha=0.806e-30), # [m^3]
                He = dict(VMR=0.15, alpha=0.204956e-30),
            )

        # Check that the perturber information is complete
        self.pressure_broadening_info = {}
        for perturber, info in perturber_info.items():

            # Add mass to the perturber_info dictionary if it doesn't exist
            mass = info.get('mass', None)
            if perturber == 'H2' and mass is None:
                
                mass = sc.m_H2 / sc.amu
            elif perturber == 'He' and mass is None:
                mass = sc.m_He / sc.amu
            elif mass is None:
                raise ValueError(f'Mass of perturber {perturber} must be provided (in amu) in the configuration file.')
            info['mass'] = mass*sc.amu # [kg]

            if 'VMR' not in info:
                raise ValueError(f'Volume mixing ratio of perturber {perturber} must be provided in the configuration file.')

            # Store the broadening parameters
            self.pressure_broadening_info[perturber] = info.copy()

            # Use gamma and temperature-exponent from dictionary (gamma*(T0/T)^n*(P/1bar))
            gamma = info.get('gamma', 0.)
            if callable(gamma):
                # Callable parameterisation (e.g. Gharib-Nezhad et al. 2024)
                self.pressure_broadening_info[perturber]['gamma'] = gamma
            else:
                self.pressure_broadening_info[perturber]['gamma'] = np.atleast_1d(gamma)
            self.pressure_broadening_info[perturber]['n'] = np.atleast_1d(info.get('n', 0.))

            # Load the broadening parameters from file, if given
            file = info.get('file', None)
            if file is None:
                continue
            broadening_params = read_fwf(file, header=None)
            diet = np.array(broadening_params[0], dtype=str)
            gamma, n = np.array(broadening_params[[1,2]].T)

            # Currently only handles these 2 ('a0', 'm0') broadening diets
            mask_diet = (diet == 'a0')
            if not mask_diet.any():
                mask_diet = (diet == 'm0')
                self.pressure_broadening_info[perturber]['diet'] = 'm0'
            
            self.pressure_broadening_info[perturber]['gamma'] = gamma[mask_diet]
            self.pressure_broadening_info[perturber]['n'] = n[mask_diet]

            if (len(gamma)==1) or broadening_params.shape[1]==3:
                # One row or no quantum-number columns, ignore quantum-number dependence
                self.pressure_broadening_info[perturber]['gamma'] = np.nanmean(
                    self.pressure_broadening_info[perturber]['gamma'], keepdims=True
                    )
                self.pressure_broadening_info[perturber]['n'] = np.nanmean(
                    self.pressure_broadening_info[perturber]['n'], keepdims=True
                    )
                continue

            # Total angular momentum quantum numbers
            self.pressure_broadening_info[perturber]['J'] = \
                np.array(broadening_params[3])[mask_diet]

    def _read_partition_function(self):
        raise NotImplementedError("This method should be implemented in the subclass.")

    def save_merged_outputs(self, **kwargs):
        """
        Merge the temporary files and save the final output. Same for all LineByLine classes.
        """
        raise NotImplementedError
    
    def plot_merged_outputs(self, **kwargs):
        """
        Plot the merged outputs. Same for all LineByLine classes.
        """
        raise NotImplementedError

    def convert_to_pRT2(self):
        raise NotImplementedError('Conversion to petitRADTRANS-v2.0 format not implemented.')

    def convert_to_pRT3(self):
        raise NotImplementedError('Conversion to petitRADTRANS-v3.0 format not implemented.')

class HITRAN(LineByLine):

    def download_data(self, config):
        """
        Download CIA data from HITRAN.
        """

        print('\nDownloading CIA data from HITRAN')

        for url in config.urls:
            file = utils.wget_if_not_exist(url, config.input_data_dir)
    
    def __init__(self, config, **kwargs):

        print('-'*60)
        print('  Line-by-line Absorption from HITRAN/HITEMP')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent LineByLine class

    def _read_from_config(self, config):
        
        # Read the common parameters
        super()._read_from_config(config)
        
        # Read the isotope information
        self.isotope_idx       = getattr(config, 'isotope_idx', 1)
        self.isotope_abundance = getattr(config, 'isotope_abundance', 1.0)
    
    def calculate_tmp_outputs(self, **kwargs):
        raise NotImplementedError

class ExoMol(LineByLine):

    def download_data(self, config):
        raise NotImplementedError
    
    def __init__(self, config, **kwargs):

        print('-'*60)
        print('  Line-by-line Absorption from ExoMol')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent LineByLine class

        raise NotImplementedError

    def calculate_tmp_outputs(self, **kwargs):
        raise NotImplementedError

class Kurucz(LineByLine):

    parent_dir = pathlib.Path(__file__).parent.resolve()
    atoms_info = read_csv(parent_dir/'atoms_info.csv', index_col=0)

    def download_data(self, config):
        raise NotImplementedError
    
    def __init__(self, config, **kwargs):

        print('-'*60)
        print('  Line-by-line Absorption from Kurucz')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent LineByLine class

        raise NotImplementedError
    
    def _read_from_config(self, config):
        
        # Read the common parameters
        super()._read_from_config(config)
        
        # If alkali, read the E_ion and Z for different vdW-broadening
        self.E_ion = getattr(config, 'E_ion', None) # [cm^-1]
        if self.E_ion is not None:
            self.E_ion *= 1e2 # [cm^-1] -> [m^-1]
        self.Z = getattr(config, 'Z', None)

        # Transition energies to ignore
        self.nu_0_to_ignore = np.atleast_1d(getattr(config, 'nu_0_to_ignore', []))
        self.nu_0_to_ignore *= 1e2*sc.c # [cm^-1] -> [s^-1]

    def calculate_tmp_outputs(self, **kwargs):
        raise NotImplementedError

class VALD(Kurucz):

    def download_data(self, config):
        raise NotImplementedError
    
    def __init__(self, config, **kwargs):

        print('-'*60)
        print('  Line-by-line Absorption from VALD')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent Kurucz class

        raise NotImplementedError
    
    def calculate_tmp_outputs(self, **kwargs):
        raise NotImplementedError