import numpy as np
from utils import sc

class LineByLine(CrossSections):
    """
    Base class for line-by-line cross-sections.
    """

    def __init__(self, config):
        # Initialise the CrossSections parent class
        super().__init__(config)

        # Mean molecular weight of the atmosphere
        self.VMR_H2 = getattr(self.config, 'VMR_H2', 0.74)
        self.VMR_He = getattr(self.config, 'VMR_He', 0.24)

        # Mean molecular weight of the atmosphere
        # TODO: make it possible for non-H2/He atmospheres
        self.mean_mass = self.VMR_H2*self.m_H2 + self.VMR_He*self.m_He

        self._load_equation_of_state()
        self._load_pressure_broadening()

        self.N_lines_in_chunk = getattr(self.config, 'N_lines_in_chunk', 10_000_000)

    def _load_equation_of_state(self):
        """
        Load the equation of state from the configuration file.
        """

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

    def _load_pressure_broadening(self):
        """
        Load the pressure broadening parameters from the configuration file.
        """
        raise NotImplementedError

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

    def download_data(self):
        raise NotImplementedError
    
    def __init__(self, config):

        print('-'*60)
        print('  Line-by-line Absorption from HITRAN/HITEMP')
        print('-'*60+'\n')

        super().__init__(config) # Initialise the parent LineByLine class

        raise NotImplementedError
    
    def calculate_tmp_outputs(self, **kwargs):
        raise NotImplementedError

class ExoMol(LineByLine):

    def download_data(self):
        raise NotImplementedError
    
    def __init__(self, config):

        print('-'*60)
        print('  Line-by-line Absorption from ExoMol')
        print('-'*60+'\n')

        super().__init__(config) # Initialise the parent LineByLine class

        raise NotImplementedError
    
    def calculate_tmp_outputs(self, **kwargs):
        raise NotImplementedError

class Kurucz(LineByLine):

    def download_data(self):
        raise NotImplementedError
    
    def __init__(self, config):

        print('-'*60)
        print('  Line-by-line Absorption from Kurucz')
        print('-'*60+'\n')

        super().__init__(config) # Initialise the parent LineByLine class

        raise NotImplementedError
    
    def calculate_tmp_outputs(self, **kwargs):
        raise NotImplementedError

class VALD(Kurucz):

    def download_data(self):
        raise NotImplementedError
    
    def __init__(self, config):

        print('-'*60)
        print('  Line-by-line Absorption from VALD')
        print('-'*60+'\n')

        super().__init__(config) # Initialise the parent Kurucz class

        raise NotImplementedError
    
    def calculate_tmp_outputs(self, **kwargs):
        raise NotImplementedError