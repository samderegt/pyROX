import numpy as np
from scipy.interpolate import interp1d

import pathlib

from utils import sc
from .cia import CIA

class CIA_Borysow(CIA):
    """
    Class for handling CIA data from Borysow.
    """
    def __init__(self, config, **kwargs):
        """
        Initialize the CIA_Borysow object.

        Parameters:
        config (object): Configuration object containing settings and file paths.
        **kwargs: Additional arguments for initialization.
        """
        
        print('-'*60)
        print('  Collision-Induced Absorption from Borysow')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent CIA class

    def _read_absorption_coefficients(self, file):
        """
        Read absorption coefficients from a Borysow CIA file.

        Parameters:
        file (str): Path to the Borysow CIA file.

        Returns:
        tuple: Temperature grid, absorption coefficients (k), and absorption coefficients (alpha).

        Raises:
        FileNotFoundError: If the specified file does not exist.
        """

        file = pathlib.Path(file)
        print(f'  Reading from \"{file}\"')

        if not file.exists():
            raise FileNotFoundError(f'File \"{file}\" not found.')

        # Read CIA data
        with open(file, 'r') as f:
            lines = f.readlines()

        # Infer column-widths
        import re
        col_0 = re.findall('\s+\S+', lines[3]) # Skip the header
        col_widths = [len(col) for col in col_0]
        
        # Infer the temperatures
        T_grid = [el.replace('K','') for el in lines[1].split()[1:]]
        T_grid = np.array(T_grid, dtype=np.float64)

        from pandas import read_fwf
        data = np.asarray(read_fwf(file, widths=col_widths, skiprows=3, header=None))
        nu_native = data[:,0] * 100.0*sc.c # [cm^-1] -> [s^-1]
        abs_coeff_alpha_native = data[:,1:] * 1e2 # [cm^-1 molecule^-2] -> [m^-1 molecule^-2]

        abs_coeff_alpha = np.zeros((len(T_grid), self.N_nu))
        for i, abs_coeff_alpha_i in enumerate(abs_coeff_alpha_native.T):

            # Remove any empty entries
            nu_i = nu_native[~np.isnan(abs_coeff_alpha_i)]
            abs_coeff_alpha_i = abs_coeff_alpha_i[~np.isnan(abs_coeff_alpha_i)]

            # Interpolate onto nu_grid
            interp_func = interp1d(
                nu_i, np.log10(abs_coeff_alpha_i), kind='linear', 
                fill_value=np.nan, bounds_error=False
                )
            abs_coeff_alpha[i] = 10**interp_func(self.nu_grid)
            abs_coeff_alpha[i] = np.nan_to_num(abs_coeff_alpha[i], nan=0.0)

        # [m^-1 molecule^-2] -> [m^5 molecule^-2]
        abs_coeff_k = abs_coeff_alpha / sc.L0**2

        return T_grid, abs_coeff_k, abs_coeff_alpha