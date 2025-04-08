import numpy as np
from scipy.interpolate import interp1d

import pathlib

from pyROX import sc
from .cia import CIA

class CIA_HITRAN(CIA):
    """
    Class for handling CIA data from HITRAN.
    """    
    def __init__(self, config, **kwargs):
        """
        Initialise the CIA_HITRAN object.

        Args:
            config (object): Configuration object containing settings and file paths.
            **kwargs: Additional arguments for initialisation.
        """

        print('-'*60)
        print('  Collision-Induced Absorption from HITRAN')
        print('-'*60+'\n')
        
        super().__init__(config, **kwargs) # Initialise the parent CIA class

    def _read_absorption_coefficients(self, file):
        """
        Reads absorption coefficients from a HITRAN CIA file.

        Args:
            file (str): Path to the HITRAN CIA file.

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

        # Select the headers
        N_native_grid, T_grid = [], []
        for line in lines:
            if len(line) < 101:
                continue

            # Header line
            N_native_grid.append(int(line[40:40+7]))
            T_grid.append(float(line[47:47+7]))

        N_native_grid = np.array(N_native_grid)
        T_grid   = np.array(T_grid)

        # Read absorption coefficients for each temperature
        abs_coeff_k = np.zeros((len(T_grid), self.N_nu))
        for i, (N_native_grid_i, T_i) in enumerate(zip(N_native_grid, T_grid)):
            
            idx_min = 1 + i*(N_native_grid_i+1)
            idx_max = idx_min + N_native_grid_i
            selected_lines = lines[idx_min:idx_max]

            # Convert wavenumber [cm^-1] to frequency [s^-1]
            nu_i = np.array([line[0:10].strip() for line in selected_lines], dtype=np.float64)
            nu_i *= 100.0*sc.c # [cm^-1] -> [s^-1]

            # Convert abs. coefficients to SI units
            abs_coeff_k_i = np.array([line[10:21].strip() for line in selected_lines], dtype=np.float64)
            abs_coeff_k_i *= (1e-2)**5 # [cm^5 molecule^-2] -> [m^5 molecule^-2]

            # Interpolate onto nu_grid
            interp_func = interp1d(
                nu_i, np.log10(abs_coeff_k_i), kind='linear', fill_value=np.nan, bounds_error=False
                )
            abs_coeff_k[i] = 10**interp_func(self.nu_grid)
            abs_coeff_k[i] = np.nan_to_num(abs_coeff_k[i], nan=0.0)

        # [m^5 molecule^-2] * [m^-6] = [m^-1 molecule^-2]
        abs_coeff_alpha = abs_coeff_k * sc.L0**2

        return T_grid, abs_coeff_k, abs_coeff_alpha