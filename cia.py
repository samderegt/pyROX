import numpy as np
from scipy.interpolate import interp1d

import pathlib

from cross_sections import CrossSections
from utils import sc
import utils

class CIA_HITRAN(CrossSections):

    def download_data(self):
        """
        Download CIA data from HITRAN.
        """

        print('\nDownloading CIA data from HITRAN')

        files = []
        for url in self.config.urls:
            file = utils.wget_if_not_exist(url, self.config.input_dir)
            files.append(file)
        
        print()
        for file in files:
            print(file)
                
    def __init__(self, config):
        
        print('\nInitialising CIA_HITRAN')
        super().__init__(config)

    def _read_abs_coeff(self, file):
        
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
        abs_coeff_k = np.zeros((len(T_grid), self.N_grid))
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

    def calculate_tmp_outputs(self, **kwargs):

        files = getattr(self.config, 'files', None)
        if files is None:
            raise ValueError('No files specified in the configuration.')
        cia_files = files.get('cia', None)
        if cia_files is None:
            raise ValueError('No CIA files specified in the configuration.')

        self.T_grid, self.abs_coeff_k, self.abs_coeff_alpha = [], [], []
        for i, (file, *masks) in enumerate(cia_files):
            
            print(f'\nReading CIA data from file \"{file}\"')

            # Compute absorption coefficients
            T_grid, abs_coeff_k, abs_coeff_alpha = self._read_abs_coeff(file)

            if len(masks) != 0:
                # Remove temperatures outside the range
                mask_T = masks[0](T_grid)
                T_grid, abs_coeff_k, abs_coeff_alpha = self.mask_arrays(
                    [T_grid, abs_coeff_k, abs_coeff_alpha], mask=mask_T, axis=0
                    )

                # Remove wavenumbers outside the range
                mask_nu = masks[1](self.nu_grid/(100.0*sc.c))
                abs_coeff_k[:,~mask_nu]     = 0.
                abs_coeff_alpha[:,~mask_nu] = 0.

            # Temporarily save the data
            tmp_output_file = str(self.tmp_output_file)
            tmp_output_file = tmp_output_file.format(pathlib.Path(file).name)

            utils.save_to_hdf5(
                tmp_output_file, 
                data={
                    'wave': self.wave_grid, 
                    'T': T_grid, 
                    'k': abs_coeff_k.T, 
                    'alpha': abs_coeff_alpha.T
                    },
                attrs={
                    'wave': {'units': 'm'}, 
                    'T': {'units': 'K'}, 
                    'k': {'units': 'm^5 molecule^-2'}, 
                    'alpha': {'units': 'm^-1 molecule^-2'}
                    }
                )

    def save_merged_outputs(self, *args, **kwargs):

        if self.final_output_file.exists():
            response = input(f"Warning: Final output file '{self.final_output_file}' already exists. Do you want to overwrite it? (yes/no): ")
            if response.lower() not in ['y', 'yes']:
                raise FileExistsError(f"Not overwriting final output file '{self.final_output_file}'.")

        # Merge the temporary files
        self.merge_tmp_outputs(keys_to_merge=['k', 'alpha'])

        # Convert to log10
        # self.merged_datasets['k'] = np.around(np.log10(self.merged_datasets['k']+1e-250), decimals=3)
        # self.merged_attrs['k']    = {'units': 'log10(m^5 molecule^-2)'}

        # self.merged_datasets['alpha'] = np.around(np.log10(self.merged_datasets['alpha']+1e-250), decimals=3)
        # self.merged_attrs['alpha']    = {'units': 'log10(m^-1 molecule^-2)'}

        # Flip arrays to be ascending in wavelength
        if np.diff(self.merged_datasets['wave'])[0] < 0:
            self.merged_datasets['wave']  = self.merged_datasets['wave'][::-1]
            self.merged_datasets['k']     = self.merged_datasets['k'][::-1]
            self.merged_datasets['alpha'] = self.merged_datasets['alpha'][::-1]

        print(self.merged_attrs)

        # Save the merged data
        utils.save_to_hdf5(
            self.final_output_file, 
            data=self.merged_datasets, 
            attrs=self.merged_attrs
            )

    def plot_merged_outputs(self, *args, **kwargs):

        # Load the merged data
        self.merged_datasets = utils.read_from_hdf5(
            self.final_output_file, keys_to_read=['wave', 'T', 'k', 'alpha']
            )

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(12,6), nrows=2, sharex=True)

        for i, T_i in enumerate(self.merged_datasets['T']):
            c_i = plt.get_cmap('RdBu_r')(T_i/self.merged_datasets['T'].max())

            ax[0].plot(1e6*self.merged_datasets['wave'], self.merged_datasets['k'][:,i], c=c_i, label=f'{T_i} K')
            ax[1].plot(1e6*self.merged_datasets['wave'], self.merged_datasets['alpha'][:,i], c=c_i)

        ax[0].set(xscale='log', yscale='log', ylabel='k [m^5 molecule^-2]')
        ax[1].set(xscale='log', yscale='log', xlabel='wave [m]', ylabel='alpha [m^-1 molecule^-2]')

        plt.savefig(self.output_data_dir / 'cia_coeffs.pdf', bbox_inches='tight')
        plt.close()
