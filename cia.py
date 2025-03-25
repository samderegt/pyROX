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

    def calculate_cross_sections(self, **kwargs):

        self.T_grid, self.abs_coeff_k, self.abs_coeff_alpha = [], [], []
        for i, (file, *masks) in enumerate(self.config.files['cia']):
            
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
                mask_nu = masks[1](self.nu_grid)
                abs_coeff_k[:,~mask_nu]     = 0.
                abs_coeff_alpha[:,~mask_nu] = 0.

            # Temporarily save the data
            tmp_output_file = str(self.tmp_output_file)
            tmp_output_file = tmp_output_file.format(pathlib.Path(file).name)

            utils.save_to_hdf5(
                tmp_output_file, 
                data={
                    'wave': sc.c/self.nu_grid, 
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

        # Combine the temporary files and save the data
        datasets, attrs = self.combine_tmp_cross_sections(keys_to_read=['k', 'alpha'])
        
        datasets['k']     = np.around(np.log10(datasets['k']+1e-250), decimals=3)
        datasets['alpha'] = np.around(np.log10(datasets['alpha']+1e-250), decimals=3)

        utils.save_to_hdf5(
            self.output_file, 
            data=datasets, 
            attrs=attrs
            )

        # Save the combined data
        print(datasets.keys(), attrs.keys())

"""
    def get_cross_sections(self, config, tmp_output_file='cia{}.hdf5', show_pbar=True, **kwargs):

        # Load CIA data from HITRAN
        file = conf.files['cia']
        with open(file, 'r') as f:
            lines = f.readlines()

        N_grid, self.T_grid = [], []
        for line in lines:
            if len(line) < 101:
                continue
            # Header line
            N_grid.append(int(line[40:40+7]))
            self.T_grid.append(float(line[47:47+7]))

        self.T_grid = np.array(self.T_grid)
        N_grid = np.array(N_grid) # TODO: use different name

        self.k_grid = np.zeros((len(self.T_grid), self.N_grid))
        for i, (N_grid_i, T_i) in enumerate(zip(N_grid, self.T_grid)):
            idx_min = 1 + i*(N_grid_i+1)
            idx_max = idx_min + N_grid_i

            nu = np.array(
                [line[0:10].strip() for line in lines[idx_min:idx_max]], dtype=np.float64
                )
            nu *= 100.0*sc.c # [cm^-1] -> [s^-1]

            k = np.array(
                [line[10:21].strip() for line in lines[idx_min:idx_max]], dtype=np.float64
                )
            k *= (1e-2)**5 # [cm^5 molecule^-2] -> [m^5 molecule^-2]

            # Interpolate onto nu_grid
            interp_func = interp1d(nu, np.log10(k), kind='linear', fill_value=np.nan, bounds_error=False)
            self.k_grid[i] = 10**interp_func(self.nu_grid)
            self.k_grid[i] = np.nan_to_num(self.k_grid[i], nan=0.0)

        # [m^5 molecule^-2] * [m^-6] = [m^-1 molecule^-2]
        self.alpha_grid = self.k_grid * sc.L0**2

        # Save cross-sections to file
        tmp_output_file = f'{conf.output_dir}/' + tmp_output_file.format('')
        self.save_cross_sections(tmp_output_file)

    def save_cross_sections(self, file):
        
        print(f'\nSaving cross-sections to file \"{file}\"')

        # Create directory if not exist
        pathlib.Path(file).parent.mkdir(parents=True, exist_ok=True)

        with h5py.File(file, 'w') as f:
            # Flip arrays to be ascending in wavelength
            wave       = sc.c / self.nu_grid[::-1]
            alpha_grid = self.alpha_grid[:,::-1]
            k_grid     = self.k_grid[:,::-1]

            f.create_dataset('T', compression='gzip', data=self.T_grid) # [K]
            f.create_dataset('wave', compression='gzip', data=wave) # [m]
            f.create_dataset(
                'alpha', compression='gzip', 
                data=np.log10(alpha_grid+1e-250) # [m^-1 molecule^-2] -> log10([m^-1 molecule^-2])
                ) 
            f.create_dataset(
                'k', compression='gzip', 
                data=np.log10(k_grid+1e-250) # [m^5 molecule^-2] -> log10([m^5 molecule^-2])
                )

        return

'''
        plt.figure(figsize=(8,6))
        for i, T_i in enumerate(T):
            if T_i%100 != 0:
                continue
            if T_i > 3000:
                continue
            plt.plot(1e6*sc.c/self.nu_grid, self.alpha_grid[i], c=plt.get_cmap('RdBu_r')(T_i/3000), label=str(T_i))
        plt.xscale('log'); plt.yscale('log')
        plt.xlim(0.5,250); plt.ylim(1e-14,1e-4)
        plt.legend(ncols=2)
        plt.savefig(conf.output_dir + 'cia.pdf')
        plt.close()

        import h5py
        hdf5_file = '/net/schenk/data2/regt/pRT3_input_data/input_data/opacities/continuum/collision_induced_absorptions/H2--H2/H2--H2-NatAbund/H2--H2-NatAbund__BoRi.R831_0.6-250mu.ciatable.petitRADTRANS.h5'
        with h5py.File(hdf5_file, 'r') as f:
            nu = f['wavenumbers'][:]
            alpha = f['alpha'][:]
            temperature = f['t'][:]

        plt.figure(figsize=(8,6))
        for i, T_i in enumerate(temperature):
            if T_i > 3000:
                continue
            plt.plot(1e4/nu, alpha[i], c=plt.get_cmap('RdBu_r')(T_i/3000), label=str(T_i))
        plt.xscale('log'); plt.yscale('log')
        plt.xlim(0.5,250); plt.ylim(1e-14,1e-4)
        plt.legend(ncols=2)
        plt.savefig(conf.output_dir + 'cia_petitRADTRANS.pdf')
        plt.close()
'''
"""