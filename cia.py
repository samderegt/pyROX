import numpy as np
from scipy.interpolate import interp1d

import datetime
import pathlib

from cross_sections import CrossSections
from utils import sc
import utils

class CIA(CrossSections):
    """
    Base class for CIA cross-sections.
    """

    def __init__(self, config):
        # Initialise the CrossSections parent class
        super().__init__(config)

    def calculate_tmp_outputs(self, **kwargs):
        """
        Calculate the CIA coefficients.
        """

        print('\nCalculating CIA coefficients')

        files = getattr(self.config, 'files', None)
        if files is None:
            raise ValueError('No files specified in the configuration.')
        cia_files = files.get('cia', None)
        if cia_files is None:
            raise ValueError('No CIA files specified in the configuration.')

        self.T_grid, self.abs_coeff_k, self.abs_coeff_alpha = [], [], []
        for i, (file, *masks) in enumerate(cia_files):

            # Compute absorption coefficients
            T_grid, abs_coeff_k, abs_coeff_alpha = self._read_absorption_coefficients(file)

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

    def save_merged_outputs(self, **kwargs):
        """
        Merge the temporary files and save the final output. Same for all CIA classes.
        """

        print('\nMerging temporary files and saving final output')

        # Ask to overwrite the final output file if it already exists
        response = ''
        if not self.final_output_file.exists():
            response = 'yes'

        while response == '':
            response = input(f'  Warning: Final output file \"{self.final_output_file}\" already exists. Do you want to overwrite it? (yes/no): ')
            if response == '':
                continue
            elif response.lower() not in ['y', 'yes']:
                raise FileExistsError(f"Not overwriting final output file '{self.final_output_file}'.")

        # Merge the temporary files
        self.merge_tmp_outputs(keys_to_merge=['k', 'alpha'])

        # Flip arrays to be ascending in wavelength
        if np.diff(self.merged_datasets['wave'])[0] < 0:
            self.merged_datasets['wave']  = self.merged_datasets['wave'][::-1]
            self.merged_datasets['k']     = self.merged_datasets['k'][::-1]
            self.merged_datasets['alpha'] = self.merged_datasets['alpha'][::-1]

        # Save the merged data
        print(f'  Saving final output to \"{self.final_output_file}\"')
        utils.save_to_hdf5(
            self.final_output_file, 
            data=self.merged_datasets, 
            attrs=self.merged_attrs
            )

    def plot_merged_outputs(self, cmap='coolwarm', xscale='log', yscale='log', xlim=None, ylim=None, **kwargs):
        """
        Plot the merged outputs. Same for all CIA classes.
        """
        
        print(f'\nPlotting CIA coefficients')

        # Load the merged data
        self.merged_datasets = utils.read_from_hdf5(
            self.final_output_file, keys_to_read=['wave', 'T', 'k', 'alpha']
            )

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9,6), nrows=2, sharex=True)

        # Plot for certain temperatures
        T_to_plot = kwargs.get('T_to_plot', self.merged_datasets['T'])
        indices_T, T_to_plot = utils.find_nearest(self.merged_datasets['T'], T_to_plot)

        indices_T = np.unique(indices_T)
        T_to_plot = np.unique(T_to_plot)

        for idx_T, T in zip(indices_T, T_to_plot):

            c = plt.get_cmap(cmap)((T-T_to_plot.min())/(T_to_plot.max()-T_to_plot.min()))

            ax[0].plot(self.merged_datasets['wave'], self.merged_datasets['k'][:,idx_T], c=c, label=f'{T:.0f} K')
            ax[1].plot(self.merged_datasets['wave'], self.merged_datasets['alpha'][:,idx_T], c=c)

        ax[0].set(xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim, ylabel='k [m^5 molecule^-2]')
        ax[1].set(xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim, xlabel='wave [m]', ylabel='alpha [m^-1 molecule^-2]')

        handles, _ = ax[0].get_legend_handles_labels()
        ncols = 1 + len(handles)//8
        ax[0].legend(loc='upper right', ncol=ncols, labelcolor='linecolor')

        plt.savefig(self.output_data_dir / 'cia_coeffs.pdf', bbox_inches='tight')
        plt.close()

    def convert_to_pRT2(self):
        raise NotImplementedError('Conversion to petitRADTRANS-v2.0 format is not implemented.')
    
    def convert_to_pRT3(self, **kwargs):

        print(f'\nConverting to petitRADTRANS-v3.0 format')

        # Load the attributes for the hdf5 file
        pRT3_metadata = getattr(self.config, 'pRT3_metadata', None)
        if pRT3_metadata is None:
            raise ValueError('No pRT3_metadata specified in the configuration.')

        # Check if required keys are in pRT3_metadata
        for key in ['DOI', 'mol_mass', 'mol_name']:
            if key in pRT3_metadata:
                continue
            raise KeyError(f"Required key '{key}' not found in pRT3_metadata.")


        data = {
            'DOI': np.atleast_1d(pRT3_metadata['DOI']),
            'Date_ID': np.atleast_1d(f'petitRADTRANS-v3_{datetime.datetime.now(datetime.timezone.utc).isoformat()}'),
            'mol_mass': np.atleast_1d(pRT3_metadata['mol_mass']),
            'mol_name': np.atleast_1d(pRT3_metadata['mol_name'])
        }

        # Attributes of all datasets
        attrs = dict(
            DOI = {'additional_description': 'None', 'long_name': 'Data object identifier linked to the data'},
            Date_ID = {'long_name': 'ISO 8601 UTC time (https://docs.python.org/3/library/datetime.html) at which the table has been created, along with the version of petitRADTRANS'},
            alpha = {'long_name': 'Table of monochromatic absorption with axes (temperature, wavenumber)', 'units': 'cm^-1'},
            mol_mass = {'long_name': 'Masses of the colliding species', 'units': 'AMU'},
            mol_name = {'long_name': 'Names of the colliding species described'},
            t = {'long_name': 'Temperature grid', 'units': 'K'},
            wavenumbers = {'long_name': 'CIA wavenumbers', 'units': 'cm^-1'},
            wlrange = {'long_name': 'Wavelength range covered', 'units': 'µm'},
            wnrange = {'long_name': 'Wavenumber range covered', 'units': 'cm^-1'}
        )
        # Add contributors if given
        contributor = kwargs.get('contributor', None)
        if contributor is not None:
            attrs['DOI']['contributor'] = contributor

        # Load the merged data
        self.merged_datasets = utils.read_from_hdf5(
            self.final_output_file, keys_to_read=['wave', 'T', 'k', 'alpha']
            )

        wave_min = 1e6*self.merged_datasets['wave'].min()
        wave_max = 1e6*self.merged_datasets['wave'].max()
        resolution = 1e6*sc.c/self.delta_nu # At 1 um

        # Fill the dictionary
        data['alpha'] = 1e-2*self.merged_datasets['alpha'][::-1,:].T # [m^-1 molecule^-2] -> [cm^-1 molecule^-2], ascending in wavenumber
        data['t'] = self.merged_datasets['T'] # [K]
        data['wavenumbers'] = 1e-2/self.merged_datasets['wave'][::-1] # [m] -> [cm^-1], ascending in wavenumber
        data['wlrange'] = [wave_min, wave_max] # [um]
        data['wnrange'] = [1e4/wave_max, 1e4/wave_min] # [cm^-1]

        # Complete the filename
        if isinstance(self, CIA_HITRAN):
            database = 'HITRAN'
        elif isinstance(self, CIA_Borysow):
            database = 'Borysow'

        pRT_file = '{}--{}-NatAbund__{}.R{:.0f}_{:.1f}-{:.0f}mu.ciatable.petitRADTRANS.h5'
        pRT_file = pRT_file.format(
            pRT3_metadata['mol_name'][0], pRT3_metadata['mol_name'][1], 
            database, resolution, wave_min, wave_max
            )
        pRT_file = self.output_data_dir / pRT_file

        # Save the datasets
        utils.save_to_hdf5(pRT_file, data=data, attrs=attrs, compression=None)

class CIA_HITRAN(CIA):

    def download_data(self):
        """
        Download CIA data from HITRAN.
        """

        print('\nDownloading CIA data from HITRAN')

        files = []
        for url in self.config.urls:
            file = utils.wget_if_not_exist(url, self.config.input_dir)
            files.append(file)
        
    def __init__(self, config):

        print('-'*60)
        print('  Collision-Induced Absorption from HITRAN')
        print('-'*60+'\n')
        
        super().__init__(config) # Initialise the parent CIA class

    def _read_absorption_coefficients(self, file):
        """
        Read absorption coefficients from a HITRAN CIA file.
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

class CIA_Borysow(CIA):
    def __init__(self, config):
        
        print('-'*60)
        print('  Collision-Induced Absorption from Borysow')
        print('-'*60+'\n')
        super().__init__(config) # Initialise the parent CIA class

    def download_data(self):
        raise NotImplementedError('Please download the data manually from https://www.astro.ku.dk/~aborysow/programs/')

    def _read_absorption_coefficients(self, file):
        """
        Read absorption coefficients from a Borysow CIA file.
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

        abs_coeff_alpha = np.zeros((len(T_grid), self.N_grid))
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
