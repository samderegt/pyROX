import numpy as np

import datetime

from pyROX import utils, sc, CrossSections

class CIA(CrossSections):
    """
    Base class for Collision-Induced Absorption (CIA) data handling.
    """
    def download_data(self, config):
        """
        Download CIA data from HITRAN or Borysow.

        Args:
            config (object): Configuration object containing URLs and input data directory.

        Raises:
            ValueError: If any file fails to download.
        """

        print('\nDownloading CIA data')

        files = []
        for url in config.urls:
            file = utils.download(url, config.input_data_dir)
            files.append(file)

        if None in files:
            raise ValueError('Failed to download all urls.')

    def __init__(self, config, **kwargs):
        """
        Initialises the CIA object.

        Args:
            config (object): Configuration object containing settings and file paths.
            **kwargs: Additional arguments for initialisation.
        """
        # Initialise the CrossSections parent class
        super().__init__(config, **kwargs)

    def calculate_temporary_outputs(self, overwrite=False, **kwargs):
        """
        Calculate the CIA coefficients and save temporary outputs.

        Args:
            overwrite (bool): Whether to overwrite existing temporary files.
            **kwargs: Additional arguments for calculation.
        """

        print('\nCalculating CIA coefficients')

        cia_files_masks = self.config.files.get('cia', None)
        if cia_files_masks is None:
            raise ValueError('No CIA files specified in the configuration.')

        cia_files, cia_masks = [], []
        for file in cia_files_masks:
            if isinstance(file, tuple):
                file, *masks = file
            else:
                masks = []
            cia_files.append(file)
            cia_masks.append(masks)

        # Check if the output files already exist
        tmp_output_files = self._check_existing_output_files(
            cia_files, overwrite_all=overwrite
            )

        self.T_grid, self.abs_coeff_k, self.abs_coeff_alpha = [], [], []
        for i, (file, masks) in enumerate(zip(cia_files, cia_masks)):

            # Compute absorption coefficients
            T_grid, abs_coeff_k, abs_coeff_alpha = self._read_absorption_coefficients(file)

            if len(masks) != 0:
                # Remove temperatures outside the range
                mask_T = masks[0](T_grid)
                T_grid, abs_coeff_k, abs_coeff_alpha = self._mask_arrays(
                    [T_grid, abs_coeff_k, abs_coeff_alpha], mask=mask_T, axis=0
                    )

                # Remove wavenumbers outside the range
                mask_nu = masks[1](self.nu_grid/(1e2*sc.c))
                abs_coeff_k[:,~mask_nu]     = 0.
                abs_coeff_alpha[:,~mask_nu] = 0.

            # Temporarily save the data
            utils.save_to_hdf5(
                tmp_output_files[i], 
                data={
                    'wave': self.wave_grid, 
                    'T': T_grid, 
                    'log10(k)': utils.log10_round(abs_coeff_k.T, decimals=3), 
                    'log10(alpha)': utils.log10_round(abs_coeff_alpha.T, decimals=3)
                    },
                attrs={
                    'wave': {'units': 'm'}, 
                    'T': {'units': 'K'}, 
                    'log10(k)': {'units': 'log10(m^5 molecule^-2)'}, 
                    'log10(alpha)': {'units': 'log10(m^-1 molecule^-2)'}
                    }
                )

    def save_combined_outputs(self, **kwargs):
        """
        Save the merged CIA outputs to a file.

        Args:
            **kwargs: Additional arguments for saving.
        """
        super().save_combined_outputs(keys_to_merge=['log10(k)','log10(alpha)'], **kwargs)

    def plot_combined_outputs(self, return_fig_ax=False, cmap='coolwarm', xscale='log', yscale='log', xlim=None, ylim=None, **kwargs):
        """
        Plot the merged CIA coefficients.

        Args:
            return_fig_ax (bool): Whether to return the figure and axes objects.
            cmap (str): Colormap for the plot.
            xscale (str): Scale for the x-axis ('linear' or 'log').
            yscale (str): Scale for the y-axis ('linear' or 'log').
            xlim (tuple, optional): Limits for the x-axis.
            ylim (tuple, optional): Limits for the y-axis.
            **kwargs: Additional arguments for plotting.
        """
        
        print(f'\nPlotting CIA coefficients')

        # Load the merged data
        self.combined_datasets = utils.read_from_hdf5(
            self.final_output_file, keys_to_read=['wave', 'T', 'log10(k)', 'log10(alpha)']
            )
        wave = self.combined_datasets['wave'] * 1e6 # [m] -> [um]
        k = 10**self.combined_datasets['log10(k)'] * (1e2)**5 # [m^5 molecule^-2] -> [cm^5 molecule^-2]
        alpha = 10**self.combined_datasets['log10(alpha)'] * (1e2)**(-2) # [m^-1 molecule^-2] -> [cm^-1 molecule^-2]

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9,6), nrows=2, sharex=True)

        # Plot for certain temperatures
        T_to_plot = kwargs.get('T_to_plot', self.combined_datasets['T'])
        indices_T, T_to_plot = utils.find_closest_indices(self.combined_datasets['T'], T_to_plot)

        indices_T = np.unique(indices_T)
        T_to_plot = np.unique(T_to_plot)

        for idx_T, T in zip(indices_T, T_to_plot):
            if len(T_to_plot)==1:
                c = plt.get_cmap(cmap)(0.4)
            else:
                c = plt.get_cmap(cmap)((T-T_to_plot.min())/(T_to_plot.max()-T_to_plot.min()))

            ax[0].plot(wave, k[:,idx_T], c=c, label=f'{T:.0f} K')
            ax[1].plot(wave, alpha[:,idx_T], c=c)

        ax[0].set(xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim, ylabel=r'k [cm$^5$ molecule$^{-2}$]')
        ax[1].set(xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim, xlabel=r'Wavelength [$\mu$m]', ylabel=r'alpha [cm$^{-1}$ molecule$^{-2}$]')

        handles, _ = ax[0].get_legend_handles_labels()
        ncols = 1 + len(handles)//8
        ax[0].legend(loc='upper right', ncol=ncols, labelcolor='linecolor')

        plt.savefig(self.output_data_dir / 'cia_coeffs.pdf', bbox_inches='tight')
        if return_fig_ax:
            return fig, ax
        plt.close()
    
    def convert_to_pRT3(self, contributor=None, **kwargs):
        """
        Convert the CIA data to petitRADTRANS v3.0 format.

        Args:
            contributor (str): Name of the contributor for these data.
            **kwargs: Additional arguments for conversion.

        Raises:
            ValueError: If required metadata is missing in the configuration.
            KeyError: If required keys are missing in the metadata.
        """

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
        if contributor is not None:
            attrs['DOI']['contributor'] = contributor

        # Load the merged data
        self.combined_datasets = utils.read_from_hdf5(
            self.final_output_file, keys_to_read=['wave', 'T', 'log10(alpha)']
            )

        wave_min = self.combined_datasets['wave'].min() / sc.micron
        wave_max = self.combined_datasets['wave'].max() / sc.micron

        resolution = self.resolution
        if not np.isnan(self.delta_nu):
            # Use the given delta_nu, resolution at 1 um
            resolution = sc.c / sc.micron / self.delta_nu
        elif not np.isnan(self.delta_wave):
            # Use the given delta_wave, resolution at 1 um
            resolution = sc.micron / self.delta_wave
        if np.isnan(resolution):
            idx = np.searchsorted(self.combined_datasets['wave'], 1e-6)
            delta_wave = np.diff(self.combined_datasets['wave'][idx-1:idx+1])
            resolution = self.combined_datasets['wave'][idx] / delta_wave

        # Fill the dictionary
        data['alpha'] = 1e-2 * 10**self.combined_datasets['log10(alpha)'][::-1,:].T # [m^-1 molecule^-2] -> [cm^-1 molecule^-2]
        data['t'] = self.combined_datasets['T'] # [K]
        data['wavenumbers'] = 1e-2/self.combined_datasets['wave'][::-1] # [m] -> [cm^-1], ascending in wavenumber
        data['wlrange'] = [wave_min, wave_max] # [um]
        data['wnrange'] = [1e4/wave_max, 1e4/wave_min] # [cm^-1]

        # Complete the filename
        pRT_file = '{}--{}-NatAbund__{}.R{:.0f}_{:.1f}-{:.0f}mu.ciatable.petitRADTRANS.h5'
        pRT_file = pRT_file.format(
            pRT3_metadata['mol_name'][0], pRT3_metadata['mol_name'][1], 
            self.database.replace('cia_','').upper(), resolution, 
            wave_min, wave_max
            )
        pRT_file = self.output_data_dir / pRT_file

        # Save the datasets
        utils.save_to_hdf5(pRT_file, data=data, attrs=attrs, compression=None)
        print(f'  Saved to \"{pRT_file}\"')