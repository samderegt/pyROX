import numpy as np
import pathlib
import itertools

from tqdm import tqdm

from pyROX import utils, sc

def load_data_object(config, **kwargs):
    """
    Load the data object based on the given configuration.

    Args:
        config (object): Configuration object containing database information.
        **kwargs: Additional arguments to pass to the data object.

    Returns:
        object: An instance of the appropriate data object based on the database type.

    Raises:
        NotImplementedError: If the specified database is not implemented.
    """
    
    database = getattr(config, 'database', '').lower()
    
    from pyROX import CIA_HITRAN, CIA_Borysow
    if database == 'cia_hitran':
        return CIA_HITRAN(config, **kwargs)
    elif database == 'cia_borysow':
        return CIA_Borysow(config, **kwargs)

    from pyROX import LBL_HITRAN, LBL_ExoMol, LBL_Kurucz
    if database == 'exomol':
        return LBL_ExoMol(config, **kwargs)
    elif database in ['hitemp', 'hitran']:
        return LBL_HITRAN(config, **kwargs)
    elif database in ['kurucz', 'vald']:
        return LBL_Kurucz(config, **kwargs)

    raise NotImplementedError(f"Database '{config.database}' not implemented.")


class CrossSections:
    """
    Base class for handling cross-sections.
    """

    @staticmethod
    def _mask_arrays(arrays, mask, **kwargs):
        """
        Apply a mask to the given arrays.

        Args:
            arrays (list): List of arrays to mask.
            mask (array-like): Boolean mask to apply.
            **kwargs: Additional arguments for numpy.compress.

        Returns:
            list: List of masked arrays.
        """
        return [np.compress(mask, array, **kwargs) for array in arrays]

    def __init__(self, config, download=False):
        """
        Initialises the CrossSections object.

        Args:
            config (object): Configuration object containing parameters.
            download (bool, optional): Whether to download data during initialisation. Defaults to False.
        """
        self.database = config.database.lower()
        if download:
            # First download the data
            self.download_data(config)
            
        # Read the common variables
        self._read_configuration_parameters(config)

    def download_data(self, *args, **kwargs):
        """
        Download data required for cross-section calculations.

        Raises:
            NotImplementedError: This method should be implemented in a subclass.
        """
        raise NotImplementedError("This method should be implemented in the subclass.")

    def calculate_temporary_outputs(self, *args, **kwargs):
        """
        Calculate temporary outputs for cross-sections.

        Raises:
            NotImplementedError: This method should be implemented in a subclass.
        """
        raise NotImplementedError("This method should be implemented in the subclass.")

    def plot_combined_outputs(self, *args, **kwargs):
        """
        Plot merged outputs for cross-sections.

        Raises:
            NotImplementedError: This method should be implemented in a subclass.
        """
        raise NotImplementedError("This method should be implemented in the subclass.")

    def save_combined_outputs(self, keys_to_merge, progress_bar=False, overwrite=False, **kwargs):
        """
        Combine temporary files and save the final output.

        Args:
            keys_to_merge (list): List of keys to combine.
            progress_bar (bool, optional): Whether to show a progress bar. Defaults to False.
            overwrite (bool, optional): Whether to overwrite the final output file if it exists. Defaults to False.
            **kwargs: Additional arguments for merging.
        """
        print('\nCombining temporary files and saving final output')

        # Check if the output file already exist
        self._check_existing_output_files(
            final_output_file=self.final_output_file, 
            overwrite_all=overwrite
            )

        # Combine the temporary files
        self.combine_temporary_outputs(
            keys_to_merge=keys_to_merge, 
            progress_bar=progress_bar,
            )

        # Flip arrays to be ascending in wavelength
        if np.diff(self.combined_datasets['wave'])[0] < 0:
            self.combined_datasets['wave']  = self.combined_datasets['wave'][::-1]
            for key in keys_to_merge:
                self.combined_datasets[key] = self.combined_datasets[key][::-1]

        # Save the merged data
        print(f'  Saving final output to \"{self.final_output_file}\"')
        utils.save_to_hdf5(
            self.final_output_file, 
            data=self.combined_datasets, 
            attrs=self.combined_attrs
            )

    def combine_temporary_outputs(self, keys_to_merge, progress_bar=False, **kwargs):
        """
        Combine temporary output-files into a single dataset.

        Args:
            keys_to_merge (list): List of keys to read from the temporary files. 
            progress_bar (bool, optional): Whether to show a progress bar. Defaults to False.
            **kwargs: Additional arguments for merging.

        Raises:
            FileNotFoundError: If no temporary cross-section files are found.
            ValueError: If the temperature grid is not found in a temporary file.
        """

        # Check if the outputs should be summed
        sum_outputs = self._check_if_sum_outputs()
        
        # Select all temporary files
        tmp_files = list(self.tmp_output_dir.glob('*.hdf5'))
        if len(tmp_files) == 0:
            raise FileNotFoundError(f'No temporary output files found in \"{self.tmp_output_dir}\".')
        # Sort by modification date
        tmp_files.sort(key=lambda x: x.stat().st_mtime)

        if sum_outputs:
            print(f'  Summing {len(tmp_files)} temporary files')
        else:
            print(f'  Merging {len(tmp_files)} temporary files into a single grid')
        print(f'  Temporary files:')
        for tmp_file in tmp_files:
            print(f'    - \"{tmp_file.name}\"')
        
        # Check compatibility of files before combining (wave_main is sorted)
        wave_main, P_main, T_main, flip_wave_axis = self._combine_PT_grids(tmp_files, sum_outputs)

        # Combine all files into a single array
        self.combined_datasets = {
            key: np.zeros((len(wave_main), len(P_main), len(T_main)), dtype=np.float64) 
            for key in keys_to_merge
            }
        # Make a nice progress bar
        pbar_kwargs = dict(
            disable=(not progress_bar), 
            bar_format='{l_bar}{bar:20}{r_bar}{bar:-20b}', 
        )
        for i, tmp_file in enumerate(tqdm(tmp_files, **pbar_kwargs)):

            datasets, attrs = utils.read_from_hdf5(
                tmp_file, keys_to_read=keys_to_merge+['P','T'], return_attrs=True
                )

            # Check if dataset has pressure-axis (line-by-line vs. absorption coeffs)
            has_pressure_axis = 'P' in datasets.keys()
            P = datasets.get('P', [0.])
            T = datasets.get('T')
            if T is None:
                raise ValueError("Temperature grid not found in temporary file.")
            
            # Iterate over PT-grid of the temporary file
            iterables = zip(
                itertools.product(P, T), itertools.product(range(len(P)),range(len(T)))
                )
            for PT, (j,k) in iterables:
                idx_P = np.argwhere(
                    np.isclose(P_main, PT[0], rtol=1e-5, atol=0.)
                    ).flatten()[0]
                idx_T = np.argwhere(
                    np.isclose(T_main, PT[1], rtol=1e-5, atol=0.)
                    ).flatten()[0]
                
                for key in keys_to_merge:

                    # Current data
                    combined_data = self.combined_datasets[key][:,idx_P,idx_T]
                    
                    if not has_pressure_axis:
                        # Absorption coefficients
                        data_to_add = datasets[key][:,k]
                    else:
                        # Line-by-line cross-sections
                        data_to_add = datasets[key][:,j,k]

                    if flip_wave_axis[i]:
                        # Flip the wavelength axis to be ascending
                        data_to_add = data_to_add[::-1]
                    
                    if sum_outputs:
                        # Sum cross-sections
                        if key.startswith('log10('):
                            # Convert to linear scale
                            data_to_add = 10**data_to_add

                            if np.all(combined_data==0.):
                                combined_data = combined_data # First file
                            else:
                                # Convert to linear scale
                                combined_data = 10**combined_data
                        
                        # Add to existing data
                        combined_data += data_to_add

                        if key.startswith('log10('):
                            # Convert back to log scale
                            combined_data = np.log10(combined_data)

                    else:
                        # Avoid summing when only one transition file was used
                        combined_data = data_to_add

                    # Store the combined data
                    self.combined_datasets[key][:,idx_P,idx_T] = combined_data

        # Add the grid-definition
        self.combined_datasets['wave'] = wave_main
        self.combined_datasets['T'] = T_main
        if has_pressure_axis:
            self.combined_datasets['P'] = P_main
        else:
            for key in keys_to_merge:
                # Remove the pressure axis
                self.combined_datasets[key] = np.squeeze(self.combined_datasets[key], axis=1)

        self.combined_attrs = attrs

    def _combine_PT_grids(self, tmp_files, sum_outputs):
        """
        Combine the PT-grids of the temporary files.

        Args:
            tmp_files (list): List of temporary files to combine.
            sum_outputs (bool): Whether to sum the outputs.

        Returns:
            tuple: A tuple containing the main wavelength grid, pressure grid, temperature grid, 
            and a list of booleans indicating if the wavelength axis was flipped.

        Raises:
            ValueError: If the PT grid is not rectangular or grids are incompatible.
        """
        print('  Combining PT-grids of temporary files')

        all_PT = []
        flip_wave_axis = []
        for i, tmp_file in enumerate(tmp_files):

            datasets = utils.read_from_hdf5(tmp_file, keys_to_read=['wave','P','T'])

            wave = datasets.get('wave')
            if wave[1] < wave[0]:
                # Flip the wavelength axis to be ascending
                wave = wave[::-1]
                flip_wave_axis.append(True)
            else:
                flip_wave_axis.append(False)
            P = datasets.get('P', [0.])
            T = datasets.get('T')

            if wave is None:
                raise ValueError("Wavelength grid not found in temporary file.")
            if T is None:
                raise ValueError("Temperature grid not found in temporary file.")

            if i == 0:
                # Initialise the main arrays
                wave_main = wave.copy()
                P_main = P.copy()
                T_main = T.copy()

            assert np.isclose(wave_main, wave, rtol=1e-5, atol=0.).all(), "Wavelength grids are not compatible."
            if sum_outputs:
                assert np.isclose(P_main, P, rtol=1e-5, atol=0.).all(), "Pressure grids are not compatible."
                assert np.isclose(T_main, T, rtol=1e-5, atol=0.).all(), "Temperature grids are not compatible."

                # Make an exact copy to avoid extending the arrays
                P = P_main.copy()
                T = T_main.copy()

            # Add the PT combination
            for PT in itertools.product(P, T):
                all_PT.append(PT)

        all_PT = np.array(all_PT)

        # Check if the PT grid is rectangular
        P_main = np.unique(all_PT[:,0])
        for i, P in enumerate(P_main):
            T = np.unique(all_PT[all_PT[:,0]==P,1])
            if i == 0:
                T_main = T.copy()

            assert(T_main == T).all(), "PT-grid is not rectangular."

        P_main = np.sort(P_main)
        T_main = np.sort(T_main)

        return wave_main, P_main, T_main, flip_wave_axis
        
    def _check_if_sum_outputs(self):
        """
        Check if the cross-sections should be summed.

        Returns:
            bool: True if the outputs should be summed, False otherwise.
        """
        transitions_files = self.config.files.get('transitions', [])
        if not isinstance(transitions_files, (list, tuple, np.ndarray)):
            return False
        if len(transitions_files) == 0:
            return False
        return True

    def _check_existing_output_files(self, input_files=[''], overwrite_all=False, final_output_file=None):
        """
        Check if the output files already exist.

        Args:
            input_files (list): List of input file paths.
            overwrite_all (bool): Whether to overwrite all existing files.
            final_output_file (str): Expected output file name.

        Returns:
            list: List of output file paths.
        """
        output_files = []
        for i, input_file in enumerate(input_files):
            # Check if the transition file exists
            input_file = pathlib.Path(input_file)

            if final_output_file is None:
                # Is a temporary output file, create the filename based on the input
                output_file = pathlib.Path(
                    str(self.tmp_output_basename).replace('.hdf5', f'_{input_file.stem}.hdf5')
                )
            else:
                # Use the given, final output filename
                output_file = pathlib.Path(final_output_file)

            # Check if the output file already exists
            if output_file.exists() and not overwrite_all:
                # Ask the user if they want to overwrite the file
                response = ''
                while response not in ['y', 'yes', 'n', 'no', 'all']:
                    response = input(f'  Warning: Output file \"{output_file}\" already exists. Overwrite? (yes/no/all): ')
                    response = response.strip().lower()
                    if response in ['no', 'n']:
                        raise FileExistsError(f'Not overwriting existing file: \"{output_file}\".')
                    elif response in ['yes', 'y']:
                        break
                    elif response == 'all':
                        overwrite_all = True
                        break
                    else:
                        print('  Invalid input. Please enter \"yes\", \"no\", or \"all\".')
                        continue
            output_files.append(output_file)

        input_files = [pathlib.Path(input_file) for input_file in input_files]
        return input_files, output_files

    def _read_configuration_parameters(self, config):
        """
        Read parameters from the configuration object.

        Args:
            config (object): Configuration object containing parameters.

        Raises:
            ValueError: If no input-data files are specified in the configuration.
        """
        print('\nReading parameters from the configuration file')
        self.config = config
        utils.warn_about_units(self.config)
        
        self.database = getattr(self.config, 'database', None)
        self.species  = getattr(self.config, 'species', None)

        # Wavenumber/wavelength grid
        self._setup_nu_grid()
        
        # Input/output directories
        self.input_data_dir = getattr(self.config, 'input_data_dir', f'./{self.species}/input_data/')
        self.input_data_dir = pathlib.Path(self.input_data_dir).resolve()
        self.input_data_dir.mkdir(parents=True, exist_ok=True)

        self.output_data_dir = getattr(self.config, 'output_data_dir', f'./{self.species}/')
        self.output_data_dir = pathlib.Path(self.output_data_dir).resolve()
        self.output_data_dir.mkdir(parents=True, exist_ok=True)

        # Files to read
        files = getattr(self.config, 'files', {})
        if len(files) == 0:
            raise ValueError('No input-data files specified in the configuration.')

        # Temporary output files
        self.tmp_output_dir = self.output_data_dir / 'tmp'

        default = 'xsec.hdf5'
        if self.database.startswith('cia'):
            default = 'cia_coeffs.hdf5'
        self.tmp_output_basename = getattr(self.config, 'tmp_output_basename', default)
        self.tmp_output_basename = (self.tmp_output_dir / self.tmp_output_basename).resolve()
        # Ensure the prefix ends with '.hdf5'
        if not self.tmp_output_basename.suffix == '.hdf5':
            self.tmp_output_basename = self.tmp_output_basename.with_suffix('.hdf5')

        # Final output file
        default = default.format(self.species)
        self.final_output_file = getattr(self.config, 'final_output_file', default)
        self.final_output_file = (self.output_data_dir / self.final_output_file).resolve()

    def _setup_nu_grid(self):
        """
        Configure the wavenumber grid based on the configuration parameters.
        """
        # Read the necessary parameters
        self.wave_file = getattr(self.config, 'wave_file', '')
        self.wave_file = pathlib.Path(self.wave_file).resolve()

        self.resolution = getattr(self.config, 'resolution', np.nan)

        self.delta_nu   = getattr(self.config, 'delta_nu', np.nan) # [cm^-1]
        self.delta_nu   = self.delta_nu * 1e2*sc.c # [s^-1]
        self.delta_wave = getattr(self.config, 'self.delta_wave', np.nan) # [um]
        self.delta_wave = self.delta_wave * sc.micron # [um] -> [m]
        
        self.wave_min = getattr(self.config, 'wave_min', 1.0/3.0) # [um]
        self.wave_min = self.wave_min * sc.micron # [um] -> [m]
        self.wave_max = getattr(self.config, 'wave_max', 250.0) # [um]
        self.wave_max = self.wave_max * sc.micron # [um] -> [m]

        self.nu_min = getattr(self.config, 'nu_min', 1e-2/self.wave_max) * (1e2*sc.c) # [m] -> [s^-1]
        self.nu_max = getattr(self.config, 'nu_max', 1e-2/self.wave_min) * (1e2*sc.c) # [m] -> [s^-1]

        self.adaptive_nu_grid = False
        if self.wave_file.is_file():
            # Use a custom wavelength grid from a file
            self.wave_file = pathlib.Path(self.wave_file).resolve()
            self.wave_grid = np.loadtxt(self.wave_file)

        elif not np.isnan(self.resolution):
            # Use a fixed resolution
            self.wave_grid = utils.fixed_resolution_wavelength_grid(
                self.wave_min, self.wave_max, resolution=self.resolution
                ) # [m]

        elif not np.isnan(self.delta_nu):
            # Use an equal wavenumber-spacing
            self.N_nu = int((self.nu_max-self.nu_min)/self.delta_nu) + 1
            self.delta_nu = (self.nu_max-self.nu_min) / (self.N_nu-1)
            
            self.nu_grid   = np.linspace(self.nu_min, self.nu_max, num=self.N_nu, endpoint=True)
            self.wave_grid = sc.c/self.nu_grid # [s^-1] -> [m]

            # Decrease the wavenumber-grid resolution for high pressures
            self.adaptive_nu_grid = getattr(self.config, 'adaptive_nu_grid', False)

        elif not np.isnan(self.delta_wave):
            # Use an equal wavelength-spacing
            self.N_nu  = int((self.wave_max-self.wave_min)/self.delta_wave) + 1
            self.delta_wave = (self.wave_max-self.wave_min) / (self.N_nu-1)
            
            self.wave_grid = np.linspace(self.wave_min, self.wave_max, num=self.N_nu, endpoint=True)

        else:
            raise ValueError('No wavelength/wavenumber grid can be made with the configuration parameters.')

        self.nu_grid = sc.c/self.wave_grid # [m] -> [s^-1]

        self.wave_min = np.min(self.wave_grid)
        self.wave_max = np.max(self.wave_grid)
        self.nu_min   = np.min(self.nu_grid)
        self.nu_max   = np.max(self.nu_grid)
        self.N_nu = len(self.wave_grid)

        # Sort by wavenumber
        idx = np.argsort(self.nu_grid) 
        self.wave_grid = self.wave_grid[idx]
        self.nu_grid   = self.nu_grid[idx]

        print('\nWavelength-grid:')
        print(f'  Wavelength: {self.wave_min/sc.micron:.2f} - {self.wave_max/sc.micron:.0f} um')
        print(f'  Wavenumber: {self.nu_min/(1e2*sc.c):.0f} - {self.nu_max/(1e2*sc.c):.0f} cm^-1')
        if not np.isnan(self.resolution):
            print(f'  Fixed resolution: {self.resolution:.0e}')
        elif not np.isnan(self.delta_nu):
            print(f'  Fixed wavenumber-spacing: {self.delta_nu/(1e2*sc.c):.3f} cm^-1')
        elif not np.isnan(self.delta_wave):
            print(f'  Fixed wavelength-spacing: {self.delta_wave/sc.micron:.3f} um')
        print(f'  Number of grid points: {self.N_nu}')

        print(f'  Adaptive grid: {self.adaptive_nu_grid}')

    def _setup_coarse_nu_grid(self, adaptive_delta_nu):
        """
        Configure a coarse wavenumber grid for adaptive calculations.

        Args:
            adaptive_delta_nu (float): Desired resolution for the coarse grid.

        Returns:
            numpy.ndarray: Coarse wavenumber grid.
        """
        # Use the original grid
        if not self.adaptive_nu_grid:
            return self.nu_grid
        if self.delta_nu > adaptive_delta_nu:
            return self.nu_grid
        
        # Decrease number of points in wavenumber grid
        inflation_factor = adaptive_delta_nu / self.delta_nu
        coarse_N_nu = int(self.N_nu/inflation_factor)

        # Expand wavenumber grid slightly
        delta_nu_ends = adaptive_delta_nu * (inflation_factor-1)/2
        coarse_nu_grid = np.linspace(
            self.nu_min-delta_nu_ends, self.nu_max+delta_nu_ends, 
            num=coarse_N_nu, endpoint=True
        )
        coarse_delta_nu = coarse_nu_grid[1] - coarse_nu_grid[0]

        return coarse_nu_grid