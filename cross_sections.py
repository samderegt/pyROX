import numpy as np
import pathlib
import itertools

from utils import sc
import utils

def load_data_object(config, **kwargs):
    """
    Load the data object from the given configuration.
    """
    
    database = getattr(config, 'database', '').lower()

    import cia, line_by_line
    if database == 'cia_hitran':
        return cia.CIA_HITRAN(config, **kwargs)
    elif database == 'cia_borysow':
        return cia.CIA_Borysow(config, **kwargs)

    if database == 'exomol':
        return line_by_line.ExoMol(config, **kwargs)
    elif database in ['hitemp', 'hitran']:
        return line_by_line.HITRAN(config, **kwargs)
    elif database in ['kurucz', 'vald']:
        return line_by_line.Kurucz(config, **kwargs)

    raise NotImplementedError(f"Database '{config.database}' not implemented.")


class CrossSections:
    """
    Base class for cross-sections.
    """

    @staticmethod
    def mask_arrays(arrays, mask, **kwargs):
        """Apply a mask to the given arrays."""
        return [np.compress(mask, array, **kwargs) for array in arrays]

    def __init__(self, config, download=False):

        self.database = config.database.lower()
        if download:
            # First download the data
            self.download_data(config)
            
        # Read the common variables
        self._read_from_config(config)

        # Get the wavenumber grid
        self._configure_nu_grid()

    def download_data(self, *args, **kwargs):
        raise NotImplementedError("This method should be implemented in the subclass.")

    def calculate_tmp_outputs(self, *args, **kwargs):
        raise NotImplementedError("This method should be implemented in the subclass.")

    def plot_merged_outputs(self, *args, **kwargs):
        raise NotImplementedError("This method should be implemented in the subclass.")

    def save_merged_outputs(self, keys_to_merge=['k','alpha'], **kwargs):
        """
        Merge the temporary files and save the final output. 
        """
        print('\nMerging temporary files and saving final output')

        # Ask to overwrite the final output file if it already exists
        response = ''
        if not self.final_output_file.exists():
            response = 'yes'

        while response not in ['yes', 'no', 'y', 'n']:
            response = input(f'  Warning: Final output file \"{self.final_output_file}\" already exists. Overwrite? (yes/no): ')
            response = response.strip().lower()
            if response in ['no', 'n']:
                raise FileExistsError(f'Not overwriting existing file \"{self.final_output_file}\".')
            elif response in ['yes', 'y']:
                break
            else:
                print('  Invalid input. Please enter \"yes\" or \"no\".')
                continue

        # Merge the temporary files
        self.merge_tmp_outputs(keys_to_merge=keys_to_merge)

        # Flip arrays to be ascending in wavelength
        if np.diff(self.merged_datasets['wave'])[0] < 0:
            self.merged_datasets['wave']  = self.merged_datasets['wave'][::-1]
            for key in keys_to_merge:
                self.merged_datasets[key] = self.merged_datasets[key][::-1]

        # Save the merged data
        print(f'  Saving final output to \"{self.final_output_file}\"')
        utils.save_to_hdf5(
            self.final_output_file, 
            data=self.merged_datasets, 
            attrs=self.merged_attrs
            )

    def merge_tmp_outputs(self, keys_to_merge=['xsec'], **kwargs):
        """
        Combines temporary cross-section files into a single dataset.

        This method reads temporary HDF5 files containing cross-section data,
        combines them into a single dataset, and returns the combined data along
        with the associated attributes.

        Parameters:
        -----------
        keys_to_merge : list of str, optional
            List of keys to read from the temporary files. Default is ['xsec'].
        **kwargs : dict
            Additional keyword arguments.

        Raises:
        -------
        FileNotFoundError
            If no temporary cross-section files are found.
        ValueError
            If the temperature grid is not found in a temporary file.
        """

        # Check if the outputs should be summed
        sum_outputs = self._check_if_sum_outputs()
        
        # Select all temporary files
        tmp_files = list(self.tmp_output_dir.glob('*.hdf5'))
        if len(tmp_files) == 0:
            raise FileNotFoundError(f'No temporary output files found in \"{self.tmp_output_dir}\".')
        # Sort by modification date
        tmp_files.sort(key=lambda x: x.stat().st_mtime)
        
        # Check compatibility of files before combining
        wave_main, P_main, T_main = self._merge_PT_grids(tmp_files, sum_outputs)

        # Combine all files into a single array
        self.merged_datasets = {
            key: np.zeros((len(wave_main), len(P_main), len(T_main)), dtype=np.float64) 
            for key in keys_to_merge
            }
        for i, tmp_file in enumerate(tmp_files):

            datasets, attrs = utils.read_from_hdf5(
                tmp_file, keys_to_read=keys_to_merge+['P','T','wave'], return_attrs=True
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
                idx_P = np.argwhere((P_main == PT[0])).flatten()[0]
                idx_T = np.argwhere((T_main == PT[1])).flatten()[0]

                for key in keys_to_merge:
                    
                    if not has_pressure_axis:
                        # Absorption coefficients
                        dataset_to_add = datasets[key][:,k]
                    else:
                        # Line-by-line cross-sections
                        dataset_to_add = datasets[key][:,j,k]
                    
                    if sum_outputs:
                        # Sum cross-sections
                        self.merged_datasets[key][:,idx_P,idx_T] += dataset_to_add
                    else:
                        # Avoid summing when only one transition file was used
                        self.merged_datasets[key][:,idx_P,idx_T] = dataset_to_add

        # Add the grid-definition
        self.merged_datasets['wave'] = wave_main
        self.merged_datasets['T'] = T_main
        if has_pressure_axis:
            self.merged_datasets['P'] = P_main
        else:
            for key in keys_to_merge:
                # Remove the pressure axis
                self.merged_datasets[key] = np.squeeze(self.merged_datasets[key], axis=1)

        self.merged_attrs = attrs

    def _merge_PT_grids(self, tmp_files, sum_outputs):
        """
        Combine the PT-grids of the temporary files.

        Parameters:
        tmp_files (list): List of temporary files to combine.
        sum_outputs (bool): Whether to sum the outputs.

        Returns:
        tuple: A tuple containing the main wavelength grid, pressure grid, and temperature grid.
        """
        all_PT = []
        for i, tmp_file in enumerate(tmp_files):

            datasets = utils.read_from_hdf5(tmp_file, keys_to_read=['wave','P','T'])

            wave = datasets.get('wave')
            P = datasets.get('P', [0.])
            T = datasets.get('T')

            if wave is None:
                raise ValueError("Wavelength grid not found in temporary file.")
            if T is None:
                raise ValueError("Temperature grid not found in temporary file.")

            if i == 0:
                # Initialize the main arrays
                wave_main = wave.copy()
                P_main = P.copy()
                T_main = T.copy()

            assert (wave_main == wave).all(), "Wavelength grids are not compatible."
            if sum_outputs:
                assert (P_main == P).all(), "Pressure grids are not compatible."
                assert (T_main == T).all(), "Temperature grids are not compatible."

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

        return wave_main, P_main, T_main
        
    def _check_if_sum_outputs(self):
        """
        Check if the cross-sections should be summed.
        """
        transitions_files = self.config.files.get('transitions', [])
        if not isinstance(transitions_files, (list, tuple, np.ndarray)):
            return False
        if len(transitions_files) == 0:
            return False
        return True

    def _read_from_config(self, config):
        """
        Read the configuration file.

        Parameters:
        config (dict): Configuration dictionary.
        """
        print('\nReading parameters from the configuration file')
        self.config = config
        utils.units_warning(self.config)
        
        self.database = getattr(self.config, 'database', None)
        self.species  = getattr(self.config, 'species', None)

        # Wavenumber/wavelength grid
        self.wave_min = getattr(self.config, 'wave_min', 1.0/3.0) # [um]
        self.wave_max = getattr(self.config, 'wave_max', 250.0)
        self.delta_nu = getattr(self.config, 'delta_nu', 0.01)    # [cm^-1]
        
        # Convert to SI units
        self.wave_min *= sc.micron # [um] -> [m]
        self.wave_max *= sc.micron
        self.delta_nu *= 1e2*sc.c # Change to frequency [cm^-1] -> [s^-1]
        
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
            raise ValueError(f'Temporary output \"{self.tmp_output_basename}\" should end with \".hdf5\".')

        # Final output file
        default = default.format(self.species)
        self.final_output_file = getattr(self.config, 'final_output_file', default)
        self.final_output_file = (self.output_data_dir / self.final_output_file).resolve()

    def _configure_nu_grid(self):
        """
        Configure the wavenumber grid.
        """
        self.nu_min = sc.c/self.wave_max # [m] -> [s^-1]
        self.nu_max = sc.c/self.wave_min # [m] -> [s^-1]

        # Number of grid points
        self.N_nu = int((self.nu_max-self.nu_min)/self.delta_nu) + 1

        # Not exact value of delta_nu given above, but done to keep final lambda values fixed
        self.delta_nu = (self.nu_max-self.nu_min) / (self.N_nu-1)
        self.nu_grid  = np.linspace(self.nu_min, self.nu_max, num=self.N_nu, endpoint=True)
        
        self.wave_grid = sc.c/self.nu_grid # [s^-1] -> [m]

        print('\nWavelength-grid:')
        print(f'  Wavelength: {self.wave_min/sc.micron:.2f} - {self.wave_max/sc.micron:.0f} um')
        print(f'  Wavenumber: {self.nu_min/(1e2*sc.c):.0f} - {self.nu_max/(1e2*sc.c):.0f} cm^-1')
        print(f'  Delta nu:   {self.delta_nu/(1e2*sc.c):.3f} cm^-1')
        print(f'  Number of grid points: {self.N_nu}')

    def _configure_coarse_nu_grid(self, adaptive_delta_nu):
        
        # Use the original grid
        if not self.adaptive_nu_grid:
            return self.nu_grid, self.delta_nu
        if self.delta_nu > adaptive_delta_nu:
            return self.nu_grid, self.delta_nu
        
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

        return coarse_nu_grid, coarse_delta_nu