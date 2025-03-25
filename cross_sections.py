import numpy as np
import pathlib
import itertools

from utils import sc
import utils

def load_data_object(config):
    """
    Load the data object from the given configuration.
    """
    
    database = getattr(config, 'database', '').lower()

    import cia
    if database == 'cia_hitran':
        return cia.CIA_HITRAN(config)
    elif database == 'cia_borysov':
        return cia.CIA_Borysov(config)

    import lines
    if database == 'exomol':
        return lines.ExoMol(config)
    elif database in ['hitemp', 'hitran']:
        return lines.HITEMP(config)
    elif database in ['vald', 'kurucz']:
        return lines.VALD_Kurucz(config)

    raise NotImplementedError(f"Database '{config.database}' not implemented.")


class CrossSections:
    """
    Base class for cross-sections.
    """

    @staticmethod
    def mask_arrays(arrays, mask, **kwargs):
        """Apply a mask to the given arrays."""
        return [np.compress(mask, array, **kwargs) for array in arrays]

    def __init__(self, config):

        # Read the common variables
        self._read_from_config(config)
        
        # Get the wavenumber grid
        self.nu_grid, self.delta_nu, self.nu_min, self.nu_max, self.N_grid = \
            utils.get_nu_grid(self.wave_min, self.wave_max, self.delta_nu)

    def download_data(self, *args, **kwargs):
        raise NotImplementedError("This method should be implemented in the subclass.")

    def calculate_cross_sections(self, *args, **kwargs):
        raise NotImplementedError("This method should be implemented in the subclass.")

    def combine_tmp_cross_sections(self, keys_to_read=['xsec'], **kwargs):

        # Check if the cross-sections should be summed
        sum_cross_sections = self._check_if_sum_cross_sections()
        
        print('\nCombining temporary cross-sections')
        
        # Select all temporary files
        tmp_files = list(self.tmp_output_dir.glob('*.hdf5'))
        if len(tmp_files) == 0:
            raise FileNotFoundError("No temporary cross-section files found.")
        # Sort by modification date
        tmp_files.sort(key=lambda x: x.stat().st_mtime)
        
        # Check compatibility of files before combining
        wave_main, P_main, T_main = self._combine_PT_grid(tmp_files, sum_cross_sections)

        # Combine all files into a single array
        datasets_main = {
            key: np.zeros((len(wave_main), len(P_main), len(T_main)), dtype=np.float64) 
            for key in keys_to_read
            }
        for i, tmp_file in enumerate(tmp_files):

            datasets, attrs = utils.read_from_hdf5(
                tmp_file, keys_to_read=keys_to_read+['P','T','wave'], return_attrs=True
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

                for key in keys_to_read:
                    
                    if not has_pressure_axis:
                        # Absorption coefficients
                        dataset_to_add = datasets[key][:,k]
                    else:
                        # Line-by-line cross-sections
                        dataset_to_add = datasets[key][:,j,k]
                    
                    if sum_cross_sections:
                        # Sum cross-sections
                        datasets_main[key][:,idx_P,idx_T] += dataset_to_add
                    else:
                        # Avoid summing when only one transition file was used
                        datasets_main[key][:,idx_P,idx_T] = dataset_to_add

        # Add the grid-definition
        datasets_main['wave'] = wave_main
        datasets_main['T'] = T_main
        if has_pressure_axis:
            datasets_main['P'] = P_main

        return datasets_main, attrs

    def _combine_PT_grid(self, tmp_files, sum_cross_sections):
        """
        Combine the PT-grid of the temporary files.
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
            if sum_cross_sections:
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
        
    def _check_if_sum_cross_sections(self):
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
        """Read the configuration file."""

        self.config = config
        utils.units_warning(self.config)
        
        self.database = getattr(self.config, 'database', None)
        self.species  = getattr(self.config, 'species', None)

        # Wavenumber/wavelength grid
        self.wave_min = getattr(self.config, 'wave_min', 1.0/3.0) # [um]
        self.wave_max = getattr(self.config, 'wave_max', 250.0)
        self.delta_nu = getattr(self.config, 'delta_nu', 0.01)    # [cm^-1]
        
        # Convert to SI units
        self.wave_min *= 1e-6 # [um] -> [m]
        self.wave_max *= 1e-6
        self.delta_nu *= 100.0*sc.c # Change to frequency [cm^-1] -> [s^-1]
        
        # Input/output directories
        self.input_data_dir = getattr(self.config, 'input_data_dir', f'./{self.species}/input_data/')
        self.input_data_dir = pathlib.Path(self.input_data_dir).resolve()
        self.input_data_dir.mkdir(parents=True, exist_ok=True)

        self.output_data_dir = getattr(self.config, 'output_data_dir', f'./{self.species}/')
        self.output_data_dir = pathlib.Path(self.output_data_dir).resolve()
        self.output_data_dir.mkdir(parents=True, exist_ok=True)

        # Temporary output files
        self.tmp_output_dir = self.output_data_dir / 'tmp'
        self.tmp_output_file = getattr(self.config, 'tmp_output_file', 'xsec_{}.hdf5')
        self.tmp_output_file = (self.tmp_output_dir / self.tmp_output_file).resolve()

        # Final output file
        self.final_output_file = getattr(self.config, 'final_output_file', f'{self.species}.hdf5')
        self.final_output_file = (self.output_data_dir / self.final_output_file).resolve()
        if self.final_output_file.exists():
            response = input(f"Warning: Final output file '{self.final_output_file}' already exists. Do you want to overwrite it? (yes/no): ")
            if response.lower() != 'yes':
                raise FileExistsError(f"Final output file '{self.final_output_file}' already exists and will not be overwritten.")

        assert hasattr(self.config, 'files'), 'No files specified'
