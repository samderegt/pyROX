import numpy as np
from pandas import read_fwf

import pathlib
import re

from pyROX import utils, sc
from .lbl import LineByLine

class LBL_ExoMol(LineByLine):
    """
    Class for handling line-by-line cross-sections from ExoMol data.
    """

    def download_data(self, config):
        """
        Downloads data from ExoMol.

        Args:
            config (object): Configuration object containing parameters.
        """
        print('\nDownloading data from ExoMol')

        files = []
        file_def_json = None
        for url in config.urls:
            file = utils.download(url, config.input_data_dir)
            files.append(file)

            if url.endswith('.def.json'):
                # Check if the url is a definition file
                url_base = url.replace('.def.json','')
                file_def_json = file
        
        # Download partition-function and states files
        file_partition = utils.download(f'{url_base}.pf', config.input_data_dir)
        file_states = utils.download(f'{url_base}.states.bz2', config.input_data_dir)

        # Read definition file
        import json
        with open(file_def_json, 'r') as f:
            file_def_dict = json.load(f)

        # Download transitions files
        file_transitions = [utils.download(f'{url_base}.trans.bz2', config.input_data_dir)]
        if None in file_transitions:
            # Split up into multiple files
            wavenumbers = np.linspace(
                0, file_def_dict['dataset']['transitions']['max_wavenumber'], 
                file_def_dict['dataset']['transitions']['number_of_transition_files'], 
                endpoint=False, dtype=int
            )
            file_transitions = []
            for i in range(len(wavenumbers)-1):
                nu_min_i, nu_max_i = wavenumbers[i], wavenumbers[i+1]
                file_transitions_i = utils.download(
                    f'{url_base}__{nu_min_i:05d}-{nu_max_i:05d}.trans.bz2', config.input_data_dir
                )
                file_transitions.append(file_transitions_i)

        if None in [*files, file_partition, file_states, *file_transitions]:
            raise ValueError('Failed to download all urls.')
    
    def __init__(self, config, **kwargs):
        """
        Initialises the ExoMol object.

        Args:
            config (object): Configuration object containing parameters.
            **kwargs: Additional arguments for initialisation.
        """
        print('-'*60)
        print('  Line-by-line Absorption from ExoMol')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent LineByLine class

    def _read_configuration_parameters(self, config):
        """
        Reads parameters specific to ExoMol calculations from the configuration.

        Args:
            config (object): Configuration object containing parameters.
        """
        # Read the common parameters
        super()._read_configuration_parameters(config)

        # Copy so that we can modify the broadening parameters per transition
        self.pressure_broadening_info_copy = self.pressure_broadening_info.copy()

    def _read_broadening_per_transition(self, J_l, J_u, chunk_size=100_000):
        """
        Reads broadening parameters for each transition.

        Args:
            J_l (array): Lower state rotational quantum numbers.
            J_u (array): Upper state rotational quantum numbers.
            chunk_size (int): Size of chunks to process at a time.
        """
        for perturber, info in self.pressure_broadening_info_copy.items():
            # Get the broadening parameters
            gamma = info['gamma']
            n     = info['n']

            if callable(gamma):
                # User-provided function
                self.pressure_broadening_info[perturber]['gamma'] = gamma(J_l)
            else:
                # Used as default if no quantum-number match is found
                self.pressure_broadening_info[perturber]['gamma'] = \
                    np.nanmean(gamma)*np.ones_like(J_l) # Extend to all lines

            if callable(n):
                # User-provided function
                self.pressure_broadening_info[perturber]['n'] = n(J_l)
            else:
                self.pressure_broadening_info[perturber]['n'] = \
                    np.nanmean(n)*np.ones_like(J_l)

            if callable(gamma) or callable(n):
                continue

            J_in_table = info.get('J')
            diet = info.get('diet')
            if (J_in_table is None) or (diet not in ['a0', 'm0']):
                # No quantum-number dependence
                continue
            
            # Read in chunks to avoid memory overload
            for idx_l in range(0, len(J_l), chunk_size):
                idx_h = min([idx_l+chunk_size, len(J_l)])

                J_l_to_match = J_l[idx_l:idx_h]
                
                if diet == 'm0':
                    # Check if transition is in R-branch (i.e. lower J quantum 
                    # number is +1 higher than upper state).
                    # In that case, 4th column in .broad is |m|=J_l+1.
                    delta_J = J_l[idx_l:idx_h] - J_u[idx_l:idx_h]
                    J_l_to_match[(delta_J==+1)] += 1

                # Indices in .broad table corresponding to each transition
                indices_J = np.argwhere(J_in_table[None,:]==J_l_to_match[:,None])

                # Update the gamma and n values
                self.pressure_broadening_info[perturber]['gamma'][idx_l+indices_J[:,0]] = \
                    np.array(gamma)[indices_J[:,1]]
                self.pressure_broadening_info[perturber]['n'][idx_l+indices_J[:,0]] = \
                    np.array(n)[indices_J[:,1]]

    def _read_states(self):
        """
        Reads the states from the states file.
        """
        states_file = self.config.files.get('states', None)
        if states_file is None:
            raise ValueError('No states file specified in the configuration.')
        states_file = pathlib.Path(states_file)

        print(f'  Reading states from \"{states_file}\"')

        # How to handle bz2-compression
        compression = str(states_file.suffix).replace('.','')
        if compression == 'bz2':
            import bz2
            with bz2.open(states_file, 'rb') as f:
                first_line = f.readline()
        else:
            with open(states_file, 'r') as f:
                first_line = f.readline()

        # Infer column-widths
        col_widths = [len(col) for col in re.findall('\s+\S+', str(first_line))]
        
        # Load states (ID, E, g, J)
        states = read_fwf(
            states_file, widths=col_widths[:4], header=None, compression=compression
            )
        states = np.array(states)

        # Sort the states by their IDs
        states = states[np.argsort(states[:,0])]

        self.states_ID = states[:,0].astype(int)
        self.states_E  = states[:,1].astype(float) * sc.h*(1e2*sc.c) # [cm^-1] -> [J]
        self.states_g  = states[:,2].astype(int)
        self.states_J  = states[:,3].astype(float)

        ID_diff = np.diff(self.states_ID)
        if np.any(ID_diff == 0):
            raise ValueError(f'Duplicate state IDs found in states file (ID: {self.states_ID[ID_diff==0]}).')
        if np.any(ID_diff > 1):
            raise ValueError(f'Some state IDs are skipped in states file (ID: {self.states_ID[ID_diff>1]}).')

    def process_transitions(self, input_file, **kwargs):
        """
        Reads transitions from the input file and compute cross-sections.

        Args:
            input_file (str): Path to the input file.
            **kwargs: Additional arguments.
        """
        self._read_states()

        print(f'  Reading transitions from \"{input_file}\"')
        i = 0
        state_ID_up = []; state_ID_low = []; A = []
        is_end_of_file = False

        # Read the transitions file in chunks to prevent memory overloads
        import bz2
        with bz2.open(input_file) as f:

            while not is_end_of_file:
                # Read the next line
                line = f.readline()

                is_end_of_file = (not line)
                is_chunk = (i%self.N_lines_in_chunk == 0)

                if i == 0:
                    # Replace any tabs with spaces
                    line = line.replace(b'\t', b' ')

                    # Infer column-widths from first line
                    col_widths = [len(col) for col in re.findall(b'\s+\S+', line)]
                    col_indices = np.cumsum(col_widths)
                
                elif is_chunk or is_end_of_file:

                    if len(A) == 0:
                        break

                    # Compute when N_lines_in_chunk is reached
                    print(f'  Number of lines: {len(A)}')

                    idx_l = np.searchsorted(self.states_ID, np.array(state_ID_low, dtype=int))
                    idx_u = np.searchsorted(self.states_ID, np.array(state_ID_up, dtype=int))
                    A = np.array(A, dtype=float)

                    E_up  = self.states_E[idx_u] # [J]
                    E_low = self.states_E[idx_l]
                    nu_0  = np.abs(E_up - E_low) / sc.h # [J] -> [s^-1]

                    g_up = self.states_g[idx_u] # State degeneracy

                    J_l = self.states_J[idx_l] # Rotational quantum number
                    J_u = self.states_J[idx_u]

                    # Remove any nu = 0 transitions
                    A, E_low, nu_0, g_up, J_l, J_u = self._mask_arrays(
                        [A, E_low, nu_0, g_up, J_l, J_u], mask=(nu_0 > 0)
                    )

                    # Line-strengths at reference temperature
                    S_0 = (
                        A*g_up / (8*np.pi*(nu_0/sc.c)**2) *
                        np.exp(-E_low/(sc.k*self.T_0)) / self.q_0 * 
                        (1-np.exp(-sc.h*nu_0/(sc.k*self.T_0)))
                    ) # [s^-1/(molec. m^-2)]

                    # Sort by wavenumber
                    idx_sort = np.argsort(nu_0)
                    nu_0  = nu_0[idx_sort]
                    S_0   = S_0[idx_sort]
                    E_low = E_low[idx_sort]
                    A     = A[idx_sort]

                    J_l = J_l[idx_sort]
                    J_u = J_u[idx_sort]
                    
                    # Get J-specific broadening parameters
                    self._read_broadening_per_transition(J_l, J_u)

                    # Compute the cross-sections, looping over the PT-grid
                    self.iterate_over_PT_grid(
                        function=self.calculate_cross_sections,
                        nu_0=nu_0, S_0=S_0, E_low=E_low, A=A, **kwargs
                    )

                    # Reset
                    state_ID_up = []; state_ID_low = []; A = []

                if not is_end_of_file:
                    # Access info on upper and lower states
                    state_ID_up.append(line[0:col_indices[0]])
                    state_ID_low.append(line[col_indices[0]:col_indices[1]])
                    A.append(
                        line[col_indices[1]:col_indices[2]] # Einstein A-coefficient [s^-1]
                        )

                i += 1