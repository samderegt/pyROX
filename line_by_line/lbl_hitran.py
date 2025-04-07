import numpy as np
from pandas import read_fwf
import pathlib

from utils import sc
import utils
from .lbl import LineByLine

class HITRAN(LineByLine):

    def download_data(self, config):
        """
        Download data from HITRAN.

        Parameters:
        config (object): Configuration object containing parameters.
        """
        print('\nDownloading data from HITRAN')

        files = []
        for url in config.urls:
            file = utils.download(url, config.input_data_dir)
            files.append(file)

        if None in files:
            raise ValueError('Failed to download all urls.')
    
    def __init__(self, config, **kwargs):
        """
        Initialize the HITRAN object.

        Parameters:
        config (object): Configuration object containing parameters.
        **kwargs: Additional arguments for initialization.
        """
        print('-'*60)
        print('  Line-by-line Absorption from HITRAN/HITEMP')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent LineByLine class

    def _read_configuration_parameters(self, config):
        """
        Read parameters specific to HITRAN calculations from the configuration.

        Parameters:
        config (object): Configuration object containing parameters.
        """
        # Read the common parameters
        super()._read_configuration_parameters(config)
        
        # Read the isotope information
        self.isotope_idx       = getattr(config, 'isotope_idx', 1)
        self.isotope_abundance = getattr(config, 'isotope_abundance', 1.0)
        print(f'\nIsotope index: {self.isotope_idx}, with abundance {self.isotope_abundance:.2e}')
        if self.isotope_abundance == 1.0:
            utils.warnings.warn(
                'HITRAN line-strengths are scaled by the terrestrial isotope abundance. A value of 1.0 may be incorrect.'
            )

        # Remove any quantum-number dependency
        for perturber, info in self.pressure_broadening_info.items():
            self.pressure_broadening_info[perturber]['gamma'] = np.nanmean(info['gamma'])
            self.pressure_broadening_info[perturber]['n']     = np.nanmean(info['n'])
            # Remove 'diet' and 'J' keys if they exist
            self.pressure_broadening_info[perturber].pop('diet', None)
            self.pressure_broadening_info[perturber].pop('J', None)

    def _read_transitions(self, input_file, **kwargs):
        """
        Read transitions from the input file and compute cross-sections.

        Parameters:
        input_file (str): Path to the input file.
        **kwargs: Additional arguments.
        """
        input_file = pathlib.Path(input_file)
        print(f'  Reading transitions from \"{input_file}\"')
        
        # How to handle bz2-compression
        compression = str(input_file.suffix).replace('.','')
        if compression != 'bz2':
            compression = 'infer' # Likely decompressed

        # Read the transitions file in chunks to prevent memory overloads
        transitions_in_chunks = read_fwf(
            input_file, 
            widths=(2,1,12,10,10,5,5,10,4,8), 
            header=None, 
            chunksize=self.N_lines_in_chunk, 
            compression=compression, 
            )
        for transitions in transitions_in_chunks:

            isotope_indices = np.array(transitions.iloc[:,1])
            transitions = np.array(transitions)

            if self.isotope_idx is not None:
                # Select the isotope index
                transitions = transitions[np.isin(isotope_indices.astype(str), list('0123456789'))]
                transitions = transitions[transitions[:,1].astype(int)==self.isotope_idx]

            # Unit conversion
            nu_0  = transitions[:,2].astype(float) * 1e2*sc.c        # [cm^-1] -> [s^-1]
            E_low = transitions[:,7].astype(float) * sc.h*(1e2*sc.c) # [cm^-1] -> [J]
            A = transitions[:,4].astype(float) # [s^-1]

            # [cm^-1/(molec. cm^-2)] -> [s^-1/(molec. m^-2)]
            S_0 = transitions[:,3].astype(float) * (1e2*sc.c) * 1e-4
            S_0 /= self.isotope_abundance # Remove terrestrial abundance ratio

            # Sort by wavenumber
            idx_sort = np.argsort(nu_0)
            nu_0  = nu_0[idx_sort]
            S_0   = S_0[idx_sort]
            E_low = E_low[idx_sort]
            A     = A[idx_sort]

            # Compute the cross-sections, looping over the PT-grid
            print(f'  Number of lines: {len(nu_0)}')
            self.iterate_over_PT_grid(
                function=self.calculate_cross_sections,
                nu_0=nu_0, S_0=S_0, E_low=E_low, A=A, **kwargs
            )