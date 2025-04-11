import numpy as np
from pandas import read_fwf, read_csv
from scipy.interpolate import interp1d

import pathlib

from pyROX import utils, sc
from .lbl import LineByLine

class LBL_Kurucz(LineByLine):
    """
    Class for handling line-by-line cross-sections from Kurucz data.
    """

    parent_dir = pathlib.Path(__file__).parent.resolve()
    parent_dir = parent_dir.parent / 'data'
    atoms_info = read_csv(parent_dir/'atoms_info.csv', index_col=0)

    def download_data(self, config):
        """
        Download data from Kurucz.

        Args:
            config (object): Configuration object containing parameters.
        """
        # Download states from NIST (used for partition function)
        element = config.species
        element = element.capitalize()
        if element not in self.atoms_info.index:
            raise ValueError(f'Element {element} not found in atoms_info.csv.')
        
        url = (
            'https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&'+
            f'spectrum={element}+I&units=0&format=3&output=0&page_size=15&'+
            'multiplet_ordered=0&level_out=on&g_out=on&temp=&submit=Retrieve+Data'
        )
        files = getattr(config, 'files', {})
        file_states = files.get('states', 'NIST_states.tsv')
        file_states = pathlib.Path(file_states).name
        file_states = utils.download(
            url, out_dir=config.input_data_dir, out_name=file_states
            )

        if self.database == 'vald':
            raise NotImplementedError('Please download VALD transition-data manually.')

        # Download Kurucz data
        atomic_number = self.atoms_info.loc[element, 'number']
        ionisation_state = getattr(config, 'ionisation_state', 0)
        atom_id = f'{atomic_number:02d}{ionisation_state:02d}'

        for extension in ['all','pos']:
            # Try different extensions
            url = f'http://kurucz.harvard.edu/atoms/{atom_id}/gf{atom_id}.{extension}'
            file_transitions = utils.download(url, out_dir=config.input_data_dir)
            if file_transitions is not None:
                break
        
        if None in [file_states, file_transitions]:
            raise ValueError('Failed to download all urls.')
    
    def __init__(self, config, **kwargs):
        """
        Initialises the Kurucz object.

        Args:
            config (object): Configuration object containing parameters.
            **kwargs: Additional arguments for initialisation.
        """
        print('-'*60)
        print('  Line-by-line Absorption from VALD/Kurucz')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent LineByLine class
    
    def _read_configuration_parameters(self, config):
        """
        Read parameters specific to Kurucz calculations from the configuration.

        Args:
            config (object): Configuration object containing parameters.
        """
        # Read the common parameters
        super()._read_configuration_parameters(config)

        element = config.species.capitalize()
        self.is_alkali = element in ['Li', 'Na', 'K', 'Rb', 'Cs', 'Fr']
        
        # If alkali, read the E_ion and Z for different vdW-broadening
        self.E_ion = getattr(config, 'E_ion', None) # [cm^-1]
        if self.E_ion is not None:
            self.E_ion *= (1e2*sc.c) * sc.h # [cm^-1] -> [J]
        self.Z = getattr(config, 'Z', 1.0) # Charge of absorber

        # Impact width and shift (Allard et al. 2023)
        self.impact_width_info = self._read_impact_info(
            impact_info=getattr(config, 'impact_width_info', {})
            )
        self.impact_shift_info = self._read_impact_info(
            impact_info=getattr(config, 'impact_shift_info', {})
            )        

    def _read_impact_info(self, impact_info={}):
        """
        Read impact width/shift information from the configuration.

        Args:
            impact_info (dict): Dictionary containing impact information.

        Returns:
            dict: Updated impact information.
        """
        if len(impact_info) == 0:
            # No impact width info provided
            return impact_info

        for perturber, info in impact_info.items():

            # Convert to SI units
            nu_0 = np.atleast_1d(info.get('nu_0', np.nan))
            impact_info[perturber]['nu_0'] = nu_0 * (1e2*sc.c) # [cm^-1] -> [s^-1]

            A = np.atleast_1d(info.get('A', np.nan))
            impact_info[perturber]['A'] = A * (1e2*sc.c) # [cm^-1](-ish) -> [s^-1]

            impact_info[perturber]['b'] = np.atleast_1d(info.get('b', np.nan))
            
            number_density_reference = info.get('number_density_reference', 1e20)
            number_density_reference *= (1e-2)**3 # [cm^-3] -> [m^-3]
            impact_info[perturber]['number_density_reference'] = number_density_reference

            # Read from pressure-broadening info alternatively
            VMR = info.get('VMR', self.pressure_broadening_info[perturber]['VMR'])
            impact_info[perturber]['VMR'] = VMR

        return impact_info

    def _read_partition_function(self, T_grid=np.arange(1,7001+1e-6,1)):
        """
        Read the partition function from the configuration file.

        Args:
            T_grid (array): Temperature grid for partition function calculation.
        """
        file = self.config.files.get('states', None)
        if file is None:
            raise ValueError('States file must be provided in the configuration file.')
        # Load the states from NIST
        states = read_csv(file, sep='\t', engine='python', header=0, index_col=False)

        g = np.array(states['g'])
        E = np.array(states['Level (cm-1)'])

        # (Higher) Ions beyond this index
        idx_u = np.min(np.argwhere(np.isnan(g)))
        self.E_ion = E[idx_u] * sc.h*(1e2*sc.c) # [cm^-1] -> [J]

        g = g[:idx_u]
        E = E[:idx_u]

        # Partition function, sum over states, keep temperature-axis
        partition_function = np.array([np.sum(g*np.exp(-sc.c2*1e2*E/T_i)) for T_i in T_grid])        

        # Make interpolation function, extrapolate outside temperature-range
        interpolation_function = interp1d(
            x=T_grid, y=partition_function, kind='linear', fill_value='extrapolate'
            )
        self.calculate_partition_function = interpolation_function

    def _read_kurucz_transitions(self, input_file):
        """
        Read transitions from a Kurucz input file.

        Args:
            input_file (str): Path to the input file.

        Returns:
            tuple: Arrays of transition parameters (nu_0, E_low, A, g_up, gamma_vdW, gamma_N).
        """
        # Load all transitions at once
        transitions = read_fwf(
            input_file, widths=(11,7,6,12,5,1,10,12,5,1,10,6,6,6), header=None, 
            )
        transitions = np.array(transitions)

        # Oscillator strength
        gf = 10**transitions[:,1].astype(float)

        # Calculate the statistical weights as 2J+1
        g_1 = 2*transitions[:,4].astype(float) + 1
        g_2 = 2*transitions[:,8].astype(float) + 1

        # Transition energies
        E_1 = np.abs(transitions[:,3].astype(float)) * 1e2*sc.c # [cm^-1] -> [s^-1]
        E_2 = np.abs(transitions[:,7].astype(float)) * 1e2*sc.c
        nu_0 = np.abs(E_1 - E_2) # [s^-1]

        E_low = np.minimum(E_1, E_2) # Ensure lower state has lower energy
        E_low *= sc.h # [s^-1] -> [J]

        # Choose statistical weights
        g_up = np.array([
            g_1[i] if E_low[i]!=E_1[i] else g_2[i] for i in range(len(E_low))
        ])
        g_low = np.array([
            g_1[i] if E_low[i]==E_1[i] else g_2[i] for i in range(len(E_low))
        ])

        # Einstein A-coefficient [s^-1]
        A = 2*np.pi*sc.e**2 / (sc.epsilon_0*sc.m_e*sc.c**3) * nu_0**2 * gf / g_up

        # Damping constants given
        gamma_vdW     = np.nan*np.ones_like(nu_0)
        log_gamma_vdW = transitions[:,13].astype(float)
        idx = (log_gamma_vdW != 0.)
        gamma_vdW[idx] = 10**log_gamma_vdW[idx]/(4*np.pi) * (1e-2)**3 # [s^-1 m^3]

        gamma_N     = np.nan*np.ones_like(nu_0)
        log_gamma_N = transitions[:,11].astype(float)
        idx = (log_gamma_N != 0.)
        gamma_N[idx] = 10**log_gamma_N[idx]/(4*np.pi) # [s^-1]
        
        return nu_0, E_low, A, g_up, gamma_vdW, gamma_N

    def _read_vald_transitions(self, input_file):
        """
        Read transitions from a VALD input file.

        Args:
            input_file (str): Path to the input file.

        Raises:
            NotImplementedError: If the method is not implemented.
        """
        raise NotImplementedError('VALD transitions are not implemented yet.')

    def compute_vdw_broadening(self, P, T, E_low, nu_0):
        """
        Calculate Van der Waals broadening.

        Args:
            P (float): Pressure in Pa.
            T (float): Temperature in Kelvin.
            E_low (array): Lower state energies in Joules.
            nu_0 (array): Transition frequencies in s^-1.

        Returns:
            array: Van der Waals broadening in s^-1.
        """
        # Get number density from equation-of-state or ideal-gas
        number_density = self.calculate_number_density(P, T) # [m^-3]

        reduced_mass_H = (sc.m_H*self.mass) / (sc.m_H + self.mass) # [kg]

        idx = np.isfinite(self.gamma_vdW_in_table)
        gamma_vdW = np.zeros_like(self.gamma_vdW_in_table)
        for perturber, info in self.pressure_broadening_info.items():
            # Get the broadening parameters
            number_density_perturber = number_density * info['VMR'] # Perturber density [m^-3]

            mass = info['mass'] # [kg]
            reduced_mass = mass*self.mass / (mass + self.mass)

            alpha = info.get('alpha', None)
            C = info.get('C', None)
            if C is None and alpha is not None:
                C = (reduced_mass_H/reduced_mass)**(0.3) * (alpha/sc.alpha_H)**(0.4)
            
            if C is not None:
                # Use the C-coefficient (Eq. 19 Sharp & Burrows 2007)
                gamma_vdW[idx] += (
                    self.gamma_vdW_in_table[idx] * (T/10e3)**(0.3) * number_density_perturber * C
                )

            if not self.is_alkali or np.all(idx):
                # Only replace empty values for alkali atoms
                continue
            E_up = E_low + nu_0*sc.h # [J]

            # For alkali atoms (Schweitzer et al. 1996) (Eq. A10 Lacy & Burrows 2024)
            C_6 = (
                1.01e-44 * (alpha/sc.alpha_H) * (self.Z+1)**2 *
                np.abs((sc.E_H/(self.E_ion-E_low))**2 - (sc.E_H/(self.E_ion-E_up))**2)
            ) # vdW interaction constant [m^6 s^-1]
            gamma_vdW[~idx] += (
                1.664461/2 * (sc.k*10e3/reduced_mass)**(0.3) * C_6[~idx]**(0.4) * 
                (T/10e3)**(0.3) * number_density_perturber
            ) # [s^-1]

        # Impact width (Allard et al. 2023)
        gamma_vdW = self._apply_impact_Allard_ea_2023(P, T, nu_0, gamma=gamma_vdW)

        return gamma_vdW

    def compute_natural_broadening(self, A):
        """
        Calculate natural broadening.

        Args:
            A (float): Einstein A-coefficient in s^-1.

        Returns:
            float: Natural broadening in s^-1.
        """
        # Natural broadening from Einstein A-coefficient
        gamma_N = super().compute_natural_broadening(A)

        # Use experimental data if available
        idx = np.isfinite(self.gamma_N_in_table)
        gamma_N[idx] = self.gamma_N_in_table[idx]
        return gamma_N

    def pressure_shift(self, P, T, nu_0, delta=None):
        """
        Apply pressure shift to the transition frequency.

        Args:
            P (float): Pressure in Pa.
            T (float): Temperature in Kelvin.
            nu_0 (array): Transition frequencies in s^-1.
            delta (float, optional): Pressure shift coefficient.

        Returns:
            array: Pressure-shifted frequencies.
        """
        # Use the default pressure-shift from the parent class
        nu_0 = super().pressure_shift(P, T, nu_0, delta=delta)
        
        # Apply impact shift (Allard et al. 2023)
        nu_0 = self._apply_impact_Allard_ea_2023(P, T, nu_0)
        return nu_0

    def _apply_impact_Allard_ea_2023(self, P, T, nu_0, gamma=None):
        """
        Apply impact width/shift (Allard et al. 2023).

        Args:
            P (float): Pressure in Pa.
            T (float): Temperature in Kelvin.
            nu_0 (array): Transition frequencies in s^-1.
            gamma (array, optional): Line widths.

        Returns:
            array: Modified transition frequencies or line widths.
        """
        # Copy the nu_0 so that we can add a shift
        nu_0_static = nu_0.copy()

        # Get number density from equation-of-state or ideal-gas
        number_density = self.calculate_number_density(P, T) # [m^-3]

        impact_info = self.impact_width_info # Impact width
        if gamma is None:
            impact_info = self.impact_shift_info # Shift

        for perturber, info in impact_info.items():

            # Get the width/shift parameters
            number_density_perturber = number_density * info['VMR'] # [m^-3]
            number_density_reference = info['number_density_reference']
            A = info['A'] # [s^-1](-ish)
            b = info['b']

            for nu_0_i in info['nu_0']:
                # Check where to replace
                idx = np.isclose(nu_0_static, nu_0_i, atol=1e-8)
                if np.sum(idx) != 1:
                    # No match found, or multiple matches
                    continue

                if gamma is not None:
                    # Width
                    gamma[idx] = (
                        A * T**b * number_density_perturber / number_density_reference
                    ) # [s^-1]
                    continue
                
                # Shift
                nu_0[idx] += (
                    A * T**b * number_density_perturber / number_density_reference
                ) # [s^-1]

        if gamma is not None:
            return gamma
    
        return nu_0

    def process_transitions(self, input_file, **kwargs):
        """
        Read transitions from the input file and compute cross-sections.

        Args:
            input_file (str): Path to the input file.
            **kwargs: Additional arguments.
        """
        input_file = pathlib.Path(input_file)
        print(f'  Reading transitions from \"{input_file}\"')

        if self.database == 'vald':
            nu_0, E_low, A, g_up, gamma_vdW, gamma_N = \
                self._read_vald_transitions(input_file)
        elif self.database == 'kurucz':
            nu_0, E_low, A, g_up, gamma_vdW, gamma_N = \
                self._read_kurucz_transitions(input_file)

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

        self.gamma_vdW_in_table = gamma_vdW[idx_sort]
        self.gamma_N_in_table   = gamma_N[idx_sort]

        # Compute the cross-sections, looping over the PT-grid
        print(f'  Number of lines loaded: {len(A)}')
        self.iterate_over_PT_grid(
            function=self.calculate_cross_sections,
            nu_0=nu_0, S_0=S_0, E_low=E_low, A=A, **kwargs
        )