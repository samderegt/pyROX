import numpy as np
from pandas import read_fwf
from scipy.interpolate import interp1d
from scipy.special import wofz

import pathlib
import datetime

from tqdm import tqdm

from pyROX import utils, sc, CrossSections

class LineProfileHelper:
    """
    Helper class for line profile calculations.
    """

    def compute_line_strength(self, T, S_0, E_low, nu_0):
        """
        Calculates the line strength for a given temperature.

        Args:
            T (float): Temperature in Kelvin.
            S_0 (float): Line strength at reference temperature.
            E_low (float): Lower state energy in Joules.
            nu_0 (float): Transition frequency in s^-1.

        Returns:
            float: Line strength in s^-1/(molecule m^-2).
        """
        # Partition function
        q = self.calculate_partition_function(T)

        # Gordon et al. (2017) (E_low: [J]; nu_0: [s^-1])
        S = (
            S_0 * (self.q_0/q) * np.exp(E_low/sc.k*(1/self.T_0-1/T)) *
            (1-np.exp(-sc.h*nu_0/(sc.k*T))) / (1-np.exp(-sc.h*nu_0/(sc.k*self.T_0)))
        )
        return S # [s^-1/(molecule m^-2)]

    def normalise_wing_cutoff(self, S, cutoff_distance, gamma_L):
        """
        Normalises the line strength to account for wing cutoff.

        Args:
            S (array): Line strengths.
            cutoff_distance (float): Distance for wing cutoff.
            gamma_L (array): Lorentzian widths.

        Returns:
            array: Normalised line strengths.
        """
        # Eq. A6 Lacy & Burrows (2023)
        return S / ((2/np.pi)*np.arctan(cutoff_distance/gamma_L))

    def pressure_shift(self, P, T, nu_0, delta=None):
        """
        Applies pressure shift to the transition frequency.

        Args:
            P (float): Pressure in Pa.
            T (float): Temperature in Kelvin.
            nu_0 (array): Transition frequencies in s^-1.
            delta (float, optional): Pressure shift coefficient.

        Returns:
            array: Pressure-shifted frequencies.
        """
        if delta is None:
            return nu_0 # No pressure shift

        # Calculate the pressure shift
        return nu_0 + delta*(P/sc.atm) # [s^-1]

    def compute_vdw_broadening(self, P, T, E_low, **kwargs):
        """
        Calculates Van der Waals broadening.

        Args:
            P (float): Pressure in Pa.
            T (float): Temperature in Kelvin.
            E_low (float): Lower state energy in Joules.
            **kwargs: Additional parameters.

        Returns:
            float: Van der Waals broadening in s^-1.
        """
        # Van der Waals broadening
        gamma_vdW = np.zeros_like(E_low)
        for perturber, info in self.pressure_broadening_info.items():
            # Get the broadening parameters
            VMR = info['VMR']
            gamma = info['gamma'] # [s^-1]
            n = info['n']

            # Calculate the broadening parameter
            gamma_vdW += gamma * (self.T_0/T)**n * (P/sc.atm) * VMR # [s^-1]

        return gamma_vdW

    def compute_natural_broadening(self, A):
        """
        Calculates natural broadening.

        Args:
            A (float): Einstein A-coefficient in s^-1.

        Returns:
            float: Natural broadening in s^-1.
        """
        return A / (4*np.pi) # [s^-1]

    def compute_doppler_broadening(self, T, nu_0):
        """
        Calculates Doppler broadening.

        Args:
            T (float): Temperature in Kelvin.
            nu_0 (float): Transition frequency in s^-1.

        Returns:
            float: Doppler broadening in s^-1.
        """
        return np.sqrt(2*sc.k*T/self.mass) * nu_0/sc.c # [s^-1]

    def compute_lorentz_width(self, gamma_vdW, gamma_N):
        """
        Calculates Lorentzian width.

        Args:
            gamma_vdW (float): Van der Waals broadening in s^-1.
            gamma_N (float): Natural broadening in s^-1.

        Returns:
            float: Lorentzian width in s^-1.
        """
        return gamma_vdW + gamma_N # [s^-1]

    def compute_voigt_width(self, gamma_G, gamma_L):
        """
        Calculates Voigt width.

        Args:
            gamma_G (float): Gaussian width in s^-1.
            gamma_L (float): Lorentzian width in s^-1.

        Returns:
            float: Voigt width in s^-1.
        """
        return 0.5346*gamma_L+np.sqrt(0.2166*gamma_L**2 + gamma_G**2) # [s^-1]

    def _bin_transitions_for_local_cutoff(self, nu_0):
        """
        Bins the transitions to the nearest sub-sampled grid point.

        Args:
            nu_0 (array): Transition frequencies in s^-1.

        Returns:
            array: Indices of transitions per bin.
        """

        # Mask the grid to reduce array sizes
        idx_l = np.searchsorted(self.nu_grid, nu_0.min()) - 2
        idx_l = np.maximum(idx_l, 0)
        idx_h = np.searchsorted(self.nu_grid, nu_0.max()) + 2
        idx_h = np.minimum(idx_h, len(self.nu_grid))

        # Divide the grid into sub-sampled bins
        nu_grid_fine = self.nu_grid[idx_l:idx_h]
        delta_nu_fine = (
            np.linspace(0., 1., self.sub_sampling_local_cutoff, endpoint=False)[None,:] * 
            np.diff(nu_grid_fine)[:,None]
            )
        nu_grid_fine = np.sort((nu_grid_fine[:-1,None] + delta_nu_fine).flatten())
        
        # Bin the lines to the nearest (fine) grid point
        nu_bin = np.searchsorted(nu_grid_fine, nu_0)
        nu_bin = np.maximum(nu_bin, 0) # First and last bin
        nu_bin = np.minimum(nu_bin, nu_grid_fine.size-1)

        # Upper and lower indices of lines within bins
        _, nu_bin_idx = np.unique(nu_bin, return_index=True)
        nu_bin_idx    = np.append(nu_bin_idx, len(nu_bin))
        
        return nu_bin_idx

    def apply_local_cutoff(self, S, nu_0, factor):
        """
        Applies local cutoff to line strengths. Retain lines with 
        a certain fraction of cumulative strength.

        Args:
            S (array): Line strengths.
            nu_0 (array): Transition frequencies in s^-1.
            factor (float): Fraction of cumulative strength to retain.

        Returns:
            array: Modified line strengths.
        """
        if factor == 0.:
            return S

        # Bin the lines to the nearest sub-sampled grid point
        nu_bin_idx = self._bin_transitions_for_local_cutoff(nu_0)
         
        for k in range(len(nu_bin_idx)-1):
            idx_l = nu_bin_idx[k]
            idx_h = nu_bin_idx[k+1]
            if idx_h-idx_l < 2:
                continue # One/no lines in bin

            # All lines in this bin
            S_range = S[idx_l:idx_h]

            # Sort by line-strength
            idx_sort = np.argsort(S_range)
            S_sort   = S_range[idx_sort]
            S_cumsum = np.cumsum(S_sort)

            # Lines contributing less than 'factor' to total strength
            S_cutoff = S_cumsum[-1] * factor
            i_search = np.searchsorted(S_cumsum, S_cutoff)
            if i_search == 0:
                continue # No lines are weak enough

            # Add weaker line-strengths to strongest line
            S_cumsum_weak = S_cumsum[i_search-1]
            S_range[idx_sort[-1]] += S_cumsum_weak

            # Ignore weak lines
            S_above_cutoff = S_sort[i_search]
            S_range[S_range<S_above_cutoff] = 0.

            S[idx_l:idx_h] = S_range

        return S
    
    def calculate_line_profiles(
            self, 
            nu_0, 
            S, 
            gamma_L, 
            gamma_G, 
            nu_grid, 
            wing_cutoff_distance, 
            delta_nu, 
            ):
        """
        Calculates line profiles in chunks to optimise speed.

        Args:
            nu_0 (array): Transition frequencies in s^-1.
            S (array): Line strengths.
            gamma_L (array): Lorentzian widths.
            gamma_G (array): Gaussian widths.
            nu_grid (array): Wavenumber grid.
            wing_cutoff_distance (float): Wing cutoff distance in s^-1.
            delta_nu (float): Wavenumber grid spacing in s^-1.

        Returns:
            array: Opacity cross-section on the grid.
        """
        if not np.isnan(self.delta_nu):
            # Fixed wavenumber spacing
            return self._calculate_line_profiles_in_chunks_fixed_delta_nu(
                nu_0, S, gamma_L, gamma_G, nu_grid, wing_cutoff_distance, delta_nu
                )
        
        # Other wavenumber grids
        return self._calculate_line_profiles_in_chunks(
            nu_0, S, gamma_L, gamma_G, nu_grid, wing_cutoff_distance, delta_nu
            )

    def _calculate_line_profiles_one_by_one(
            self,
            nu_0, 
            strength, 
            y, 
            gamma_G,
            nu_grid,
            slice_nu_grid,
            ):
        """
        Calculates line profiles one by one.

        Args:
            nu_0 (array): Transition frequencies in s^-1.
            strength (array): Line strengths.
            y (array): Complex part of Faddeeva coordinate.
            gamma_G (array): Gaussian widths.
            nu_grid (array): Wavenumber grid.
            slice_nu_grid (array): Slices of the wavenumber grid.

        Returns:
            array: Opacity cross-section on the grid.
        """

        sigma = np.zeros_like(nu_grid)
        for nu_0_i, strength_i, y_i, gamma_G_i, slice_nu_grid_i in zip(nu_0, strength, y, gamma_G, slice_nu_grid):
            # Coordinates for the Faddeeva function
            z = (nu_grid[slice_nu_grid_i]-nu_0_i)/gamma_G_i + y_i*1j

            # Faddeeva function
            R_w = np.real(wofz(z))

            # Add line to total cross-section
            sigma[slice_nu_grid_i] += strength_i*R_w

        return sigma
    
    def _calculate_line_profiles_in_chunks(
            self, 
            nu_0, 
            S, 
            gamma_L, 
            gamma_G, 
            nu_grid, 
            wing_cutoff_distance, 
            delta_nu,
            N_lines_in_chunk=1,
            ):
        """
        Calculates line profiles in chunks to optimise speed.

        Args:
            nu_0 (array): Transition frequencies in s^-1.
            S (array): Line strengths.
            gamma_L (array): Lorentzian widths.
            gamma_G (array): Gaussian widths.
            nu_grid (array): Wavenumber grid.
            wing_cutoff_distance (float): Wing cutoff distance in s^-1.
            delta_nu (float): Wavenumber grid spacing in s^-1.
            N_lines_in_chunk (int): Number of lines to process in each chunk.

        Returns:
            array: Opacity cross-section on the grid.
        """

        # Only consider number of lines at a time
        N_chunks = int(np.ceil(len(S)/N_lines_in_chunk))
        N_nu_grid = len(nu_grid)

        # Left and right indices of the line profile within the grid
        idx_nu_grid_l = np.searchsorted(nu_grid, nu_0-wing_cutoff_distance) - 1
        idx_nu_grid_l = np.maximum(0, idx_nu_grid_l)
        idx_nu_grid_h = np.searchsorted(nu_grid, nu_0+wing_cutoff_distance)
        idx_nu_grid_h = np.minimum(N_nu_grid, idx_nu_grid_h)
        slice_nu_grid = [slice(idx_nu_grid_l[i], idx_nu_grid_h[i]) for i in range(len(nu_0))]

        # Complex part of Faddeeva coordinate
        y = gamma_L / gamma_G
        strength = S / (gamma_G*np.sqrt(np.pi))

        if N_lines_in_chunk == 1:
            # Only one line per chunk
            return self._calculate_line_profiles_one_by_one(nu_0, strength, y, gamma_G, nu_grid, slice_nu_grid)

        sigma = np.zeros_like(nu_grid)
        for ch in range(N_chunks):
        
            # Upper and lower indices of lines in current chunk
            idx_chunk_l = int(ch*N_lines_in_chunk)
            idx_chunk_h = idx_chunk_l + N_lines_in_chunk
            idx_chunk_h = np.minimum(idx_chunk_h, len(S)) # At last chunk
            slice_chunk = slice(idx_chunk_l, idx_chunk_h)

            # Lines in current chunk
            nu_0_chunk    = nu_0[slice_chunk]
            y_chunk       = y[slice_chunk]
            gamma_G_chunk = gamma_G[slice_chunk]
            strength_chunk = strength[slice_chunk]
            
            # Wavenumber-grid slices of lines in current chunk
            slice_nu_grid_chunk = slice_nu_grid[slice_chunk]

            # Coordinates for the Faddeeva function
            z = [(nu_grid[sl_nu_grid_i]-nu_0_chunk[i])/gamma_G_chunk[i] + y_chunk[i]*1j
                for i, sl_nu_grid_i in enumerate(slice_nu_grid_chunk)]
            
            # Indices to separate back into lines
            len_profile = [len(z_i) for z_i in z]
            idx_profile = [0, *np.cumsum(len_profile)]

            # Collapse into a 1D array
            z = np.concatenate(z)
            R_w = np.real(wofz(z))

            for i, slice_nu_grid_i in enumerate(slice_nu_grid_chunk):

                idx_profile_l = idx_profile[i]
                idx_profile_h = idx_profile[i+1]

                # Add line to total cross-section
                sigma[slice_nu_grid_i] += strength_chunk[i]*R_w[idx_profile_l:idx_profile_h]
            
        return sigma

    def _calculate_line_profiles_in_chunks_fixed_delta_nu(
            self, 
            nu_0, 
            S, 
            gamma_L, 
            gamma_G, 
            nu_grid, 
            wing_cutoff_distance, 
            delta_nu, 
            N_lines_in_chunk=1,
            ):
        """
        Calculates line profiles for a fixed wavenumber spacing. 
        This method is more efficient than calculate_line_profiles_in_chunks().

        Args:
            nu_0 (array): Transition frequencies in s^-1.
            S (array): Line strengths.
            gamma_L (array): Lorentzian widths.
            gamma_G (array): Gaussian widths.
            nu_grid (array): Wavenumber grid.
            wing_cutoff_distance (float): Wing cutoff distance in s^-1.
            delta_nu (float): Wavenumber grid spacing in s^-1.
            N_lines_in_chunk (int): Number of lines to process in each chunk.
        
        Returns:
            array: Opacity cross-section on the grid.
        """
        
        # Only consider number of lines at a time
        N_chunks = int(np.ceil(len(S)/N_lines_in_chunk))
        N_nu_grid = len(nu_grid)

        # Construct the temporary line-profile grid        
        nu_line = np.arange(0, wing_cutoff_distance+delta_nu, delta_nu)
        nu_line = np.concatenate([-nu_line[1:][::-1], nu_line])
        nu_line = np.repeat(nu_line[None,:], N_lines_in_chunk, axis=0)
        
        # Indices to insert the lines back into nu_grid
        idx_nu_grid = np.searchsorted(nu_grid, nu_0) - 1
        dnu_0 = nu_0 - nu_grid[idx_nu_grid] # Separation between grid-point and nu_0

        # Left and right indices of the line profile within the grid
        idx_nu_grid_l = np.maximum(0, idx_nu_grid-nu_line.shape[1]//2)
        idx_nu_grid_h = np.minimum(N_nu_grid, idx_nu_grid+nu_line.shape[1]//2+1)
        slice_nu_grid = [slice(idx_nu_grid_l[i], idx_nu_grid_h[i]) for i in range(len(nu_0))]
        
        # Left and right indices within the temporary grid
        idx_nu_line_l = nu_line.shape[1]//2 - (idx_nu_grid-idx_nu_grid_l)
        idx_nu_line_l = np.maximum(0, idx_nu_line_l)
        idx_nu_line_h = nu_line.shape[1]//2 + (idx_nu_grid_h-idx_nu_grid)
        idx_nu_line_h = np.minimum(nu_line.shape[1], idx_nu_line_h)
        slice_nu_line = [slice(idx_nu_line_l[i], idx_nu_line_h[i]) for i in range(len(nu_0))]

        # Reshape before the loop
        nu_0 = nu_0[:,None]
        dnu_0 = dnu_0[:,None]
        S = S[:,None]
        gamma_L = gamma_L[:,None]
        gamma_G = gamma_G[:,None]

        # Complex part of Faddeeva coordinate
        y = gamma_L / gamma_G
        strength = S / (gamma_G*np.sqrt(np.pi))

        sigma = np.zeros_like(nu_grid)
        for ch in range(N_chunks):

            # Upper and lower indices of lines in current chunk
            idx_chunk_l = int(ch*N_lines_in_chunk)
            idx_chunk_h = idx_chunk_l + N_lines_in_chunk
            if idx_chunk_h > len(S):
                idx_chunk_h = len(S) # Last chunk
                nu_line = nu_line[:idx_chunk_h-idx_chunk_l]

            slice_chunk = slice(idx_chunk_l, idx_chunk_h)

            # Lines in current chunk
            dnu_0_chunk   = dnu_0[slice_chunk]
            y_chunk       = y[slice_chunk]
            gamma_G_chunk = gamma_G[slice_chunk]
            strength_chunk = strength[slice_chunk]

            # Faddeeva function on multiple lines at once
            z = (nu_line-dnu_0_chunk)/gamma_G_chunk + y_chunk*1j
            V = strength_chunk * np.real(wofz(z))

            # Wavenumber-grid slices of lines in current chunk
            slice_nu_grid_chunk = slice_nu_grid[slice_chunk]
            slice_nu_line_chunk = slice_nu_line[slice_chunk]

            for i, (sl_nu_grid_i, sl_nu_line_i, V_i) in enumerate(zip(slice_nu_grid_chunk, slice_nu_line_chunk, V)):
                # Add line to total cross-section
                sigma[sl_nu_grid_i] += V_i[sl_nu_line_i]

        return sigma

    
class LineByLine(CrossSections, LineProfileHelper):
    """
    Base class for line-by-line cross-sections.
    """

    def __init__(self, config, **kwargs):
        """
        Initialises the LineByLine object.

        Args:
            config (object): Configuration object containing parameters.
            **kwargs: Additional arguments for initialisation.
        """
        # Initialise the CrossSections parent class
        super().__init__(config, **kwargs)

    def _read_configuration_parameters(self, config):
        """
        Reads parameters specific to line-by-line calculations from the configuration.

        Args:
            config (object): Configuration object containing parameters.
        """
        # Read the common parameters
        super()._read_configuration_parameters(config)
        
        # Read the parameters for line-by-line calculations
        self._read_mass()
        self._read_pressure_broadening_info()
        self._read_equation_of_state()
        self._read_partition_function()

        self.N_lines_in_chunk = getattr(self.config, 'N_lines_in_chunk', 10_000_000)
        
        # Reference temperature and partition function
        self.T_0 = getattr(self.config, 'T_0', 296.)
        self.q_0 = self.calculate_partition_function(self.T_0)

        # Transition energies to ignore
        self.nu_0_to_ignore = np.atleast_1d(getattr(config, 'nu_0_to_ignore', []))
        self.nu_0_to_ignore *= 1e2*sc.c # [cm^-1] -> [s^-1]

        # Read the cutoff parameters
        wing_cutoff = getattr(self.config, 'wing_cutoff', None) # [cm^-1]
        if wing_cutoff is None:
            # Use Gharib-Nezhad et al. (2024) as default
            wing_cutoff = lambda gamma_V, P: 25 if P<=200 else 100

        # Convert input [s^-1]->[cm^-1] and [Pa]->[bar], and output [cm^-1]->[s^-1]
        self.wing_cutoff = lambda gamma_V, P: wing_cutoff(gamma_V/(1e2*sc.c), P/sc.bar) * 1e2*sc.c
        
        # Maximum separation, in case of pressure-dependent cutoff
        self.wing_cutoff_max = getattr(self.config, 'wing_cutoff_max', np.inf) # [cm^-1]
        self.wing_cutoff_max *= 1e2*sc.c # [cm^-1] -> [s^-1]

        # Line-strength cutoffs
        self.global_cutoff = getattr(self.config, 'global_cutoff', 0.) # [cm^1 molecule^-1]
        self.global_cutoff *= 1e-2 # [m^1 molecule^-1]

        self.local_cutoff = getattr(self.config, 'local_cutoff', 0.)  # [fraction of cumulative]
        # Finer resolution to determine which lines to keep
        self.sub_sampling_local_cutoff = getattr(self.config, 'sub_sampling_local_cutoff', 5)

        # (P,T)-grid to calculate cross-sections on
        self.P_grid = np.atleast_1d(self.config.P_grid) * 1e5 # [bar] -> [Pa]
        self.T_grid = np.atleast_1d(self.config.T_grid)
        self.N_PT = len(self.P_grid) * len(self.T_grid)

        print('\nPT-grid:')
        print(f'  P: {self.P_grid/1e5} bar')
        print(f'  T: {self.T_grid} K')

    def _read_equation_of_state(self):
        """
        Reads the equation of state for the gas mixture.
        """
        EOS_table = getattr(self.config, 'EOS_table', None)
        if EOS_table is None:
            # Assume ideal gas
            self.calculate_number_density = lambda P, T: P/(sc.k*T)
            return

        # Load equation-of-state table and convert into SI
        from pandas import read_csv
        EOS_table = read_csv(EOS_table, delim_whitespace=True, skiprows=1)
        T = np.unique(EOS_table['T[K]'])
        P = np.unique(EOS_table['P[dyne_cm-2]'])
        P *= 1e-1 # [dyne cm^-2] -> [Pa]            

        mass_density_grid = EOS_table.pivot(
            index='P[dyne_cm-2]', columns='T[K]', values='rho[g_cm-3]'
            ).values
        number_density_grid = mass_density_grid*1e3 / self.mean_mass # [g cm^-3] -> [m^-3]
        
        # Create an interpolation function for the density grid
        from scipy.interpolate import RegularGridInterpolator
        interp_func = RegularGridInterpolator(
            points=(P, T), values=number_density_grid, method='cubic', 
            bounds_error=False, fill_value=None
            )

        # Define the number density function
        self.calculate_number_density = lambda P, T: interp_func([P, T])

    def _read_mass(self):
        """
        Reads the mass of the species from the configuration or database.
        """
        self.mass = getattr(self.config, 'mass', None)
        if (self.mass is None) and (self.database in ['hitran', 'exomol']):
            raise ValueError('Mass of the species must be provided in the configuration file.')

        elif (self.mass is None) and (self.database in ['vald', 'kurucz']):
            # Atomic species, read from table
            self.mass = self.atoms_info.loc[self.species.capitalize(), 'mass']
        self.mass *= sc.amu # [kg]

    def _read_pressure_broadening_info(self):
        """
        Reads the pressure broadening parameters from the configuration.
        """

        perturber_info = getattr(self.config, 'perturber_info', None)
        if perturber_info is None:
            if self.database in ['hitran','exomol']:
                raise ValueError('Perturber information must be provided in the configuration file.')

            # For atomic lines, use polarisabilities (Eq. 23, Sharp & Burrows 2007)
            perturber_info = dict(
                H2 = dict(VMR=0.85, alpha=0.806e-30), # [m^3]
                He = dict(VMR=0.15, alpha=0.204956e-30),
            )

        # Check that the perturber information is complete
        self.pressure_broadening_info = {}
        for perturber, info in perturber_info.items():

            # Add mass to the perturber_info dictionary if it doesn't exist
            mass = info.get('mass', None)
            
            if perturber == 'H2' and mass is None:
                mass = sc.m_H2 / sc.amu
            elif perturber == 'He' and mass is None:
                mass = sc.m_He / sc.amu
            elif mass is None:
                raise ValueError(f'Mass of perturber {perturber} must be provided (in amu) in the configuration file.')
            info['mass'] = mass*sc.amu # [kg]

            if 'VMR' not in info:
                raise ValueError(f'Volume mixing ratio of perturber {perturber} must be provided in the configuration file.')

            # Store the broadening parameters
            self.pressure_broadening_info[perturber] = info.copy()

            # Use gamma and temperature-exponent from dictionary (gamma*(T0/T)^n*(P/1bar))
            gamma = info.get('gamma', 0.)
            self.pressure_broadening_info[perturber]['method'] = 'constant(s) from dictionary'
            if callable(gamma):
                # Callable parameterisation (e.g. Gharib-Nezhad et al. 2024)
                self.pressure_broadening_info[perturber]['gamma'] = gamma
                self.pressure_broadening_info[perturber]['method'] = 'function'
            else:
                self.pressure_broadening_info[perturber]['gamma'] = np.atleast_1d(gamma)
            self.pressure_broadening_info[perturber]['n'] = np.atleast_1d(info.get('n', 0.))

            # Load the broadening parameters from file, if given
            file = info.get('file', None)
            if file is None:
                continue
            file = pathlib.Path(file)
            broadening_params = read_fwf(file, header=None)
            diet = np.array(broadening_params[0], dtype=str)
            gamma, n = np.array(broadening_params[[1,2]].T)
            self.pressure_broadening_info[perturber]['method'] = f'\"{file}\"'

            # Currently only handles these 2 ('a0', 'm0') broadening diets
            mask_diet = (diet == 'a0')
            if not mask_diet.any():
                mask_diet = (diet == 'm0')
                self.pressure_broadening_info[perturber]['diet'] = 'm0'
            
            self.pressure_broadening_info[perturber]['gamma'] = gamma[mask_diet]
            self.pressure_broadening_info[perturber]['n'] = n[mask_diet]

            if (len(gamma)==1) or broadening_params.shape[1]==3:
                # One row or no quantum-number columns, ignore quantum-number dependence
                self.pressure_broadening_info[perturber]['gamma'] = np.nanmean(
                    self.pressure_broadening_info[perturber]['gamma'], keepdims=True
                    )
                self.pressure_broadening_info[perturber]['n'] = np.nanmean(
                    self.pressure_broadening_info[perturber]['n'], keepdims=True
                    )
                continue

            # Total angular momentum quantum numbers
            self.pressure_broadening_info[perturber]['J'] = \
                np.array(broadening_params[3])[mask_diet]

        print('\nPressure broadening info:')
        self.mean_mass, VMR_total = 0., 0.
        for perturber, info in self.pressure_broadening_info.items():
            VMR, mass, method = info['VMR'], info['mass'], info['method']
            print(f'  - {perturber}: VMR={VMR:.2f}, mass={mass/sc.amu:.2f} amu | {method}')
            self.mean_mass += info['VMR']*info['mass'] # [kg]
            VMR_total += info['VMR']

            # Convert gamma to SI units [cm^-1] -> [s^-1]
            if method != 'function':
                self.pressure_broadening_info[perturber]['gamma'] = info['gamma'] * 1e2*sc.c

        if VMR_total > 1.0:
            raise ValueError('Total volume mixing ratio of perturbers exceeds 1.0.')
        if VMR_total < 1.0:
            utils.warnings.warn('Total volume mixing ratio of perturbers is less than 1.0.')

        print(f'  Mean molecular weight of perturbers: {self.mean_mass/sc.amu:.2f} amu')

    def _read_partition_function(self):
        """
        Reads the partition function from the configuration file.
        """
        file = self.config.files.get('partition_function', None)
        if file is None:
            raise ValueError('Partition function file must be provided in the configuration file.')
        partition_function = np.loadtxt(file)

        # Make interpolation function, extrapolate outside temperature-range
        interpolation_function = interp1d(
            x=partition_function[:,0], y=partition_function[:,1], 
            kind='linear', fill_value='extrapolate'
            )
        self.calculate_partition_function = interpolation_function

    def iterate_over_PT_grid(self, function, progress_bar=True, **kwargs):
        """
        Iterates over the pressure-temperature grid and apply a function.

        Args:
            function (callable): Function to apply at each grid point.
            progress_bar (bool): Whether to show a progress bar.
            **kwargs: Additional arguments for the function.
        """
        # Make a nice progress bar
        pbar_kwargs = dict(
            total=self.N_PT, disable=(not progress_bar), 
            bar_format='{l_bar}{bar:20}{r_bar}{bar:-20b}', 
        )
        with tqdm(**pbar_kwargs) as pbar:

            # Loop over all PT-points
            for idx_P, P in enumerate(self.P_grid):
                for idx_T, T in enumerate(self.T_grid):

                    function(P, T, **kwargs)
                    postfix = {
                        'P': f'{P*1e-5:.0e} bar', 'T': f'{T:.0f} K',
                    }
                    if hasattr(self, 'N_lines_computed'):
                        postfix['N_lines_computed'] = '{}'.format(self.N_lines_computed)
                    pbar.set_postfix(**postfix, refresh=False)
                    pbar.update(1)
    
    def calculate_cross_sections(self, P, T, nu_0, S_0, E_low, A, delta=None, **kwargs):
        """
        Computes the cross-sections for given parameters.

        Args:
            P (float): Pressure in Pa.
            T (float): Temperature in Kelvin.
            nu_0 (array): Transition frequencies in s^-1.
            S_0 (array): Line strengths at reference temperature.
            E_low (array): Lower state energies in Joules.
            A (array): Einstein A-coefficients in s^-1.
            delta (array, optional): Pressure shift coefficients.
            **kwargs: Additional arguments.
        """
        self.N_lines_computed = 0

        # Get the line-widths
        gamma_N   = self.compute_natural_broadening(A) # Lorentzian components
        gamma_vdW = self.compute_vdw_broadening(P, T, E_low=E_low, nu_0=nu_0)
        gamma_L   = self.compute_lorentz_width(gamma_vdW, gamma_N)
        gamma_G   = self.compute_doppler_broadening(T, nu_0) # Gaussian component

        for nu_0_i in self.nu_0_to_ignore:
            # Ignore lines with a wavenumber close to the one to ignore
            mask = np.isclose(nu_0, nu_0_i, atol=1e-8)
            nu_0[mask] = np.nan

        # Apply the pressure shift
        nu_0 = self.pressure_shift(P, T, nu_0, delta=delta)

        # Select only the lines within the wavelength range
        nu_0, S_0, E_low, gamma_L, gamma_G = self._mask_arrays(
            [nu_0, S_0, E_low, gamma_L, gamma_G], mask=(nu_0>self.nu_min) & (nu_0<self.nu_max)
            )
        if len(S_0) == 0:
            return # No more lines

        # Get the line-strengths
        S = self.compute_line_strength(T, S_0, E_low, nu_0) # [s^-1/(molecule m^-2)]

        # Apply local and global line-strength cutoffs
        S = self.apply_local_cutoff(S, nu_0, factor=self.local_cutoff)
        S_min = max([0., self.global_cutoff])
        nu_0, S, gamma_L, gamma_G = self._mask_arrays(
            [nu_0, S, gamma_L, gamma_G], mask=(S>S_min)
        )
        if len(S) == 0:
            return # No more lines
        self.N_lines_computed = len(S)
        
        # Change to a coarse grid if lines are substantially broadened
        gamma_V = self.compute_voigt_width(gamma_G, gamma_L) # Voigt width
        nu_grid_to_use = self._setup_coarse_nu_grid(adaptive_delta_nu=np.mean(gamma_V)/6)
        delta_nu_to_use = np.diff(nu_grid_to_use[:2]) # [s^-1]
        
        # Wing cutoff from given lambda-function
        wing_cutoff_distance = self.wing_cutoff(np.mean(gamma_V), P) # [s^-1]
        wing_cutoff_distance = np.minimum(wing_cutoff_distance, self.wing_cutoff_max) # [s^-1]

        # Account for the lost line-strength
        S = self.normalise_wing_cutoff(S, wing_cutoff_distance, gamma_L)

        # Compute the line-profiles in chunks
        sigma = self.calculate_line_profiles(
            nu_0=nu_0,
            S=S,
            gamma_L=gamma_L,
            gamma_G=gamma_G,
            nu_grid=nu_grid_to_use,
            wing_cutoff_distance=wing_cutoff_distance,
            delta_nu=delta_nu_to_use,
        )

        if len(nu_grid_to_use) != len(self.nu_grid):
            # Interpolate to the original grid
            sigma = np.interp(self.nu_grid, nu_grid_to_use, sigma)

        # Add to the total array
        idx_P = np.searchsorted(self.P_grid, P)
        idx_T = np.searchsorted(self.T_grid, T)
        self.sigma[:,idx_P,idx_T] += sigma

    def calculate_temporary_outputs(self, overwrite=False, save_in_one_file=False, files_range=None, **kwargs):
        """
        Calculates temporary outputs for cross-sections.

        Args:
            overwrite (bool): Whether to overwrite existing files.
            save_in_one_file (bool): Whether to save all outputs in one file.
            files_range (tuple, optional): Range of files to process.
            **kwargs: Additional arguments.
        """
        print('\nCalculating cross-sections')

        transitions_files = self.config.files.get('transitions', None)
        if transitions_files is None:
            raise ValueError('No transitions files specified in the configuration.')
        transitions_files = np.atleast_1d(transitions_files)

        if files_range is not None:
            transitions_files = transitions_files[files_range[0]:files_range[1]]

        # Check if the output files already exist
        transitions_files, tmp_output_files = self._check_existing_output_files(
            transitions_files, overwrite_all=overwrite
            )

        self.sigma = np.zeros(
            (len(self.nu_grid), len(self.P_grid), len(self.T_grid)), dtype=float
            )
        for input_file, tmp_output_file in zip(transitions_files, tmp_output_files):
            
            # Compute the cross-sections
            self.process_transitions(input_file, **kwargs)

            if np.all(self.sigma == 0.):
                print(f'  No lines calculated, no need to save \"{tmp_output_file}\"')
                continue

            if save_in_one_file and (tmp_output_file != tmp_output_files[-1]):
                continue

            print(f'  Saving temporary cross-sections to \"{tmp_output_file}\"')
            self.sigma[(self.sigma == 0.)] = 1e-250

            # Temporarily save the data
            utils.save_to_hdf5(
                tmp_output_file, 
                data={
                    'wave': self.wave_grid,
                    'P': self.P_grid,
                    'T': self.T_grid,
                    'log10(xsec)': utils.log10_round(self.sigma, decimals=3), # Round to save memory
                },
                attrs={
                    'wave': {'units': 'm'},
                    'P': {'units': 'Pa'},
                    'T': {'units': 'K'},
                    'log10(xsec)': {'units': 'log10(m^2 molecule^-1)'},
                }
            )

    def save_combined_outputs(self, **kwargs):
        """
        Saves the merged outputs to a file.

        Args:
            **kwargs: Additional arguments for saving.
        """
        super().save_combined_outputs(keys_to_merge=['log10(xsec)'], **kwargs)

    def plot_combined_outputs(self, fig=None, ax=None, return_fig_ax=False, cmaps=['coolwarm','viridis'], xscale='log', yscale='log', xlim=None, ylim=None, ls='-', lw=0.7, **kwargs):
        """
        Plots the merged outputs.

        Args:
            fig (matplotlib.figure.Figure): Figure object for plotting.
            ax (matplotlib.axes.Axes): Axes object for plotting.
            return_fig_ax (bool): Whether to return the figure and axes.
            cmaps (list): List of colormaps for plotting.
            xscale (str): Scale for the x-axis.
            yscale (str): Scale for the y-axis.
            xlim (tuple, optional): Limits for the x-axis.
            ylim (tuple, optional): Limits for the y-axis.
            ls (str): Line style for plotting.
            **kwargs: Additional arguments for plotting.
        """
        
        import matplotlib.pyplot as plt
        if (fig is None) or (ax is None):
            # Make a new figure
            fig, ax = plt.subplots(figsize=(9,6), nrows=2, sharex=True, sharey=True)
            
            print('\nPlotting cross-sections')

        self.combined_datasets = utils.read_from_hdf5(
            self.final_output_file, keys_to_read=['wave', 'P', 'T', 'log10(xsec)']
        )
        wave = self.combined_datasets['wave'] * 1e6 # [m] -> [um]
        xsec = 10**self.combined_datasets['log10(xsec)'] * (1e2)**2
        P_grid = self.combined_datasets['P'] / sc.bar # [Pa] -> [bar]

        # Avoid plotting the whole dataset
        if xlim is None:
            xlim = (wave.min(), wave.max())
        xsec = xsec[(wave>=xlim[0]) & (wave<=xlim[1])]
        wave = wave[(wave>=xlim[0]) & (wave<=xlim[1])]

        # Plot for certain temperatures
        T_to_plot = kwargs.get('T_to_plot', self.combined_datasets['T'])
        indices_T, T_to_plot = utils.find_closest_indices(self.combined_datasets['T'], T_to_plot)
        indices_T = np.unique(indices_T)
        T_to_plot = np.unique(T_to_plot)

        idx_P = np.searchsorted(P_grid, 1) # 1 bar
        idx_P = np.minimum(idx_P, len(P_grid)-1)

        for idx_T, T in zip(indices_T, T_to_plot):
            idx_c = 0.3
            if len(T_to_plot)>1:
                idx_c = (T-T_to_plot.min())/(T_to_plot.max()-T_to_plot.min())

            cmap = cmaps[0]
            if isinstance(cmap, str):
                cmap = plt.get_cmap(cmap)

            ax[0].plot(wave, xsec[:,idx_P,idx_T], c=cmap(idx_c), lw=lw, ls=ls, label=f'T={T:.0f} K')

        handles, _ = ax[0].get_legend_handles_labels()
        ncols = 1 + len(handles)//8
        ax[0].legend(
            loc='upper right', ncol=ncols, labelcolor='linecolor', handlelength=0.5, 
            title=f'P={P_grid[idx_P]:.0e} bar',
            )

        # Plot for certain pressures
        P_to_plot = kwargs.get('P_to_plot', P_grid) # [bar]
        indices_P, P_to_plot = utils.find_closest_indices(P_grid, P_to_plot)
        indices_P = np.unique(indices_P)
        P_to_plot = np.unique(P_to_plot)

        idx_T = np.searchsorted(self.combined_datasets['T'], 1000.)
        idx_T = np.minimum(idx_T, len(self.combined_datasets['T'])-1)

        for idx_P, P in zip(indices_P, P_to_plot):   
            idx_c = 0.5
            if len(P_to_plot)>1:
                idx_c = np.log10(P/P_to_plot.min())/np.log10(P_to_plot.max()/P_to_plot.min())
            
            cmap = cmaps[1]
            if isinstance(cmap, str):
                cmap = plt.get_cmap(cmap)

            ax[1].plot(wave, xsec[:,idx_P,idx_T], c=cmap(idx_c), lw=lw, ls=ls, label=f'P={P:.0e} bar')

        handles, _ = ax[1].get_legend_handles_labels()
        ncols = 1 + len(handles)//8
        ax[1].legend(
            loc='upper right', ncol=ncols, labelcolor='linecolor', handlelength=0.5, 
            title=f'T={self.combined_datasets["T"][idx_T]:.0f} K',
            )

        if ylim is None:
            xsec[xsec<=1e-150] = np.nan
            ylim = (np.nanpercentile(xsec, 3), np.nanmax(xsec)*10)

        ax[0].set(xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim, ylabel=r'xsec [cm$^2$ molecule$^{-1}$]')
        ax[1].set(xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim, xlabel=r'Wavelength [$\mu$m]', ylabel=r'xsec [cm$^2$ molecule$^{-1}$]')

        plt.savefig(self.output_data_dir / 'xsec.pdf', bbox_inches='tight')
        if return_fig_ax:
            return fig, ax
        plt.close()

    def convert_to_pRT3(self, contributor=None, **kwargs):
        """
        Converts the cross-sections to petitRADTRANS v3.0 format.

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
        for key in ['DOI', 'mol_mass', 'mol_name', 'isotopologue_id']:
            if key in pRT3_metadata:
                continue
            if key == 'mol_mass':
                # If mass is not given (for atoms)
                pRT3_metadata[key] = self.mass / sc.amu
                continue
            # Key not given
            raise KeyError(f"Required key '{key}' not found in pRT3_metadata.")

        data = {
            'DOI': np.atleast_1d(pRT3_metadata['DOI']),
            'Date_ID': np.atleast_1d(f'petitRADTRANS-v3_{datetime.datetime.now(datetime.timezone.utc).isoformat()}'),
            'mol_mass': np.atleast_1d(pRT3_metadata['mol_mass']),
            'mol_name': np.atleast_1d(pRT3_metadata['mol_name']), 
            'temperature_grid_type': np.atleast_1d('regular'),
        }

        # Attributes of all datasets
        attrs = dict(
            DOI = {'additional_description': 'Calculated with pyROX', 'long_name': 'Data object identifier linked to the data'},
            Date_ID = {'long_name': 'ISO 8601 UTC time (https://docs.python.org/3/library/datetime.html) at which the table has been created, along with the version of petitRADTRANS'},
            xsecarr = {'long_name': 'Table of the cross-sections with axes (pressure, temperature, wavenumber)', 'units': 'cm^2/molecule'},
            mol_mass = {'long_name': 'Mass of the species', 'units': 'AMU'},
            mol_name = {'long_name': 'Name of the species described'},
            p = {'long_name': 'Pressure grid', 'units': 'bar'},
            t = {'long_name': 'Temperature grid', 'units': 'K'},
            temperature_grid_type = {'long_name': 'Whether the temperature grid is "regular" (same temperatures for all pressures) or "pressure-dependent"'},
            bin_edges = {'long_name': 'Wavenumber grid', 'units': 'cm^-1'},
            wlrange = {'long_name': 'Wavelength range covered', 'units': 'µm'},
            wnrange = {'long_name': 'Wavenumber range covered', 'units': 'cm^-1'}
        )
        # Add contributors if given
        if contributor is not None:
            attrs['DOI']['contributor'] = contributor

        # Load the merged data
        self.combined_datasets = utils.read_from_hdf5(
            self.final_output_file, keys_to_read=['wave', 'P', 'T', 'log10(xsec)']
        )
        xsec = 10**self.combined_datasets['log10(xsec)'] * (1e2)**2 # [m^2 molecule^-1] -> [cm^2 molecule^-1]
        xsec = np.moveaxis(xsec, 0, -1) # (wave, P, T) -> (P, T, wave)

        # Interpolate onto pRT's wavelength grid
        pRT_wave = utils.fixed_resolution_wavelength_grid(0.11, 250, resolution=1e6)
        pRT_wave *= sc.micron # [um] -> [m]
        interp_func = interp1d(
            x=self.combined_datasets['wave'], y=np.log10(xsec), kind='linear', 
            bounds_error=False, fill_value=-250, axis=-1
        )
        xsec = 10**interp_func(pRT_wave) # [cm^2 molecule^-1]

        wave_min = max(self.combined_datasets['wave'].min(), pRT_wave.min()) / sc.micron # [m] -> [um]
        wave_max = min(self.combined_datasets['wave'].max(), pRT_wave.max()) / sc.micron
        mask_wave = (pRT_wave>=wave_min*sc.micron) & (pRT_wave<=wave_max*sc.micron)
        xsec = xsec[:,:,mask_wave] # Crop to the range of valid wavelengths
        pRT_wave = pRT_wave[mask_wave]

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
        data['xsecarr'] = xsec[:,:,::-1] # Ascending in wavenumber
        data['p'] = self.combined_datasets['P'] / sc.bar  # [Pa] -> [bar]
        data['t'] = self.combined_datasets['T']  # [K]
        data['bin_edges'] = 1e-2 / pRT_wave[::-1]  # [m] -> [cm^-1], ascending in wavenumber
        data['wlrange'] = [wave_min, wave_max]  # [µm]
        data['wnrange'] = [1e4 / wave_max, 1e4 / wave_min]  # [cm^-1]

        # Complete the filename
        isotopologue_id = pRT3_metadata['isotopologue_id']
        if isinstance(isotopologue_id, dict):
            isotopologue_id = '-'.join([
                f'{mass_number}{element}' for element, mass_number in isotopologue_id.items()
                ])
        
        linelist = pRT3_metadata.get('linelist', self.database.capitalize())
        if linelist in ['Hitran', 'Hitemp']:
            linelist = linelist.upper()

        pRT_file = '{}__{}.R{:.0e}_{:.1f}-{:.1f}mu.xsec.petitRADTRANS.h5'
        pRT_file = pRT_file.format(
            isotopologue_id, linelist, resolution, wave_min, wave_max
        )
        pRT_file = self.output_data_dir / pRT_file

        # Save the datasets
        utils.save_to_hdf5(pRT_file, data=data, attrs=attrs, compression=None)
        print(f'  Saved to \"{pRT_file}\"')