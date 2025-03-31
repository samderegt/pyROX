import numpy as np
from pandas import read_fwf, read_csv
from scipy.interpolate import interp1d
from scipy.special import wofz

import pathlib
import warnings
import re

from tqdm import tqdm

from cross_sections import CrossSections
from utils import sc
import utils

class LineProfileHelper:

    def line_strength(self, T, S_0, E_low, nu_0):

        # Partition function
        q = self.calculate_partition_function(T)

        # Gordon et al. (2017) (E_low: [J]; nu_0: [s^-1])
        S = (
            S_0 * (self.q_0/q) * np.exp(E_low/sc.k*(1/self.T_0-1/T)) *
            (1-np.exp(-sc.h*nu_0/(sc.k*T))) / (1-np.exp(-sc.h*nu_0/(sc.k*self.T_0)))
        )
        return S

    def normalise_wing_cutoff(self, S, cutoff_distance, gamma_L):
        # Eq. A6 Lacy & Burrows (2023)
        return S / ((2/np.pi)*np.arctan(cutoff_distance/gamma_L))

    def gamma_vdW(self, P, T):

        # TODO: Implement the atomic broadening
        # ...

        # Van der Waals broadening
        gamma_vdW = 0
        for perturber, info in self.pressure_broadening_info.items():
            # Get the broadening parameters
            VMR = info['VMR']
            gamma = info['gamma'] # [s^-1]
            n = info['n']

            # Calculate the broadening parameter
            gamma_vdW += gamma * (T/self.T_0)**n * (P/sc.atm) # [s^-1]

        return gamma_vdW

    def gamma_N(self, nu_0, log_gamma_N=None):

        # Natural broadening
        gamma_N = 0.222 * (nu_0/(1e2*sc.c))**2 / (4*np.pi*sc.c) # [s^-1]

        if log_gamma_N is not None:
            # Natural damping constant is given (for atoms)
            mask_valid = (log_gamma_N != 0)
            gamma_N[mask_valid] = 10**log_gamma_N[mask_valid] / (4*np.pi) # [s^-1]

        return gamma_N

    def gamma_G(self, T, nu_0):
        # Doppler broadening
        return np.sqrt(2*sc.k*T/self.mass) * nu_0/sc.c # [s^-1]

    def gamma_L(self, gamma_vdW, gamma_N):
        # Lorentz width
        return gamma_vdW + gamma_N # [s^-1]

    def gamma_V(self, gamma_G, gamma_L):
        # Voigt width
        return 0.5346*gamma_L+np.sqrt(0.2166*gamma_L**2 + gamma_G**2) # [s^-1]

    def apply_local_cutoff(self, S, nu_0, factor):

        # Round to zero-th decimal
        nu_bin = np.around((nu_0-self.nu_min)/self.delta_nu_local_cutoff).astype(int)

        # Upper and lower indices of lines within bins
        _, nu_bin_idx = np.unique(nu_bin, return_index=True)
        nu_bin_idx    = np.append(nu_bin_idx, len(nu_bin))

        for k in range(len(nu_bin_idx)-1):
            # Cumulative sum of lines in bin
            S_range = S[nu_bin_idx[k]:nu_bin_idx[k+1]]
            if len(S_range) < 2:
                continue # One/no lines in bin
            idx_sort = np.argsort(S_range)
            S_sort   = S_range[idx_sort]
            S_cumsum = np.cumsum(S_sort)

            # Lines contributing less than 'factor' to total strength
            i_search = np.searchsorted(S_cumsum, factor*S_cumsum[-1])
            S_cutoff = S_sort[i_search]

            # Add weak line-strengths to strongest line
            S_sum_others = S_cumsum[i_search-1]
            S_range[idx_sort[-1]] += S_sum_others

            # Ignore weak lines
            S_range[S_range<S_cutoff] = 0.

            S[nu_bin_idx[k]:nu_bin_idx[k+1]] = S_range

        return S

    def calculate_line_profiles_in_chunks(
            self, 
            nu_0, 
            S, 
            gamma_L, 
            gamma_G, 
            nu_grid, 
            delta_nu, 
            wing_cutoff_distance, 
            N_lines_in_chunk=200, 
            ):

        # Indices where lines should be inserted
        idx_to_insert = np.searchsorted(nu_grid, nu_0) - 1

        # Wing-length in number of grid points
        wing_length = int(np.around(wing_cutoff_distance/delta_nu))

        # Array of wavenumbers from line-center
        nu_line = np.linspace(
            -wing_length*delta_nu, wing_length*delta_nu, 2*wing_length+1, endpoint=True
            )
        nu_line = nu_line[None,:]

        N_nu_grid = len(nu_grid)
        N_nu_line = nu_line.shape[1]

        # Relative width of Lorentzian vs. Gaussian
        a = gamma_L / gamma_G # Gandhi et al. (2020)

        # Only consider number of lines at a time
        N_chunks = int(np.ceil(len(S)/N_lines_in_chunk))

        sigma = np.zeros_like(nu_grid)
        for ch in range(N_chunks):
            
            # Upper and lower indices of lines in current chunk
            idx_chunk_l = int(ch*N_lines_in_chunk)
            idx_chunk_h = idx_chunk_l + N_lines_in_chunk
            idx_chunk_h = np.minimum(idx_chunk_h, len(S)) # At last chunk

            # Indices of nu_grid_coarse to insert current lines
            idx_to_insert_chunk = idx_to_insert[idx_chunk_l:idx_chunk_h]

            # Lines in current chunk | (N_lines,1)
            nu_0_chunk    = nu_0[idx_chunk_l:idx_chunk_h,None]
            gamma_G_chunk = gamma_G[idx_chunk_l:idx_chunk_h,None]
            a_chunk = a[idx_chunk_l:idx_chunk_h,None]
            S_chunk = S[idx_chunk_l:idx_chunk_h,None]

            # Correct for coarse grid | (N_lines,1)
            nu_grid_chunk = nu_grid[idx_to_insert_chunk,None]

            # Eq. 10 (Gandhi et al. 2020) | (N_lines,N_wave[cut])
            u = ((nu_line+nu_grid_chunk) - nu_0_chunk) / gamma_G_chunk

            # (Scaled) Faddeeva function for Voigt profiles | (N_lines,N_wave[cut])
            sigma_chunk = S_chunk * np.real(wofz(u+a_chunk*1j)) / (gamma_G_chunk*np.sqrt(np.pi))

            # Upper and lower index of these lines in sigma_coarse
            idx_sigma_l = np.maximum(0, idx_to_insert_chunk-wing_length)
            idx_sigma_h = np.minimum(N_nu_grid, idx_to_insert_chunk+wing_length+1)

            # Upper and lower index of these lines in sigma_ch
            idx_sigma_chunk_l = np.maximum(0, wing_length-idx_to_insert_chunk)
            idx_sigma_chunk_h = np.minimum(N_nu_line, wing_length-idx_to_insert_chunk+N_nu_grid)

            # Loop over each line profile
            for i, sigma_i in enumerate(sigma_chunk):
                # Add line to total cross-section
                sigma[idx_sigma_l[i]:idx_sigma_h[i]] += sigma_i[idx_sigma_chunk_l[i]:idx_sigma_chunk_h[i]]
                
        return sigma
    
class LineByLine(CrossSections, LineProfileHelper):
    """
    Base class for line-by-line cross-sections.
    """

    def __init__(self, config, **kwargs):

        # Initialise the CrossSections parent class
        super().__init__(config, **kwargs)

    def _read_from_config(self, config):
        
        # Read the common parameters
        super()._read_from_config(config)
        
        # Read the parameters for line-by-line calculations
        self._read_mass()
        self._read_pressure_broadening_info()
        self._read_equation_of_state()
        self._read_partition_function()

        self.N_lines_in_chunk = getattr(self.config, 'N_lines_in_chunk', 10_000_000)
        
        # Reference temperature and partition function
        self.T_0 = getattr(self.config, 'T_0', 296.)
        self.q_0 = self.calculate_partition_function(self.T_0)

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

        self.local_cutoff = getattr(self.config, 'local_cutoff', None)  # [fraction of cumulative]
        # Finer resolution to determine which lines to keep
        self.delta_nu_local_cutoff = getattr(self.config, 'delta_nu_local_cutoff', 1e-3) # [cm^-1]
        self.delta_nu_local_cutoff *= 1e2*sc.c # [cm^-1] -> [s^-1]

        # Decrease the wavenumber-grid resolution for high pressures
        #self._set_nu_grid(self.config)
        self.adaptive_nu_grid = getattr(self.config, 'adaptive_nu_grid', False)

        # (P,T)-grid to calculate cross-sections on
        self.P_grid = np.atleast_1d(self.config.P_grid) * 1e5 # [bar] -> [Pa]
        self.T_grid = np.atleast_1d(self.config.T_grid)
        self.N_PT = len(self.P_grid) * len(self.T_grid)

        print('\nPT-grid:')
        print(f'  P: {self.P_grid/1e5} bar')
        print(f'  T: {self.T_grid} K')

    def _read_equation_of_state(self):
        """
        Read the equation of state.
        """
        EOS_table = getattr(self.config, 'EOS_table', None)
        if EOS_table is None:
            # Assume ideal gas
            self.get_number_density = lambda P, T: P/(sc.k*T)
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
        self.get_number_density = lambda P, T: interp_func([P, T])

    def _read_mass(self):
        """
        Read the mass of the species.
        """
        self.mass = getattr(self.config, 'mass', None)
        if (self.mass is None) and (self.database in ['hitran', 'exomol']):
            raise ValueError('Mass of the species must be provided in the configuration file.')
        elif (self.mass is None) and (self.database in ['vald', 'kurucz']):
            # Atomic species, read from table
            self.mass = self.atoms_info.loc[self.species, 'mass']
        self.mass *= sc.amu # [kg]

    def _read_pressure_broadening_info(self):
        """
        Read the pressure broadening parameters from the configuration file.
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

        for perturber, info in self.pressure_broadening_info.items():
            self.pressure_broadening_info[perturber]['gamma'] = \
                np.array(info['gamma'], dtype=float)

        print(f'  Pressure broadening info:')
        self.mean_mass, VMR_total = 0., 0.
        for perturber, info in self.pressure_broadening_info.items():
            VMR, mass, method = info['VMR'], info['mass'], info['method']
            print(f'    - {perturber}: VMR={VMR:.2f}, mass={mass/sc.amu:.2f} amu | {method}')
            self.mean_mass += info['VMR']*info['mass'] # [kg]
            VMR_total += info['VMR']

            # Convert gamma to SI units [cm^-1] -> [s^-1]
            if method != 'function':
                self.pressure_broadening_info[perturber]['gamma'] = info['gamma'] * 1e2*sc.c
        
        if VMR_total > 1.0:
            raise ValueError('Total volume mixing ratio of perturbers exceeds 1.0.')
        if VMR_total < 1.0:
            warnings.warn('Total volume mixing ratio of perturbers is less than 1.0.')

        print(f'  Mean molecular weight of perturbers: {self.mean_mass/sc.amu:.2f} amu')

    def _read_partition_function(self):
        """
        Read the partition function from the configuration file.
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

    def _check_if_output_exists(self, input_files, overwrite_all=False):
        """
        Check if the output files already exist.
        """
        output_files = []
        for i, input_file in enumerate(input_files):
            # Check if the transition file exists
            input_file = pathlib.Path(input_file)
            
            # Check if the temporary output file already exists
            output_file = pathlib.Path(
                str(self.tmp_output_basename).replace('.hdf5', f'_{input_file.stem}.hdf5')
            )

            if output_file.exists() and not overwrite_all:
                # Ask the user if they want to overwrite the file
                response = ''
                while response not in ['y', 'yes', 'n', 'no', 'all']:
                    response = input(f'  Warning: Temporary output file \"{output_file}\" already exists. Overwrite? (yes/no/all): ')
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
        return output_files
        
    def loop_over_PT_grid(self, function, show_progress_bar=True, **kwargs):
        
        # Make a nice progress bar
        pbar_kwargs = dict(
            total=self.N_PT, disable=(not show_progress_bar), 
            bar_format='{l_bar}{bar:20}{r_bar}{bar:-20b}', 
        )
        with tqdm(**pbar_kwargs) as pbar:

            # Loop over all PT-points
            for idx_P, P in enumerate(self.P_grid):
                for idx_T, T in enumerate(self.T_grid):

                    pbar.set_postfix(P='{:.0e} bar'.format(P*1e-5), T='{:.0f} K'.format(T), refresh=False)
                    function(P, T, **kwargs)
                    pbar.update(1)
    
    def calculate_cross_sections(self, P, T, nu_0, S_0, E_low, delta_P=0., **kwargs):
        """
        Compute the cross-sections.
        """
        
        # Get the line-widths
        gamma_N   = self.gamma_N(nu_0) # Lorentzian components
        gamma_vdW = self.gamma_vdW(P, T)
        gamma_L = self.gamma_L(gamma_vdW, gamma_N)
        gamma_G = self.gamma_G(T, nu_0) # Gaussian component

        # Select only the lines within the wavelength range
        nu_0, S_0, E_low, gamma_N, gamma_vdW, gamma_L, gamma_G = self.mask_arrays(
            [nu_0, S_0, E_low, gamma_N, gamma_vdW, gamma_L, gamma_G], 
            mask=(nu_0>self.nu_min) & (nu_0<self.nu_max)
            )
        if len(S_0) == 0:
            return # No more lines

        # Get the line-strengths
        S = self.line_strength(T, S_0, E_low, nu_0)

        # Apply local and global line-strength cutoffs
        S = self.apply_local_cutoff(S, nu_0, self.local_cutoff)
        S_min = max([0., self.global_cutoff])
        nu_0, S, gamma_L, gamma_G = self.mask_arrays(
            [nu_0, S, gamma_L, gamma_G], mask=(S>S_min)
        )
        if len(S) == 0:
            return # No more lines
        
        # Change to a coarse grid if lines are substantially broadened
        gamma_V = self.gamma_V(gamma_G, gamma_L) # Voigt width
        nu_grid_to_use, delta_nu_to_use = \
            self._configure_coarse_nu_grid(adaptive_delta_nu=np.mean(gamma_V)/6)

        # Indices where lines should be inserted
        idx_to_insert = np.searchsorted(nu_grid_to_use, nu_0) - 1
        
        # Wing cutoff from given lambda-function
        wing_cutoff_distance = self.wing_cutoff(np.mean(gamma_V), P) # [s^-1]
        wing_cutoff_distance = np.minimum(wing_cutoff_distance, self.wing_cutoff_max) # [s^-1]

        # Account for the lost line-strength
        S = self.normalise_wing_cutoff(S, wing_cutoff_distance, gamma_L)

        # Compute the line-profiles in chunks
        sigma = self.calculate_line_profiles_in_chunks(
            nu_0=nu_0, 
            S=S, 
            gamma_L=gamma_L, 
            gamma_G=gamma_G, 
            nu_grid=nu_grid_to_use, 
            delta_nu=delta_nu_to_use, 
            wing_cutoff_distance=wing_cutoff_distance, 
        )
        if len(nu_grid_to_use) != len(self.nu_grid):
            # Interpolate to the original grid
            sigma = np.interp(self.nu_grid, nu_grid_to_use, sigma)

        # Add to the total array
        idx_P = np.searchsorted(self.P_grid, P)
        idx_T = np.searchsorted(self.T_grid, T)
        self.sigma[:,idx_P,idx_T] += sigma

    def calculate_tmp_outputs(self, overwrite=False, **kwargs):
        print('\nCalculating cross-sections')

        transitions_files = self.config.files.get('transitions', None)
        if transitions_files is None:
            raise ValueError('No transitions files specified in the configuration.')
        transitions_files = np.atleast_1d(transitions_files)

        # Check if the output files already exist
        tmp_output_files = self._check_if_output_exists(
            transitions_files, overwrite_all=overwrite
            )

        for input_file, tmp_output_file in zip(transitions_files, tmp_output_files):
            
            # Compute the cross-sections
            self.sigma = np.zeros(
                (len(self.nu_grid), len(self.P_grid), len(self.T_grid)), dtype=float
                )
            self._read_transitions_in_chunks(input_file, tmp_output_file, **kwargs)

            if np.all(self.sigma == 0.):
                continue # No lines in this file, no need to save

            self.sigma[(self.sigma == 0.)] = 1e-250

            # Temporarily save the data
            utils.save_to_hdf5(
                tmp_output_file, 
                data={
                    'wave': self.wave_grid,
                    'P': self.P_grid,
                    'T': self.T_grid,
                    'xsec': self.sigma, 
                },
                attrs={
                    'wave': {'units': 'm'},
                    'P': {'units': 'Pa'},
                    'T': {'units': 'K'},
                    'xsec': {'units': 'm^2 molecule^-1'},
                }
            )

    def save_merged_outputs(self, **kwargs):
        super().save_merged_outputs(keys_to_merge=['xsec'], **kwargs)

    def plot_merged_outputs(self, cmaps=['coolwarm','viridis'], xscale='log', yscale='log', xlim=None, ylim=None, **kwargs):
        """
        Plot the merged outputs. Same for all LineByLine classes.
        """
        
        print('\nPlotting cross-sections')

        self.merged_datasets = utils.read_from_hdf5(
            self.final_output_file, keys_to_read=['wave', 'P', 'T', 'xsec']
        )
        #self.merged_datasets['xsec'] *= 1e4 # [m^2] -> [cm^2]

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9,6), nrows=2, sharex=True, sharey=True)

        # Plot for certain temperatures
        T_to_plot = kwargs.get('T_to_plot', self.merged_datasets['T'])
        indices_T, T_to_plot = utils.find_nearest(self.merged_datasets['T'], T_to_plot)
        indices_T = np.unique(indices_T)
        T_to_plot = np.unique(T_to_plot)

        #idx_P = np.searchsorted(self.merged_datasets['P'], 1e5) # 1 bar
        idx_P = np.searchsorted(self.merged_datasets['P'], 1e4) # 1 bar

        for idx_T, T in zip(indices_T, T_to_plot):
            if len(T_to_plot)==1:
                c = plt.get_cmap(cmaps[0])(0.4)
            else:
                c = plt.get_cmap(cmaps[0])((T-T_to_plot.min())/(T_to_plot.max()-T_to_plot.min()))

            xsec = self.merged_datasets['xsec'][:,idx_P,idx_T]
            ax[0].plot(self.merged_datasets['wave'], xsec, c=c, lw=1, label=f'T={T:.0f} K')

        # Plot for certain pressures
        P_to_plot = kwargs.get('P_to_plot', self.merged_datasets['P'])
        indices_P, P_to_plot = utils.find_nearest(self.merged_datasets['P'], P_to_plot)
        indices_P = np.unique(indices_P)
        P_to_plot = np.unique(P_to_plot)

        idx_T = np.searchsorted(self.merged_datasets['T'], 1000.)

        for idx_P, P in zip(indices_P, P_to_plot):
            if len(P_to_plot)==1:
                c = plt.get_cmap(cmaps[1])(0.5)
            else:
                c = plt.get_cmap(cmaps[1])(np.log10(P/P_to_plot.min())/np.log10(P_to_plot.max()/P_to_plot.min()))

            xsec = self.merged_datasets['xsec'][:,idx_P,idx_T]
            ax[1].plot(self.merged_datasets['wave'], xsec, c=c, lw=1, label=f'P={P/sc.bar:.0e} bar')

        if ylim is None:
            xsec = self.merged_datasets['xsec'][self.merged_datasets['xsec']>1e-250]
            ylim = (np.nanpercentile(xsec, 3), np.nanmax(xsec)*10)

        ax[0].set(xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim, ylabel='xsec [m^2 molecule^-1]')
        ax[1].set(xscale=xscale, yscale=yscale, xlim=xlim, ylim=ylim, xlabel='wavelength [m]', ylabel='xsec [m^2 molecule^-1]')
        
        handles, _ = ax[0].get_legend_handles_labels()
        ncols = 1 + len(handles)//8
        ax[0].legend(loc='upper right', ncol=ncols, labelcolor='linecolor')
        ax[1].legend(loc='upper right', ncol=ncols, labelcolor='linecolor')

        plt.savefig(self.output_data_dir / 'xsec.pdf', bbox_inches='tight')
        plt.close()

    def convert_to_pRT2(self):
        raise NotImplementedError('Conversion to petitRADTRANS-v2.0 format not implemented.')

    def convert_to_pRT3(self):
        raise NotImplementedError('Conversion to petitRADTRANS-v3.0 format not implemented.')

class HITRAN(LineByLine):

    def download_data(self, config):
        """
        Download data from HITRAN.
        """

        print('\nDownloading data from HITRAN')

        files = []
        for url in config.urls:
            file = utils.wget_if_not_exist(url, config.input_data_dir)
            files.append(file)

        if None in files:
            raise ValueError('Failed to download all urls.')
    
    def __init__(self, config, **kwargs):

        print('-'*60)
        print('  Line-by-line Absorption from HITRAN/HITEMP')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent LineByLine class

    def _read_from_config(self, config):
        
        # Read the common parameters
        super()._read_from_config(config)
        
        # Read the isotope information
        self.isotope_idx       = getattr(config, 'isotope_idx', 1)
        self.isotope_abundance = getattr(config, 'isotope_abundance', 1.0)
        print(f'\nIsotope index: {self.isotope_idx}, with abundance {self.isotope_abundance:.2e}')

        # Remove any quantum-number dependency
        for perturber, info in self.pressure_broadening_info.items():
            self.pressure_broadening_info[perturber]['gamma'] = np.nanmean(info['gamma'])
            self.pressure_broadening_info[perturber]['n']     = np.nanmean(info['n'])
            # Remove 'diet' and 'J' keys if they exist
            self.pressure_broadening_info[perturber].pop('diet', None)
            self.pressure_broadening_info[perturber].pop('J', None)

    def _read_transitions_in_chunks(self, input_file, tmp_output_file, **kwargs):
        """
        Compute the cross-sections.
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
            
            # Sort lines by wavenumber
            transitions = transitions[np.argsort(transitions[:,2])]

            # Unit conversion
            nu_0  = transitions[:,2].astype(float) * 1e2*sc.c        # [cm^-1] -> [s^-1]
            E_low = transitions[:,7].astype(float) * sc.h*(1e2*sc.c) # [cm^-1] -> [J]

            # [cm^-1/(molec. cm^-2)] -> [s^-1/(molec. m^-2)]
            S_0 = transitions[:,3].astype(float) * (1e2*sc.c) * 1e-4
            S_0 /= self.isotope_abundance # Remove terrestrial abundance ratio

            idx = np.argsort(nu_0)
            nu_0 = nu_0[idx]
            S_0 = S_0[idx]
            E_l = E_l[idx]
            
            # Compute the cross-sections, looping over the PT-grid
            print(f'  Number of lines: {len(nu_0)}')
            self.loop_over_PT_grid(
                function=self.calculate_cross_sections,
                nu_0=nu_0, S_0=S_0, E_low=E_low, 
                **kwargs
            )

class ExoMol(LineByLine):

    def download_data(self, config):
        """
        Download data from ExoMol.
        """

        print('\nDownloading data from ExoMol')

        files = []
        file_def_json = None
        for url in config.urls:
            file = utils.wget_if_not_exist(url, config.input_data_dir)
            files.append(file)

            if url.endswith('.def.json'):
                # Check if the url is a definition file
                url_base = url.replace('.def.json','')
                file_def_json = file
        
        # Download partition-function and states files
        file_partition = utils.wget_if_not_exist(f'{url_base}.pf', config.input_data_dir)
        file_states = utils.wget_if_not_exist(f'{url_base}.states.bz2', config.input_data_dir)

        # Read definition file
        import json
        with open(file_def_json, 'r') as f:
            file_def_dict = json.load(f)

        # Download transitions files
        file_transitions = [utils.wget_if_not_exist(f'{url_base}.trans.bz2', config.input_data_dir)]
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
                file_transitions_i = utils.wget_if_not_exist(
                    f'{url_base}__{nu_min_i:05d}-{nu_max_i:05d}.trans.bz2', config.input_data_dir
                )
                file_transitions.append(file_transitions_i)

        if None in [*files, file_partition, file_states, *file_transitions]:
            raise ValueError('Failed to download all urls.')
    
    def __init__(self, config, **kwargs):

        print('-'*60)
        print('  Line-by-line Absorption from ExoMol')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent LineByLine class

    def _read_from_config(self, config):

        # Read the common parameters
        super()._read_from_config(config)

        # Copy so that we can modify the broadening parameters per transition
        self.pressure_broadening_info_copy = self.pressure_broadening_info.copy()

    def _read_broadening_per_transition(self, J_l, J_u):

        for perturber, info in self.pressure_broadening_info_copy.items():
            # Get the broadening parameters
            gamma = info['gamma']
            n     = info['n']
            #label = info.get('label', None)

            if callable(gamma):
                # User-provided function
                self.pressure_broadening_info[perturber]['gamma'] = gamma(J_l)
            else: 
                self.pressure_broadening_info[perturber]['gamma'] = \
                    np.nanmean(gamma)*np.ones_like(J_l)

            if callable(n):
                # User-provided function
                self.pressure_broadening_info[perturber]['n'] = n(J_l)
            else:
                self.pressure_broadening_info[perturber]['n'] = \
                    np.nanmean(n)*np.ones_like(J_l)

            if callable(gamma) or callable(n):
                continue

            J_in_table = info.get('J')
            label = info.get('label')
            if (J_in_table is None) or (label not in ['a0', 'm0']):
                # No quantum-number dependence
                continue
            
            # Read in chunks to avoid memory overload
            for idx_l in range(0, len(J_l), 100_000):
                idx_h = min([idx_l+100_000, len(J_l)])

                J_l_to_match = J_l[idx_l:idx_h]
                
                if label == 'm0':
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

    def _read_transitions_in_chunks(self, input_file, tmp_output_file, **kwargs):

        self._read_states()

        print(f'  Reading transitions from \"{input_file}\"')
        i = 0
        state_ID_u = []
        state_ID_l = []
        A = []

        # Read the transitions file in chunks to prevent memory overloads
        import bz2
        with bz2.open(input_file) as f:

            while True:
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

                    idx_l = np.searchsorted(self.states_ID, np.array(state_ID_l, dtype=int))
                    idx_u = np.searchsorted(self.states_ID, np.array(state_ID_u, dtype=int))
                    A = np.array(A, dtype=float)

                    E_u = self.states_E[idx_u] # [J]
                    E_l = self.states_E[idx_l]
                    nu_0 = np.abs(E_u - E_l) / sc.h # [J] -> [s^-1]

                    g_u = self.states_g[idx_u] # State degeneracy

                    J_l = self.states_J[idx_l] # Rotational quantum number
                    J_u = self.states_J[idx_u]

                    # Remove any nu = 0 transitions
                    A, E_l, nu_0, g_u, J_l, J_u = self.mask_arrays(
                        [A, E_l, nu_0, g_u, J_l, J_u], mask=(nu_0 > 0)
                    )

                    # Line-strengths at reference temperature
                    S_0 = (
                        A*g_u / (8*np.pi*(nu_0/sc.c)**2) *
                        np.exp(-E_l/(sc.k*self.T_0)) / self.q_0 * 
                        (1-np.exp(-sc.h*nu_0/(sc.k*self.T_0)))
                    ) # [s^-1/(molec. m^-2)]

                    idx = np.argsort(nu_0)
                    nu_0 = nu_0[idx]
                    S_0 = S_0[idx]
                    E_l = E_l[idx]
                    J_l = J_l[idx]
                    J_u = J_u[idx]
                    
                    # Get J-specific broadening parameters
                    self._read_broadening_per_transition(J_l, J_u)

                    # Compute the cross-sections, looping over the PT-grid
                    self.loop_over_PT_grid(
                        function=self.calculate_cross_sections,
                        nu_0=nu_0, S_0=S_0, E_low=E_l, **kwargs
                        )

                    # Reset
                    state_ID_u = []
                    state_ID_l = []
                    A = []

                if is_end_of_file:
                    break

                # Access info on upper and lower states
                state_ID_u.append(line[0:col_indices[0]])
                state_ID_l.append(line[col_indices[0]:col_indices[1]])
                A.append(
                    line[col_indices[1]:col_indices[2]] # Einstein A-coefficient [s^-1]
                    )

                i += 1

class Kurucz(LineByLine):

    parent_dir = pathlib.Path(__file__).parent.resolve()
    atoms_info = read_csv(parent_dir/'atoms_info.csv', index_col=0)

    def download_data(self, config):
        
        # TODO: Download states from NIST
        # ...

        if self.database == 'vald':
            raise NotImplementedError('Please download VALD transition-data manually.')

        # TODO: Download Kurucz data
        # ...

        raise NotImplementedError
    
    def __init__(self, config, **kwargs):

        print('-'*60)
        print('  Line-by-line Absorption from VALD/Kurucz')
        print('-'*60+'\n')

        super().__init__(config, **kwargs) # Initialise the parent LineByLine class

        raise NotImplementedError
    
    def _read_from_config(self, config):
        
        # Read the common parameters
        super()._read_from_config(config)
        
        # If alkali, read the E_ion and Z for different vdW-broadening
        self.E_ion = getattr(config, 'E_ion', None) # [cm^-1]
        if self.E_ion is not None:
            self.E_ion *= 1e2 # [cm^-1] -> [m^-1]
        self.Z = getattr(config, 'Z', None)

        # Transition energies to ignore
        self.nu_0_to_ignore = np.atleast_1d(getattr(config, 'nu_0_to_ignore', []))
        self.nu_0_to_ignore *= 1e2*sc.c # [cm^-1] -> [s^-1]

    def _read_partition_function(self, T_grid=np.arange(1,5001+1e-6,1)):
        """
        Read the partition function from the configuration file.
        """
        raise NotImplementedError('Partition function not implemented for Kurucz data.')

        file = self.config.files.get('partition_function', None)
        if file is None:
            raise ValueError('Partition function file must be provided in the configuration file.')
        # Load the states from NIST
        partition_function = read_csv(file, sep='\t', engine='python', header=0)

        g = np.array(partition_function['g'])
        E = np.array(partition_function['Level (cm-1)'])

        # (Higher) Ions beyond this index
        idx_u = np.min(np.argwhere(np.isnan(g)))
        g = g[:idx_u]
        E = E[:idx_u]

        # Partition function, sum over states, keep temperature-axis
        partition_function = np.array([np.sum(g*np.exp(-sc.c2*1e2*E/T_i)) for T_i in T_grid])

        # Make interpolation function, extrapolate outside temperature-range
        interpolation_function = interp1d(
            x=T_grid, y=partition_function, kind='linear', fill_value='extrapolate'
            )
        self.calculate_partition_function = interpolation_function

    def _read_transitions_in_chunks(self, input_file, tmp_output_file, **kwargs):
        raise NotImplementedError
