import numpy as np

import warnings
import wget
import pathlib
import h5py
import time
import datetime

import scipy.constants as sc
sc.c2 = 1.438777e-2 # [m K]
sc.L0 = sc.physical_constants['Loschmidt constant (273.15 K, 101.325 kPa)'][0] # [m^-3]
sc.amu = sc.physical_constants['atomic mass constant'][0] # [kg]

sc.E_H = sc.physical_constants['Rydberg constant times hc in J'][0] # [J]

sc.m_H = 1.00784*sc.amu # [kg]
sc.m_H2 = 2.01588*sc.amu  # [kg]
sc.m_He = 4.002602*sc.amu # [kg]
sc.alpha_H = 0.666793e-30 # Polarisability [m^3]

def download(url, out_dir, out_name=None):
    """
    Download a file from a URL if it does not already exist.

    Parameters:
    url (str): URL of the file to download.
    out_dir (str): Output directory.
    out_name (str): Output file name.

    Returns:
    str: Path to the downloaded file.
    """

    # Ensure the output directory exists
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    if out_name is None:
        out_name = url.split('/')[-1]
    out_name = pathlib.Path(out_dir) / out_name

    if pathlib.Path(out_name).is_file():
        print(f'  File \"{out_name}\" already exists, skipping download')
        return str(out_name)
    else:
        print(f'  Downloading \"{url}\"')
    
    # Download and rename
    try:
        tmp_file = wget.download(url, out=str(out_dir))
    except Exception as e:
        warnings.warn(f'Failed to download \"{url}\": {e}')

        try:
            # If wget fails, try using requests
            import requests
            tmp_file = str(out_dir / 'tmp.out')

            response = requests.get(url)
            if response.status_code == 200:
                with open(tmp_file, 'wb') as f:
                    f.write(response.content)
            else:
                raise Exception(f'{response.status_code} - {response.reason}')

        except Exception as e:
            warnings.warn(f'Failed to download \"{url}\": {e}')
            return
    print()
    tmp_file = pathlib.Path(tmp_file)

    out_name = tmp_file.rename(out_name)
    return str(out_name)

def save_to_hdf5(file, data, attrs, compression='gzip', **kwargs):
    """
    Save data to an HDF5 file.

    Parameters:
    file (str): Path to the output file.
    data (dict): Dictionary containing the data to save.
    attrs (dict): Dictionary containing the attributes to save.
    """
    # Make sure the output directory exists
    pathlib.Path(file).parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(file, 'w') as f:
        for key, value in data.items():

            value = np.atleast_1d(value)
            if value.dtype.kind in {'U', 'S'}:
                value = value.astype('S')  # Convert unicode to bytes
                
            dat_i = f.create_dataset(
                name=key, data=value, compression=compression, **kwargs
                )

            attrs_i = attrs.get(key, None)
            if attrs_i is None:
                continue
            for key_j, value_j in dict(attrs_i).items():
                dat_i.attrs[key_j] = value_j

def read_from_hdf5(file, keys_to_read, return_attrs=False):
    """
    Read data from an HDF5 file.

    Parameters:
    file (str): Path to the input file.

    Returns:
    dict: Dictionary containing the data.
    """

    datasets, datasets_attrs = {}, {}
    with h5py.File(file, 'r') as f:
        for key_i in keys_to_read:
            if key_i not in f.keys():
                continue
            datasets[key_i] = f[key_i][:]
            datasets_attrs[key_i] = dict(f[key_i].attrs)

    if return_attrs:
        return datasets, datasets_attrs

    return datasets

def log10_round(array, decimals=3):
    """
    Compute the base-10 logarithm of an array and round the result to save memory.

    Parameters:
    array (array-like): Input array for which the logarithm is computed.
    decimals (int): Number of decimal places to round to (default is 3).

    Returns:
    numpy.ndarray: Array with the base-10 logarithm values rounded to the specified decimals.
    """
    return np.around(np.log10(array), decimals=decimals)

def display_welcome_message():
    """
    Display a welcome message.
    """
    print('\n'+'='*80)
    print('  Welcome to pyROX: Rapid Opacity X-sections for Python')
    print('='*80+'\n')

    return time.time()

def display_finish_message(time_start):
    """
    Display a finish message and the elapsed time.

    Parameters:
    time_start (float): Start time of the process.
    """
    time_finish = time.time()
    time_elapsed = time_finish - time_start

    print('\nTime elapsed: {}'.format(str(datetime.timedelta(seconds=time_elapsed))))
    print('='*80+'\n')

def update_config_with_args(config=None, **kwargs):
    """
    Update the configuration object with command-line arguments.

    Parameters:
    config (object): Configuration object to update.
    kwargs (dict): Keyword arguments representing the parameters to update.

    Returns:
    object: Updated configuration object.
    """
    print('\nUpdating configuration with new parameters')

    if config is None:
        class Config:
            pass
        config = Config()

    for key, value in kwargs.items():
        if value is None:
            continue  # Parameter not given

        new_value = value
        if isinstance(value, str):
            new_value = f'\"{new_value}\"'

        # Different warning messages
        if hasattr(config, key):
            old_value = getattr(config, key)
            warnings.warn(f'Overwriting parameter \"{key}\" from {old_value} to {new_value}.')
        else:
            warnings.warn(f'Adding parameter \"{key}\" as {new_value}.')

        # Update or add the parameter
        setattr(config, key, value)
    print()

    return config

def warn_about_units(config):
    """
    Display a warning message about expected units for specific parameters.

    Parameters:
    config (object): Configuration object containing parameter definitions.
    """
    default_units = {
        'mass': 'amu',
        'P_grid': 'bar',
        'T_grid': 'K',
        'wave_min': 'um',
        'wave_max': 'um',
        'delta_wave': 'um',
        'delta_nu': 'cm^-1',
        'wave_file': 'um',
        'wing_cutoff': 'cm^-1',
        'global_cutoff': 'cm^1 molecule^-1',
    }

    keys_units_to_warn = []
    for key in dir(config):
        unit = default_units.get(key, None)
        if unit is None:
            continue

        keys_units_to_warn.append((key, unit))

    warnings.warn('Please make sure that the following parameters are given in the expected units:')
    for key, unit in keys_units_to_warn:
        print(f'  - {key} [{unit}]')
    print()

def find_closest_indices(a, b):
    """
    Find the indices of the closest elements in arrays a and b.
    """
    a_is_array = isinstance(a, (list, tuple, np.ndarray))
    b_is_array = isinstance(b, (list, tuple, np.ndarray))

    if a_is_array and b_is_array:
        idx = np.abs(np.asarray(a)[:,None] - np.asarray(b)[None,:]).argmin(axis=0)
    else:
        idx = np.abs(a - b).argmin()
    return idx, a[idx]

def fixed_resolution_wavelength_grid(wave_min, wave_max, resolution):
    """
    Adopted from petitRADTRANS. 
    Return a fixed resolution wavelength grid.

    Parameters:
    wave_min (float): Minimum wavelength.
    wave_max (float): Maximum wavelength.
    resolution (float): Desired resolution.

    Returns:
    numpy.ndarray: Wavelength grid.
    """

    # Get maximum space length (much higher than required)
    size_max = int(np.ceil((wave_max-wave_min) / (wave_min/resolution)))

    # Start generating space
    samples = [wave_min]
    i = 0

    inverse_resolution = 1 / resolution
    for i in range(size_max):
        samples.append(samples[-1] * np.exp(inverse_resolution))

        if samples[-1] >= wave_max:
            break

    if i == size_max - 1 and samples[-1] < wave_max:
        raise ValueError(f"maximum size ({size_max}) reached before reaching stop ({samples[-1]} < {stop})")

    return np.array(samples)


class Broaden_Gharib_Nezhad_ea_2021:
    """
    Broadening parameterisation from Gharib-Nezhad et al. (2021).
    """

    def pade_equation(self, J, a, b):
        """
        Pade approximation for the broadening coefficient.

        Parameters:
        J (float): Rotational quantum number.
        a (array-like): Coefficients for the numerator.
        b (array-like): Coefficients for the denominator.

        Returns:
        float: Broadening coefficient.
        """
        numerator = a[0] + a[1]*J + a[2]*J**2 + a[3]*J**3
        denominator = 1 + b[0]*J + b[1]*J**2 + b[2]*J**3 + b[3]*J**4
        
        return numerator / denominator * (1e2*sc.c) # [cm^-1] -> [s^-1]

    def __init__(self, species='AlH'):
        """
        Initialise the broadening parameters for a given species.

        Parameters:
        species (str): Name of the species.
        """
        
        if species in ['AlH']:
            self.a_H2 = [+7.6101e-02, -4.3376e-02, +1.9967e-02, +2.4755e-03]
            self.b_H2 = [-5.6857e-01, +2.7436e-01, +3.6216e-02, +1.5350e-05]
            self.a_He = [+4.8630e-02, +2.1731e+03, -2.5351e+02, +3.8607e+01]
            self.b_He = [+4.4644e+04, -4.4438e+03, +6.9659e+02, +4.7331e+00]

        elif species in ['CaH', 'MgH']:
            self.a_H2 = [+8.4022e-02, -8.2171e+03, +4.6171e+02, -7.9708e+00]
            self.b_H2 = [-9.7733e+04, -1.4141e+03, +2.0290e+02, -1.2797e+01]
            self.a_He = [+4.8000e-02, +7.1656e+02, -3.9616e+01, +6.7367e-01]
            self.b_He = [+1.4992e+04, +1.2361e+02, -1.4988e+01, +1.5056e+00]

        elif species in ['CrH', 'FeH', 'TiH']:
            self.a_H2 = [+7.0910e-02, -6.5083e+04, +2.5980e+03, -3.3292e+01]
            self.b_H2 = [-9.0722e+05, -4.3668e+03, +6.1772e+02, -2.4038e+01]
            self.a_He = [+4.2546e-02, -3.0981e+04, +1.2367e+03, -1.5848e+01]
            self.b_He = [-7.1977e+05, -3.4645e+03, +4.9008e+02, -1.9071e+01]

        elif species in ['SiO']:
            self.a_H2 = [+4.7273e-02, -2.7597e+04, +1.1016e+03, -1.4117e+01]
            self.b_H2 = [-5.7703e+05, -2.7774e+03, +3.9289e+02, -1.5289e+01]
            self.a_He = [+2.8364e-02, -6.7705e+03, +2.7027e+02, -3.4634e+00]
            self.b_He = [-2.3594e+05, -1.1357e+03, +1.6065e+02, -6.2516e+00]

        elif species in ['TiO', 'VO']:
            self.a_H2 = [+1.0000e-01, -2.4549e+05, +8.7760e+03, -8.7104e+01]
            self.b_H2 = [-2.3874e+06, +1.6350e+04, +1.7569e+03, -4.1520e+01]
            self.a_He = [+4.0000e-02, -2.8682e+04, +1.0254e+03, -1.0177e+01]
            self.b_He = [-6.9735e+05, +4.7758e+03, +5.1317e+02, -1.2128e+01]

        else:
            raise ValueError(f'Species \"{species}\" not recognised.')
    
    def gamma_H2(self, J):
        """
        Calculate the broadening coefficient for H2.

        Parameters:
        J (float): Rotational quantum number.

        Returns:
        float: Broadening coefficient.
        """
        return self.pade_equation(J, a=self.a_H2, b=self.b_H2) # [s^-1]
    
    def gamma_He(self, J):
        """
        Calculate the broadening coefficient for He.

        Parameters:
        J (float): Rotational quantum number.

        Returns:
        float: Broadening coefficient.
        """
        return self.pade_equation(J, a=self.a_He, b=self.b_He) # [s^-1]
