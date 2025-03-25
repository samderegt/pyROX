from scipy.interpolate import interp1d
import scipy.constants as sc
sc.L0 = sc.physical_constants['Loschmidt constant (273.15 K, 101.325 kPa)'][0]

import matplotlib.pyplot as plt
import numpy as np

import pathlib
import h5py

class CIA_HITRAN:

    def download_data(self, conf):
        raise NotImplementedError

        # Download CIA data from HITRAN
        # ...

    def __init__(self, conf):

        self._set_nu_grid(conf)

        # Load CIA data from HITRAN
        # ...
        pass

    def _set_nu_grid(self, conf):

        # Pre-computed wavelength-grid
        wave_file = conf.files.get('wavelength')
        if wave_file is not None:
            raise NotImplementedError('Custom wavelength-grids are not yet implemented.')
        
        # Create new wavenumber grid
        self.delta_nu = conf.delta_nu # [cm^-1]
        self.delta_nu *= 100.0*sc.c # [cm^-1] -> [s^-1]

        self.wave_min = conf.wave_min*1e-6 # [um] -> [m]
        self.wave_max = conf.wave_max*1e-6 # [um] -> [m]

        self.nu_min = sc.c/self.wave_max # [m] -> [s^-1]
        self.nu_max = sc.c/self.wave_min # [m] -> [s^-1]

        # Number of grid points
        self.N_grid   = int((self.nu_max-self.nu_min)/self.delta_nu) + 1

        # Not exact value of delta_nu given above, but done to keep final lambda values fixed
        self.delta_nu = (self.nu_max-self.nu_min) / (self.N_grid-1)
        self.nu_grid  = np.linspace(
            self.nu_min, self.nu_max, num=self.N_grid, endpoint=True
        )

    def get_cross_sections(self, conf, tmp_output_file='cia{}.hdf5', show_pbar=True, **kwargs):

        # Load CIA data from HITRAN
        file = conf.files['cia']
        with open(file, 'r') as f:
            lines = f.readlines()

        N_grid, self.T_grid = [], []
        for line in lines:
            if len(line) < 101:
                continue
            # Header line
            N_grid.append(int(line[40:40+7]))
            self.T_grid.append(float(line[47:47+7]))

        self.T_grid = np.array(self.T_grid)
        N_grid = np.array(N_grid) # TODO: use different name

        self.k_grid = np.zeros((len(self.T_grid), self.N_grid))
        for i, (N_grid_i, T_i) in enumerate(zip(N_grid, self.T_grid)):
            idx_min = 1 + i*(N_grid_i+1)
            idx_max = idx_min + N_grid_i

            nu = np.array(
                [line[0:10].strip() for line in lines[idx_min:idx_max]], dtype=np.float64
                )
            nu *= 100.0*sc.c # [cm^-1] -> [s^-1]

            k = np.array(
                [line[10:21].strip() for line in lines[idx_min:idx_max]], dtype=np.float64
                )
            k *= (1e-2)**5 # [cm^5 molecule^-2] -> [m^5 molecule^-2]

            # Interpolate onto nu_grid
            interp_func = interp1d(nu, np.log10(k), kind='linear', fill_value=np.nan, bounds_error=False)
            self.k_grid[i] = 10**interp_func(self.nu_grid)
            self.k_grid[i] = np.nan_to_num(self.k_grid[i], nan=0.0)

        # [m^5 molecule^-2] * [m^-6] = [m^-1 molecule^-2]
        self.alpha_grid = self.k_grid * sc.L0**2

        # Save cross-sections to file
        tmp_output_file = f'{conf.output_dir}/' + tmp_output_file.format('')
        self.save_cross_sections(tmp_output_file)

    def save_cross_sections(self, file):
        
        print(f'\nSaving cross-sections to file \"{file}\"')

        # Create directory if not exist
        pathlib.Path(file).parent.mkdir(parents=True, exist_ok=True)

        with h5py.File(file, 'w') as f:
            # Flip arrays to be ascending in wavelength
            wave       = sc.c / self.nu_grid[::-1]
            alpha_grid = self.alpha_grid[:,::-1]
            k_grid     = self.k_grid[:,::-1]

            f.create_dataset('T', compression='gzip', data=self.T_grid) # [K]
            f.create_dataset('wave', compression='gzip', data=wave) # [m]
            f.create_dataset(
                'alpha', compression='gzip', 
                data=np.log10(alpha_grid+1e-250) # [m^-1 molecule^-2] -> log10([m^-1 molecule^-2])
                ) 
            f.create_dataset(
                'k', compression='gzip', 
                data=np.log10(k_grid+1e-250) # [m^5 molecule^-2] -> log10([m^5 molecule^-2])
                )

        return

'''
        plt.figure(figsize=(8,6))
        for i, T_i in enumerate(T):
            if T_i%100 != 0:
                continue
            if T_i > 3000:
                continue
            plt.plot(1e6*sc.c/self.nu_grid, self.alpha_grid[i], c=plt.get_cmap('RdBu_r')(T_i/3000), label=str(T_i))
        plt.xscale('log'); plt.yscale('log')
        plt.xlim(0.5,250); plt.ylim(1e-14,1e-4)
        plt.legend(ncols=2)
        plt.savefig(conf.output_dir + 'cia.pdf')
        plt.close()

        import h5py
        hdf5_file = '/net/schenk/data2/regt/pRT3_input_data/input_data/opacities/continuum/collision_induced_absorptions/H2--H2/H2--H2-NatAbund/H2--H2-NatAbund__BoRi.R831_0.6-250mu.ciatable.petitRADTRANS.h5'
        with h5py.File(hdf5_file, 'r') as f:
            nu = f['wavenumbers'][:]
            alpha = f['alpha'][:]
            temperature = f['t'][:]

        plt.figure(figsize=(8,6))
        for i, T_i in enumerate(temperature):
            if T_i > 3000:
                continue
            plt.plot(1e4/nu, alpha[i], c=plt.get_cmap('RdBu_r')(T_i/3000), label=str(T_i))
        plt.xscale('log'); plt.yscale('log')
        plt.xlim(0.5,250); plt.ylim(1e-14,1e-4)
        plt.legend(ncols=2)
        plt.savefig(conf.output_dir + 'cia_petitRADTRANS.pdf')
        plt.close()
'''