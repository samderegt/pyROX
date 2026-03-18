import numpy as np
import pathlib

from pyROX import utils, cross_sections
import sys
import shutil

import h5py

# Add the root directory to the Python path
parent_dir = pathlib.Path(__file__).resolve().parent
sys.path.append(str(parent_dir.parent))

reference_data_dir = parent_dir / 'reference_data'
reference_data_dir.mkdir(parents=True, exist_ok=True)

test_output_data_dir = parent_dir / 'test_output_data'

def test_lbl_exomol_alh():
    """
    Test the execution of the ExoMol AlH example.
    """
    import examples.exomol_alh.exomol_alh as config

    # Remove previous in- and output files if they exist
    if pathlib.Path(config.input_data_dir).is_dir():
        shutil.rmtree(config.input_data_dir, ignore_errors=True)
    if test_output_data_dir.is_dir():
        shutil.rmtree(test_output_data_dir, ignore_errors=True)

    # Reduce computational demand
    kwargs = dict(
        config=config, 
        output_data_dir=test_output_data_dir,
        P_grid=[1e0,1e2], T_grid=[500.,2000.], # [bar], [K]
        nu_min=1e4-200, nu_max=1e4, delta_nu=0.1, 
    )
    config = utils.update_config_with_args(**kwargs)

    # Run the calculation steps
    data = cross_sections.load_data_object(config, download=True)
    data.calculate_temporary_outputs(overwrite=True)
    data.save_combined_outputs(overwrite=True)
    data.plot_combined_outputs()
    data.convert_to_pRT3()

    # Check consistency with reference data
    compare_outputs(data.final_output_file, reference_data_dir/'exomol_alh_xsec.hdf5')
    compare_outputs(data.pRT_file, reference_data_dir/'exomol_alh_pRT3.hdf5')

def test_lbl_hitemp_co():
    """
    Test the execution of the HITEMP CO example.
    """
    import examples.hitemp_co.hitemp_co as config

    # Remove previous in- and output files if they exist
    if pathlib.Path(config.input_data_dir).is_dir():
        shutil.rmtree(config.input_data_dir, ignore_errors=True)
    if test_output_data_dir.is_dir():
        shutil.rmtree(test_output_data_dir, ignore_errors=True)

    # Avoid password prompt from HITRAN
    config.urls.remove('https://hitran.org/files/HITEMP/bzip2format/05_HITEMP2019.par.bz2')
    
    # Create a dummy transitions file
    test_output_data_dir.mkdir(parents=True, exist_ok=True)
    transitions_file = test_output_data_dir / pathlib.Path(config.files['transitions']).stem
    print(transitions_file)
    transition = ' 51  100.0000001.691E-031 1.284E-09.05550.061   50.00000.72-.024493             11              0                    R 13      447664 5 8 2 2 1 7    29.0   27.0'
    with open(transitions_file, 'w') as f:
        f.write(transition + '\n')
    config.files['transitions'] = str(transitions_file)

    # Reduce computational demand
    kwargs = dict(
        config=config, 
        output_data_dir=test_output_data_dir,
        P_grid=[1e0,1e2], T_grid=[500.,2000.], # [bar], [K]
        delta_nu=np.nan, wave_min=90, wave_max=120, resolution=1e5, # Test fixed resolution
    )
    config = utils.update_config_with_args(**kwargs)

    # Run the calculation steps
    data = cross_sections.load_data_object(config, download=True)
    data.calculate_temporary_outputs(overwrite=True)
    data.save_combined_outputs(overwrite=True)
    # data.plot_combined_outputs()

    # Check consistency with reference data
    compare_outputs(data.final_output_file, reference_data_dir/'hitemp_co_xsec.hdf5')

def test_lbl_kurucz_fe():
    """
    Test the execution of the Kurucz Fe example.
    """
    import examples.kurucz_fe.kurucz_fe as config

    # Remove previous in- and output files if they exist
    if pathlib.Path(config.input_data_dir).is_dir():
        shutil.rmtree(config.input_data_dir, ignore_errors=True)
    if test_output_data_dir.is_dir():
        shutil.rmtree(test_output_data_dir, ignore_errors=True)

    # Reduce computational demand
    kwargs = dict(
        config=config, 
        output_data_dir=test_output_data_dir,
        P_grid=[1e0,1e2], T_grid=[500.,2000.], # [bar], [K]
        nu_min=1e4-200, nu_max=1e4, delta_nu=0.1, 
    )
    config = utils.update_config_with_args(**kwargs)

    # Run the calculation steps
    data = cross_sections.load_data_object(config, download=True)
    data.calculate_temporary_outputs(overwrite=True)
    data.save_combined_outputs(overwrite=True)
    # data.plot_combined_outputs()

    # Check consistency with reference data
    compare_outputs(data.final_output_file, reference_data_dir/'kurucz_fe_xsec.hdf5')

def test_cia_hitran_h2_h2():
    """
    Test the execution of the HITRAN H2-H2 CIA example.
    """
    import examples.cia_hitran_h2_h2.cia_hitran_h2_h2 as config

    # Remove previous in- and output files if they exist
    if pathlib.Path(config.input_data_dir).is_dir():
        shutil.rmtree(config.input_data_dir, ignore_errors=True)
    if test_output_data_dir.is_dir():
        shutil.rmtree(test_output_data_dir, ignore_errors=True)

    # Reduce computational demand
    kwargs = dict(
        config=config, 
        wave_min=1.0, wave_max=2.0, 
        output_data_dir=test_output_data_dir,
    )
    config = utils.update_config_with_args(**kwargs)

    # Run the calculation steps
    data = cross_sections.load_data_object(config, download=True)
    data.calculate_temporary_outputs(overwrite=True)
    data.save_combined_outputs(overwrite=True)
    data.plot_combined_outputs(T_to_plot=[500.,1000.,1500.,2000.,2500.,3000.])
    data.convert_to_pRT3()

    # Check consistency with reference data
    compare_outputs(data.final_output_file, reference_data_dir/'cia_hitran_h2_h2_xsec.hdf5')
    compare_outputs(data.pRT_file, reference_data_dir/'cia_hitran_h2_h2_pRT3.hdf5')

def compare_outputs(new_file, ref_file):
    """
    Compare the calculated output with the reference data.

    Args:
        new_file (pathlib.Path): Path to the new .h5 file.
        ref_file (pathlib.Path): Path to the reference .h5 file.
    """

    if not ref_file.exists():
        response = input(f'Reference output file \"{ref_file}\" does not exist. Create it from the new output \"{new_file}\"? (y/n) ')
        if response.lower() not in ['y', 'yes']:
            raise FileNotFoundError(f'Reference file \"{ref_file}\" does not exist.')
        new_file.rename(ref_file)
        return
    
    # Load the calculated and reference data
    with h5py.File(new_file, 'r') as f_new, h5py.File(ref_file, 'r') as f_ref:
        assert set(f_new.keys()) == set(f_ref.keys()), (
            'The keys in the new output-file do not match the keys in the reference file.'
        )
        
        # Compare all datasets
        for key in f_new.keys():
            ds_new = f_new[key][:]
            ds_ref = f_ref[key][:]

            if (ds_new.dtype not in [float, int]) or (ds_ref.dtype not in [float, int]):
                continue
            assert np.allclose(f_new[key][:], f_ref[key][:], rtol=1e-5, atol=1e-5), (
                f'Key \"{key}\" does not match between the new output and the reference data.'
            )

if __name__ == '__main__':
    # Run the tests
    test_lbl_exomol_alh()
    test_lbl_hitemp_co()
    test_lbl_kurucz_fe()
    test_cia_hitran_h2_h2()