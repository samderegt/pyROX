import numpy as np
import pathlib

from pyROX import utils, cross_sections
import shutil
import sys

# Add the root directory to the Python path
parent_dir = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(parent_dir.parent))

reference_data_dir = parent_dir / 'reference_data'
reference_data_dir.mkdir(parents=True, exist_ok=True)

input_data_dir  = parent_dir / 'input_data'
output_data_dir = parent_dir / 'output_data'

def test_lbl_kurucz():
    """
    Test the Kurucz line profile calculation.
    """
    # Filenames of synthetic data
    input_data_dir.mkdir(parents=True, exist_ok=True)
    states_file      = input_data_dir / 'kurucz_states.txt'
    transitions_file = input_data_dir / 'kurucz_transitions.txt'

    # Create the synthetic states file
    states_data = [
        '2\t\"\"\t\"0.0000\"\t\"\"',
        '2\t\"\"\t\"12985.185724\"\t\"\"',
        '4\t\"\"\t\"13042.896027\"\t\"\"',
        '2\t\"\"\t\"21026.551\"\t\"\"',
        '6\t\"\"\t\"21534.680\"\t\"\"',
        '4\t\"\"\t\"21536.988\"\t\"\"',
        '\t\"\"\t\"35009.8140\"\t\"\"'
        
    ]
    with open(states_file, 'w') as f:
        f.write("g\tPrefix\tLevel (cm-1)\tSuffix\n")
        f.write("\n".join(states_data) + "\n")

    # Create the synthetic transitions file
    transition = '   100.0000 -2.503 19.00      50.000  0.5 4s 2S          100.000  1.5 3d3D4s 2D   9.41 -5.87 -7.72K12  0 0  0 0.000  0 0.000                     2002  801     0'
    with open(transitions_file, 'w') as f:
        f.write(transition + '\n')

    # Import the example configuration
    import examples.kurucz_k.kurucz_k as config
    files = dict(states=states_file, transitions=transitions_file)
    config = setup_lbl_config(config, files)

    # Compare the line profiles to a reference file
    compare_outputs(config, reference_data_dir/'kurucz_xsec.txt')
   
def test_lbl_hitran():
    """
    Test the HITRAN line profile calculation.
    """
    # Filenames of synthetic data
    input_data_dir.mkdir(parents=True, exist_ok=True)
    transitions_file = input_data_dir / 'hitran_transitions.txt'
    partition_function_file = input_data_dir / 'hitran_pf.txt'

    # Create the synthetic transition file
    transition = ' 51  100.0000001.691E-031 1.284E-09.05550.061   50.00000.72-.024493             11              0                    R 13      447664 5 8 2 2 1 7    29.0   27.0'
    with open(transitions_file, 'w') as f:
        f.write(transition + '\n')

    # Create the synthetic partition function
    T = np.arange(0, 3000., 1)
    with open(partition_function_file, 'w') as f:
        for t in T:
            f.write(f'{t} {2*t}\n')
    
    # Import the example configuration
    import examples.hitemp_co.hitemp_co as config
    files = dict(
        transitions=transitions_file, 
        partition_function=partition_function_file
    )
    config = setup_lbl_config(config, files)

    # Compare the line profiles to a reference file
    compare_outputs(config, reference_data_dir/'hitran_xsec.txt')

def test_lbl_exomol():
    """
    Test the ExoMol line profile calculation.
    """
    # Filenames of synthetic data
    input_data_dir.mkdir(parents=True, exist_ok=True)
    states_file      = input_data_dir / 'exomol_states.txt'
    transitions_file = input_data_dir / 'exomol_transitions.txt'
    partition_function_file = input_data_dir / 'exomol_pf.txt'

    # Create the synthetic states file
    upper_state = '           1    150.000000     36       1'
    lower_state = '           2     50.000000     12       0'
    with open(states_file, 'w') as f:
        f.write(upper_state + '\n')
        f.write(lower_state + '\n')

    # Create the synthetic transitions file
    transition = '           1            2  3.5960E-05'
    with open(transitions_file, 'w') as f:
        f.write(transition + '\n')

    # Create the synthetic partition function
    T = np.arange(0, 3000., 1)
    with open(partition_function_file, 'w') as f:
        for t in T:
            f.write(f'{t} {2*t}\n')

    # Import the example configuration
    import examples.exomol_alh.exomol_alh as config
    files = dict(
        states=states_file, 
        transitions=transitions_file, 
        partition_function=partition_function_file
    )
    config = setup_lbl_config(config, files, resolution=1e4) # Test fixed resolution

    # Compare the line profiles to a reference file
    compare_outputs(config, reference_data_dir/'exomol_xsec.txt')

def setup_lbl_config(config, files, resolution=np.nan):
    """
    Set up the configuration for the line profile calculation.

    Args:
        config (module): Configuration module containing the parameters.
        files (dict): Dictionary containing the paths to the input files.
        resolution (float): Resolution of the line profile calculation.

    Returns:
        module: Updated configuration module.
    """
    kwargs = dict(
        config=config, 
        output_data_dir=output_data_dir, 
        input_data_dir=input_data_dir, 
        files=files, 
        P_grid=np.logspace(-4, 2, 5), T_grid=2000., # [bar], [K]
        nu_min=1, nu_max=200, wave_file='', adaptive_nu=False, 
        delta_nu=0.1, delta_wave=np.nan, resolution=resolution, # delta_nu is used
        perturber_info={'H2':{'VMR':1.0,'gamma':0.07,'n':0.5}}
    )

    if not np.isnan(resolution):
        kwargs['resolution'] = resolution
        kwargs['wave_min'] = 90
        kwargs['wave_max'] = 120
        kwargs['delta_nu'] = np.nan
        
    # Update the configuration object
    config = utils.update_config_with_args(**kwargs)
    return config


def test_cia_hitran():
    """
    Test the CIA calculation using HITRAN data.
    """
    # Filenames of synthetic data
    input_data_dir.mkdir(parents=True, exist_ok=True)
    cia_file = input_data_dir / 'cia_hitran.txt'

    # Create the synthetic CIA file
    cia_data = [
        '               H2-H2    20.000   200.000      5 2000.0 8.788E-45 -.999                             6', 
        '   20.0000  2.668E-47', 
        '   21.0000  2.931E-47', 
        '   22.0000  3.203E-47', 
        '   23.0000  3.485E-47', 
        '  200.0000  2.074E-47', 
    ]
    with open(cia_file, 'w') as f:
        f.write("\n".join(cia_data) + "\n")
    
    # Import the example configuration
    import examples.cia_hitran_h2_h2.cia_hitran_h2_h2 as config
    files = dict(cia=[cia_file])
    config = setup_cia_config(config, files)

    # Compare the CIA data to a reference file
    compare_outputs(
        config, reference_data_dir/'cia_hitran_alpha.txt', 
        key_to_compare='log10(alpha)'
        )

def setup_cia_config(config, files):
    """
    Set up the configuration for the CIA calculation.

    Args:
        config (module): Configuration module containing the parameters.
        files (dict): Dictionary containing the paths to the input files.
    
    Returns:
        module: Updated configuration module.
    """
    # Update the configuration object
    config = utils.update_config_with_args(
        config=config, 
        output_data_dir=output_data_dir, 
        input_data_dir=input_data_dir, 
        files=files, 
        T_grid=2000., # [K]
        nu_min=1, nu_max=200, delta_nu=1.0, adaptive_nu=False,
        wave_file='', delta_wave=np.nan, resolution=np.nan, # delta_nu is used
    )
    return config


def compare_outputs(config, reference_file, key_to_compare='log10(xsec)'):
    """
    Compare the calculated output with the reference data.

    Args:
        config (module): Configuration module containing the parameters.
        reference_file (pathlib.Path): Path to the reference file.
        key_to_compare (str): Key to compare in the calculated data.
    """
    # Run the calculation steps
    data = cross_sections.load_data_object(config, download=False)
    data.calculate_temporary_outputs(overwrite=True)
    data.save_combined_outputs(overwrite=True)
    data.plot_combined_outputs()

    # Clean up the input_data and output_data directories
    shutil.rmtree(input_data_dir, ignore_errors=True)
    shutil.rmtree(output_data_dir, ignore_errors=True)

    # Compare the calculated data with the reference data
    calculated_output = data.combined_datasets[key_to_compare].flatten()
    if not reference_file.exists():
        print('Reference file does not exist. Creating it.')
        np.savetxt(reference_file, calculated_output)
        return
    reference_output = np.loadtxt(reference_file)

    assert np.allclose(calculated_output, reference_output, rtol=1e-5, atol=1e-5), (
        'Calculated output does not match the reference data.'
    )

if __name__ == '__main__':
    # Run the tests
    test_lbl_kurucz()
    test_lbl_hitran()
    test_lbl_exomol()

    test_cia_hitran()