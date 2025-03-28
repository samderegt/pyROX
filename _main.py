import argparse
import cross_sections, utils

if __name__ == '__main__':
    
    utils.print_welcome_message()

    # Instantiate the parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'config_file', type=str, help='Configuration file (e.g. path/to/config.py)'
        )
    parser.add_argument(
        '--download', '-d', action='store_true', help='Download data from the given database.'
        )
    parser.add_argument(
        '--calculate', '-c', action='store_true', 
        help='Calculate the outputs (e.g. cross-sections) and save to temporary files.'
        )
    parser.add_argument(
        '--save', '-s', action='store_true', help='Merge the outputs and save to a final file.'
        )
    parser.add_argument(
        '--plot', action='store_true', help='Plot the merged outputs.'
        )

    parser.add_argument(
        '--convert_to_pRT2', action='store_true', help='Convert to petitRADTRANS v2 format.'
        )
    parser.add_argument(
        '--convert_to_pRT3', action='store_true', help='Convert to petitRADTRANS v3 format.'
        )

    parser.add_argument(
        '--show_progress_bar', action='store_true', default=False, 
        help='Show progress bar during calculations.'
        )

    # Overwrite some parameters from the configuration file
    parser.add_argument(
        '--tmp_output_basename', type=str, default=None, 
        help='Basename of temporary output file (e.g. xsec_new.hdf5).'
        )

    # Optionally, overwrite the P and T values
    parser.add_argument(
        '--P_grid', type=float, nargs='+', default=None, 
        help='Specify the pressure grid (e.g. 0.1 1 10).'
        )
    parser.add_argument(
        '--T_grid', type=float, nargs='+', default=None, 
        help='Specify the temperature grid (e.g. 100 200 500 1000).'
        )

    parser.add_argument(
        '--transition_files_range', type=int, nargs=2, default=None, 
        help='Specify the range of transition-files to read (e.g. 0 10).'
        )

    args = parser.parse_args()

    # Import input file as 'conf'
    config_string = str(args.config_file).replace('.py', '').replace('/', '.')
    config = __import__(config_string, fromlist=[''])

    # Overwrite some configuration parameters with command line arguments
    config = utils.add_to_config(
        config, 
        tmp_output_basename=args.tmp_output_basename, 
        P_grid=args.P_grid, 
        T_grid=args.T_grid, 
        )

    # Perform the requested operations
    data = cross_sections.load_data_object(config, download=args.download)
    if args.calculate:
        data.calculate_tmp_outputs(
            show_progress_bar=args.show_progress_bar, 
            transition_files_range=args.transition_files_range
            )
    if args.save:
        data.save_merged_outputs()
    if args.plot:
        data.plot_merged_outputs(T_to_plot=[100,200,500,1000,2000,3000])

    # Optional conversions to petitRADTRANS format
    if args.convert_to_pRT2:
        data.convert_to_pRT2()
    if args.convert_to_pRT3:
        data.convert_to_pRT3()

    print('\n'+'='*80+'\n')