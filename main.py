import argparse
import cross_sections, utils

if __name__ == '__main__':
    
    time_start = utils.display_welcome_message()

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
        '--save', '-s', action='store_true', help='Combine the outputs and save to a final file.'
        )
    parser.add_argument(
        '--plot', '-p', action='store_true', help='Plot the merged outputs.'
        )
    parser.add_argument(
        '--overwrite', '-o', action='store_true', help='Overwrite the existing output files.'
        )

    parser.add_argument(
        '--convert_to_pRT3', action='store_true', help='Convert to petitRADTRANS v3 format.'
        )

    parser.add_argument(
        '--progress_bar', '-pbar', action='store_true', default=False, 
        help='Show progress bar during calculations.'
        )

    # Overwrite some parameters from the configuration file
    parser.add_argument(
        '--tmp_output_basename', '-out', type=str, default=None, 
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
        '--files_range', type=int, nargs=2, default=None, 
        help='Specify the range of transition-files to read (e.g. 0 10).'
        )

    args = parser.parse_args()

    # Import input file as 'conf'
    config_string = str(args.config_file).replace('.py', '').replace('/', '.')
    config = __import__(config_string, fromlist=[''])

    # Overwrite some configuration parameters with command line arguments
    config = utils.update_config_with_args(
        config, 
        tmp_output_basename=args.tmp_output_basename, 
        P_grid=args.P_grid, 
        T_grid=args.T_grid, 
        )

    # Perform the requested operations
    data = cross_sections.load_data_object(config, download=args.download)
    if args.calculate:
        data.calculate_temporary_outputs(
            progress_bar=args.progress_bar, 
            files_range=args.files_range, 
            overwrite=args.overwrite,
            )
    if args.save:
        data.save_combined_outputs(
            overwrite=args.overwrite, 
            )
    if args.plot:
        data.plot_combined_outputs(
            T_to_plot=[100,200,500,1000,2000,3000], 
            P_to_plot=[1e3,1e5,1e7], # 1e-2, 1, 100 bar
            #xscale='linear', xlim=(0.3, 5)
            )

    # Optional conversions to petitRADTRANS format
    if args.convert_to_pRT3:
        data.convert_to_pRT3()

    utils.display_finish_message(time_start)