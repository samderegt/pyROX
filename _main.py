import argparse
import cross_sections, utils

if __name__ == '__main__':

    # Instantiate the parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'input_file', type=str, help='Name of input file (e.g. input.example_ExoMol)', 
        )
    parser.add_argument('--download', '-d', action='store_true')
    parser.add_argument('--cross_sections', '-cs', action='store_true')
    parser.add_argument('--save', '-s', action='store_true')
    parser.add_argument('--plot', action='store_true')

    parser.add_argument('--convert_to_pRT2', action='store_true')
    parser.add_argument('--convert_to_pRT3', action='store_true')

    args = parser.parse_args()

    # Import input file as 'conf'
    input_string = str(args.input_file).replace('.py', '').replace('/', '.')
    config = __import__(input_string, fromlist=[''])

    utils.print_welcome_message()

    # Perform the requested operations
    data = cross_sections.load_data_object(config)
    if args.download:
        data.download_data()
    if args.cross_sections:
        data.calculate_tmp_outputs()
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