import argparse
import cross_sections

if __name__ == '__main__':

    # Instantiate the parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'input_file', type=str, help='Name of input file (e.g. input.example_ExoMol)', 
        )
    parser.add_argument('--download', '-d', action='store_true')
    parser.add_argument('--cross_sections', '-cs', action='store_true')

    args = parser.parse_args()

    # Import input file as 'conf'
    input_string = str(args.input_file).replace('.py', '').replace('/', '.')
    config = __import__(input_string, fromlist=[''])


    data = cross_sections.load_data_object(config)
    if args.download:
        data.download_data()
    if args.cross_sections:
        data.calculate_cross_sections()