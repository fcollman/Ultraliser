import h5py
import argparse
import os

####################################################################################################
# @parse_command_line_arguments
####################################################################################################
def parse_command_line_arguments(arguments=None):

    # Add all the options
    description = 'Aggregating a group of spine files into a single HDF5 file to reduce the ' \
                  'overhead to store too many files on the file system'
    parser = argparse.ArgumentParser(description=description)

    arg_help = 'Input file'
    parser.add_argument('--input-file', action='store', dest='input_file', help=arg_help)

    arg_help = 'Output directory where output file is written'
    parser.add_argument('--output-directory', action='store', dest='output_directory', help=arg_help)

    # Parse the arguments
    return parser.parse_args()


####################################################################################################
# @extract_hdf5_to_files
####################################################################################################
def extract_hdf5_to_files(hdf5_filename, output_directory):

    # Ensure the output directory exists
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Open the HDF5 file in read mode
    with h5py.File(hdf5_filename, 'r') as hdf5_file:

        # Iterate over all datasets in the HDF5 file
        for dataset_name in hdf5_file:

            # Read the content of the dataset
            data = hdf5_file[dataset_name][()]

            # Determine the output file path
            output_file_path = os.path.join(output_directory, dataset_name)

            # Write the content to a new file
            with open(output_file_path, 'w') as output_file:
                output_file.write(data.decode('utf-8'))  # Decode bytes to string before writing

####################################################################################################
# @main
####################################################################################################
if __name__ == "__main__":

    args = parse_command_line_arguments()
    hdf5_filename = args.input_file

    extract_hdf5_to_files(args.input_file, args.output_directory)
    print(f"All files from {args.input_file} have been extracted to {args.output_directory}.")


