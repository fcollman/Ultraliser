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

    arg_help = 'Input directory'
    parser.add_argument('--input-directory', action='store', dest='input_directory', help=arg_help)

    arg_help = 'Output directory where output file is written'
    parser.add_argument('--output-directory', action='store', dest='output_directory', help=arg_help)

    # Parse the arguments
    return parser.parse_args()


####################################################################################################
# @aggregate_files_to_hdf5
####################################################################################################
def aggregate_files_to_hdf5(directory, hdf5_filename):

    # Open an HDF5 file in write mode
    with h5py.File(hdf5_filename, 'w') as hdf5_file:

        # Iterate over all files in the directory
        for file_name in os.listdir(directory):
            print(file_name)

            file_path = os.path.join(directory, file_name)
            if os.path.isfile(file_path):

                # Read the content of the file
                with open(file_path, 'r') as f:
                    content = f.read()

                # Create a dataset for each file, using the file name as the dataset name
                # Store the content as a variable-length string
                data_set = hdf5_file.create_dataset(file_name, data=content)


####################################################################################################
# @main
####################################################################################################
if __name__ == "__main__":

    args = parse_command_line_arguments()
    hdf5_filename = os.path.basename(os.path.normpath(args.input_directory))

    aggregate_files_to_hdf5(args.input_directory, '%s/%s.h5' % (args.output_directory, hdf5_filename))
    print(f"All files from {args.input_directory} have been aggregated into {hdf5_filename}.")


