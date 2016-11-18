#!/usr/bin/python

# Name: Ryan Hagenson
# Email: rhagenson@unomaha.edu
import os
from datetime import datetime
from getopt import GetoptError, getopt
import sys
from os import path, makedirs, walk

profile_dir = ""  # Location of the profiles directory
now = datetime.now().strftime("%d-%m-%y")  # Default run time
isoformsSubDirName = "isoforms"

# General directory tree within dataDir is:
# ./allMuts
# ./refSeq
# ./refSeq/iupredLong
# ./refSeq/iupredShort
# ./refSeq/hmmer
# ./cds
# ./allMAFs
# ./pfam30.0
# ./profiles

# ./profiles will have subdirectories based on cancer type, then by gene id
# ./profiles/isoforms/ contains cancer-independent profiles of all isoforms


def main():
    """
    A simple wrapper for all CLI options
    """
    global profile_dir, now

    # Enables command-line options via getopt and sys packages
    try:
        opts, args = getopt(sys.argv[1:],
                            'p:d:',
                            ["profDir=", "date="]
                            )
    except GetoptError as err:
        # Redirect STDERR to STDOUT (ensures screen display)
        sys.stdout = sys.stderr
        # Print help information
        print(str(err))
        # Print usage information
        # usage()
        # Exit
        sys.exit(2)

    # Configure the action of each CLI option
    # First loop for global variables with defaults
    for (opt, arg) in opts:
        if opt in ("-d", "--date"):
            now = arg

    for (opt, arg) in opts:
        if opt in ("-p", "--profDir"):  # Where profiles are found
            # Create profiles directory with now date
            profile_dir = path.join(arg, now)


#
#
# Important note: this should be looping through the R/outputs/DD-MM-YY/
# directory tree, NOT the profiles directory. Preferably the concatenated
# file will reside at the peak of that tree too, see below example:
#       .../DD-MM-YY/BRCA/BRCA_LOG.csv
#
#
def generate_cancer_logs(profiles_dir=profile_dir):
    """
    :return: A single file per cancer type of the complete LOGs
    """

    cancer_types = os.listdir(profiles_dir)
    cancer_types.remove(isoformsSubDirName)

    # for type in cancer_types:
        # CONCATLOG = open()
    # for (dirpath, dirnames, filenames) in walk(profiles_dir):

if __name__ == "__main__":
    # Run the CLI wrapper to change global variables
    main()

    # Run the program
    generate_cancer_logs(profile_dir)
