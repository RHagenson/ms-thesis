#!/usr/bin/python

# Name: Ryan Hagenson
# Email: rhagenson@unomaha.edu

from operator import itemgetter
import sys
from csv import reader, writer
from datetime import datetime
from distutils.dir_util import mkpath
from getopt import GetoptError, getopt
from multiprocessing import Pool, cpu_count
from os import path, makedirs, listdir, walk
from os.path import basename
from re import search, compile

from shutil import rmtree

dataDir = ""  # Default False, should be overwritten at CLI
allMAFsName = "allMAFs"  # The name of the allMAFs dir in dataDir
allMutsName = "allMuts"  # The name of the allMuts dir in dataDir
cdsName = "cds"  # The name of the cds dir in dataDir
pfamName = "pfam30.0"  # The name of the pfam dir in dataDir
refSeqName = "refSeq"  # The name of the refSeq dir in dataDir
profilesName = "profiles"  # The name of the final mutation profile csv's dir
isoformsSubDirName = "isoforms"

cancerTypes = ['BRCA']
now = datetime.now().strftime("%d-%m-%y")  # Default run time

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
    global dataDir, now, cancerTypes

    # Enables command-line options via getopt and sys packages
    try:
        opts, args = getopt(sys.argv[1:],
                            'd:c:',
                            ["date=", "dataDir=", "cancerTypes="]
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
    for (opt, arg) in opts:
        if opt in ("-d", "--dataDir"):  # Set high-level data directory location
            dataDir = arg

            # Create profiles directory with now date
            profile_dir = path.join(dataDir, profilesName, now)
            if not path.exists(profile_dir):
                makedirs(profile_dir)
                del profile_dir

        if opt in ("-c",  "--cancerTypes"):
            cancerTypes = arg.split(',')

        if opt == "--date":
            now = arg


def create_csv_profile((mut_file, long_short_file)):
    """
    Run by Pool.map() with data from generate_data_pairs()
    :return:
    """
    global dataDir, allMutsName, refSeqName, profilesName, isoformsSubDirName

    # Captures three groups:
    # .group(1): position number
    # .group(2): amino acid 1-letter code
    # .group(3): disorder score
    long_short_re = compile('\s+(\d+)\s+(\w+)\s+(.+)')

    # Build and open the allMuts file
    try:
        mut_file_handle = open(path.join(dataDir,
                                         allMutsName,
                                         mut_file), 'r')
    except IOError as e:
        print(str(e))  # send the error out for bug tracking

        # Return if the file does not exists,
        # therefore the datapair is invalid
        return

    # Build and open the long or short file based on extension
    try:
        if '.long' in long_short_file:
            long_short_file_handle = open(path.join(dataDir,
                                                    refSeqName,
                                                    "iupredLong",
                                                    long_short_file), 'r')
        elif '.short' in long_short_file:
            long_short_file_handle = open(path.join(dataDir,
                                                    refSeqName,
                                                    "iupredShort",
                                                    long_short_file), 'r')
        else:
            # Break if no long/short file is found
            sys.exit(2)
    except IOError as e:
        print(str(e))  # send the error out for bug tracking

        # Return if the file does not exists,
        # therefore the datapair is invalid
        return

    # Generate profiles directory tree with each cancer type and gene id
    cancer_type = search("(\w+)\_.+\.txt", mut_file).group(1)
    geneMatch = search("([\w|-]+)+\.\d+\.([long|short]+)",
                       long_short_file)
    # GENE.long or GENE.short, separates long and short at the gene level
    gene_name = ".".join(map(str, geneMatch.group(1, 2)))

    # Create by cancer-type and cancer-independent paths
    full_path = path.join(dataDir,
                          profilesName,
                          now,
                          cancer_type,
                          gene_name)

    isoform_path = path.join(dataDir,
                             profilesName,
                             now,
                             isoformsSubDirName)

    # If full_path does not exist, isoform_path should not
    if not path.exists(full_path):
        mkpath(full_path)
        mkpath(isoform_path)

    # profile_file is a profile for each isoform, dependent on cancer type
    profile_file = open(path.join(dataDir,
                                  profilesName,
                                  now,
                                  cancer_type,
                                  gene_name,
                                  long_short_file + ".prof"), 'w')
    profile_csv = writer(profile_file, delimiter='\t')

    # Extract the mutations from the allMuts file
    isoform_name, long_short = path.splitext(long_short_file)

    print "Processing " + str(long_short_file)  # Inform user what is being done

    # Sort the file based on isoform name so mutations are sequential
    mut_csv = sorted(reader(mut_file_handle,
                            delimiter='\t'),
                     key=itemgetter(0))

    # Set the mutation counts for each position, requires sorted isoform lines in mut_csv
    mutations = {}  # Should be {pos# : count}
    for row in mut_csv:
        # If we have found our mutations and there are no more for this
        # isoform break the loop, depends on isoform lines being sequential
        # which is handled by sorting the file
        if len(mutations) > 0:
            if row[0] != isoform_name:
                break

        # Found a mutation in the protein
        if row[0] == isoform_name:
            # If the pos has a mutation, iterate or initialize it
            if row[6] in mutations:
                mutations[row[5]] += 1
            else:
                mutations[row[5]] = 1

    # Gather the rest of the information from long_short file
    for line in long_short_file_handle:
        # Skip comment lines at start
        if line.startswith('#'):
            continue

        long_short_match = search(long_short_re, line)
        if long_short_match:
            pos_muts = 0
            if long_short_match.group(1) in mutations:
                pos_muts = mutations[long_short_match.group(1)]

            # Output the individual isoform results to the pertinent file: profile_csv
            # Combining these files into gene, isoform, and cancer-level is done later
            profile_csv.writerow([long_short_match.group(1),
                                  long_short_match.group(2),
                                  long_short_match.group(3),
                                  pos_muts])

    # Be sure to release the file to free resources
    profile_file.close()


def generate_data_pairs(ctype):
    """
    :arg ctype: Which cancer is currently being processed
    :type ctype: str

    Builds an iterable list of data pairs for Pool.map()
    Reads each file in data/allMuts/, for each line it determines if that
    protein has a corresponding file in iupredLong|iupredShort, if it does it
    adds a new entry in datapairs in style ['<allMuts filename>',
    '<iupredLong prop.XXX>.long']
    """
    global dataDir, allMutsName, refSeqName

    # Define the absolute path to the allMuts directory
    mut_loc = path.join(dataDir, allMutsName)

    print("Generating data pairs")

    # Define the return iterable
    # Entries should be in form [[allMutsFile, iupredFile]]
    # Each with full absolute path
    datapairs = []

    # Collect the files in allMuts with absolute pathing
    mut_filepaths = []
    for f in listdir(mut_loc):
        # Only process the file that matches ctype
        if ctype+"_" in f:
            mut_filepaths.append(path.join(mut_loc, f))

    # Process each Mut file in turn
    # If the gene/isoform combination has an iupred entry, create a profile
    # for that entry
    for mutFile in mut_filepaths:
        print "Now processing: " + mutFile

        mut_name = path.basename(mutFile)
        with open(mutFile, 'r') as FILE:
            csv = reader(FILE, delimiter='\t')
            past_isoform = ""
            for row in csv:
                protein_isoform = row[0]

                # Ensure each isoform only gets processed once by skipping
                # repetitions
                if protein_isoform == past_isoform:
                    continue
                else:
                    past_isoform = protein_isoform

                # Check if a long, short, or both files exist
                long_path = path.join(dataDir, refSeqName,
                                      "iupredLong", protein_isoform + ".long")
                short_path = path.join(dataDir, refSeqName,
                                       "iupredShort",
                                       protein_isoform + ".short")

                # Create a new datapairs entry for each file found
                if path.exists(long_path):
                    datapairs.append([mut_name, protein_isoform + ".long"])
                if path.exists(short_path):
                    datapairs.append([mut_name, protein_isoform + ".short"])

    # Once all the allMuts files have been fully processed, return the datapairs
    return datapairs


# Need to add functionality that concatenates isoforms from across cancer
# types once all types have been processed, can be done recursively by
# matching ctype isoforms to ones in profiles/now/isoforms/ or as true post
# processing once all types are done (do not concate non-isoform files)
def concatenate_isoforms(cancer_type):
    """
    :param cancer_type: The cancer type that should be walked through for
    concatenation
    :type cancer_type: str
    :return: None, outputs concatenated files within the same directory the
    individual isoform files are found and a single full file for all within
    a cancer type
    """
    cancer_dir = path.join(dataDir, profilesName, now, cancer_type)
    cancerProfile = open(path.join(cancer_dir, cancer_type + ".prof"), 'w')

    for (dirpath, dirnames, filenames) in walk(cancer_dir):
        # Open the concatenation file for writing
        with open(str(path.join(dirpath,
                                basename(dirpath) + ".prof")),
                  'w') as outfile:

            for fname in filenames:
                # Check to make sure the file has an isoform number before
                # reading it
                if search("\w+\.\d+\.\w+", fname):
                    with open(path.join(dirpath, fname), 'r') as infile:
                        isoformProfile = open(path.join(dataDir,
                                                        profilesName,
                                                        now,
                                                        isoformsSubDirName,
                                                        fname),
                                              'a')
                        for line in infile:
                            outfile.write(line)
                            cancerProfile.write(line)
                            isoformProfile.write(line)

                        # Close the isoform profile after appending all lines to it
                        isoformProfile.close()

    # Be sure to close the whole cancer profile
    cancerProfile.close()


if __name__ == "__main__":
    # Run the CLI wrapper to change global variables
    main()

    for ctype in cancerTypes:
        # Create the CANCER root or clear the CANCER root
        cancer_dir = path.join(dataDir, profilesName, now, ctype)
        if not path.exists(cancer_dir):
            makedirs(cancer_dir)
            del cancer_dir
        else:
            rmtree(cancer_dir)
            makedirs(cancer_dir)
            del cancer_dir

        # Create a Pool with a life of 100 tasks each before replacement
        if cpu_count() < 16:
            # Set processes to size cpu_count(), local workaround
            pool = Pool(maxtasksperchild=100)
        else:
            # Set processes size to 16 directly, remote workaround
            pool = Pool(maxtasksperchild=100, processes=16)

        # Runs the function once per worker on the next available pair in the
        # dataset
        pool.map(create_csv_profile, generate_data_pairs(ctype))

        # Close the Pool
        pool.close()

        # Walk through the directory for each type in cancerTypes, concatenating
        # isoform files
        concatenate_isoforms(ctype)
