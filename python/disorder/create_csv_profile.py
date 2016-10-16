#!/usr/bin/python

# Name: Ryan Hagenson
# Email: rhagenson@unomaha.edu

import sys
from csv import reader, writer
from distutils.dir_util import mkpath
from getopt import GetoptError, getopt
from multiprocessing import Pool, cpu_count
from os import path, makedirs, listdir
from re import search, compile
from shutil import rmtree
import datetime

# Global variables
dataDir = False  # Default False, should be overwritten at CLI
allMAFsName = "allMAFs"  # The name of the allMAFs dir in dataDir
allMutsName = "allMuts"  # The name of the allMuts dir in dataDir
cdsName = "cds"  # The name of the cds dir in dataDir
pfamName = "pfam30.0"  # The name of the pfam dir in dataDir
refSeqName = "refSeq"  # The name of the refSeq dir in dataDir
profilesName = "profiles"  # The name of the final mutation profile csv's dir
isoformsSubDirName = "isoforms"

cancerTypes = ['BRCA']
now = datetime.datetime.now().strftime("%d-%m-%y")

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
        opts, args = getopt(sys.argv[1:], 'd:c:', ["date=", "dataDir=", "cancerTypes="])
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
        mut_file_handle = open(path.join(dataDir, allMutsName, mut_file), 'r')
    except IOError as e:
        print(str(e))
        return  # Return if the file does not exists, therefore the datapair is invalid

    # Build and open the long or short file based on extension
    try:
        if '.long' in long_short_file:
            long_short_file_handle = open(path.join(dataDir, refSeqName,
                                             "iupredLong", long_short_file), 'r')
        elif '.short' in long_short_file:
            long_short_file_handle = open(path.join(dataDir, refSeqName,
                                             "iupredShort", long_short_file), 'r')
        else:
            # Break if a non long/short file is found
            sys.exit(2)
    except IOError as e:
        print(str(e))
        return  # Return if the file does not exists, therefore the datapair is invalid

    # Generate profiles directory tree with each cancer type and gene id
    cancer_type = search("(\w+)\_.+\.txt", mut_file).group(1)
    geneMatch = search("([\w|-]+)+\.\d+\.([long|short]+)",
                       long_short_file)
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

    # Open each of the cancer type, gene id, isoform,
    # and individual profile files
    # Done as two steps each:
    # 1) opening a file handle, and
    # 2) passing the handle to csv.writer() to allow directed closing of the
    # file for guaranteed resource release for other workers.
    #
    #
    # cancer_file is a profile for the combination of all genes and isoforms
    # found within that cancer
    cancer_file = open(path.join(dataDir,
                                 profilesName,
                                 now,
                                 cancer_type,
                                 cancer_type + ".prof"), "a")
    cancer_csv = writer(cancer_file, delimiter='\t')

    # gene_file is a profile for each gene that combines all its isoforms
    gene_file = open(path.join(dataDir,
                               profilesName,
                               now,
                               cancer_type,
                               gene_name,
                               gene_name + ".prof"), "a")
    gene_csv = writer(gene_file, delimiter='\t')

    # isoform_file is a profile for each isoform, independent of cancer type
    isoform_file = open(path.join(dataDir,
                                  profilesName,
                                  now,
                                  isoformsSubDirName,
                                  long_short_file + ".prof"), "a")
    isoform_csv = writer(isoform_file, delimiter='\t')

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
    mut_csv = reader(mut_file_handle, delimiter='\t')
    mutations = {}  # Should be {pos# : count}
    for row in mut_csv:
        # If we have found our mutations and there are no more for this
        # isoform break the loop, depends on isoform lines being sequential
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

            # Output the individual isoform results to each of the
            # pertinent files. cancer_csv, gene_csv, and isoform_csv
            # are appending/non-unique, while profile_csv is writing/unique.
            cancer_csv.writerow([long_short_match.group(1),
                                 long_short_match.group(2),
                                 long_short_match.group(3),
                                 pos_muts])
            gene_csv.writerow([long_short_match.group(1),
                               long_short_match.group(2),
                               long_short_match.group(3),
                               pos_muts])
            isoform_csv.writerow([long_short_match.group(1),
                                  long_short_match.group(2),
                                  long_short_match.group(3),
                                  pos_muts])
            profile_csv.writerow([long_short_match.group(1),
                                  long_short_match.group(2),
                                  long_short_match.group(3),
                                  pos_muts])

    # Be sure to release the files for other workers to take over control of
    cancer_file.close()
    gene_file.close()
    isoform_file.close()
    profile_file.close()


def generate_data_pairs():
    """
    Builds an iterable list of data pairs for Pool.map()
    Reads each file in data/allMuts/, for each line it determines if that
    protein has a corresponding file in iupredLong|iupredShort, if it does it
    adds a new entry in datapairs in style ['<allMuts filename>',
    '<iupredLong prop.XXX>.long']
    """
    global dataDir, allMutsName, refSeqName, cancerTypes

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
        if cancerTypes:
            for cancer in cancerTypes:
                if cancer+"_" in f:
                    mut_filepaths.append(path.join(mut_loc, f))
        else:
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


if __name__ == "__main__":
    # Run the CLI wrapper
    main()

    # Create a Pool with a life of 100 tasks each before replacement
    if cpu_count() < 16:
        pool = Pool(maxtasksperchild=100)  # Set processes to size cpu_count(), local workaround
    else:
        pool = Pool(maxtasksperchild=100, processes=16)  # Set processes size to 16 directly, remote workaround

    # Runs the function once per worker on the next available pair in the
    # dataset
    pool.map(create_csv_profile, generate_data_pairs())

    # Close the Pool
    pool.close()
