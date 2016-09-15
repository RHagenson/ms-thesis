#!/usr/bin/python

# Name: Ryan Hagenson
# Email: rhagenson@unomaha.edu

import csv
import multiprocessing
import re
import shutil
import sys
from distutils.dir_util import mkpath
from getopt import GetoptError, getopt
from os import path, makedirs, listdir, remove

# Global variables
dataDir = False  # Default False, should be overwritten at CLI
allMAFsName = "allMAFs"  # The name of the allMAFs dir in dataDir
allMutsName = "allMuts"  # The name of the allMuts dir in dataDir
cdsName = "cds"  # The name of the cds dir in dataDir
pfamName = "pfam30.0"  # The name of the pfam dir in dataDir
refSeqName = "refSeq"  # The name of the refSeq dir in dataDir
profilesName = "profiles"  # The name of the final mutation profile csv's dir


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

# ./profiles will have subdirectories based on cancer type and gene


def main():
    """
    A simple wrapper for all CLI options
    """

    # Enables command-line options via getopt and sys packages
    try:
        opts, args = getopt(sys.argv[1:], 'd:')
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
        if opt == "-d":  # Set high-level data directory location
            global dataDir
            dataDir = arg

            # Create/empty profiles directory
            profileDir = path.join(dataDir, profilesName)
            if not path.exists(profileDir):
                makedirs(profileDir)
                del profileDir
            else:
                shutil.rmtree(profileDir)
                makedirs(profileDir)
                del profileDir


def create_csv_profile((MutFile, LongShortFile)):
    """
    Run by Pool.map() with data from generate_data_pairs()
    :type MutFile: str
    :type LongShortFile: str
    :return:
    """
    global dataDir, allMAFsName, allMutsName, \
        refSeqName, cdsName, pfamName, profilesName

    LongShortRE = re.compile('\s+(\d+)\s+(\w+)\s+(.+)')

    # Generate profiles directory tree with each cancer type and gene id
    cancerType = re.search("(\w+)\_.+\.txt", MutFile).group(1)
    geneName = re.search("(\w+)\.\d+\.[long|short]+", LongShortFile).group(1)
    fullPath = path.join(dataDir,
                         profilesName,
                         cancerType,
                         geneName)
    if not path.exists(fullPath):
        mkpath(fullPath)

    # Open each of the cancer type, gene id, and individual profile files
    cancerFile = open(path.join(dataDir,
                                profilesName,
                                cancerType,
                                cancerType + ".prof"), "a")
    cancerCSV = csv.writer(cancerFile, delimiter='\t')
    geneFile = open(path.join(dataDir,
                              profilesName,
                              cancerType,
                              geneName,
                              geneName + ".prof"), "a")
    geneCSV = csv.writer(geneFile, delimiter='\t')
    profileFile = open(path.join(dataDir,
                                 profilesName,
                                 cancerType,
                                 geneName,
                                 LongShortFile + ".prof"), 'w')
    profileCSV = csv.writer(profileFile, delimiter='\t')

    # Build and open the allMuts file
    MutFileHandle = open(path.join(dataDir, allMutsName, MutFile), 'r')

    # Build and open the long or short file based on extension
    if '.long' in LongShortFile:
        LongShortFileHandle = open(path.join(dataDir, refSeqName,
                                             "iupredLong", LongShortFile), 'r')
    elif '.short' in LongShortFile:
        LongShortFileHandle = open(path.join(dataDir, refSeqName,
                                             "iupredShort", LongShortFile), 'r')
    else:
        # Break if a non LongShort file is found
        sys.exit(2)

    # Extract the mutations from the allMuts file
    IsoformName, LongShort = path.splitext(LongShortFile)
    print "Processing " + str(LongShortFile)  # Inform user what is being done
    MutCSV = csv.reader(MutFileHandle, delimiter='\t')
    mutations = {}  # Should be {pos# : count}
    for row in MutCSV:
        # If we have found our mutations and there are no more for this
        # isoform break the loop, depends on isoform lines being sequential
        if len(mutations) > 0:
            if row[0] != IsoformName:
                break

        # Found a mutation in the protein
        if row[0] == IsoformName:
            # If the pos has a mutation, iterate or initialize it
            if row[6] in mutations:
                mutations[row[5]] += 1
            else:
                mutations[row[5]] = 1

    # Gather the rest of the information from LongShort file
    for line in LongShortFileHandle:
        # Skip comment lines at start
        if line.startswith('#'):
            continue

        LongShortMatch = re.search(LongShortRE, line)
        if LongShortMatch:
            posMuts = 0
            if LongShortMatch.group(1) in mutations:
                posMuts = mutations[LongShortMatch.group(1)]

            # Output the individual isoform results to each of the three
            # pertinent files. cancerCSV and geneCSV are appending, while
            # profileCSV is writing.
            cancerCSV.writerow([LongShortMatch.group(1),
                                LongShortMatch.group(2),
                                LongShortMatch.group(3),
                                posMuts])
            geneCSV.writerow([LongShortMatch.group(1),
                              LongShortMatch.group(2),
                              LongShortMatch.group(3),
                              posMuts])
            profileCSV.writerow([LongShortMatch.group(1),
                                 LongShortMatch.group(2),
                                 LongShortMatch.group(3),
                                 posMuts])

    # Be sure to release the files for other workers to take over control of
    cancerFile.close()
    geneFile.close()
    profileFile.close()


def generate_data_pairs():
    """
    Builds an iterable list of data pairs for Pool.map()
    Reads each file in data/allMuts/, for each line it determines if that
    protein has a corresponding file in iupredLong|iupredShort, if it does it
    adds a new entry in datapairs in style ['<allMuts filename>',
    '<iupredLong prop.XXX>.long']
    :return:
    """
    global dataDir, allMAFsName, allMutsName, \
        refSeqName, cdsName, pfamName

    print "Generating data pairs"

    # Define the data file pairs, allMutsFile should be found in
    # dataDir/allMuts/, while iupredFile should be found in either
    # dataDir/refSeq/iupredLong/ or dataDir/refSeq/iupredShort/, which path
    # is used is determined by the file extension (.long or .short)
    allMutsFile = ""
    iupredFile = ""

    # Define the return iterable
    # Entries should be in form [[allMutsFile, iupredFile]]
    datapairs = []

    # Collect the files in allMuts with abs pathing
    MutFilepaths = []
    MutLoc = path.join(dataDir, allMutsName)
    for f in listdir(MutLoc):
        MutFilepaths.append(path.join(MutLoc, f))

    # Process each Mut file in turn
    for mutFile in MutFilepaths:
        print "Now processing: " + mutFile

        mutName = mutFile.lstrip(MutLoc)
        with open(mutFile, 'r') as FILE:
            CSV = csv.reader(FILE, delimiter='\t')
            past_isoform = ""
            for row in CSV:
                proteinIsoform = row[0]

                # Ensure each isoform only gets processed once by skipping
                # repetitions
                if proteinIsoform == past_isoform:
                    continue
                else:
                    past_isoform = proteinIsoform

                # Check if a file exists
                longPath = path.join(dataDir, refSeqName,
                                     "iupredLong", proteinIsoform + ".long")
                shortPath = path.join(dataDir, refSeqName,
                                      "iupredShort", proteinIsoform + ".short")

                if path.exists(longPath):
                    datapairs.append([mutName, proteinIsoform + ".long"])
                elif path.exists(shortPath):
                    datapairs.append([mutName, proteinIsoform + ".short"])

    return datapairs


if __name__ == "__main__":
    main()

    # Eventual multiprocessing code
    pool = multiprocessing.Pool(maxtasksperchild=100)  # Returns Pool of size
    # os.cpu_count()

    # Runs the function once per worker on the next available pair in the
    # dataset

    pool.map(create_csv_profile, generate_data_pairs())
