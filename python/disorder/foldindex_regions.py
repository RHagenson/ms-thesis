#!/usr/bin/python

# Name: Ryan Hagenson
# Email: rhagenson@unomaha.edu

import sys
from getopt import GetoptError, getopt
from multiprocessing import Pool, cpu_count
from os import path, makedirs, listdir
from os.path import basename
import urllib2
import xml.etree.ElementTree as ET
import csv

dataDir = "../../../disorderCancer/data/"  # Default relative path from pwd/current dir
allMAFsName = "allMAFs"  # The name of the allMAFs dir in dataDir
allMutsName = "allMuts"  # The name of the allMuts dir in dataDir
cdsName = "cds"  # The name of the cds dir in dataDir
pfamName = "pfam30.0"  # The name of the pfam dir in dataDir
refSeqName = "refSeq"  # The name of the refSeq dir in dataDir
foldindexName = "foldindex"

fasta_directory = path.join(dataDir, refSeqName)
output_directory = path.join(dataDir, refSeqName, foldindexName)

foldindex_url = "http://bioportal.weizmann.ac.il/fldbin/findex?m=xml&sq="


# General directory tree within dataDir is:
# ./allMuts
# ./refSeq
# ./refSeq/iupredLong
# ./refSeq/iupredShort
# ./refSeq/hmmer
# ./refSeq/foldindex  # Should be made by script


def main():
    """
    A simple wrapper for all CLI options
    """
    global fasta_directory, output_directory

    # Enables command-line options via getopt and sys packages
    try:
        opts, args = getopt(sys.argv[1:],
                            'd:o:',
                            ["directory=", "output="]
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
        # Reassign fasta_directory
        if opt in ("-d", "--directory"):
            fasta_directory = str(arg)

        # Reassign output_directory
        if opt in ("-o", "--output"):
            output_directory = str(arg)

    # Recursively build output_directory path
    if not path.exists(output_directory):
        makedirs(output_directory)


def create_foldindex_file((gene_w_isoform_num, fasta_sequence)):
    """
    :arg gene_w_isoform_num: string of gene name with isoform number
    :arg fasta_sequence: string of gene FASTA sequence

    :type gene_w_isoform_num: str
    :type fasta_sequence: str

    Run by Pool.map() with data from generate_pairs()
    :return: file at output_directory/<gene_w_isoform_num>.csv
    """
    global output_directory, foldindex_url

    print("Now processing: " + gene_w_isoform_num)

    try:
        # Append the sequence and submit to server
        response = urllib2.urlopen(foldindex_url + fasta_sequence)

        # Read the server response (XML string)
        page = response.read()
    except urllib2.HTTPError as err:
        # Redirect STDERR to STDOUT (ensures screen display)
        sys.stdout = sys.stderr
        # Print help information
        print(str(err))
        sys.exit(2)

    # Parse XML string
    root = ET.fromstring(page)

    # Only process if 'segments' can be found and is not empty
    if root.find("segments").findall("segment"):
        print("Segments found")

        # Define output file
        foldindex_file = open(path.join(output_directory,
                                        gene_w_isoform_num + ".csv"), "w")
        foldindex_csv = csv.writer(foldindex_file)

        # Writer header row to CSV
        foldindex_csv.writerow(['Start',
                                'End',
                                'Length',
                                'Score',
                                'STD'])

        for segment in root.find("segments").findall("segment"):
            start = segment.get('start')
            end = segment.get('end')
            length = segment.get('len')
            score = segment.get('score')
            std = segment.get('std')

            # Define order of elements to match headers
            entry = [start, end, length, score, std]

            # Write entry to file
            foldindex_csv.writerow(entry)

        foldindex_file.close()
    else:
        print("No segments found in: " + gene_w_isoform_num)


def generate_pairs(fasta_dir):
    """
    :arg fasta_dir: a string representation of where the FASTA files are
    :type str

    Builds an iterable list of data pairs for use with Pool.map()
    Reads each file in fasta_directory and generates a:
        (<GENE.ISOFORM #>, <FASTA Sequence>)
    data pair
    """

    print("Generating data pairs")

    # Define the return iterable
    # Entries should be in form [[gene_w_iso_num, FASTA_sequence]]
    # Each with full absolute path
    datapairs = []

    # Collect fasta files with absolute path into fasta_filepaths
    fasta_filepaths = []
    for f in listdir(fasta_dir):
        # Collect only fasta files
        if ".fasta" in f:
            fasta_filepaths.append(path.join(fasta_dir, f))

    # Process each FASTA file in turn
    for fasta_file in fasta_filepaths:
        gene_w_iso_num = path.splitext(basename(fasta_file))[0]
        sequence = ""  # Need to concatenate lines to build seq

        print("Now processing: " + gene_w_iso_num)

        with open(fasta_file, 'r') as FILE:
            for line in FILE:
                if ">" in line:
                    continue
                else:
                    sequence += line.strip()

        datapairs.append([gene_w_iso_num, sequence])

    return datapairs


if __name__ == "__main__":
    # Run the CLI wrapper to change global variables
    main()

    # Run process in parallel via Pool.map()
    # Create a Pool with a life of 100 tasks each before replacement
    if cpu_count() < 16:
        # Set processes to size cpu_count()-1, local workaround
        pool = Pool(maxtasksperchild=100, processes=cpu_count()-1)
    else:
        # Set processes size to 16 directly, remote workaround
        pool = Pool(maxtasksperchild=100, processes=16)

    # Run the function pipeline once per entry in datapairs
    pool.map(create_foldindex_file,
             generate_pairs(fasta_dir=fasta_directory))

    # Close the Pool
    pool.close()
