#!/usr/bin/python

################################################################################
# extractAnnotations.py                                                        #
# Author:  Dario Ghersi                                                        #
# Version: 20110202                                                            #
# Goal:    parse a Gene Ontology 'gene_association' file and prints either:    #
#          1. all the annotations (useful for computing the Information        #
#             Content of terms)                                                #
#          2. selected annotations if a file with a list of genes is supplied  #
#                                                                              #
# Usage:   ./extractAnnotations.py gene_association namespace IEA,ND,RCA,IPI   #
#                                  idType [genes list file]                    #
#          where the first argument is the annotation file (as downloaded from #
#          http://www.geneontology.org/), the second is the blacklist for the  #
#          evidence codes that should be excluded and the third optional       #
#          argument is a list of genes                                         #
# N.B.:    the format of the output is as follows:                             #
#          ID TERM1 TERM2 ...                                                  #
################################################################################

import sys

################################################################################
# CONSTANTS                                                                    #
################################################################################

pos = {"id": 1, "symbol": 2, "qualifier": 3, "term": 4, "evidence": 6,
       "namespace": 8, "name": 10}
allowedIdType = ["id", "symbol"]

################################################################################
# FUNCTIONS                                                                    #
################################################################################

def skipHeader(infile):
  """
  skip the lines beginning with a '!'
  """
  ## find where the header ends
  counter = 0
  while infile.readline().startswith("!"):
    counter += 1

  ## reposition the file iterator
  infile.seek(0)
  for i in range(0, counter):
    infile.readline()

  return infile

################################################################################

def parseGeneList(genesListFileName):
  """
  store the selected genes in a list
  """
  genesList = []

  ## parse the geneList file
  genesListFile = open(genesListFileName)
  for line in genesListFile:
    genesList.append(line[:-1])
  genesListFile.close()

  return genesList

################################################################################

def parseAnnFile(goAnnFileName, namespace, blacklist, idType):
  """
  parse the annotations file and return two dictionaries, one with the
  annotations assigned to the primary IDs and the other with the gene names
  assigned to the primary IDs
  """
  ann = {}
  dictID = {}
  goAnnFile = open(goAnnFileName, "r")
  goAnnFile = skipHeader(goAnnFile) # skip the header
  for line in goAnnFile:
    fields = line.split("\t")
    if fields[pos["namespace"]] == namespace and not fields[pos["evidence"]] in\
       blacklist and fields[pos["qualifier"]] != "NOT":
      if ann.has_key(fields[pos[idType]]):
        if not fields[pos["term"]] in ann[fields[pos[idType]]]:
          ann[fields[pos[idType]]].append(fields[pos["term"]])
      else:
        ann[fields[pos[idType]]] = [fields[pos["term"]]]
      if not dictID.has_key(fields[pos[idType]]):
        dictID[fields[pos[idType]]] = fields[pos["name"]]

  ## close the file
  goAnnFile.close()

  return [ann, dictID]

################################################################################
# MAIN PROGRAM                                                                 #
################################################################################

## parse the parameters
if len(sys.argv) < 5:
  print "Usage: ./extractAnnotations.py gene_association namespace IEA,ND,RCA,IPI idType [genes list file]"
  sys.exit(1)
goAnnFileName, namespace, blacklist, idType = sys.argv[1:5]
if len(sys.argv) == 6:
  genesListFileName = sys.argv[5]
  onlySelected = True
  genesList = parseGeneList(genesListFileName)
else:
  onlySelected = False

## make sure idType is one of ["id", "symbol"]
if not idType in allowedIdType:
  print "idType can only be one of: [" + ", ".join(allowedIdType) + "]"
  sys.exit(1)

## process the blacklist
blacklist = blacklist.split(",")

## parse the annotations file
ann, dictID = parseAnnFile(goAnnFileName, namespace, blacklist, idType)

## print the results (with the name of the gene)
for id in ann:
  if onlySelected:
    printAnn = False
    genesNames = dictID[id].split("|")
    for gene in genesNames:
      if gene in genesList:
        printAnn = gene
        break
    if printAnn:
      print gene + "\t" + "\t".join(ann[id])
  else:
    print id + "\t" + "\t".join(ann[id])
