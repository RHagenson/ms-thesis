# A functional wrapper to determine the length of a given isoform

isoform_length <- function(isoformName) {
  len <- numeric()
  
  # The assumed directory tree is based off running script in ms-thesis/R/
  assumedDataDir = "../../disorderCancer/data"
  assumedRefSeqDir <- paste(assumedDataDir, "refSeq/", sep="/")
  assumedIupredLong <- paste(assumedRefSeqDir, "iupredLong", sep="/")
  assumedIupredShort <- paste(assumedRefSeqDir, "iupredShort", sep="/")
  
  # Build location of file
  if (grepl("long", isoformName)) {
    isoformFile <- paste(assumedIupredLong, isoformName, sep="/")
  } else {
    isoformFile <- paste(assumedIupredShort, isoformName, sep="/")
  }
  
  # Read in the file
  FILE <- read.table(isoformFile)
  
  # Get the last line
  lastLine <- tail(FILE, 1)
  
  # Retrieve the length from the lastLine
  len <- lastLine$V1
  
  return(len)
}