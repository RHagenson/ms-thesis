# This function is intended to take the name of a cancer type, isoform name (with ".long" or .short"); 
# returning a data.frame of Position, AA, and Disorder

extract_disorder <- function(cancer_type, isoform_name) {
  # The assumed directory tree is based off running script in ms-thesis/R/
  assumedDataDir = "../../disorderCancer/data"
  assumedRefSeqDir <- paste(assumedDataDir, "refSeq/", sep="/")
  assumedIupredLong <- paste(assumedRefSeqDir, "iupredLong", sep="/")
  assumedIupredShort <- paste(assumedRefSeqDir, "iupredShort", sep="/")
  
  # Build location of file
  if (grepl("long", isoform_name)) {
    isoformFile <- paste(assumedIupredLong, isoform_name, sep="/")
  } else {
    isoformFile <- paste(assumedIupredShort, isoform_name, sep="/")
  }
  
  # read.table should ignore the commented header lines
  FILE.frame <- read.table(isoformFile)
  
  # Apply names to columns
  colnames(FILE.frame) <- c("Position", "AA", "Disorder")
  
  return(FILE.frame)
}