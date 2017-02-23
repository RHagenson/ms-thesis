# A functional wrapper to determine the length of a given isoform
# Example input: A1CF.004

isoform_length <- function(isoformName) {
  len <- numeric()
  
  # The assumed directory tree is based off running script in ms-thesis/R/
  assumedDataDir = "../../disorderCancer/data"
  assumedRefSeqDir <- paste(assumedDataDir, "refSeq/", sep="/")
  
  # Build location of file
  isoformFile <- paste(assumedRefSeqDir, paste0(isoformName, ".fasta"),
                       sep="/")
  
  # Read in the fasta file, removing the first line
  file_content <- readLines(isoformFile)[-1]
  
  # Retrieve the length by collapsing the character vector and counting the length
  len <- nchar(paste0(file_content, collapse = ""))
  
  return(len)
}