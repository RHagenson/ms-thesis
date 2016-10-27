# This script takes the results from every isoform within a cancer and provides the Bonferroni 
# corrected p-values. The best isoform from each gene is taken as the representative for that gene. 
# Best is defined as taking the top isoform when sorted by:
#     decreasing number of mutations, decreasing length of isoform, and increasing numerically/alphabetically
# Only one isoform from each gene is taken because Bonferroni correction assumes independent values, while
# isoforms of the same gene are dependent values.


correction <- function(logsDir, cancerType) {
  # Add cancerType onto logsDir to build the cancer-specific director
  logsDir <- paste(logsDir, cancerType, sep="/")
  
  # Remove any double slashes
  logsDir <- gsub("//", "/", logsDir)
  
  ###
  ### Begin creating long and short file vectors
  ###
  # Ensure a date is included
  if(grepl("[0-9]{2}-[0-9]{2}-[0-9]{2}", logsDir)) {
    # Do nothing for now
  } else {
    stop(paste("The logsDir provided:", logsDir, "does not appear to end with a date"), call.=FALSE)
  }
  
  # Gather all the LOG files into one vector
  logsFilesVector <- dir(path = logsDir, 
                         pattern = "LOG.csv", 
                         full.names = TRUE, 
                         recursive = TRUE)
  logsFilesVectorLong <- vector("character")
  logsFilesVectorShort <- vector("character")
  
  # Loop through and separate into long and short vector, and remove concatenated files
  for(logF in logsFilesVector) {
    # Ensure it is an isoform file LOG
    if(grepl("[0-9]{3}.[long]|[short]", logF)) {
      # Divide into long or short
      if(grepl("[0-9]{3}.[long]", logF)) {
        append(logsFilesVectorLong, logF)
      } else {
        append(logsFilesVectorShort, logF)
      }
    }
  }
  
  rm(logsFilesVector) # Clear memory
  
  ###
  ### Begin determining the top isoform in each gene
  ### Reminder: this is done by number of mutations, length, then numerically
  ###
  
  
  
}