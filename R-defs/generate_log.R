# This function takes a profile, number of samples to take, location of profile, where figs should be placed
# where output LOGs should be placed, and a pValCut off and generates a LOG and  PDF

source("R-defs/build_plot.R")  # Generate a plot if significant

generate_log <- function(filename, number, profileDir, figsDir, outputDir, pValCut=0.05, plot=FALSE) {
  # Remove trailing '/' because paste() below adds it when creating path
  # profilesDir <- sub(pattern = "/$", replacement = "", profilesDir)
  
  # Generate log location and file
  logTree <- sub(pattern = "^.+/profiles/", replacement = "", profileDir)
  CSVPath <- paste(outputDir, 
                   logTree, 
                   sub(pattern="*.prof", replacement = "", filename),
                   sep="/")
  if(! dir.exists(CSVPath)) {
    dir.create(CSVPath, recursive = TRUE)
  }
  CSV <- "LOG.csv"
  
  # Generate PDF plot location
  figsTree <- sub(pattern = "^.+/profiles/", replacement = "", profileDir)
  pdfPath <- paste(figsDir, figsTree, sep = "/")
  
  # The parameters that will eventual change
  filename <- filename # The profile being processed, e.g. "MUC16.001.long.prof"
  
  N = as.numeric(number) # 100000  # The number of samples to take
  
  profileDir = profileDir  # The location of where the filename is found
                          # Ex ~/Thesis/data/profiles/isoforms/
  figsDir = figsDir  # The location of where figures should be written, do not ending forget '/'
  
  # Read in the source profile file
  path <- paste(profileDir, filename, sep = "/")
  
  # Try to read in the profile file, if it does not exist
  FILE <- tryCatch(read.delim(path, header = FALSE), error=function(e) NULL)

  # Create a vector to hold the sum of random mutations
  normalVector = vector(mode = "double")
  
  # Determine how many mutations are present and the observed disorder score of mutations
  realLevel = as.numeric(sum(FILE$V3 * FILE$V4))
  numMutations = as.numeric(sum(FILE$V4))
  
  for (i in 1:N) {
    # Print a helpful message to the user for what is being done periodically
    if ((i %% 1000) == 0) {
      print(paste("On sample number", as.character(i), "for", filename))
    }
    
    normalVector <-
      append(normalVector, 
             round(sum(sample(FILE$V3, replace = TRUE, size = numMutations)), 
                   digits = 4))
  }
  
  # Determine average disorder score from normal curve
  avgDisorder = round(sum(normalVector) / length(normalVector), digits=3)
  
  # Generate a data frame with values and freq as percent
  frame <-
    as.data.frame(table(normalVector) / length(normalVector) * 100)

  # Calculate the empirical pValue and directionality
  lessPValue <-  sum(realLevel <= normalVector) / length(normalVector)
  morePValue <- sum(realLevel >= normalVector) / length(normalVector)
  if (morePValue < lessPValue) {
    pValue <- morePValue
    pValDirection <- "+"  # The observed is above the expected average
  } else {
    pValue <- lessPValue
    pValDirection <- "-"  # The observed is below the expected average
  }
  
  # Add column names
  colnames(frame) <- c("TotalDisOrderScore", "PercentFreq")
  
  # Coerce factor from use of table() into numeric
  frame$TotalDisOrderScore <-
    as.numeric(levels(frame$TotalDisOrderScore))[frame$TotalDisOrderScore]
  
  # No matter the pValue add the isoform to the results CSV
  # Columns should be in the following order
  #   1. isoform name
  #   2. observed disorder score
  #   3. average random disorder score
  #   4. total number of mutations
  #   5. empirical p-value
  #   6. Direction of p-value, '+' meaning the real, observed level is above the average disorder
  write.table(x=data.frame(sub(".prof", "", filename), 
                           realLevel, 
                           avgDisorder, 
                           numMutations, 
                           pValue,
                           pValDirection), 
              file=paste(CSVPath, CSV, sep="/"), append = TRUE, row.names = FALSE,
              quote = FALSE, sep=",", col.names = FALSE)
  
  # Produce a PDF plot if the pValue is below the threshold and a plot is desired
  if (pValue <= pValCut) {
    if (plot) {
      build_plot(pdfPath = pdfPath, normalVector = normalVector, filename = filename,
                 pValue = pValue, pValDirection = pValDirection, avgDisorder = avgDisorder,
                 realLevel = realLevel, numMutations = numMutations)
    }
  }
}
