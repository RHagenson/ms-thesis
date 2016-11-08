# This script takes the results from every isoform within a cancer and provides the Bonferroni
# corrected p-values. The best isoform from each gene is taken as the representative for that gene.
# Best is defined as taking the top isoform when sorted by:
#     decreasing number of mutations, increasing length of isoform, and progressing numerically/alphabetically
# Only one isoform from each gene is taken because Bonferroni correction assumes independent values, while
# isoforms of the same gene are dependent values.

source("R-defs/isoform_length.R")

correction <- function(date, cancerType) {
  # Select global variables
  pValCutoff = 0.1
  preferredDirection = "+" # Should correspond to which pValues represent higher disorder than average in LOG.csv files
  correctionMethod = "fdr" # Options are below:
  # p.adjust.methods
  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  #   "fdr", "none")
  
  # Add cancerType onto logsDir to build the cancer-specific directory
  logsDir <- paste("./outputs", date, cancerType, sep = "/")
  
  # Remove any double slashes
  logsDir <- gsub(pattern = "//", replacement = "/", x = logsDir)
  
  ###
  ### Begin creating long and short file vectors
  ###
  # Ensure a date is included in correct format
  if (grepl("[0-9]{2}-[0-9]{2}-[0-9]{2}", logsDir)) {
    # Do nothing for now
  } else {
    stop(paste(
      "The date provided:", date, "does not appear to be in form DD-MM-YY"
    ), call. = FALSE)
  }
  
  # Gather all the directories in the cancer type subdirectory (all gene directories)
  logsFilesVector <- list.dirs(path = logsDir,
                               full.names = TRUE,
                               recursive = FALSE)
  
  ###
  ### Begin determining the top isoform in each gene
  ### Reminder: this is done by number of mutations, length, then numerically/alphabetically
  ###
  selectIsoformsLong <-
    data.frame()  # The data.frame for the "best" isoforms found in long files
  selectIsoformsShort <-
    data.frame()  # The data.frame for the "best" isoforms found in short files
  
  # Explore each directory and select the top isoform from each for both long and short files
  # These should be the same isoform, but they are processed separately.
  for (directory in logsFilesVector) {
    longFiles <- vector("character")
    shortFiles <- vector("character")
    longDataFrame <- data.frame()
    shortDataFrame <- data.frame()
    
    # Get only the individual isoform files
    for (file in list.files(
      directory, pattern = "LOG.csv", full.names = TRUE, recursive = TRUE
    )) {
      if (grepl("long", file)) {
        longFiles <- c(longFiles, file)
      } else if (grepl("short", file)) {
        shortFiles <- c(shortFiles, file)
      }
    } # End files for loop
    
    # Create a dataframe for easy mutliple column sorting, can be combined into one by isoform name only
    # absent of .long or .short file name
    longDataFrame <- as.data.frame(longFiles)
    shortDataFrame <- as.data.frame(shortFiles)
    
    # Find the top longFile and add its p-value to the best vector
    if (length(longFiles) > 0) {
      # Set number of mutations column via entire longFiles vector
      longDataFrame$numMuts <-
        apply(longDataFrame, 1, function(row)
          read.table(as.character(row["longFiles"]), sep = ",")$V4)
      
      # Set the pValue of each entry in the longFiles vector
      longDataFrame$pVal <-
        apply(longDataFrame, 1, function(row)
          read.table(as.character(row["longFiles"]), sep = ",")$V5)
      
      # Set the isoName by separating it from the full file path
      longDataFrame$isoName <-
        apply(longDataFrame, 1, function(row)
          basename(dirname(as.character(row["longFiles"]))))
      
      # Set length of each isofom by running isoform_length()
      longDataFrame$isoLength <-
        apply(longDataFrame, 1, function(row)
          isoform_length(row["isoName"]))
      
      # Set direction of p-value
      longDataFrame$pValDir <-
        apply(longDataFrame, 1, function(row)
          read.table(as.character(row["longFiles"]), sep = ",")$V6)
      
      # Gives the preferred direction precedence in sorting if it is present
      if (preferredDirection %in% longDataFrame$pValDir) {
        longDataFrame$pValDir <-
          relevel(longDataFrame$pValDir, preferredDirection)
      }
      
      # Add top entry/row to selectIsoformsLong data.frame
      selectIsoformsLong <- rbind(selectIsoformsLong,
                                  longDataFrame[with(longDataFrame,
                                                     order(-numMuts, isoLength, isoName)), ][1, ])
      
    } # End finding the top longFile and adding it p-value to the best vector
    
    # Find the top shortFile and add its p-value to the best vector
    if (length(shortFiles) > 0) {
      # Set number of mutations column via entire shortFiles vector
      shortDataFrame$numMuts <-
        apply(shortDataFrame, 1, function(row)
          read.table(as.character(row["shortFiles"]), sep = ",")$V4)
      
      # Set the pValue of each entry in the shortFiles vector
      shortDataFrame$pVal <-
        apply(shortDataFrame, 1, function(row)
          read.table(as.character(row["shortFiles"]), sep = ",")$V5)
      
      # Set the isoName by separating it from the full file path
      shortDataFrame$isoName <-
        apply(shortDataFrame, 1, function(row)
          basename(dirname(as.character(row["shortFiles"]))))
      
      # Set length of each isofom by running isoform_length()
      shortDataFrame$isoLength <-
        apply(shortDataFrame, 1, function(row)
          isoform_length(row["isoName"]))
      
      # Set direction of p-value
      shortDataFrame$pValDir <-
        apply(shortDataFrame, 1, function(row)
          read.table(as.character(row["shortFiles"]), sep = ",")$V6)
      
      # Gives '+' direction precedence in sorting if it is present
      if (preferredDirection %in% shortDataFrame$pValDir) {
        shortDataFrame$pValDir <-
          relevel(shortDataFrame$pValDir, preferredDirection)
      }
      
      # Add top entry/row to selectIsoformsShort data.frame
      selectIsoformsShort <- rbind(selectIsoformsShort,
                                   shortDataFrame[with(shortDataFrame,
                                                       order(-numMuts, isoLength, isoName)),][1,])
    } # End finding the top shortFile and adding it p-value to the best vector
  } # End directory for loop
  
  # Output to the user
  print("Computing the adjusted long p-values")
  selectIsoformsLong$pValAdj <- p.adjust(apply(selectIsoformsLong,
                                            1,
                                            function(row)
                                              if (row['pValDir'] != preferredDirection) {
                                                1 - as.numeric(row['pVal'])
                                              }else{
                                                as.numeric(row['pVal'])
                                              }),
                                      method = correctionMethod)
  
  print("Computing the adjusted short p-values")
  selectIsoformsShort$pValAdj <- p.adjust(apply(selectIsoformsShort,
                                             1,
                                             function(row)
                                               if (row['pValDir'] != preferredDirection) {
                                                 1 - as.numeric(row['pVal'])
                                               }else{
                                                 as.numeric(row['pVal'])
                                               }),
                                       method = correctionMethod)
  
  
  ### This function should return the two data.frames in a list/vector/etc
  ### Once that is done, the print/output below should be extracted to the
  ### top-level script directory alongside a log of each data.frame by type
  ### (long v short) and by method used to generate the corrections, hence
  ### the function should return the method used as well or have the method be
  ### a required parameter.
  
  return(list("long" = selectIsoformsLong, 
              "short" = selectIsoformsShort))
  
  # Create a directory for the signficant results
#   outDir <- paste("./outputs", date, "p-adjusted", cancerType, sep = "/")
#   dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
#   
#   # Output the significant results
#   # Write the two data.frames to file
# #   write.csv(selectIsoformsLong, file = paste(outDir, "long_log.csv", sep = "/"))
# #   write.csv(selectIsoformsShort, file = paste(outDir, "short_log.csv", sep = "/"))
#   
#   # Redirect output to file and STOUT
#   sink(paste(outDir, "adjResults.txt", sep = "/"), split = TRUE)
#   
#   # Output the significant results for long files
#   print(paste("The number of significant long values is:", as.character(sum(selectIsoformsLong$pValAdj < pValCutoff))))
#   print("Value(s):")
#   print(as.character(selectIsoformsLong$pVal[which(selectIsoformsLong$pValAdj < pValCutoff)]))
#   print("Related isoform:")
#   print(as.character(selectIsoformsLong$isoName[which(selectIsoformsLong$pValAdj < pValCutoff)]))
#   print("Direction of significance:")
#   print(as.character(selectIsoformsLong$pValDir[which(selectIsoformsLong$pValAdj < pValCutoff)]))
#   
#   # Output the significant results for short files
#   print(paste("The number of significant short values is:", as.character(sum(selectIsoformsShort$pValAdj < pValCutoff))))
#   print("Value(s):")
#   print(as.character(selectIsoformsShort$pVal[which(selectIsoformsShort$pValAdj < pValCutoff)]))
#   print("Related isoform:")
#   print(as.character(selectIsoformsShort$isoName[which(selectIsoformsShort$pValAdj < pValCutoff)]))
#   print("Direction of significance:")
#   print(as.character(selectIsoformsShort$pValDir[which(selectIsoformsShort$pValAdj < pValCutoff)]))
}