# This script takes the results from every isoform within a cancer and provides the Bonferroni
# corrected p-values. The best isoform from each gene is taken as the representative for that gene.
# Best is defined as taking the top isoform when sorted by:
#     decreasing number of mutations, increasing length of isoform, and progressing numerically/alphabetically
# Only one isoform from each gene is taken because Bonferroni correction assumes independent values, while
# isoforms of the same gene are dependent values.

source("R-defs/isoform_length.R")

correction <- function(date,
                       cancerType,
                       method = c("holm",
                                  "hochberg",
                                  "hommel",
                                  "bonferroni",
                                  "BH",
                                  "BY",
                                  "fdr",
                                  "none"),
                       preferredDirection = "+") {
  # Select global variables
  pValCutoff = 0.1
  # Should correspond to which pValues represent higher disorder than average in LOG.csv files
  preferredDirection = preferredDirection
  
  if (length(method) > 1) {
    correctionMethod = "fdr"
  } else {
    correctionMethod = method
  }
  
  # Add cancerType onto logsDir to build the cancer-specific directory
  logsDir <- paste("./outputs", date, cancerType, sep = "/")
  
  # Remove any double slashes
  logsDir <- gsub(pattern = "//",
                  replacement = "/",
                  x = logsDir)
  
  ###
  ### Begin creating long and short file vectors
  ###
  # Ensure a date is included in correct format
  if (grepl("[0-9]{2}-[0-9]{2}-[0-9]{2}", logsDir)) {
    # Do nothing for now
  } else {
    stop(paste(
      "The date provided:",
      date,
      "does not appear to be in form DD-MM-YY"
    ),
    call. = FALSE)
  }
  
  # Gather all the directories in the cancer type subdirectory (all gene directories)
  logsFilesVector <- list.dirs(path = logsDir,
                               full.names = TRUE,
                               recursive = FALSE)
  
  # Approximate how many files will be the result in the end, assume most genes have at least 1 isoform
  SELECT_LEN <- length(list.files(logsDir)) / 2
  
  # Track empty rows
  L_SELECT_ROW <-
    1 # Used to track where the next "empty" row is in the select longs
  S_SELECT_ROW <-
    1 # Used to track where the next "empty" row is in the select shorts
  
  ###
  ### Begin determining the top isoform in each gene
  ### Reminder: this is done by number of mutations, length, then numerically/alphabetically
  ###
  # The data.frame for the "best" isoforms found in long files
  selectIsoformsLong <- data.frame(
    longFiles = character(SELECT_LEN),
    numMuts = numeric(SELECT_LEN),
    pVal = numeric(SELECT_LEN),
    isoName = character(SELECT_LEN),
    isoLength = numeric(SELECT_LEN),
    pValDir = character(SELECT_LEN),
    pValAdj = numeric(SELECT_LEN),
    stringsAsFactors=FALSE
  )
  
  # The data.frame for the "best" isoforms found in short files
  selectIsoformsShort <-
    data.frame(
      shortFiles = character(SELECT_LEN),
      numMuts = numeric(SELECT_LEN),
      pVal = numeric(SELECT_LEN),
      isoName = character(SELECT_LEN),
      isoLength = numeric(SELECT_LEN),
      pValDir = character(SELECT_LEN),
      pValAdj = numeric(SELECT_LEN),
      stringsAsFactors=FALSE
    )
  
  selectLongCache <- paste(tempdir(), 
                           paste(cancerType, ".long.tmp", sep=""), sep="/")
  selectShortCache <- paste(tempdir(), 
                            paste(cancerType, ".short.tmp", sep=""), sep="/")
  
  # Load in the past select isoforms if it has been cached
  if (file.exists(selectLongCache) &
      file.exists(selectShortCache)) {
    # Eventually this cache should just be the long/short files with full abs path 
    # then this part of the if/else will fill in the remaining values as they would be done below
    # For now it is the entire selectIsoforms data.frame object written and read as a csv file
    selectIsoformsLong <- read.table(selectLongCache)
    selectIsoformsShort <- read.table(selectShortCache)
  } else {
    # Explore each directory and select the top isoform from each for both long and short files
    # These should be the same isoform, but they are processed separately.
    for (directory in logsFilesVector) {
      longFiles <- vector("character")
      shortFiles <- vector("character")
      longDataFrame <- data.frame()
      shortDataFrame <- data.frame()
      
      # Get only the individual isoform files
      for (file in list.files(
        directory,
        pattern = "LOG.csv",
        full.names = TRUE,
        recursive = TRUE
      )) {
        if (grepl("\\.[0-9]{3}\\.long", file)) {
          longFiles <- c(longFiles, file)
        } else if (grepl("\\.[0-9]{3}\\.short", file)) {
          shortFiles <- c(shortFiles, file)
        }
      } # End files for loop
      
      # Create a dataframe for easy mutliple column sorting
      longDataFrame <- as.data.frame(
        longFiles,
        numMuts = numeric(length(longFiles)),
        pVal = numeric(length(longFiles)),
        isoName = character(length(longFiles)),
        isoLength = numeric(length(longFiles)),
        pValDir = character(length(longFiles)),
        pValAdj = numeric(length(longFiles)),
        stringsAsFactors=FALSE
      )
      shortDataFrame <- as.data.frame(
        shortFiles,
        numMuts = numeric(length(shortFiles)),
        pVal = numeric(length(shortFiles)),
        isoName = character(length(shortFiles)),
        isoLength = numeric(length(shortFiles)),
        pValDir = character(length(shortFiles)),
        pValAdj = numeric(length(longFiles)),
        stringsAsFactors=FALSE
      )
      
      # Find the top longFile and add its p-value to the best vector
      if (length(longFiles) > 0) {
        # Read files into memory
        HANDLE <- apply(longDataFrame,
                        1,
                        function(row)
                          read.table(as.character(row["longFiles"]),
                                     sep = ","))
        
        # Set number of mutations column via entire longFiles vector
        longDataFrame$numMuts <-
          lapply(HANDLE, function(row)
            as.numeric(as.character(row$V4)))
        
        # Set the pValue of each entry in the longFiles vector
        longDataFrame$pVal <-
          lapply(HANDLE, function(row)
            as.numeric(as.character(row$V5)))
        
        # Set the isoName by separating it from the full file path, .long is maintained for passing into
        # isoform_length()
        longDataFrame$isoName <- lapply(longFiles,
                                        function(full)
                                          basename(dirname(full)))
        
        
        # Set length of each isofom by running isoform_length()
        longDataFrame$isoLength <- apply(longDataFrame,
                                         1,
                                         function(row)
                                           as.numeric(as.character(isoform_length(row["isoName"]))))
        
        # Set direction of p-value
        longDataFrame$pValDir <-
          lapply(HANDLE, function(row)
            as.character(row$V6))
        
        # Gives the preferred direction precedence in sorting if it is present
        if (preferredDirection %in% longDataFrame$pValDir) {
          longDataFrame$pValDir <- relevel(factor(longDataFrame$pValDir,
                                                  levels = c("+", "-")),
                                           preferredDirection)
        }
        
        # Add top entry/row to selectIsoformsLong data.frame
        longDataFrame <- as.data.frame(lapply(longDataFrame, unlist), stringsAsFactors=FALSE)
        SELECT_ROW <- longDataFrame[with(longDataFrame,
                                         order(-numMuts, isoLength, isoName)),][1,]
        
        # Individually file values
        selectIsoformsLong[L_SELECT_ROW,]$longFiles <- SELECT_ROW$longFiles
        selectIsoformsLong[L_SELECT_ROW,]$numMuts <- SELECT_ROW$numMuts
        selectIsoformsLong[L_SELECT_ROW,]$pVal <- SELECT_ROW$pVal
        selectIsoformsLong[L_SELECT_ROW,]$isoName <- SELECT_ROW$isoName
        selectIsoformsLong[L_SELECT_ROW,]$isoLength<- SELECT_ROW$isoLength
        selectIsoformsLong[L_SELECT_ROW,]$pValDir <- SELECT_ROW$pValDir
        
        # Iterate row in question
        L_SELECT_ROW <- L_SELECT_ROW + 1
      } # End finding the top longFile and adding it p-value to the best vector
      
      # Find the top shortFile and add its p-value to the best vector
      if (length(shortFiles) > 0) {
        # Read files into memory
        HANDLE <- apply(shortDataFrame,
                        1,
                        function(row)
                          read.table(as.character(row["shortFiles"]),
                                     sep = ","))
        
        # Set number of mutations column via entire shortFiles vector
        shortDataFrame$numMuts <- lapply(HANDLE, function(row)
          row$V4)
        
        # Set the pValue of each entry in the shortFiles vector
        shortDataFrame$pVal <- lapply(HANDLE, function(row)
          row$V5)
        
        # Set the isoName by separating it from the full file path
        shortDataFrame$isoName <- lapply(shortFiles,
                                         function(full)
                                           basename(dirname(full)))
        
        # Set length of each isofom by running isoform_length()
        shortDataFrame$isoLength <- apply(shortDataFrame,
                                          1,
                                          function(row)
                                            isoform_length(row["isoName"]))
        
        # Set direction of p-value
        shortDataFrame$pValDir <-
          lapply(HANDLE, function(row)
            as.character(row$V6))
        
        # Gives '+' direction precedence in sorting if it is present
        if (preferredDirection %in% shortDataFrame$pValDir) {
          shortDataFrame$pValDir <- relevel(factor(shortDataFrame$pValDir,
                                                   levels = c("+", "-")),
                                            preferredDirection)
        }
        
        # Add top entry/row to selectIsoformsLong data.frame
        shortDataFrame <- as.data.frame(lapply(shortDataFrame, unlist), stringsAsFactors=FALSE)
        SELECT_ROW <- shortDataFrame[with(shortDataFrame,
                                          order(-numMuts, isoLength, isoName)), ][1, ]
        
        # Individually file values
        selectIsoformsShort[S_SELECT_ROW,]$longFiles <- SELECT_ROW$longFiles
        selectIsoformsShort[S_SELECT_ROW,]$numMuts <- SELECT_ROW$numMuts
        selectIsoformsShort[S_SELECT_ROW,]$pVal <- SELECT_ROW$pVal
        selectIsoformsShort[S_SELECT_ROW,]$isoName <- SELECT_ROW$isoName
        selectIsoformsShort[S_SELECT_ROW,]$isoLength<- SELECT_ROW$isoLength
        selectIsoformsShort[S_SELECT_ROW,]$pValDir <- SELECT_ROW$pValDir
          
        # Iterate row in question
        S_SELECT_ROW <- S_SELECT_ROW + 1
      } # End finding the top shortFile and adding it p-value to the best vector
    } # End directory for loop
    
    # Remove empty rows from select data.frames
    selectIsoformsLong <- na.omit(selectIsoformsLong)
    selectIsoformsShort <- na.omit(selectIsoformsShort)
  } # End else for nonexistent cache 
  
  # Cache the previous results
  write.csv(file = selectLongCache, selectIsoformsLong)
  write.csv(file = selectShortCache, selectIsoformsShort)
  
  # Output to the user
  print("Computing the adjusted long p-values")
  selectIsoformsLong$pValAdj <- p.adjust(apply(selectIsoformsLong,
                                               1,
                                               function(row)
                                                 as.numeric(as.character(row['pVal']))),
                                         method = correctionMethod)
  
  print("Computing the adjusted short p-values")
  selectIsoformsShort$pValAdj <- p.adjust(apply(selectIsoformsShort,
                                                1,
                                                function(row)
                                                  as.numeric(as.character(row['pVal']))),
                                          method = correctionMethod)
  
  
  return(list("long" = selectIsoformsLong,
              "short" = selectIsoformsShort))
}