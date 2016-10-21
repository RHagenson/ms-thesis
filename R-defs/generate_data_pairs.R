# This function takes the path to where cancer profiles tree is and returns values for lapply
# into build_plot <- function(...) {...} found in R-defs/

generate_data_pairs <- function(profilesDir, number=1000000, figsDir="figs/", 
                                outputDir="outputs/", pValCut=0.05, plot = FALSE) {
  # Inform user what is being done
  print(paste("Generating data pairs within:", as.character(profilesDir)))
  
  # Remove trailing '/' because dir() below adds it
  profilesDir <- sub(pattern = "/$", replacement = "", profilesDir)
  
  # Gather all the .prof files in the tree with absolute path names
  profFilesVector <- dir(path = profilesDir, 
                     pattern = "*.prof", 
                     full.names = TRUE, 
                     recursive = TRUE)
  
  # Define each of the lappy options
  filenameVector <- basename(profFilesVector)
  numberVector <- rep.int(number, length(filenameVector))
  profileDirVector <- dirname(profFilesVector)
  figsDirVector <- rep(figsDir, length(filenameVector))
  outputDirVector <- rep(outputDir, length(filenameVector))
  pValCutVector <- rep.int(pValCut, length(filenameVector))
  plotVector <- rep.int(plot, length(filenameVector))
  
  # Create data.frame from vectors
  return.frame <- cbind.data.frame(filenameVector, 
                                   numberVector, 
                                   profileDirVector, 
                                   figsDirVector, 
                                   outputDirVector,
                                   pValCutVector,
                                   plotVector)
  
  return(return.frame)
}
