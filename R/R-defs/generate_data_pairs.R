# This function takes the path to where cancer profiles tree is and returns values for lapply
# into build.plot <- function(filename, number, profileDir, figsDir, pValCut=0.05) {...}
# Where filename = ABCA12.001.long.prof, number = 100000 or more, 
# profileDir = ~/Thesis/data/profiles/LUSC/ABCA12/, figsDir = "figs/"

generate_data_pairs <- function(profilesDir, number=1000000, figsDir="figs/", pValCut=0.05) {
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
  pValCutVector <- rep.int(pValCut, length(filenameVector))
  
  # Create data.frame from vectors
#   return.matrix <- matrix(rbind(filenameVector, 
#                                  numberVector, 
#                                  profileDirVector, 
#                                  figsDirVector, 
#                                  pValCutVector), ncol = 5, byrow = TRUE)
  return.frame <- cbind.data.frame(filenameVector, 
                                   numberVector, 
                                   profileDirVector, 
                                   figsDirVector, 
                                   pValCutVector)
  
  return(return.frame)
  # print(return.matrix)
}
