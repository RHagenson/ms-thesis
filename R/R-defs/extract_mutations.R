# This function is intended to take the name of a cancer type, gene, and isoform number 
# returning a list of mutations that occur in that type in the form A251G
# The protein=T switch allows extractions of either the protein mutation or genetic mutation

extract_mutations <- function(cancer_type, gene, isoform_num, protein=T) {
  # The assumed directory tree is based off running script in ms-thesis/R/
  assumedDataDir = "../../disorderCancer/data"
  assumedAllMutsDir <- paste(assumedDataDir, "allMuts/", sep="/")
  allMutsNaming <- "_mut.txt"
  
  # Build location of file
  cancer_file <- paste(assumedAllMutsDir, paste0(cancer_type, allMutsNaming), sep="/")

  # read.table should ignore the commented header lines
  FILE.frame <- read.table(cancer_file, header = FALSE)
  
  # Map names to columns based on protein switch
  Gene_w_isoform <- 1
  Pos <- if(protein) 6 else 3
  Start_code <- if (protein) 7 else 4
  End_code <- if (protein) 8 else 5
  
  # Format isoform number
  if (length(isoform_num) == 1) {
    isoform_num = paste0("00", isoform_num)
  } else if (length(isoform_num) == 2) {
    isoform_num = paste0("0", isoform_num)
  }

  
  # Find indexes of matches
  matches <- which(FILE.frame[,Gene_w_isoform] == paste0(gene, ".", isoform_num))
  
  # Use index of matches to extract start, pos, and end values
  start_code_vector <- FILE.frame[matches, Start_code]
  position_vector <- FILE.frame[matches, Pos]
  end_code_vector <- FILE.frame[matches, End_code]
  
  output_list <- vector(mode = "character", 
                        length=length(start_code_vector))
  
  # Format output from parallel vectors
  for (index in 1:length(start_code_vector)) {
    output_list[index] = paste0(start_code_vector[index],
                                position_vector[index],
                                end_code_vector[index]
                                )
  }
  
  return(output_list)
}