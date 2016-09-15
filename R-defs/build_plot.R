# This script takes a disorder mutation profile created by parsing TCGA data
# and generates a .pdf with the normal disorder plot marked by the observed disorder score
# only plots with significant P values are output.

build_plot <- function(filename, number, profileDir, figsDir, pValCut=0.05) {
  # Remove trailing '/' because paste() below adds it when creating path
  # profilesDir <- sub(pattern = "/$", replacement = "", profilesDir)
  
  CSV <- "outputs/isoformPvaluesWithTests.csv"
  
  # The parameters that will eventual change
  filename <-
    filename # "MUC16.001.long.prof"  # The profile being processed
  
  N = as.numeric(number) # 100000  # The number of samples to take
  
  profileDir = profileDir  # The location of where the filename is found
  figsDir = figsDir  # The location of where figures should be written, do not ending forget '/'
  
  # Read in the source file
  path <- paste(profileDir, filename, sep = "/")
  
  FILE <- read.delim(path, header = FALSE)
  
  # Create a vector to hold the sum of random mutations
  normalVector = vector(mode = "double")
  
  # Determine how many mutations are present
  mutationNum <- sum(FILE$V4)
  realLevel = as.numeric(sum(FILE$V3 * FILE$V4))
  numMutations = as.numeric(sum(FILE$V4))
  
  for (i in 1:N) {
    # Print a helpful message to the user for what is being done
    if ((i %% 1000) == 0) {
      print(paste("On sample number", as.character(i), "for", filename))
    }
    
    normalVector <-
      append(normalVector, 
             round(sum(sample(FILE$V3, replace = TRUE, size = mutationNum)), 
                   digits = 3))
  }
  
  # Determine average disorder score from normal curve
  avgDisorder = round(sum(normalVector) / length(normalVector), digits=5)
  
  # Generate a data frame with values and freq as percent
  frame <-
    as.data.frame(table(normalVector) / length(normalVector) * 100)

  # Calculate the empirical p-value, using min to find whether it deviates high or lower than average 
  pValue <- min((sum(realLevel < normalVector) / length(normalVector)), 
                (sum(realLevel > normalVector) / length(normalVector)))
  
  # Free resources now that the vector has served its purpose
  # rm(normalVector)
  
  # Add column names
  colnames(frame) <- c("TotalDisOrderScore", "PercentFreq")
  
  # Coerce factor from use of table() into numeric
  frame$TotalDisOrderScore <-
    as.numeric(levels(frame$TotalDisOrderScore))[frame$TotalDisOrderScore]
  
  # No matter the pValue add the isoform to the results CSV
  # Columns should be in the following order
  #   1. isoform name
  #   3. observed disorder score
  #   4. average random disorder score
  #   5. total number of mutations
  #   6. empirical p-value
  write.table(x=data.frame(filename, realLevel, avgDisorder, numMutations, pValue), 
              file=CSV, append = TRUE, row.names = FALSE,
              quote = FALSE, sep=",", col.names = FALSE)
  
  # if (pValue <= pValCut) {
  if (TRUE) {
    # Open the output pdf for writing, with naming based on test type
    pdf(paste(figsDir, sub(".prof", "", filename), ".pdf", sep = ""))
    
    # Plot
    plot(density(normalVector), xlab = "TotalDisorderScore", ylab = "Percent Frequency", type = "p", main =
           sub(".prof", "", filename), pch = 20
    )
    
    # Function brought to you by: http://eranraviv.com/adding-text-to-r-plot/
    Corner_text <- function(text, location="topright"){
      legend(location,legend=text, bty ="n", pch=NA) 
    }
    Corner_text(text = paste("p-value: ", pValue, "\n",
                             "Average score: ", avgDisorder, "\n",
                             "Observed score: ", realLevel, "\n",
                             "Number of mutations: ", numMutations, sep=""))
    
    # Add the real value to the plot
    # points(x=c(realLevel), y=c(0.09), pch=25, col=20)
    # Add a light-blue vertical line at real value
    abline(v = realLevel, col = 21)
    
    # Close the graphic device to save to file
    dev.off()
  }
}