# This script takes a disorder mutation profile created by parsing TCGA data
# and generates a .pdf with the normal disorder plot marked by the observed disorder score
# only plots with significant P values are output.

source("R-defs/corner_text.R")  # Adds text in corner of plot

build_plot <- function(pdfPath, normalVector, filename, pValue, 
                       pValDirection, avgDisorder, realLevel, numMutations) {
  if(! dir.exists(pdfPath)){
    dir.create(pdfPath, recursive = TRUE)  # Create path if it does not exist
  }
  
  # Open the output pdf for writing, with naming based on where the profile is from
  pdfName = sub(".prof", ".pdf", sub(".short", "", sub(".long", "", filename)))
  pdf(paste(pdfPath, pdfName, sep = "/"))
  
  # Plot
  plot(density(normalVector), 
       xlab = "Total Disorder Score", 
       ylab = "Probability Density", 
       type = "p", 
       main = sub(".prof", "", filename), 
       pch = 20
  )
  
  corner_text(text = paste("p-value: ", paste(pValDirection, pValue, sep=""), "\n",
                           "Average score: ", avgDisorder, "\n",
                           "Observed score: ", realLevel, "\n",
                           "Number of mutations: ", numMutations, sep=""))
  
  # Add a light-blue vertical line at real value
  abline(v = realLevel, col = 21)
  
  # Close the graphic device to save to file
  dev.off()
  
}