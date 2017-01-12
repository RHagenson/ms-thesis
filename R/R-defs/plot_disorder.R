# This function creates and entropy plot of a given isoform (with .long/.short file extension) within a given cancer type

# Currently broken, repeats legend until off the page

source("R-defs/extract_disorder.R")
source("R-defs/corner_text.R")
library(entropy)
library(ggplot2)

plot_disorder <- function(cancer_type, isoform_name) {
  # The assumed directory tree is based off running script in ms-thesis/R/
  assumedDataDir = "../../disorderCancer/data"
  assumedRefSeqDir <- paste(assumedDataDir, "refSeq/", sep="/")
  assumedIupredLong <- paste(assumedRefSeqDir, "iupredLong", sep="/")
  assumedIupredShort <- paste(assumedRefSeqDir, "iupredShort", sep="/")
  
  # The root of where disorder PDF plots should reside, within figs directory by cancer type
  pdfRoot <- paste("./figs/disorder")
  pdfPath <- paste(pdfRoot, cancer_type, sep="/")
  
  # Build location of file
  if (grepl("long", isoform_name)) {
    isoformFile <- paste(assumedIupredLong, isoform_name, sep="/")
  } else {
    isoformFile <- paste(assumedIupredShort, isoform_name, sep="/")
  }
  
  # Extract the data.frame of Position, AA, Disorder
  iso_disorder <- extract_disorder(cancer_type = cancer_type, isoform_name = isoform_name)
  
  # Create the path for PDF if it does not exist
  if(! dir.exists(pdfPath)){
    dir.create(pdfPath, recursive = TRUE)
  }
  
  # Open the output pdf for writing, with naming based on where the profile is from
  pdfName = paste(isoform_name, ".pdf", sep="")
  #pdf(paste(pdfPath, pdfName, sep = "/"))
  
  # plot(x = iso_disorder$Position,
  #      y = iso_disorder$Disorder,
  #      xlab = "Position Number", 
  #      ylab = "Disorder Score", 
  #      type = "p", 
  #      main = as.character(isoform_name), 
  #      pch = 20
  # )
  
  # Calculate plot descriptors
  entropy <- round(entropy(iso_disorder$Disorder, frequency(iso_disorder$Disorder)), digits = 4)
  avg_disorder <- round(ave(iso_disorder$Disorder), digits = 4)
  
  # Add description to plot
  #corner_text(c(paste("Entropy: ", entropy, "\n",
  #                  "Average disorder: ", avg_disorder, "\n", sep = "")))
  
  # Generate plot
  ggplot(iso_disorder, aes(x = Position,
                           y = Disorder)) +
    geom_point() +
    ggtitle(isoform_name) +
    xlab("Position Number") +
    ylab("Disorder Score")
    # geom_label(vjust="", 
    #            label = paste("Entropy: ", entropy, "\n",
    #                          "Average disorder: ", avg_disorder, 
    #                          sep = ""))
    # 
  # Save file
  ggsave(paste(pdfPath, pdfName, sep = "/"), device = "pdf")
  
  # Close the graphic device to save to file
  #dev.off()
}