# Build a lollipop plot of the observed mutations using ggplot2

library(ggplot2)
source("R-defs/isoform_length.R")

lollipop_plot <- function(isoform_name, cancer_file, legend=F) {
  len <- isoform_length(isoform_name)
  # Define the data.frame used for plotting
  df_c <- data.frame(Pos = 1:len)
  
  cancer <- read.delim(cancer_file, header=F)
  
  df_c$Muts <- sapply(df_c$Pos, function(pos) {
    sub_cancer <- cancer[which(cancer$V1 == isoform_name),]
    sum(sub_cancer$V6 == pos)
  })
  
  plot <- ggplot(df_c, aes(Pos, Muts)) +
    geom_point() +
    geom_segment(aes(x = Pos, y = 0, xend = Pos, yend = Muts)) +
    geom_hline(yintercept = 0, lty = 2) +
    xlab("Residue Position") + 
    ylab("Number of Mutations") + 
    ggtitle(as.character(isoform_name)) +
    xlim(0, len)
  
  if(legend) {
    plot <- plot + aes(colour = df_c$Muts) + theme(legend.title = element_blank())
  }
  
  return(plot)
}