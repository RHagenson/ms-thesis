# Build a lollipop plot of the observed mutations using ggplot2

library(ggplot2)
source("R-defs/isoform_length.R")

lollipop_plot <- function(isoform_name, cancer_file) {
  # Define the data.frame used for plotting
  df_c <- data.frame(Pos = 1:isoform_length(isoform_name))
  
  cancer <- read.delim(cancer_file, header=F)
  
  df_c$Muts <- sapply(df_c$Pos, function(pos) {
    sub_cancer <- cancer[which(cancer$V1 == isoform_name),]
    sum(sub_cancer$V6 == pos)
  })
  
  ggplot(df_c, aes(Pos, Muts)) +
    geom_point() +
    geom_segment(aes(x = Pos, y = 0, xend = Pos, yend = Muts)) +
    geom_hline(yintercept = 0, lty = 2)
}