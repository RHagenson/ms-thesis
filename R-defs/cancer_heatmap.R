# A functional wrapper to create the heatmap of all significant isoforms following
# p.adjust()

library(stats)

cancer_heatmap <- function(something) {
  # Prepare a data.frame() to hold the isoform adjusted p-value for each cancer type
  # In the frame: cancerType, isoName1, isoName2, isoName3, etc
  # The cancerType should then end up on the y-axis, and the different isoNames
  # should then end up on the x-axis.
  #
  # There is also the possibility of creating a heatmap of the obs/avg disorder ratio
  # only choosing to look at those that remained significant in at least one cancer
  # Thus, cancerType, obs/avg of cancerType with isoform1, 
  # obs/avg of cancerType with isoform2...and so on.
  
  frame <- data.frame()
  
  # Heatmap requires a data.matrix not a data.frame
  matrix <- data.matrix(frame)
  
  # Better heatmap col option: brewer.pal(11, "RdYlBu")
  # Must then add library(RColorBrewer) to header
  
  nba_heatmap <- heatmap(matrix, 
                         Rowv=NA, 
                         Colv=NA, 
                         col = cm.colors(256), 
                         scale="column", 
                         margins=c(5,10))
  
}