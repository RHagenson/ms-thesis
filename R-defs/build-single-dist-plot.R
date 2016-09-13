#!/usr/bin/env Rscript
# This script takes a disorder mutation profile created by parsing TCGA data
# and generates a .pdf with the normal disorder plot marked by the observed disorder score
# only plots with significant P values are output.

build.plot <- function(filename, number, profileDir, figsDir, pValCut=0.05, test="two") {
  # p-value test options "less", "greater", and "two"
  
  CSV <- "outputs/isoformPvaluesWithTests.csv"
  
  # The parameters that will eventual change
  filename <-
    filename # "MUC16.001.long"  # The profile being processed
  N = as.numeric(number) # 100000  # The number of samples to take
  
  profileDir = profileDir  # The location of where the filename is found
  figsDir = figsDir  # The location of where figures should be written, do not ending forget '/'
  
  # Read in the source file
  path <- paste(profileDir, filename, sep = "")
  
  FILE <- read.delim(path, header = FALSE)
  
  # Create a vector to hold the sum of random mutations
  normalVector = vector(mode = "double")
  
  # Determine how many mutations are present
  mutationNum <- sum(FILE$V4)
  realLevel = as.numeric(sum(FILE$V3 * FILE$V4))
  
  for (i in 1:N) {
    # Print a helpful message to the user for what is being done
    if ((i %% 1000) == 0) {
      print(paste("On sample number", as.character(i), "for", filename))
    }
    
    normalVector <-
      append(normalVector, round(sum(
        sample(FILE$V3, replace = TRUE, size = mutationNum)
      ), digits = 3))
  }
  
  # Generate a data frame with values and freq as percent
  frame <-
    as.data.frame(table(normalVector) / length(normalVector) * 100)

  # Calculate the p-value based on given test option  
  switch(test,
         less={pValue <- sum(realLevel < normalVector) / length(normalVector)},
         greater={pValue <- sum(realLevel > normalVector) / length(normalVector)},
         two={t <- (mean(normalVector)-realLevel)/(sd(normalVector)/sqrt(length(normalVector)))
              pValue <- 2*pt(-abs(t), df=length(normalVector)-1)},
         stop("Not a valid test option, use 'less', 'greater', or 'two"))

    # Two tail test, Method 2
    #     a <- realLevel
    #     s <- sd(normalVector)
    #     n <- length(normalVector)
    #     xbar <- mean(normalVector)
    #     z <- (xbar-a)/(s/sqrt(n))
    #     pValue = 2*pnorm(-abs(z))
  
  
  # Free resources now that the vector has served its purpose
  # rm(normalVector)
  
  # Add column names
  colnames(frame) <- c("TotalDisOrderScore", "PercentFreq")
  
  # Coerce factor from use of table() into numeric
  frame$TotalDisOrderScore <-
    as.numeric(levels(frame$TotalDisOrderScore))[frame$TotalDisOrderScore]
  
  # No matter the pValue add the isoform to the results CSV
  write.table(x=data.frame(filename,pValue,test), file=CSV, 
                         append = TRUE, row.names = FALSE,
                         quote = FALSE, sep=",", col.names = FALSE)
  
  if (pValue <= pValCut) {
    # Open the output pdf for writing, with naming based on test type
    pdf(paste(figsDir, filename, ".", test, ".pdf", sep = ""))
    
    # Plot
    plot(density(normalVector), xlab = "TotalDisorderScore", ylab = "Percent Frequency", type = "p", main =
        filename, pch = 20
    )
    
    # Function brought to you by: http://eranraviv.com/adding-text-to-r-plot/
    Corner_text <- function(text, location="topright"){
      legend(location,legend=text, bty ="n", pch=NA) 
    }
    Corner_text(text = paste("p-value:", pValue))
    
    # Add the real value to the plot
    # points(x=c(realLevel), y=c(0.09), pch=25, col=20)
    # Add a light-blue vertical line at real value
    abline(v = realLevel, col = 21)
    
    # Close the graphic device to save to file
    dev.off()
  }
}