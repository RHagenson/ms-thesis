#!/usr/bin/env Rscript
# This script takes a disorder mutation profile created by parsing TCGA data
# and generates a .pdf with the normal disorder plot marked by the observed disorder score

args <- commandArgs(trailingOnly = TRUE)
#print(args)

# The parameters that will eventual change
filename <- args[1] # "MUC16.001.long"  # The profile being processed
N = as.numeric(args[2]) # 100000  # The number of samples to take

profileDir = "~/Thesis/disorderCancer/data/profiles/"
figsDir = "figs/"  # The location of where figures should be written, do not ending forget '/'

# Read in the source file
path <- paste(profileDir, filename, sep="")

FILE <- read.delim(path, header=FALSE)

# Create a vector to hold the sum of random mutations
normalVector = vector(mode="double")

# Determine how many mutations are present
mutationNum <- sum(FILE$V4)
realLevel = as.numeric(sum(FILE$V3*FILE$V4))

for (i in 1:N) {
  normalVector <- append(normalVector, round(sum(sample(FILE$V3, replace=TRUE, size=mutationNum)), digits = 3))
}

# Generate a data frame with values and freq
frame <- as.data.frame(table(normalVector)/length(normalVector)*100)

# Add column names
colnames(frame) <- c("TotalDisOrderScore", "PercentFreq")

# Coerce factor from use of table() into numeric
frame$TotalDisOrderScore <- as.numeric(levels(frame$TotalDisOrderScore))[frame$TotalDisOrderScore]

# Open the output pdf for writing
pdf(paste(figsDir, filename, ".pdf", sep=""))

# Plot
plot(frame, xlab = "TotalDisorderScore", ylab = "Percent Frequency", type = "p", main=filename, pch=20)

# Add the real value to the plot
# points(x=c(realLevel), y=c(0.09), pch=25, col=20)
# Add a light-blue vertical line at real value
abline(v=realLevel, col=21)

# Close the graphic device to save to file
dev.off()

#
# Clean up memory
#
# Individual removes
rm(N)
rm(i)
rm(mutationNum)
rm(frame)

# Clear any remaining variables
rm(list = c(ls()))