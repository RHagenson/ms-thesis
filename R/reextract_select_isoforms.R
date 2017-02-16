# The code that was used to extract the previously selected isoforms

place_files <- "../R/processing/selectedIsoforms"

dir.create(place_files)

# Create a vector of the previously selected isoforms
files_list <- list.files("../R/outputs/25-10-16/p-adjusted/", 
                         pattern = "*.csv", 
                         recursive = T, 
                         full.names = T)

# Get the previously selected isoforms with all their excess
selectIsoforms <- lapply(files_list, 
                         read.csv, 
                         header=T, 
                         stringsAsFactors=F)

# Rename the entries to match that of the files_list
names(selectIsoforms) <- lapply(files_list, 
                                function(x) {paste(tail(strsplit(x, 
                                                                 split="/")[[1]],
                                                        n=2), 
                                                   collapse=".")})
# Remove my very, very extra-terrible naming scheme where
# I still have the file name attached AND the file name is bad
names(selectIsoforms) <- gsub("_log.csv", 
                              "", 
                              names(selectIsoforms))

# Where the magic happens
# Write to file by a mess of calls
for (item in seq_along(selectIsoforms)) {
filename <- paste0(place_files,"/",
                   attributes(selectIsoforms)$name[item],
                   ".select")
data <- selectIsoforms[[item]]$isoName
data <- gsub(".short", "", data)
data <- gsub(".long", "", data)
write(data, file=filename)
}
