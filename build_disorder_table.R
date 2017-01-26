# A wrapper of the extract_disorder function and the assumed directory structure to generate a table of:
# cancer type, isoform(.long/.short), disorder

source("R-defs/extract_disorder.R")
library("optparse")
library("parallel")

# Define the CLI arguments
option_list = list(
  make_option(c("-d", "--date"), 
              type="character", 
              default=format(Sys.Date(), format="%d-%m-%y"), 
              help="Date in profiles to process [form = DD-MM-YY]", 
              metavar="character")
);

# Parse the arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Make date required CLI arguments
if (is.null(opt$date)){
  print_help(opt_parser)
  stop("-d/--date is a required argument.n", call.=FALSE)
}

### Variables
# now_date <- opt$date
now_date <- "27-10-16"

# The assumed directory tree is based off running script in ms-thesis/R/
assumedDataDir = "../../disorderCancer/data"
assumedRefSeqDir <- paste(assumedDataDir, "refSeq/", sep="/")
assumedIupredLong <- paste(assumedRefSeqDir, "iupredLong", sep="/")
assumedIupredShort <- paste(assumedRefSeqDir, "iupredShort", sep="/")
assumedProfiles <- paste(assumedDataDir, "profiles/", sep="/")

profilesDirectory <- paste(assumedProfiles, now_date, sep="/")

# Determine how large the vectors need to be
FILES <- list.files(profilesDirectory, 
                    recursive = TRUE,
                    full.names = TRUE,
                    pattern = "*.prof")
vector_length <- length(FILES)

# Initialize two parallel vectors
cancer_types = vector(mode="character", length = vector_length)
isoforms = vector(mode="character", length = vector_length)

INDEX <- 1

pattern <- paste("^.+", now_date, "(.+)", ".+$", sep="/")

# Fill the parallel cancer_types and isoforms vectors
for(file in FILES) {
  
  split <- unlist(strsplit(file, split = "/"))
  cancer_types[INDEX] <- split[length(split)-2]
  
  isoforms[INDEX] <- sub(pattern=".prof", replacement = "", basename(file))
  
  INDEX <- INDEX + 1
}

data <- mcmapply(
  extract_disorder,
  cancer_type = cancer_types, 
  isoform_name = isoforms,
  mc.cores = detectCores() - 1
)

