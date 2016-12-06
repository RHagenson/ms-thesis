# Function to add text to the topright corner of a plot by default

# Function brought to you by: http://eranraviv.com/adding-text-to-r-plot/
corner_text <- function(text, location="topright"){
  legend(location, legend=text, bty ="n", pch=NA) 
}