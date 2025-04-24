
# Title:

# Purpose:

# Project name: STox-Pmicrops-2025

# Author: A. Duque
# Contact details: duque.andrea2000@gmail.com

# Date script created: Thu Apr 24 11:22:37 2025 --------------------------------
# Date last modified: Thu Apr 24 11:22:37 2025 ---------------------------------






# package dependencies

## ???


## BiocManager
library(tximport) 

## tidyverse





# ==============================================================================
# TASK 1 - Blablabla blablabla blabla
# ==============================================================================


#' Load BLASTp output data
#' 
#'  @param filename Path to the BLASTp CSV file.
#'  @return A data.frame with standardized column names. 
#' 
#'  @example
#'  data <- load_blastp_data("./Blastp/ExtremeOceans_Blastp_SwissProt.csv")


load_blastp_data <- function(filename){
  
  # Read data (table format)
  data <- read.delim(filename, sep = ";", header = FALSE)
  
  # Assign column names | Note.: Could change according to the input
  colnames(data) <- c(
    "qseqid", "sseqid","pident","length","mismatch",
    "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore","sacc","qcovs","stitle"
    )
  
  return(data)
}





















