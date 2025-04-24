
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







# Constants and Global settings

sample_name <- c( "F14.Sp_TR","F15.Su_TR","F16.Su_RA",
                  "F17.Sp_RF","F18.Sp_TR","F19.Sp_RA",
                  "F2.Su_RF" ,"F2.Su_TR" ,"F21.Sp_RA",
                  "F21.Sp_TR","F22.Sp_RF","F27.Sp_RA",
                  "F28.Su_RF","F30.Su_TR","F35.Sp_RF",
                  "F35.Su_RF","F38.Su_RA","F4.Su_RA" ) 





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
#'  
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




#' Add "GeneID" column to Blastp data
#' 
#' Add's "GeneID" column to Blastp data by splitting "qseqid" at "i"
#' @param blast_data Data frame with "qseqid" column
#' @return Data frame witha new "GeneID" column
#' 
Add_GeneID_column <- function(blast_data){
  
  blast_data$GeneID <- sapply(blast_data$qseqid,
                              function(x) strsplit(x, "i", fixed = TRUE)[[1]][1])
  
  return(blast_data)
}


















