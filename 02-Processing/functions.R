
# Title:

# Purpose:

# Project name: STox-Pmicrops-2025

# Author: A. Duque
# Contact details: duque.andrea2000@gmail.com

# Date script created: Thu Apr 24 11:22:37 2025 --------------------------------
# Date last modified: Thu Apr 24 11:22:37 2025 ---------------------------------






# package dependencies

## ???
library(openxlsx)

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





#' Select the best hits in Blast
#' 
#' Explicar o que é que acontece aqui!!!
#' 
#' @param blast_data data frame with Blastp results (must include "GeneID" column)
#' @return 
#'  
select_best_hits <- function(blast_data){
  
  # Add GeneID column 
  blast_data <- Add_GeneID_column(blast_data)
  
  # Order: lowest evalue, highest pident
  blast_data <- blast_data[order(blast_data$evalue, -blast_data$pident, decreasing = TRUE),]
  
  # Best hits (first row)
  best_hits <- blast_data[match(unique(blast_data$GeneID), blast_data$GeneID), ]
  
  return(best_hits)
  
}





#' Save a dataframe or a list of dataframes to an Excel file
#' 
#' @param data A dataframe or a named list of dataframes
#' @param output_file The filename for the Excel file
#' @return Saves the data to an Excel file
#' 
#' @example 
#' save_to_excel(data, "/03-Output/01-DEG-Analysis/fileName.xlsx")
#' save_to_excel(data, "fileName.xlsx")
#' 
save_to_excel <- function(data, output_file) {
  
  work_book <- createWorkbook()
  
  if (is.data.frame(data)) {
    
    sheet_name <- deparse(substitute(data))  # Get the variable name
    addWorksheet(work_book, sheet_name)
    writeData(work_book, sheet_name, data)
    
  } else if (is.list(data) && all(sapply(data, is.data.frame))) {
    
    for (name in names(data)) {
      addWorksheet(work_book, name)
      writeData(work_book, name, data[[name]])
    }
  }
  
  saveWorkbook(work_book, file = output_file, overwrite = TRUE)
  cat("Excel created successfully:", output_file, "\n")
}





