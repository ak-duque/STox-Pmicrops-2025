
# Title:

# Purpose:

# Project name: STox-Pmicrops-2025

# Author: A. Duque
# Contact details: duque.andrea2000@gmail.com

# Date script created: Thu Apr 24 11:22:37 2025 --------------------------------
# Date last modified: Thu Apr 24 11:22:37 2025 ---------------------------------






# package dependencies

## CRAN
library(showtext)
library(openxlsx)
library(stringr)

## BiocManager
library(tximport)
library(ggrepel)
library(PCAtools) # require package: ggrepel
library(limma)
library(edgeR) # require package: limma
library(UniprotR)
library(multiGSEA)
library(org.Dr.eg.db)


## tidyverse
library(tidyr)
library(dplyr)

## Plot
library(ggplot2)
library(ggVennDiagram)


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
Add_GeneID_column <- function(blast_data) {
  
  # Debugging: 
  # Example of correct qseqid: TRINITY_DN0_c10_g1_i1
  not_trinity <- blast_data$qseqid[!grepl("TRINITY_D", blast_data$qseqid)]
  if (length(not_trinity) > 0) {
    #cat("qseqid sem 'TRINITY_D':\n")
    #print(not_trinity)
  } else {
    #cat("Todos os qseqid têm 'TRINITY_D'.\n")
  }
  
  
  blast_data <- blast_data[grepl("TRINITY_D", blast_data$qseqid), ] 
 
  blast_data$GeneID <- sapply(
    blast_data$qseqid,
    function(x) strsplit(x, "i", fixed = TRUE)[[1]][1]
  )
  
  return(blast_data)
}








#' Select the best hits in Blast
#' 
#' Diferentes trinityID's podem dar o mesmo GeneID. Portanto vamos ter várias
#' entradas com o mesmo GeneID. 
#' 
#' @param blast_data data frame with Blastp results (must include "GeneID" column)
#' @return 
#'  
select_best_hits <- function(blast_data){
  # version 1
  # Add GeneID column 
  blast_data <- Add_GeneID_column(blast_data) # Vai haver GeneID repetidos!
  
  # Debugging step
  # Example: TRINITY_DN0_c10_g1_i1
  #example <- blast_data %>% 
  #  filter(GeneID == "TRINITY_DN100212_c0_g1_") %>%
  #  select(c("qseqid","GeneID","evalue","pident"))
  #View(example)
  
  # Order: lowest evalue and highest pident
  blast_data <- blast_data[order(blast_data$evalue, -blast_data$pident),]
  
  # Best hits (first row)
  best_hits <- blast_data[match(unique(blast_data$GeneID), blast_data$GeneID), ]
  
  # Debugging step
  # Example: TRINITY_DN0_c10_g1_i1
  #example_best_hit <- best_hits %>% 
  #  filter(GeneID == "TRINITY_DN100212_c0_g1_") %>%
  #  select(c("qseqid","GeneID","evalue","pident"))
  #View(example_best_hit)
  
  

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
      df <- data[[name]]
      if (nrow(df) > 0) {
        addWorksheet(work_book, name)
        writeData(work_book, name, df)
      }
    }
  }
  
  saveWorkbook(work_book, file = output_file, overwrite = TRUE)
  cat("Excel created successfully:", output_file, "\n")
}





#' Read excel file with multiple sheets
#' 
#' @param path the directory of the file
#' @return list of dataframes(sheets)
#' @references https://www.geeksforgeeks.org/how-to-read-a-xlsx-file-with-multiple-sheets-in-r/
#' 
read_excel <- function(path){
  
  # getting data from sheets
  sheets <- openxlsx::getSheetNames(path)
  dataframe <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=path)
  
  # assigning names to dataframe
  names(dataframe) <- sheets
  
  return(dataframe)
}












#' Process Blastp and Kallisto(counts) merged result
#' 
#' @param BpK_data
#' @return 
process_BpK_data <- function(BpK_data){
  
  
  # Test
  #BpK_data <- BpK_swissprot
  #dim(BpK_data)
  
  BpK_data <- BpK_data %>%
    mutate(GeneNameID = str_extract(sseqid, "(?<=\\|)[^|]+(?=_)"),
           # sp|P79762|ZP3_CHICK -> ZP3 
           Accession = str_extract(sacc, "(?<=\\|)[^\\|]+(?=\\|)"),
           # sp|P79762|ZP3_CHICK -> P79762 
           OS = str_extract(stitle, "(?<=OS=)[^ ]+ [^ ]+")
           # sp|Q90508|VIT1_FUNHE Vitellogenin-1 OS=Fundulus heteroclitus OX=8078 GN=vtg1 PE=1 SV=2 -> Fundulus heteroclitus
           )
  
  
  
  # Why we need to order again? 
  # When we add "GeneNameID" we are going to found duplicated values
  #table(BpK_data$GeneNameID)
  
  
  # Debugging step
  # Example: ABCA1 (SwissProt)
  #example <- BpK_data %>% 
  #  filter(GeneNameID== "ABCA1") %>%
  #  select(c("GeneNameID","evalue","pident"))
  #View(example)
  
  
  
  
  # Order: lowest evalue and highest pident
  BpK_data <- BpK_data[order(BpK_data$evalue, -BpK_data$pident), ]
  
  # Select the best hit (1st row)
  BpK_data <- BpK_data[match(unique(BpK_data$GeneNameID), BpK_data$GeneNameID), ]
  #dim(BpK_data)
  
  
  
  
  # Let's confirm we select the right row
  # Debugging step
  # Example: ABCA1  (SwissProt)
  #example1 <- BpK_data %>% 
  #  filter(GeneNameID== "ABCA1") %>%
  #  select(c("GeneNameID","evalue","pident"))
  #View(example1)
  
  
  
  
  # We also need to confirm if the string extraction went well 
  #sum(is.na(BpK_data)) # 3 NA           
  #colSums(is.na(BpK_data)) # Column "OS"
  
  # Further investigate why it went wrong, to know if we can fix it 
  #View(BpK_data[is.na(BpK_data$OS), ])
  
  # the extraction fail because stitle was incomplete (e.g., sp|Q3SYR3|ABEC2_BOVIN Probable C- )
  # how to fix? (solution specific to my data results/problems) 
  
  # BpK_swissprot -> missing: HUMAN (Homo sapiens) | BOVIN (Bos taurus)  | MOUSE (Mus musculus)
  # BpK_zebrafish -> missing: DANRE (Danio rerio) 
  # BpK_plat -> all ok!
  
  
  if (anyNA(BpK_data)) {
 
    BpK_data <- BpK_data %>%
      # ifelse(test, yes, no)

      # test
      mutate(OS = ifelse(is.na(OS),
                         
                         # yes -> Se OS for NA - executa case_when()
                         case_when(
                           grepl("HUMAN", sseqid) ~ "Homo sapiens",
                           grepl("BOVIN", sseqid) ~ "Bos taurus",
                           grepl("MOUSE", sseqid) ~ "Mus musculus",
                           grepl("DANRE", sseqid) ~ "Danio rerio"),
                         
                         # no -> Senão mantém o valor original de OS
                         OS))
  }
  
  
  #sum(is.na(BpK_data))
  #head(BpK_data$OS, 20)
  # confirm if it is well done:
  example3 <- BpK_data %>% 
    filter(Accession %in% c("Q99J72","P14060", "Q3SYR3")) %>%
    select(OS, stitle)
  #View(example3)
  

  return(BpK_data) 
}







#' Title
#' 
#' @param value
#' @param metric
#' @return
#' 
#' @example  
evaluate_sequence_quality <- function(value, metric) {
  
  # By annotation
  if (metric == "pident") {
    # pident -> percentage of identical positions
    case_when(
      value < 30 ~ "Bad",
      value >= 30 & value <= 70 ~ "Good",
      value > 70 ~ "Excellent"
    )
  } else if (metric == "evalue") {
    # evalue -> expect value
    # The lower the E value, the more significant the score and the alignment.
    case_when(
      value > 1e-10 ~ "Bad",
      value <= 1e-70 ~ "Excellent",
      TRUE ~ "Good"
    )
    
  # By expression  
  } else if (metric == "logFC") {
    case_when(
      # Genes -> overexpressed (sobrexpresso) ou underexpressed (subexpresso)
      # Proteinas -> upregulated (super-regulado) ou downregulated (sub-regulado)
      # so (although not entirely correct):
      # under -> Down
      # over -> Up
      
      value < -1 ~ "Underexpressed",
      value > 1 ~ "Overexpressed",
      TRUE ~ "No Significant"
    )
  } else if (metric == "FDR") {
    case_when(
      value < 0.01 ~ "Excellent",
      value <= 0.05 ~ "Good",
      TRUE ~ "Bad"
    )
  }
}




#' Title
#' 
#' @param reference_data
#' @return 
#' 
format_table <- function(table) {
  paste(names(table), table, sep=": ", collapse=" | ")
}





#' Title
#' 
#' @param reference_data
#' @param plat_data
#' @return 
#' 
#' @example 
#' 
identify_putative_plat_genes <- function(reference_data, plat_data){
  
  
  # Test
  #reference_data <- Data_SP
  #plat_data <- Data_P
  
  
  matched_GeneNameID <- merge(reference_data, plat_data, by = "GeneNameID",
                              suffixes = c(".Ref", ".Plat")) 
                              # R -> reference
                              # P -> plat
  
  #View(matched_GeneNameID)
  
  
  
  if(nrow(matched_GeneNameID)>0){

    matched_GeneNameID_processed <- matched_GeneNameID %>%
      mutate(
        pident_category.Ref = evaluate_sequence_quality(pident.Ref, "pident"),
        evalue_category.Ref = evaluate_sequence_quality(evalue.Ref, "evalue"),
        pident_category.Plat = evaluate_sequence_quality(pident.Plat, "pident"),
        evalue_category.Plat = evaluate_sequence_quality(evalue.Plat, "evalue")
      ) %>%
      select(c(GeneNameID,
               pident.Ref,
               pident_category.Ref,
               evalue.Ref,
               evalue_category.Ref,
               pident.Plat,
               pident_category.Plat,
               evalue.Plat,
               evalue_category.Plat, 
               OS.Ref,OS.Plat,
               GeneID.Ref, GeneID.Plat))
   
    
    #View(matched_GeneNameID_processed)
     

    
    putative_plat_data <- matched_GeneNameID_processed %>%
      # Old strategy (in case professor want's the previous one)
      #filter(
      #  evalue.Plat < evalue.Ref, #só com e-value -> 36
      #  OS.Ref == OS.Plat #evalue e mesmo OS -> 32
      #)
      # New strategy (Sat May  3 18:51:55 2025)
      filter(
        pident_category.Ref == "Excellent",
        evalue_category.Ref == "Excellent",
        pident_category.Plat == "Excellent",
        evalue_category.Plat == "Excellent",
        
        OS.Ref == OS.Plat
      )
    
    
    #View(putative_plat_data)
    
    # To confirm if the filtering was well done  
    #View(putative_plat_data[, c("evalue.Plat", "evalue.Ref","OS.Ref","OS.Plat")]) # -> 32 GeneNameID
    #putative_plat_data%>% count(OS.Ref)
    
    
    # --------------------------------------------------------------------------
    # Print summary ------------------------------------------------------------
    summary <- paste(
      "--- SUMMARY ---",
      paste("Total GeneNameID in reference data:", nrow(reference_data)),
      paste("Total GeneNameID in plat data:", nrow(plat_data)),
      paste("GeneNameID in common:", nrow(matched_GeneNameID)),
      "",
      "--- Categories ---",
      paste("pident (Ref):", format_table(table(matched_GeneNameID_processed$pident_category.Ref))),
      paste("evalue (Ref):", format_table(table(matched_GeneNameID_processed$evalue_category.Ref))),
      paste("pident (Plat):", format_table(table(matched_GeneNameID_processed$pident_category.Plat))),
      paste("evalue (Plat):", format_table(table(matched_GeneNameID_processed$evalue_category.Plat))),
      "",
      "--- OS in common (examples) ---",
      paste("Common GeneNameID with same OS:", 
            sum(matched_GeneNameID$OS.Ref == matched_GeneNameID$OS.Plat)),
      paste("Most common organism in OS.Ref:", 
            names(which.max(table(matched_GeneNameID$OS.Ref)))),
      paste("Most common organism in OS.Plat:", 
            names(which.max(table(matched_GeneNameID$OS.Plat)))),
      "",
      "--- Putative parasite genes ---",
      "",
      "criteria:",
      "e-value == Excellent (evalue <= 1e-70)",
      "pident == Excellent (pident >= 70)",
      "Equal organism",
      "",
      paste("We identify a total of", nrow(putative_plat_data), "putative parasite genes"),
      paste(putative_plat_data$GeneNameID, collapse = ", "),
      "",
      paste("pident (Ref): Min:", min(putative_plat_data$pident.Ref), " | Max:", max(putative_plat_data$pident.Ref)),
      paste("evalue (Ref): Min:", min(putative_plat_data$evalue.Ref), " | Max:", max(putative_plat_data$evalue.Ref)),
      paste("pident (Plat): Min:", min(putative_plat_data$pident.Plat), " | Max:", max(putative_plat_data$pident.Plat)),
      paste("evalue (Plat): Min:", min(putative_plat_data$evalue.Plat), " | Max:", max(putative_plat_data$evalue.Plat)),
      "---------------",
      "",
      sep = "\n"
    )
    
    # Print summary
    cat(summary, sep = "\n")
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    
    
  } else {
    cat("Coun't identify putative platyhelminthes genes")
  }
  return(putative_plat_data)
}











#' Title
#' 
#' @param 
#' @param 
#' @return 
#' 
#' @example 
#' 
analyse_overlap <- function(data1, data2, name1, name2, column_name = "GeneID"){
  
  #test
  #data1 <- BpK_swissprot
  #data2 <- BpK_zebrafish
  #name1 <- "swissprot"
  #name2 <- "zebrafish"
  #column_name <- "GeneID"
  
  
  dim(data1) # 37824
  value1 <- unique(data1[[column_name]])
  num1 <- length(value1) # 37824
  num1
  head(value1)
  
  dim(data2) # 17932
  value2 <- unique(data2[[column_name]])
  num2 <- length(value2) # 17932
  num2
  
  # find common values (e.g., GenesID | GeneNameID)
  common <- intersect(value1, value2)
  num_common <- length(common)
  num_common
  
  percentage_in_1 <- 100*(num_common/num1)
  # i.e.,
  # se de 17932 ----- 17855 em comum
  # então   x   ----- 100
  # x = Cerca de x% dos valores de data1 também estão em data2
  
  # Q: qual porcentagem dos valores de data2 aparecem em data1 ???
  percentage_in_2 <- 100*(num_common/num2)
  # R: cerca de x% dos valores data2 também estão em data1
  
  
  
  summary <- paste(
    "--- SUMMARY ---",
    paste("Overlap analysis:", name1, "vs", name2),
    paste(name1, ":", num1),
    paste(name2, ":", num2),
    paste("Common:", num_common),
    paste("% of", name1, "in", name2, ":", round(percentage_in_1,2)),
    paste("% of", name2, "in", name1, ":", round(percentage_in_2,2)),
    sep = "\n"
  )
  
  # Print summary
  cat(summary)
  
  
  # Data frames for common and unique values
  
  # common
  common_df <- merge(data1, data2, by = column_name, suffixes = c(".1", ".2"))
  # unique to data1
  unique_to_1 <- anti_join(data1, data2, by = column_name)
  #  unique to data2
  unique_to_2 <- anti_join(data2, data1, by = column_name)
  #all rows from data2 whose value in column_name does not appear in data1
  
  return(list(common_values = common, common_df = common_df, 
              unique_to_1= unique_to_1,
              unique_to_2 = unique_to_2))
  
}





#' Title
#' 
#' @param 
#' @param 
#' @return 
#' 
#' @example 
#' 
analyse_taxonomic_representation <- function(data){
  
  taxa_occurrance <- list()
  
  if(is.list(data)){
    
    for (df_name in names(data)) {
      
      df <- data[[df_name]]  
      
      occurrance_df <- df %>%
        count(OS, name = "Count", sort = TRUE) %>%
        as.data.frame()
      
      num_os <- nrow(occurrance_df)
      
     
      top_taxa <- head(occurrance_df, 5)
      
      
      total_entries <- nrow(df)
      
     
      top_taxon <- top_taxa$OS[1]
      top_taxon_count <- top_taxa$Count[1]
      top_taxon_percentage <- round(100 * top_taxon_count / total_entries, 2)
      
      
      summary <- paste(
        "--- SUMMARY ---",
        paste0("Total entries: ", total_entries),
        paste0("Number of unique taxa (OS): ", num_os),
        paste0("Most represented taxon: ", top_taxon, " (", top_taxon_count, " entries, ", top_taxon_percentage, "%)"),
        "Top 5 taxa by representation:",
        paste0(
          apply(top_taxa, 1, function(row) {
            sprintf("  - %s: %s entries (%.2f%%)", 
                    row["OS"], 
                    row["Count"], 
                    100 * as.numeric(row["Count"]) / total_entries)
          }),
          collapse = "\n"
        ),
        sep = "\n"
      )
      
      cat(summary, "\n\n")
      
      taxa_occurrance[[df_name]] <- occurrance_df
    }
  }
  
  return(taxa_occurrance)
}

  













#' Title
#' 
#' 
#' 
#' @param count_data A dataframe containing count data with gene identifiers
#' @param model_type Statistical model ("QLF" or "LRT")
#' @param use_filter Logical for applying low expression filtering
#' @param use_lfc Logical for applying log-fold change threshold
#' 
#' @return
#' 
#' @references https://doi.org/10.12688/f1000research.8987.2
#'             edgeRUsersGuide()
#'                
#' 
#' @example 
#' 
perform_DEG <-function(data, model_type, use_filter = FALSE, use_lfc =FALSE){
  
  data <- data

  # constants
  p_value <- 0.05 
  log_fold_change <- 1.5
  adjust_method <- "fdr"
 
  
  experimental_groups <- c(
    "Sp_TR","Su_TR","Su_RA",
    "Sp_RF","Sp_TR","Sp_RA",
    "Su_RF" ,"Su_TR" ,"Sp_RA",
    "Sp_TR","Sp_RF","Sp_RA",
    "Su_RF","Su_TR","Sp_RF",
    "Su_RF","Su_RA","Su_RA")
   
  
  # data preparation -----------------------------------------------------------

  
  dge <- DGEList(counts = data[ ,17:34], 
                 group = experimental_groups, 
                 genes = data[ ,1])
  
  
  dge$samples
  
  
  # pre-processing -------------------------------------------------------------
  
  
  ##  Filtering to remove low counts
  if(use_filter){
 
    keep <- filterByExpr(dge, group = experimental_groups)
            # https://rdrr.io/bioc/edgeR/man/filterByExpr.html 
            # (BIG ACHO) same as:
            # rowSum(cpm(deg) >  10/min.library size in millions) >= num_replicates
    #table(keep)
    dge <- dge[keep, ,keep.lib.sizes = FALSE]
    
  }
  
  
  
  ##  Normalization for composition bias
      #  Normalization by trimmed mean of M values (TMM)
      #  eliminate composition biases between libraries
  dge <- calcNormFactors(dge)
  dge$samples
  
  
  
  # Model setup ----------------------------------------------------------------
  
  design <- model.matrix(~0+experimental_groups, data = dge$samples)
  colnames(design) <- levels(dge$samples$group)
  design
  
  
  # Dispersion estimation ------------------------------------------------------
  
  dge <- estimateDisp(dge, design)
  #plotBCV(dge)
  
  
  # Model fitting --------------------------------------------------------------
  
  # Model: Quasi-likelihood negative binomial generalized log-linear model
  # Fit  : quasi-likelihood F-test
  if (model_type == "QLF"){
    fit <- glmQLFit(dge, design, robust = TRUE)
    head(fit$coefficients)
    #plotQLDisp(fit)
  
  # Model: Negative binomial generalized log-linear model 
  # Fit  : likelihood ratio test  
  } else {
    fit <- glmFit(dge, design)
    head(fit$coefficients)
  }
  
  
  
  # Testing for differential expression ----------------------------------------
  
  
  result_list <-list()
  decideTest_list <- list()
  
  
  ## (2 seasons + 3 sites)
  ## seasons: spring | summer
  ## sites: Ria de Aveiro | Troia | Ria Formosa
  
  ## single contrasts: 
  ## season vs season in site
  ## site vs site in season
  
  # Define contrast pairs-------
  my_contrasts <- makeContrasts(
    Sp.RAvsTR = Sp_RA - Sp_TR, 
    Sp.RAvsRF = Sp_RA - Sp_RF,
    Sp.TRvsRF = Sp_TR - Sp_RF,
    
    Su.RAvsTR = Su_RA - Su_TR,
    Su.RAvsRF = Su_RA - Su_RF,
    Su.TRvsRF = Su_TR - Su_RF,
    
    RA.SpvsSu = Sp_RA - Su_RA,
    TR.SpvsSu = Sp_TR - Su_TR,
    RF.SpvsSu = Sp_RF - Su_RF,
    levels = design
  )
  # ----------------------------
  
  contrast_name <- colnames(my_contrasts)

  for(name in contrast_name){  
    
    cat("Performing contrast: ", name, "\n")
    
    if(model_type == "QLF"){
      result <- glmQLFTest(fit, contrast = my_contrasts[,name])
    } else {
      result <- glmRT(fit, contrast = my_contrasts[,name])
    }
    
    
    
    # Assessing Results  -------------------------------------------------------
    
    if(use_lfc){
      # Identify which genes are significantly differentially expressed for each
      # contrast from a fit object containing p-values and test statistics.
      summary_result <- summary(decideTests(result,
                                            p.value = p_value,
                                            lfc = log_fold_change,
                                            adjust.method = adjust_method))
      
      decideTest_result <- as.data.frame(decideTests(result,
                                                     p.value = p_value,
                                                     lfc = log_fold_change,
                                                     adjust.method = adjust_method))
    }else{
      summary_result <- summary(decideTests(result,
                                            p.value = p_value,
                                            adjust.method = adjust_method))
      
      decideTest_result <- as.data.frame(decideTests(result,
                                                     p.value = p_value,
                                                     adjust.method = adjust_method))
    }
    
    
    # Extracts the most differentially expressed genes 
    DEGs <- topTags(result,
                    n = nrow(result$table),
                    p.value = p_value,
                    adjust.method = adjust_method)$table
    
    
    
    

    
    # Processing Results  -----------------------------------------------------
    
    
    if(is.null(DEGs)){
      cat("No significant DEG's found", "\n")
      result_list[[name]] <- data.frame()
      
      
      decideTest_list[[name]] <- data.frame()
      decideTest <- cbind(result, decideTest_result)
      decideTest_with_data <- merge(decideTest, data, by.x = "genes", by.y="GeneID")
      decideTest_list[[name]] <- decideTest_with_data
      next
    }
    
    
    
    merge_degs_with_data <- merge(DEGs, data, by.x = "genes", by.y="GeneID")
    merge_degs_with_data <- merge_degs_with_data %>% arrange(desc(logFC))
    result_list[[name]] <- merge_degs_with_data
    
    
    decideTest <- cbind(result, decideTest_result)
    decideTest_with_data <- merge(decideTest, data, by.x = "genes", by.y="GeneID")
    decideTest_list[[name]] <- decideTest_with_data
    
    print(summary_result)
  }
  
  
  # ----------------------------------------------------------------------------
  # Print parameter summary ----------------------------------------------------
  
  # INCOMPLETO !!
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  
  return(list(DEG_result = result_list, decideTest_result = decideTest_list))
  
}








#' Title
#' 
#' 
#' 
#' @param data A dataframe containing count data with gene identifiers
#' @param model_type Statistical model ("QLF" or "LRT")
#' @param use_filter Logical for applying low expression filtering
#' @param use_lfc Logical for applying log-fold change threshold
#' 
#' @return
#' 
#' @example 
#' 
run_and_save_DEG <- function(data, blast_siglas, model_type, 
                             use_filter=TRUE,use_lfc = FALSE){
  
  
  # --- run deg ---
  results <- perform_DEG(data, 
                         model_type,
                         use_filter = use_filter,
                         use_lfc = use_lfc)
  
  deg <- results[["DEG_result"]]
  dt <- results[["decideTest_result"]]
  
  
  # --- save deg ---
  out_dir <- "./03-Output/01-DEG-Analysis/DEG-results"
  
  
  filter_label <- ifelse(use_filter, # test (use_filter = TRUE or FALSE)
                         "with-filterByExpr", # yes (use_filter = TRUE)
                         "without-filterByExpr") # no (use_filter = FALSE)
  
  save_to_excel(deg, file.path(out_dir, filter_label,
                               paste0("deg-", blast_siglas, ".xlsx")))
  
  save_to_excel(dt, file.path(out_dir, filter_label,
                               paste0("decideTest-", blast_siglas, ".xlsx")))
  
  return(list(DEG_result=deg, decideTest_result= dt))              
  
}








#' Title
#' 
#' 
#' 
#' @param 
#' 
#' @return
#' 
#' @example 
#' 
get_top10_degs <- function(result_list, use_filter = TRUE, blast_siglas){
  
  
  top10 <- list()
  
  for (contraste_name in names(result_list)) {
    
    contrast <- result_list[[contraste_name]]
    
    if (nrow(contrast) > 0) {
      
      cat("\nPerforming contrast: ", contraste_name, "\n")
      
      ## Top10 ? smallest p-value| FDR (based on statistical testing)
      
      # Top 10 over-expressed|upregulated (logFC > 0)
      top10_over <- contrast %>%
        filter(logFC > 0) %>%
        # arrange -> menor -> maior
        # arrange (desc) -> maior -> menor
        arrange(desc(PValue)) %>%
        head(10) # os maiores
      
      # Top 10 under-expressed|downregulated (logFC < 0)
      top10_under <- contrast %>%
        filter(logFC < 0) %>%
        # menor -> maior
        arrange(PValue) %>%
        head(10) # os menores 
      
      
      # Top 10 global (independent of expression direction)
      top10_global <- contrast %>%
        arrange(PValue) %>%
        head(10)
      
      min_fdr_gene <- top10_global %>%
        filter(FDR == min(FDR)) %>%
        pull(GeneNameID)
      
      max_fdr_gene <- top10_global %>%
        filter(FDR == max(FDR)) %>% # devolve a linha completa
        pull(GeneNameID) # devolve apenas o valor da coluna
      
      
      # Dataframe with the top 10 under- plus the top 10 overexpressed genes
      top10_over_top10_under <- rbind(top10_over, top10_under)
      top10[[contraste_name]] <- top10_over_top10_under
      
      
      
      # Print parameter summary ----------------------------------------------------
      summary_top10 <- paste(
        "--- SUMMARY ---",
        "",
        "--- total counts ---",
        paste("Total DEGs:", nrow(contrast)),
        paste("Total over:", nrow(contrast%>%filter(logFC>0))),
        paste("Total under:", nrow(contrast%>%filter(logFC<0))),
        "",
        "--- over-expressed (logFC > 0) ---",
        paste(top10_over$GeneNameID, collapse = " |"),
        "",
        "--- under-expressed (logFC < 0) ---",
        paste(top10_under$GeneNameID, collapse = " |"),
        "",
        "--- Global Top 10 (Independent of logFC) ---",
        paste(top10_global$GeneNameID, collapse = " |"),
        "",
        "--- FDR Summary ---",
        paste("FDR (Over): Min:", min(top10_over$FDR), " | Max:", max(top10_over$FDR)),
        paste("FDR (Under): Min:", min(top10_under$FDR), " | Max:", max(top10_under$FDR)),
        paste("FDR (Global): Min:", min(top10_global$FDR), " | Max:", max(top10_global$FDR)),
        # Nota.: The min /or max GeneNameID, coud be several, there is different
        #         geneNameID with the same FDR
        "-----------------------------------------------------------------------",
        "",
        sep = "\n"
      )
      # Print summary
      cat(summary_top10, sep ="\n\n")
      
      # -------------------------------------------------------------------------
      
      

    }

  }
  
  
  # --- save deg ---
  out_dir <- "./03-Output/01-DEG-Analysis/DEG-results"
  
  filter_label <- ifelse(use_filter,"with-filterByExpr","without-filterByExpr")
  
  save_to_excel(top10, file.path(out_dir, filter_label, paste0("top10-", blast_siglas, ".xlsx")))
  
  
  return(top10)
}  










#' Title
#' 
#' 
#' 
#' @param 
#' 
#' @return
#' 
#' @example 
#' 
get_GeneNameID_contrasts <- function(deg_list) {
  
  # Remove empty data frames from the list
  deg_list <- deg_list[sapply(deg_list, function(df) nrow(df) > 0)]
  
  # For each contrast, extract GeneNameID and rename the column to the contrast name
  gene_dfs <- lapply(names(deg_list), function(contrast) {
    deg_list[[contrast]] %>%
      select(GeneNameID) %>%
      rename(!!contrast := GeneNameID)
  })
  
  # Combine all the data frames into one (uma em cima da outra)
  contrasts_GeneNameID <- bind_rows(gene_dfs)
  
  # Convert to long format and summarize which contrasts each GeneNameID appears in
  result <- contrasts_GeneNameID %>%
    pivot_longer(
      cols = everything(),
      names_to = "contrast",
      values_to = "GeneNameID",
      values_drop_na = TRUE
    ) %>%
    group_by(GeneNameID) %>%
    summarize(
      contrast = toString(unique(contrast)),
      .groups = "drop"
    )
  
  
  
  # MEGA IMCOMPLETE!! I need to think better about this one...
  b <- a %>% count(contrast) %>% arrange(desc(n))
  
  
  return(result)
}
# mudar o nome, que não se percebe bem! li outra vez e já não sabia o que era








#' Title
#' 
#' 
#' 
#' @param 
#' 
#' @return
#' 
#' @note Pequeno ERRO!! not_in_common não se percebe bem, deveria mudar para algo
#' que não me esqueça o que é suposto ser...tipo "unique_to"
#' 
#' @example 
#' 
common_deg_SP_Z <- function(deg_SP, deg_Z) {
  
  
  merged_list <- list()  
  not_in_common_list_SP <- list()
  not_in_common_list_Z <- list()
  
  
  for (contrast_name_Z in names(deg_Z)) {
    for (contrast_name_SP in names(deg_SP)) {
      if (contrast_name_SP == contrast_name_Z) {
        
        contrast_z  <- deg_Z[[contrast_name_Z]]
        contrast_sp <- deg_SP[[contrast_name_SP]]
        
        # Check if either data frame is empty
        if (nrow(contrast_z) == 0 || nrow(contrast_sp) == 0) {
          next  
        }
        
        common_degs <- merge(contrast_sp, contrast_z, by = "GeneNameID",
                             suffixes = c(".SP", ".Z"))
        
        common_degs  <- common_degs  %>%
          select(GeneNameID,
                 logFC.SP,logFC.Z, 
                 PValue.SP, PValue.Z,
                 FDR.SP, FDR.Z,
                 pident.SP, pident.Z,
                 evalue.SP, evalue.Z,
                 OS.SP,OS.Z) %>%
          arrange(logFC.SP) %>% #(menor -> maior)
          mutate(
            pident_category.SP = evaluate_sequence_quality(pident.SP, "pident"),
            pident_category.Z = evaluate_sequence_quality(pident.Z, "pident"),
            
            evalue_category.SP = evaluate_sequence_quality(evalue.SP, "evalue"),
            evalue_category.Z = evaluate_sequence_quality(evalue.Z, "evalue"),
            
            FDR_category.SP = evaluate_sequence_quality(FDR.SP, "FDR"),
            FDR_category.Z = evaluate_sequence_quality(FDR.Z, "FDR"))
        
        
        # Genes NOT in common --------------------------------------------------
        not_in_common_SP <- anti_join(contrast_sp, contrast_z, by = "GeneNameID")
        not_in_common_Z  <- anti_join(contrast_z, contrast_sp, by = "GeneNameID")
        
        
        not_in_common_SP <- not_in_common_SP %>%
          select(GeneNameID,logFC, PValue, FDR, pident, evalue, OS) %>%
          arrange(logFC) %>%
          mutate(
            pident_category = evaluate_sequence_quality(pident, "pident"),
            evalue_category = evaluate_sequence_quality(evalue, "evalue"),
            FDR_category = evaluate_sequence_quality(FDR, "FDR"))
        
        not_in_common_Z <- not_in_common_Z %>%
          select(GeneNameID,logFC, PValue, FDR, pident, evalue, OS) %>%
          arrange(logFC) %>%
          mutate(
            pident_category = evaluate_sequence_quality(pident, "pident"),
            evalue_category = evaluate_sequence_quality(evalue, "evalue"),
            FDR_category = evaluate_sequence_quality(FDR, "FDR"))
        # ----------------------------------------------------------------------
        
        # Save the data frame in the list
        merged_list[[contrast_name_SP]] <- common_degs 
        not_in_common_list_SP[[contrast_name_SP]] <- not_in_common_SP
        not_in_common_list_Z[[contrast_name_Z]]   <- not_in_common_Z
        
        
        # Mini summary
        cat("\nContrast:", contrast_name_SP)
        cat("\ndeg swiss-prot:", nrow(contrast_sp))
        cat("\ndeg zebrafish:", nrow(contrast_z))
        cat("\nIn common: ", nrow(common_degs))
        cat("\nunder: ",nrow(common_degs %>% filter(logFC.SP < 0)))
        cat("\nover: ", nrow(common_degs %>% filter(logFC.SP > 0)),"\n")
        
      }
    }
  }
  
  
  return(list(
    common = merged_list,
    unique_SP = not_in_common_list_SP,
    unique_Z = not_in_common_list_Z
  ))
  
}








#' Get uniprot information
#' 
#'  
#' 
#' @param 
#' 
#' @return
#' 
#' @example 
#'
get_uniprot_taxainfo <- function(deg_list, specific_accession = NULL) {
  
  uniprot_info <- list()

  for (contrast_name in names(deg_list)) {
    
    contrast <- deg_list[[contrast_name]]
    
    if (nrow(contrast) > 0) {
      
      # Decide which accession(s) to use
      Accession <- if (!is.null(specific_accession)) specific_accession else contrast$Accession
      
      ## 1. Get Taxonomy Information
        TaxaObj <- GetNamesTaxa(Accession) 
        merge_data <- merge(contrast, TaxaObj, by.x="Accession", by.y="Entry", all.x = TRUE)
  
      uniprot_info[[contrast_name]] <- merge_data
    }
  }
  
  return(uniprot_info)
}



  
  




#'  
#' 
#' @param 
#' 
#' @return
#' 
#' @example 
#'
process_STRINGdata <- function(uniprot_taxaobj){
  
  preprocessing <- function(df){
    df %>%
      select("Gene.Names..primary.", "logFC")%>%
      arrange(logFC)}
  
  STRING_data <- lapply(uniprot_taxaobj, preprocessing)
  
  
  return(STRING_data)
}


#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
process_GSEAdata <- function(uniprot_taxaobj){
  
  preprocessing <- function(df){
    df %>%
      select("Accession", "logFC", "PValue")%>%
      arrange(PValue)}
  
  GSEA_data <- lapply(uniprot_taxaobj, preprocessing)
  
  
  return(GSEA_data)
}









#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
process_EnrichmentScores <- function(enrichment_scores){
  
  # convert to a data frame
  ES <- as.data.frame(enrichment_scores)

  # leadingEdge -> [1] "A8E7C5" "Q6PC64" -> list 
  # collapse leadingEdge list to string (to facilitate data manipulation)
  
  ES$leadingEdge <- sapply(ES$leadingEdge,
                           function(x) paste(x, collapse = ";"))
  
  # out: leadingEdge -> "A8E7C5;Q6PC64" -> string
  
  # select only the significant and order
  ES <- ES %>% 
    filter(padj < 0.05) %>%
    arrange(padj)
  
  # remove pathway prefix if present
  pathway_prefix <- "^\\(KEGG\\)"
  if ("pathway" %in% colnames(ES)){
    ES$pathway <- sub(pathway_prefix, "", ES$pathway)
  }
  
  return(ES)
}








#'
#' @param 
#' 
#' @return
#' 
#' @references https://www.bioconductor.org/packages/devel/bioc/vignettes/multiGSEA/inst/doc/multiGSEA.html
#' 
#' @example 
#'
perform_GSEA <- function(GSEA_data){
  
  databases <- c("kegg") # options: kegg, reactome, GO
  
  pathways <- getMultiOmicsFeatures(
    dbs = databases,
    layer = "transcriptome",
    returnTranscriptome = "UNIPROT", # options: SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ
    useLocal = FALSE,
    organism = "drerio"
  )
  
  pathways_short <- lapply(names(pathways), function(name){
    head(pathways[[name]], 2)
  })
  names(pathways_short) <- names(pathways)
  #pathways_short
  
  
  # ES - [E]nrichemnt [S]core
  raw_ES <- list()
  processed_ES <- list()
  
  
  for (contraste_name in names(GSEA_data)) {
    
    
    # test
    #contrast <- GSEA_data[["Sp.RAvsTR"]]
    
    contrast <- GSEA_data[[contraste_name]]
    
    if (nrow(contrast) > 0) {
      
      
      omics_data <- initOmicsDataStructure(layer = c("transcriptome"))
      layers <- names(omics_data)
      print(layers)
      
      # ranks
      omics_data$transcriptome <- rankFeatures(
        contrast$logFC,
        contrast$PValue
      )
      
      #print(head(omics_data$transcriptome))
      
      names(omics_data$transcriptome) <- contrast$Accession
      omics_data$transcriptome <- sort(omics_data$transcriptome)
      
      #print(head(omics_data$transcriptome))
      
      set.seed(42)
      
      enrichment_scores <- multiGSEA(pathways, 
                                     omics_data) # ranks
      
      #print(head(enrichment_scores))
      #View(enrichment_scores[["transcriptome"]])
      
      
      # Save the results
      
      
      if(is.null(enrichment_scores)){
        raw_ES[[contraste_name]] <- data.frame()
        processed_ES[[contraste_name]] <- data.frame()
        next
      }
      
      
      raw_ES[[contraste_name]] <- enrichment_scores
      
      processed_enrichment_scores <- process_EnrichmentScores(enrichment_scores$transcriptome)
      processed_ES[[contraste_name]] <- processed_enrichment_scores
      
    }
  
  }
  
  return(list(ES_raw_result = raw_ES, ES_processed_result= processed_ES))
}


























# ==============================================================================
# Visualizations 
# ==============================================================================

# references:
# -> https://ggplot2-book.org/layers.html
# -> https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/#figure-sizing
# -> https://showteeth.github.io/ggpie/articles/ggpie_manual.html

# ==============================================================================


list.files("C:/Windows/Fonts", pattern = "arial", ignore.case = TRUE)
font_add("Arial",
         regular = "C:/Windows/Fonts/arial.ttf",
         bold = "C:/Windows/Fonts/arialbd.ttf",
         italic = "C:/Windows/Fonts/ariali.ttf",
         bolditalic = "C:/Windows/Fonts/arialbi.ttf")

showtext_auto()


my_col <- c("Sp_RA" = "#E41A1C", "Sp_TR" = "#377EB8", "Sp_RF" = "#4DAF4A",
            "Su_RA" = "#984EA3", "Su_TR" = "#FF7F00", "Su_RF" = "#A65628")









# ==============================================================================
# Visualizations - Venn Diagram
# ==============================================================================

# references:
# -> https://r-graph-gallery.com/14-venn-diagramm
# -> https://r-charts.com/part-whole/ggvenndiagram/ 
# -> https://venn.bio-spring.top/using-ggvenndiagram
# ?ggVennDiagram
# ==============================================================================

#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
plot_vennDiagram <- function(vennData, save_plot = FALSE, output_file = NULL){
  
    
  venn_plot <- ggVennDiagram(
    vennData,                           # Data (list(set1,set2,set3))
    label_alpha = 0,                    # Remove the background from region labels
    edge_size = 0.8,                    # Set thickness of the set borders
    edge_lty = "solid",                 # Set -> set edges | Options: "dashed", "solid
    
    set_color = c("#0072b2", "#e69f00", "#56b4e9"), # Set the color of set borders and category names
    set_size = 6,                       # Set font size for category names
    
    
    label = "both",                     # Options: "both", , "count", "percent", "none"
    label_percent_digit = 1,            # Set decimal digits for percent labels
    # label_color = "firebrick",        # Set color for region label text
    label_size =  4.5                   # Set font size for region labels
  ) +
    
    # Expand x-axis to allow space for long set labels
    scale_x_continuous(
      expand = expansion(mult = .2)     # Add space (20%) to both sides of the x-axis
    ) +
    
    theme(
      text = element_text(family = "Arial", size = 6),
      legend.position = "none"
    )
  
  
  if(save_plot){
    if(is.null(output_file)) stop("Please provide output_file when save_plot = TRUE")
    print(venn_plot)
    ggsave(output_file, plot = venn_plot, width = 8.9, height = 8.9, units = "cm")
    
  } else {
    print(venn_plot)
  }
}

















# ==============================================================================
# Visualizations - Bar plot
# ==============================================================================



#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
process_barplot_data <- function(data, blast_name, metric) {
  
  data <- data %>%
    select(pident, evalue) %>%
    mutate(blast = blast_name)
  
  
  if (metric == "pident") {
    data <- data %>%
      mutate(category = evaluate_sequence_quality(pident, "pident"))
  } else if (metric == "evalue") {
    data <- data %>%
      mutate(category = evaluate_sequence_quality(evalue, "evalue"))
  } else {
    stop("Use 'pident' or 'evalue'.")
  }
  
  return(data)
}



#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
prepare_summary <- function(data_list){
  
  #test
  #data_list <- list("Swiss-prot" = barData_sp,
  #                    "Zebra-fish" = barData_z,
  #                    "Platyhelminth" = barData_p)
  
  all_data <- bind_rows(data_list)
  #View(all_data)
  
  print(table(all_data$category))
  
  summary_table <- all_data %>%
    group_by(blast, category) %>%
    summarise(count = n(),
              .groups = "drop") %>%
    group_by(blast) %>%
    mutate(percentage = count/sum(count)*100)
  
  #View(summary_table)
  
  return(summary_table)
  
}



#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
plot_barplot <- function(barData, save_plot = TRUE, output_file){
  
  barData$category <- factor(barData$category, levels = c("Bad", "Good", "Excellent"))
  
  bar_plot <- ggplot(
    barData,
    aes(x = category, y = percentage, fill = blast)
  ) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    geom_text(
      aes(label = sprintf("%.1f%%", percentage)),
      position = position_dodge(width = 0.8),
      vjust = -0.8,
      size = 2.0,
      family = "Arial"
      ) +
    labs(
      x = NULL, y = NULL,
      fill = "Database"
    ) +
    scale_y_continuous(
      labels = scales::percent_format(scale = 1),
      expand = expansion(mult = c(0, 0.15))
    ) +
    scale_fill_manual(
      values = c("platyhelminthes" = "#0072b2",
                 "swissprot" = "#e69f00",
                 "zebrafish" = "#56b4e9")
    ) +
    theme_minimal(base_family = "Arial") +
    theme(
      text = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 6),
      legend.position = "bottom",
      legend.background = element_rect(fill = "white", color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(20,20,20,20)
    )
    
  
  print(bar_plot)
  
  
  if(save_plot){
    if(is.null(output_file)) stop("Please provide output_file when save_plot = TRUE")
    ggsave(
      filename = output_file,
      plot = bar_plot,
      width = 9,
      height = 9,
      units = "cm",
      dpi = 600
    )
  }
  
}






# ==============================================================================
# Visualizations - Pie chart
# ==============================================================================




# missing! to be done....



# ==============================================================================
# Visualizations - Principal component analysis (PCA)
# ==============================================================================





#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
get_common_GeneNameID <- function(x, y){
  
  matched_GeneNameID <- merge(x, y, by = "GeneNameID", suffixes = c(".x", ".y")) 
  
  return(matched_GeneNameID)
}



#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
get_unique_GeneID <- function(common_GeneNameID){
  
  GeneID <- unique(c(common_GeneNameID$GeneID.x,common_GeneNameID$GeneID.y))
  
  return(GeneID)
  
}




#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
process_PCAdata <- function(data = NULL,
                            df1 = NULL,
                            df2 = NULL,
                            GeneID_vector = NULL,
                            tsv_data){
  
  
  
  # Case 1: Use all data 
  # e.g.,
  #   - tsv$abundance
  #   - tsv$count
  if (is.null(data) && is.null(df1) && is.null(df2) && is.null(GeneID_vector)){
    return(tsv_data)
  }
  
  
  
  
  
  # Case 2: Subset by GeneID from data frame
  # e.g.,
  #   - differentially expressed genes (under-, over- , No sig)
  #   - differentially expressed genes (under-, over-)
  if (!is.null(data) && is.list(data)) {
    
    ids_list <- lapply(data, function(df){
      if(nrow(df)>0) df$genes
    })
    
    ids <- unique(unlist(ids_list))
    length(ids)
    PCAdata <- tsv_data[rownames(tsv_data) %in% ids, ]
    
    return(PCAdata)
  }
  
  
  
  # Case 3: Subset by common GeneNameID from two data frames
  # e.g.,
  #   - fish genes
  if(!is.null(df1) && !is.null(df2)){
    
    common <- merge(df1, df2, by = "GeneNameID", suffixes = c(".x", ".y"))
    ids <- unique(c(common$GeneID.x, common$GeneID.y))
    PCAdata <- tsv_data[rownames(tsv_data) %in% ids, ]
    
    return(PCAdata)
  }
  
  
  
  
  
  # Case 4: Subset by provided GeneID vector
  #   - putative plat genes
  if(!is.null(GeneID_vector)){
    
    PCAdata <- tsv_data[rownames(tsv_data) %in% GeneID_vector, ]
    return(PCAdata)
  }
  
  
  return(NULL)
  
}
  
  


#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
runPCA <- function(PCAdata, add_pseudocount=TRUE,removeVar = NULL){
  
  # inspired by : PraticalTutorial04
  
  # pseudocount : we add 1 TPM to each gene to avoid infinite values after log
  # removeVar : Remove this % of variables based on low variance
  #
  
  # Prepare gene expression matrix
  ## Log transform the data
  if (add_pseudocount){
    logTPMs <- log2(PCAdata + 1)
  } else {
    logTPMs <- log2(PCAdata)
  }
  
  
  ## Remove duplicated genes
  uniqueGenes <- unique(rownames(logTPMs))
  logTPMs <- logTPMs[uniqueGenes, ]
  
  
  
  # Prepare metadata with sample type
  sampleTypes <- gsub("^[^\\.]+\\.", # pattern
                      "", # replacement
                      colnames(logTPMs) # x
  )
  
  metaData <- data.frame(sampleTypes)
  rownames(metaData) <- colnames(logTPMs)
  
  
  # Run pca 
  pca_result <- pca(mat = logTPMs,
                    metadata = metaData,
                    removeVar = removeVar) 
  
  
  return(list(pca_result = pca_result, logTPMs = logTPMs, metadata = metaData))
}


#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
plot_pca <- function(PCAdata, save_plot = TRUE, output_file){
  
  
  #test
  PCAdata <- pca_result$pca_result
  
  #View(PCAdata)
  
  # Prepare data to plot
  rotated_data <- as.data.frame(PCAdata$rotated)
  #View(rotated_data)
  
  rotated_data$sampleTypes <- PCAdata$metadata$sampleTypes
  #View(rotated_data)
  
  var_pc1 <- round(PCAdata$variance[1], 2)
  #View(var_pc1)
  var_pc2 <- round(PCAdata$variance[2], 2)
  # -------------------------------------------------------
  
  
  # Plot
  pca_plot <- ggplot(
    rotated_data,
    aes(x = PC1, y = PC2)
  ) +
    ggforce::geom_mark_ellipse(
      aes(fill = sampleTypes, color = sampleTypes),
      alpha = 0.2,
      size = 1
    ) +
    geom_point(
      aes(fill = sampleTypes),
      size = 4,
      shape = 21,
      stroke = 1.2,
      color = "black"
    ) +
    scale_fill_manual(
      values = my_col
    ) +
    scale_color_manual(
      values = my_col
    ) +
    labs(
      x = paste0("PC1 (", var_pc1, "% variation)"),
      y = paste0("PC2 (", var_pc2, "% variation)")
    ) +
    theme_bw(
      base_family = "Arial",
      base_size = 14
    ) +
    theme(
      text = element_text(family = "Arial", size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      legend.position = "right",
      legend.spacing.y = unit(0.8, "cm"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      panel.border = element_rect(linewidth = 1.2),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray90", linewidth = 0.3, linetype = "dashed")
      
    ) + 
    guides(
      fill = guide_legend(nrow = 2, byrow = TRUE),
      color = guide_legend(nrow = 2, byrow = TRUE)
    )
  
  
  
  print(pca_plot)
  
  
  
  if(save_plot){
    if(is.null(output_file)) stop("Please provide output_file when save_plot = TRUE")
    ggsave(
      filename = output_file,
      plot = pca_plot,
      width = 20,
      height = 20,
      units = "cm",
      dpi = 600
    )
  }
  
  
  
}





# Tue May 27 10:36:48 2025 ------------------------------
# Fail: loading plot is a mess, is not entirely working !! 

#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
process_PCAloadings <- function(pca_result, reference_data){
  
  PCAloadings <- pca_result$pca_result$loadings
  PCAloadings$GeneID <- rownames(PCAloadings)
  #View(PCAloadings)
  #dim(PCAloadings) # 2975
  
  #reference_data <- Data_SP
  GeneNameID <- reference_data %>% 
    select(GeneNameID, GeneID) 
  #View(GeneNameID)
  #dim(GeneNameID) # 15947
  
  # Note: I coud also use dyplr::left.join()
  loadings <- merge(GeneNameID, PCAloadings, by = "GeneID", all.y = TRUE)
  #View(loadings) 
  #dim(loadings) # 2975
  
  #print(sum(is.na(loadings$GeneNameID))) # 249
  
  # Forma 1:
  # trocar NA pelo GeneID
  loadings$GeneNameID[is.na(loadings$GeneNameID)] <- loadings$GeneID[is.na(loadings$GeneNameID)]
  
  # Forma 2:
  #loadings <- loadings %>%
  #  mutate(
  #    GeneNameID = if_else(is.na(GeneNameID), GeneID, GeneNameID))
  
  #View(loadings) 
  #dim(loadings) # 2975
  
  #print(sum(is.na(loadings$GeneNameID)))
  
  rownames(loadings) <- loadings$GeneNameID
  #View(loadings) 
  
  processed_PCAloadings <- loadings %>%
    select(-c(GeneNameID, GeneID))
  
  
  pca_result$pca_result$loadings <- processed_PCAloadings
  
    
  return(pca_result)
  
}

#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
plot_pcaloadings <- function(PCAdata, save_plot = FAlSE, output_file){
  
  
  
  
  loadings_plot <- plotloadings(pcaobj = pca_res$pca_result,
                                components = c("PC1","PC2","PC3"),
                                rangeRetain = 0.01, #default
                                shape = 21,
                                shapeSizeRange = c(10,10),
                                legendPosition = "top",
                                legendLabSize = 15,
                                legendIconSize = 3,
                                labSize = 5,
                                typeConnectors = "open",
                                endsConnectors = "last")
  
  print(loadings_plot)
  
  if(save_plot){
    if(is.null(output_file)) stop("Please provide output_file when save_plot = TRUE")
    ggsave(
      filename = output_file,
      plot = loadings_plot,
      width = 20,
      height = 20,
      units = "cm",
      dpi = 600
    )
  }
  
  
  
}
# -------------------------------------------------------


















# ==============================================================================
# Visualizations - Volcano plot
# ==============================================================================

# references:
# -> https://biostatsquid.com/volcano-plots-r-tutorial/
# ==============================================================================

# Note.: I need to confirm if the top5,  is the logFC, or the -log10 P-Value....
# I used the logFC instead of the -log10 P-Value, maybe is not the correct way
# I also, according to the thesis examples, need to change -log10 P-value to -log10 FDR 


#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
process_volcano_data <- function(dt){
  
  
  dt <- deg_SP_with_filter$decideTest_result$Sp.RAvsTR
  
  volcano_data <- dt %>%
    mutate(contrast = .[[6]]) %>%  
    select(GeneNameID, logFC, PValue, contrast) %>%
    mutate(
      logPValue = -log10(PValue),
      Expression = case_when(
        contrast == 1 & logFC > 0 ~ "Overexpressed",
        contrast == -1 & logFC < 0 ~ "Underexpressed",
        TRUE ~ "Not Significant"  
      )
    )
  
  
  # Fail: acho que há qualquer erro aqui, na seleção dos top...
  top5_over <- volcano_data %>%
    filter(Expression == "Overexpressed") %>%
    arrange(desc(PValue)) %>%
    slice_head(n = 5)
  
  top5_under <- volcano_data %>%
    filter(Expression == "Underexpressed") %>%
    arrange(PValue) %>%
    slice_head(n = 5)
  
  highlighted <- bind_rows(top5_under, top5_over)
  
  return(list(highlighted = highlighted, volcano_data = volcano_data))
}


#'
#' @param 
#' 
#' @return
#' 
#' @example 
#'
plot_volcanoPlot <- function(volcano_data, highlighted, save_plot = TRUE, output_file = NULL){
  # volcano_data: Full volcano data
  # highlighted: Subset for labeling
  
  volcano_plot <- ggplot(
    volcano_data, 
    aes(x = logFC, y = logPValue, color = Expression)
  ) +
    geom_point(
      alpha = 0.8, 
      size = 4.0
    ) +
    geom_hline(
      yintercept = -log10(0.05), 
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = c(-1, 1), 
      linetype = "dashed"
    ) +
    scale_color_manual(
      values = c("Overexpressed" = "tomato",
                 "Underexpressed" = "steelblue",
                 "Not Significant" = "gray")
    ) +
    labs(
      x = bquote(Log[2] ~ "fold change"),
      y = bquote(-Log[10] ~ "P-Value"),
      color = "Gene Expression"
    ) +
    theme_bw(base_size = 14) +
    theme(
      text = element_text(family = "Arial"),
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      panel.border = element_blank()
    ) +
    geom_label_repel(
      data = highlighted, 
      aes(label = GeneNameID),  
      size = 6.0, 
      max.overlaps = 100,
      fontface = "bold", 
      color = "black",
      box.padding = 1, 
      point.padding = 0.3,
      segment.color = "black", 
      segment.size = 1.0
    )
  
  print(volcano_plot)
  
  if(save_plot){
    if(is.null(output_file)) stop("Please provide output_file when save_plot = TRUE")
    ggsave(
      filename = output_file,
      plot = volcano_plot,
      width = 20,
      height = 20,
      units = "cm",
      dpi = 600
    )
  }
}









# ==============================================================================
# Visualizations - Heatmap 
# ==============================================================================




# missing! to be done....















