
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
library(dplyr)
library(ggplot2)
library(stringr)

## BiocManager
library(tximport) 
library(limma)
library(edgeR) # require package: limma


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
  example <- blast_data %>% 
    filter(GeneID == "TRINITY_DN100212_c0_g1_") %>%
    select(c("qseqid","GeneID","evalue","pident"))
  #View(example)
  
  # Order: lowest evalue and highest pident
  blast_data <- blast_data[order(blast_data$evalue, -blast_data$pident),]
  
  # Best hits (first row)
  best_hits <- blast_data[match(unique(blast_data$GeneID), blast_data$GeneID), ]
  
  # Debugging step
  # Example: TRINITY_DN0_c10_g1_i1
  example_best_hit <- best_hits %>% 
    filter(GeneID == "TRINITY_DN100212_c0_g1_") %>%
    select(c("qseqid","GeneID","evalue","pident"))
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
      addWorksheet(work_book, name)
      writeData(work_book, name, data[[name]])
    }
  }
  
  saveWorkbook(work_book, file = output_file, overwrite = TRUE)
  cat("Excel created successfully:", output_file, "\n")
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
  example <- BpK_data %>% 
    filter(GeneNameID== "ABCA1") %>%
    select(c("GeneNameID","evalue","pident"))
  #View(example)
  
  
  
  
  # Order: lowest evalue and highest pident
  BpK_data <- BpK_data[order(BpK_data$evalue, -BpK_data$pident), ]
  
  # Select the best hit (1st row)
  BpK_data <- BpK_data[match(unique(BpK_data$GeneNameID), BpK_data$GeneNameID), ]
  #dim(BpK_data)
  
  
  
  
  # Let's confirm we select the right row
  # Debugging step
  # Example: ABCA1  (SwissProt)
  example1 <- BpK_data %>% 
    filter(GeneNameID== "ABCA1") %>%
    select(c("GeneNameID","evalue","pident"))
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
      value > 1e-10 ~ "Good",
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
  } else if (metric == "pvalue") {
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
    table(keep)
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
  
  
  if (model_type == "QLF"){
    fit <- glmQLFit(dge, design, robust = TRUE)
    head(fit$coefficients)
    plotQLDisp(fit)
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
    
    
    
    
    if(is.null(DEGs)){
      cat("No significant DEG's found", "\n")
      result_list[[name]] <- data.frame()
      decideTest_list[[name]] <- data.frame()
      next
    }
    
    
    # Processing Results  -----------------------------------------------------
    
    merge_degs_with_data <- merge(DEGs, data, by.x = "genes", by.y="GeneID")
    result_list[[name]] <- merge_degs_with_data
    
    decideTest_result <- cbind(result, decideTest_result)
    merge_decideTest_with_data <- merge(decideTest_result, data, by.x = "genes", by.y="GeneID")
    decideTest_list[[name]] <- merge_decideTest_with_data
    
    print(summary_result)
  }
  
  
  # ----------------------------------------------------------------------------
  # Print parameter summary ----------------------------------------------------
  
  # INCOMPLETO !!
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  
  return(list(DEG_result = result_list, decideTest_result=decideTest_list))
  
}














