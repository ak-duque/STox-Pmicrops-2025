
# Title:

# Purpose:

# Project name: STox-Pmicrops-2025

# Author: A. Duque
# Contact details: duque.andrea2000@gmail.com

# Date script created: Thu Apr 24 11:22:37 2025 ------------------------------
# Date last modified: Thu Apr 24 11:22:37 2025 ------------------------------






# 1. Load Blastp files

swissprot <- load_blastp_data("01-Input/Blastp/ExtremeOceans_Blastp_SwissProt.csv")
dim(swissprot) # 161.919 15

zebrafish <- load_blastp_data("01-Input/Blastp/ExtremeOceans_Blastp_Zebrafish.csv")
dim(zebrafish) # 78.410 15

plat <- load_blastp_data("01-Input/Blastp/ExtremeOceans_Blastp_Plat.csv")
dim(plat) # 13.905 15






# 2. Load Kallisto files

## File paths to access each sample's quantification data
abundance_files <- file.path("./01-Input/Kallisto", sample_name, "abundance.tsv")
### output e.g.,
###
###    [1] "./01-Input/Kallisto/F14.Sp_TR/abundance.tsv"
###    [2] "./01-Input/Kallisto/F15.Su_TR/abundance.tsv"
###    [3] "./01-Input/Kallisto/F16.Su_RA/abundance.tsv"

names(abundance_files) <- sample_name


## Import kallisto files (1st time)

abundance_tsv <- tximport(
  abundance_files,
  type = "kallisto",
  txOut = TRUE,
  countsFromAbundance = "lengthScaledTPM"
  )


## Create tx2gene dataframe
### But why? 
###          Transcripts need to be associated with gene IDs for gene-level
###          summarization. If that information is present in the files, we can
###          skip this step. For Salmon, Sailfish, and kallisto the files only
###          provide the transcript ID. We first make a data.frame called tx2gene
###          with two columns: 1) transcript ID and 2) gene ID.
###          Since the transcript ID must be the same one used in the abundance
###          files, to know the transcript ID we need to import first, and only
###          then create the tx2gene dataframe.

###  Process Summary:
###          Import the abundance files first - to extract the correct transcript IDs.
###          Create tx2gene dataframe - pairing transcript IDs to their gene IDs.
###          Re-import abundance files properly - using the mapping in tx2gene.

tx2gene_df <-unlist(sapply(rownames(abundance_tsv$counts), 
                            function(x) unlist(strsplit(x, "i", fixed = TRUE))[1]))
### logic e.g.,
###
### TranscriptID - TRINITY_DN168445_c0_g1_i1 (i.e., [c]luster| [g]ene| [i]soform])
###
### sapply(rownames(abundance_tsv$counts) - For every transcript ID, the code:
### function(x) unlist(strsplit(x, "i", fixed = TRUE))[1])) - splits the word at "i"
###                                                           [TRINITY_DN168445_c0_g1_] + [i1]
###                                                           and select the first "[1]" part
### GenesID - TRINITY_DN168445_c0_g1



### data frame => | TranscriptID | GeneID |
tx2gene_df <-data.frame(rownames(abundance_tsv$counts), tx2gene_df)
colnames(tx2gene_df)<-c("TranscriptID", "GeneID")
print(head(tx2gene_df))


## Import Kallisto files (2nd time) but know with tx2gene @param
abundance_tsv <- tximport(
  abundance_files,
  type = "kallisto",
  tx2gene = tx2gene_df,
  countsFromAbundance = "lengthScaledTPM"
  )


## Explore the object created with tximport
names(abundance_tsv)
View(head(abundance_tsv$abundance)) # matrix containing TPM's
View(head(abundance_tsv$counts)) # matrix containing read counts
View(head(abundance_tsv$length))


save(abundance_tsv, file = "03-Output/01-DEG-Analysis/processed-data/abundance_tsv")






# 3. Merge Blastp(1.) + Kallisto(2.) counts

## Preparing data to merge | Purpose .: Be easy to merge

### Kallisto data - Create a new data frame that contains a "GeneID" column 
counts <- abundance_tsv$counts
counts <- cbind(as.data.frame(row.names(counts)), counts)
colnames(counts)[1]<-"GeneID"


### Blastp data - Add "GeneID" column + clean Blastp data
swissprot <- select_best_hits(swissprot)
dim(swissprot) # 37.826 16

zebrafish <- select_best_hits(zebrafish)
dim(zebrafish) # 17.933 16 

plat <- select_best_hits(plat)
dim(plat) # 3.743 16


## Merge (Blastp + kallisto counts) | [Bp]Blastp + [k]kallisto - BpK 
swissprot <- merge(swissprot, counts, by = "GeneID")
dim(swissprot) # 37.824 34

zebrafish <- merge(zebrafish, counts, by = "GeneID")
dim(zebrafish) # 17.932 34

plat <- merge(plat, counts, by = "GeneID")
dim(plat) # 3.743 34

## Save files
save_to_excel(swissprot, "03-Output/01-DEG-Analysis/processed-data/BpK-Swissprot.xlsx")
save_to_excel(zebrafish, "03-Output/01-DEG-Analysis/processed-data/BpK-Zebrafish.xlsx")
save_to_excel(plat, "03-Output/01-DEG-Analysis/processed-data/BpK-Plat.xlsx")
