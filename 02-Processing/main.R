
# Title:

# Purpose:

# Project name: STox-Pmicrops-2025

# Author: A. Duque
# Contact details: duque.andrea2000@gmail.com

# Date script created: Thu Apr 24 11:22:37 2025 ------------------------------
# Date last modified: Thu Apr 24 11:22:37 2025 ------------------------------




# ==============================================================================
# TASK 1 - Blablabla blablabla blabla
# ==============================================================================


# 1. Load Blastp files

swissprot <- load_blastp_data("01-Input/Blastp/ExtremeOceans_Blastp_SwissProt.csv")
dim(swissprot) # 161.919 15

zebrafish <- load_blastp_data("01-Input/Blastp/ExtremeOceans_Blastp_Zebrafish.csv")
dim(zebrafish) # 78.410 15

plat <- load_blastp_data("01-Input/Blastp/ExtremeOceans_Blastp_Plat.csv")
dim(plat) # 13.905 15


# Check for missing values
sum(is.na(swissprot))           # Total count of missing values
colSums(is.na(swissprot))       # Count per column
anyNA(swissprot)                # TRUE if any missing values exist




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

save(abundance_tsv, file = "03-Output/01-DEG-Analysis/analysis-ready-data/abundance_tsv")






# 3. Merge Blastp(1.) + Kallisto(2.) counts

## Preparing data to merge | Purpose .: Be easy to merge

### Kallisto data - Create a new data frame that contains a "GeneID" column 
counts <- abundance_tsv$counts
counts <- cbind(as.data.frame(row.names(counts)), counts)
colnames(counts)[1]<-"GeneID"


### Blastp data - Add "GeneID" column + clean Blastp data

swissprot <- select_best_hits(swissprot)
dim(swissprot) # 37.824 16

zebrafish <- select_best_hits(zebrafish)
dim(zebrafish) # 17.932 16 

plat <- select_best_hits(plat)
dim(plat) # 3.743 16


## Merge (Blastp + kallisto counts) | [Bp]Blastp + [k]kallisto - BpK 
BpK_swissprot <- merge(swissprot, counts, by = "GeneID")
dim(BpK_swissprot) # 37.824 34

BpK_zebrafish <- merge(zebrafish, counts, by = "GeneID")
dim(BpK_zebrafish) # 17.932 34

BpK_plat <- merge(plat, counts, by = "GeneID")
dim(BpK_plat) # 3.743 34




## Save files
save_to_excel(BpK_swissprot, "03-Output/01-DEG-Analysis/preprocessed-data/BpK-Swissprot.xlsx")
save_to_excel(BpK_zebrafish, "03-Output/01-DEG-Analysis/preprocessed-data/BpK-Zebrafish.xlsx")
save_to_excel(BpK_plat, "03-Output/01-DEG-Analysis/preprocessed-data/BpK-Plat.xlsx")




# ==============================================================================
# TASK 2 - Blablabla blablabla blabla
# ==============================================================================

# Data needed:
BpK_swissprot <- read.xlsx("03-Output/01-DEG-Analysis/preprocessed-data/BpK-Swissprot.xlsx")
BpK_zebrafish <- read.xlsx("03-Output/01-DEG-Analysis/preprocessed-data/BpK-Zebrafish.xlsx")
BpK_plat <- read.xlsx("03-Output/01-DEG-Analysis/preprocessed-data/BpK-Plat.xlsx")

load("03-Output/01-DEG-Analysis/analysis-ready-data/abundance_tsv")






# Prepare data to plot --------------


S <- BpK_swissprot %>%
  select(c("pident","evalue")) %>%
  mutate(blast = "Swissprot") # create a new column named blast
  
Z <- BpK_zebrafish %>%
  select(c("pident","evalue")) %>%
  mutate(blast = "Zebrafish")

P <- BpK_plat %>%
  select(c("pident","evalue")) %>%
  mutate(blast = "Platyhelminthes")

all_blast <- bind_rows(S,Z,P)
head(all_blast)
View(all_blast)


all_blast <- all_blast %>%
  mutate(
    category = case_when(
      pident < 30 ~ "Bad",
      pident >= 30 & pident <= 70 ~ "Good",
      pident > 70 ~ "Excellent"
    )
  )

head(all_blast)
View(all_blast)

print(table(all_blast$category))

summary <- all_blast %>%
  count(blast, category)
print(summary)




# Data prepared --------------




Data_SP <- process_BpK_data(BpK_swissprot) #15947
save_to_excel(Data_SP, "./03-Output/01-DEG-Analysis/analysis-ready-data/BpK-data-processed-SP.xlsx")

Data_Z <- process_BpK_data(BpK_zebrafish) #2975
save_to_excel(Data_Z, "./03-Output/01-DEG-Analysis/analysis-ready-data/BpK-data-processed-Z.xlsx")

Data_P <- process_BpK_data(BpK_plat) #145
save_to_excel(Data_P, "./03-Output/01-DEG-Analysis/analysis-ready-data/BpK-data-processed-P.xlsx")


# Prepare data to plot --------------
S2 <- Data_SP %>%
  select(c("pident","evalue")) %>%
  mutate(blast = "Swissprot") # create a new column named blast

Z2 <- Data_Z %>%
  select(c("pident","evalue")) %>%
  mutate(blast = "Zebrafish")

P2 <- Data_P %>%
  select(c("pident","evalue")) %>%
  mutate(blast = "Platyhelminthes")

all_blast2 <- bind_rows(S2,Z2,P2)
head(all_blast2)
View(all_blast2)


all_blast2 <- all_blast2 %>%
  mutate(
    category = case_when(
      pident < 30 ~ "Bad",
      pident >= 30 & pident <= 70 ~ "Good",
      pident > 70 ~ "Excellent"
    )
  )

head(all_blast2)
View(all_blast2)

print(table(all_blast2$category))


summary2 <- all_blast2 %>%
  group_by(blast, category) %>% # agrupar
  summarise(count = n(), .groups = "drop") %>%  #.groups = "drop" serve para desagrupar
  group_by(blast) %>%
  mutate(percentage = count / sum(count) * 100)


ggplot(summary2, aes(x = category, y = percentage, fill = blast)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Percentagem das categorias pident por base de dados",
    x = "Categoria",
    y = "Percentagem (%)",
    fill = "Base de dados"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))


ggplot(summary2, aes(x = category, y = percentage, fill = blast)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(
    title = "Percentagem das categorias pident por base de dados",
    x = "Categoria",
    y = "Percentagem (%)",
    fill = "Base de dados"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))

# Data prepared --------------






# Organism part is also missing!
# Add PCA's part!! Is missing! 





# ==============================================================================
# TASK 3 - DEG analysis
# ==============================================================================
# Perform DEG analysis for swiss-prot and zebrafish datasets (with) and (without) 
# filterByExpr filter
# - explain what is done
# - explain what is done
# - explain what is done
# ------------------------------------------------------------------------------



# --- DEG Analysis for swiss-prot ----------------------------------------------
## -- filter by expression (with) ----------------------------------------------


deg_SP_with_filter <-run_and_save_DEG(Data_SP,
                                      "SP",
                                      model_type = "QLF" ,
                                      use_filter = TRUE,
                                      use_lfc = FALSE)


# Available results and contrasts
names(deg_SP_with_filter) 
names(deg_SP_with_filter$DEG_result) 
names(deg_SP_with_filter$decideTest_result) 

# Example: see specific contrast and columns
deg_SP <- deg_SP_with_filter$DEG_result
Sp.RAvsTR <- deg_SP[["Sp.RAvsTR"]] 
View(Sp.RAvsTR)

colnames(Sp.RAvsTR) 

View(Sp.RAvsTR[, c("logFC", "FDR", "evalue", "pident")]) 

Sp.RAvsTR <- Sp.RAvsTR %>%
  # from lowest FDR to highest
  arrange(FDR)  

View(Sp.RAvsTR)


# Get top 10 + summary
top10_SP_with_filter <- get_top10_degs(deg_SP, use_filter = TRUE, "SP") 



# --- DEG Analysis for swiss-prot ----------------------------------------------
## -- filter by expression (without) -------------------------------------------

deg_SP_without_filter <-run_and_save_DEG(Data_SP,
                                      "SP",
                                      model_type = "QLF",
                                      use_filter = FALSE,
                                      use_lfc = FALSE)

# Get top 10 + summary
top10_SP_without_filter <- get_top10_degs(deg_SP_without_filter$DEG_result, 
                                       use_filter = FALSE, "SP") 



# --- DEG Analysis for zebrafish -----------------------------------------------
## -- filter by expression (with) ----------------------------------------------

deg_Z_with_filter <-run_and_save_DEG(Data_Z,
                                      "Z",
                                      model_type = "QLF",
                                      use_filter = TRUE,
                                      use_lfc = FALSE)


# Get top 10 + summary
top10_Z_with_filter <- get_top10_degs(deg_Z_with_filter$DEG_result, 
                                       use_filter = TRUE, "Z") 



# --- DEG Analysis for zebrafish -----------------------------------------------
## -- filter by expression (without) -------------------------------------------

deg_Z_without_filter<-run_and_save_DEG(Data_Z,
                                       "Z",
                                       model_type = "QLF",
                                       use_filter = FALSE,
                                       use_lfc = FALSE)

# Get top 10 + summary
top10_Z_without_filter <- get_top10_degs(deg_Z_without_filter$DEG_result, 
                                      use_filter = FALSE, "Z") 
