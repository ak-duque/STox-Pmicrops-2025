
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




# Process data: 
Data_SP <- process_BpK_data(BpK_swissprot) #15947
dim(Data_SP)
#save_to_excel(Data_SP, "./03-Output/01-DEG-Analysis/analysis-ready-data/BpK-data-processed-SP.xlsx")

Data_Z <- process_BpK_data(BpK_zebrafish) #2975
dim(Data_Z)
#save_to_excel(Data_Z, "./03-Output/01-DEG-Analysis/analysis-ready-data/BpK-data-processed-Z.xlsx")

Data_P <- process_BpK_data(BpK_plat) #145
dim(Data_P)
#save_to_excel(Data_P, "./03-Output/01-DEG-Analysis/analysis-ready-data/BpK-data-processed-P.xlsx")




# ==============================================================================
# TASK 2.1 - Global annotation performance 
# ==============================================================================


# 2.1.1 e-value and pident behavior (before and after processing)


## pident raw data :

#barData_swissprot <- process_barplot_data(swissprot, "swissprot", metric = "pident")
#barData_zebrafish <- process_barplot_data(zebrafish, "zebrafish", metric = "pident")
#barData_plat <- process_barplot_data(plat, "platyhelminthes", metric = "pident")

#data_list <- list("Swiss-prot" = barData_swissprot, 
#                  "Zebrafish" = barData_zebrafish, 
#                  "Platyhelminth" = barData_plat)

#barData <- prepare_summary(data_list)

## Bar plot
#plot_barplot(barData, save_plot = TRUE, 
#             output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/barPlot-pident-add-geneID.pdf" )


#plot_barplot(barData, save_plot = TRUE, 
#             output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/barPlot-pident-raw-data.pdf" )




## Pident before processing :

# Prepare each data set
barData_BpKsp <- process_barplot_data(BpK_swissprot, "swissprot", metric = "pident")
barData_BpKz <- process_barplot_data(BpK_zebrafish, "zebrafish", metric = "pident")
barData_BpKp <- process_barplot_data(BpK_plat, "platyhelminthes", metric = "pident")

data_list <- list("Swiss-prot" = barData_BpKsp, 
                  "Zebrafish" = barData_BpKz, 
                  "Platyhelminth" = barData_BpKp)

barData <- prepare_summary(data_list)

## Bar plot
plot_barplot(barData, save_plot = TRUE, 
             output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/barPlot-pident-before-processing.pdf" )


## Pident after processing :

# Prepare each data set
barData_sp <- process_barplot_data(Data_SP, "swissprot", metric = "pident")
barData_z <- process_barplot_data(Data_Z, "zebrafish", metric = "pident")
barData_p <- process_barplot_data(Data_P, "platyhelminthes", metric = "pident")

data_list <- list("Swiss-prot" = barData_sp,
                  "Zebrafish" = barData_z, 
                  "Platyhelminth" = barData_p)

barData <- prepare_summary(data_list)

## Bar plot
plot_barplot(barData, save_plot = TRUE, 
             output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/barPlot-pident-after-processing.pdf" )







## E-value raw data :

#barData_swissprot <- process_barplot_data(swissprot, "swissprot", metric = "evalue")
#barData_zebrafish <- process_barplot_data(zebrafish, "zebrafish", metric = "evalue")
#barData_plat <- process_barplot_data(plat, "platyhelminthes", metric = "evalue")

#data_list <- list("Swiss-prot" = barData_swissprot, 
#                  "Zebrafish" = barData_zebrafish, 
#                  "Platyhelminth" = barData_plat)

#barData <- prepare_summary(data_list)

## Bar plot

#plot_barplot(barData, save_plot = TRUE, 
#             output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/barPlot-evalue-add-geneID.pdf" )


#plot_barplot(barData, save_plot = TRUE, 
#             output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/barPlot-evalue-raw-data.pdf" )




## E-value before processing :

barData_BpKsp <- process_barplot_data(BpK_swissprot, "swissprot", metric = "evalue")
barData_BpKz <- process_barplot_data(BpK_zebrafish, "zebrafish", metric = "evalue")
barData_BpKp <- process_barplot_data(BpK_plat, "platyhelminthes", metric = "evalue")

data_list <- list("Swiss-prot" = barData_BpKsp, 
                  "Zebrafish" = barData_BpKz, 
                  "Platyhelminth" = barData_BpKp)

barData <- prepare_summary(data_list)

## Bar plot
plot_barplot(barData, save_plot = FALSE)

plot_barplot(barData, save_plot = TRUE, 
             output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/barPlot-evalue-before-processing.pdf" )




## E-value after processing :

barData_sp <- process_barplot_data(Data_SP, "swissprot", metric = "evalue")
barData_z <- process_barplot_data(Data_Z, "zebrafish", metric = "evalue")
barData_p <- process_barplot_data(Data_P, "platyhelminthes", metric = "evalue")

data_list <- list("Swiss-prot" = barData_sp,
                  "Zebrafish" = barData_z, 
                  "Platyhelminth" = barData_p)

barData <- prepare_summary(data_list)

## Bar plot
plot_barplot(barData, save_plot = FALSE)

plot_barplot(barData, save_plot = TRUE, 
             output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/barPlot-evalue-after-processing.pdf" )









# 2.1.2 Common GeneNameID's on the different BLAST's outputs

# Before processing
Overlap_BpKSp_and_BpKZ <- analyse_overlap(BpK_swissprot, BpK_zebrafish,"swissprot", "zebrafish")
Overlap_BpKSp_and_BpKP <- analyse_overlap(BpK_swissprot, BpK_plat,"swissprot", "plat")
Overlap_BpKZ_and_BpKP <- analyse_overlap(BpK_zebrafish, BpK_plat,"zebrafish", "plat")

# After processing 
Overlap_Sp_and_Z <- analyse_overlap(Data_SP, Data_Z,"swissprot", "zebrafish")
Overlap_Sp_and_P <- analyse_overlap(Data_SP, Data_P,"swissprot", "plat")
Overlap_Z_and_P <- analyse_overlap(Data_Z, Data_P,"zebrafish", "plat")


## Extra 
# We notice that : common GeneID != common GeneNameID 
# WHY? I don't understand why this happen...

GeneID_Overlap_Sp_and_P <- analyse_overlap(Data_SP, Data_P,"swissprot", "plat",
                                           column_name = "GeneID")
GeneNameID_Overlap_Sp_and_P <- analyse_overlap(Data_SP, Data_P,"swissprot", "plat",
                                               column_name = "GeneNameID")

## Venn Diagram
# Note.: I know it doesn't follow the code good practice rules, but is good enough...
vennData <- list("Swiss-prot" = BpK_swissprot$GeneID,
                 "Zebra-fish" = BpK_zebrafish$GeneID,
                 "Platyhelminth" = BpK_plat$GeneID)

plot_vennDiagram(vennData, save_plot=TRUE, 
                 output_file="./03-Output/01-DEG-Analysis/annotation-results/figures/vennDiagram-BpKcommonGeneID-before-processing.pdf")


# ==============================================================================
# TASK 2.2 - Taxonomic representation \ Putative parasitic genes 
# ==============================================================================




# 2.2.1 Taxonomic representation

taxData <- list(SP_Organisms = Data_SP, P_Organisms = Data_P, Z_Organisms = Data_Z)

taxResults <- analyse_taxonomic_representation(taxData)
#save_to_excel(taxResults, 
#              output_file = "./03-Output/01-DEG-Analysis/annotation-results/tables/taxonomic-representation.xlsx")




## Pie chart




# 2.2.2 Fish organisms (Is not 100% correct, but is good enough)

# fish species: ----
fish <- c('Acipenser baerii',
          'Acipenser transmontanus',
          'Anguilla anguilla',
          'Anguilla japonica',
          'Anoplopoma fimbria',
          'Argyrosomus regius',
          'Boreogadus saida',
          'Brachyopsis segaliensis',
          'Carassius auratus',
          'Catostomus commersonii',
          'Chaenocephalus aceratus',
          'Chauliodus sloani',
          'Chelon auratus',
          'Chelon ramada',
          'Chionodraco hamatus',
          'Ctenopharyngodon idella',
          'Cynoscion nebulosus',
          'Cyprinus carpio',
          'Danio rerio',
          'Devario aequipinnatus',
          'Dicentrarchus labrax',
          'Diplobatis ommata',
          'Dissostichus eleginoides',
          'Electrophorus electricus',
          'Epinephelus akaara',
          'Epinephelus coioides',
          'Esox lucius',
          'Formosania lacustris',
          'Fundulus heteroclitus',
          'Gadus morhua',
          'Galeus melastomus',
          'Gasterosteus aculeatus',
          'Gillichthys mirabilis',
          'Gillichthys seta',
          'Haplochromis burtoni',
          'Haplochromis nubilus',
          'Haplochromis xenognathus',
          'Harpagifer antarcticus',
          'Hemitripterus americanus',
          'Heterodontus francisci',
          'Heteropneustes fossilis',
          'Hippocampus comes',
          'Hippoglossus hippoglossus',
          'Ichthyomyzon unicuspis',
          'Ictalurus punctatus',
          'Katsuwonus pelamis',
          'Labeo rohita',
          'Lepomis macrochirus',
          'Leucoraja ocellata',
          'Lophius americanus',
          'Makaira nigricans',
          'Megalobrama amblycephala',
          'Meiacanthus atrodorsalis',
          'Micropogonias undulatus',
          'Micropterus salmoides',
          'Misgurnus fossilis',
          'Monopterus albus',
          'Mullus surmuletus',
          'Myoxocephalus octodecemspinosus',
          'Notemigonus crysoleucas',
          'Oncorhynchus keta',
          'Oncorhynchus kisutch',
          'Oncorhynchus masou',
          'Oncorhynchus mykiss',
          'Oncorhynchus nerka',
          'Oncorhynchus tshawytscha',
          'Oplegnathus fasciatus',
          'Opsanus tau',
          'Oreochromis mossambicus',
          'Oreochromis niloticus',
          'Oryzias javanicus',
          'Oryzias latipes',
          'Oryzias luzonensis',
          'Osmerus mordax',
          'Pagrus major',
          'Paralichthys olivaceus',
          'Paramisgurnus dabryanus',
          'Petromyzon marinus',
          'Pleuronectes platessa',
          'Poecilia reticulata',
          'Poeciliopsis lucida',
          'Pomatoschistus minutus',
          'Prionace glauca',
          'Protopterus aethiopicus',
          'Psalidodon fasciatus',
          'Pseudopleuronectes americanus',
          'Salmo salar',
          'Salmo trutta',
          'Salvelinus fontinalis',
          'Scomber japonicus',
          'Scomber scombrus',
          'Scomberomorus niphonius',
          'Scorpaena plumieri',
          'Scyliorhinus stellaris',
          'Seriola quinqueradiata',
          'Siganus canaliculatus',
          'Siniperca chuatsi',
          'Solea senegalensis',
          'Sparus aurata',
          'Squalus acanthias',
          'Synanceia horrida',
          'Synanceia verrucosa',
          'Takifugu pardalis',
          'Takifugu rubripes',
          'Tetraodon fluviatilis',
          'Tetraodon miurus',
          'Tetraodon nigroviridis',
          'Tetronarce californica',
          'Thalassophryne nattereri',
          'Thunnus obesus',
          'Torpedo marmorata',
          'Trachurus japonicus',
          'Trematomus bernacchii',
          'Trichopodus trichopterus',
          'Xiphophorus hellerii',
          'Xiphophorus maculatus')

# --------------------------

fish_SP <- Data_SP %>% filter(OS %in% fish)

fish_occurrance <- fish_SP %>%
  count(OS, name = "Count", sort = TRUE) %>%
  as.data.frame()

View(fish_occurrance)

dim(fish_SP) # 2420
View(fish_SP)

fishGeneNameID <- fish_SP$GeneNameID
fishGeneID <- fish_SP$GeneID
length(fishGeneID)

common_GeneNameID_fishSP_P <- get_common_GeneNameID(fish_SP, Data_P)  
nrow(common_GeneNameID_fishSP_P) # 9
View(common_GeneNameID_fishSP_P)

common_GeneID_fishSP_P <- get_unique_GeneID(common_GeneNameID_fishSP_P)
length(common_GeneID_fishSP_P) # 17

fish_SP <- fish_SP[!fish_SP$GeneID %in% common_GeneID_fishSP_P, ] 
dim(fish_SP) # 2410 






# 2.2.3 Putative parasitic genes

platGenes <- identify_putative_plat_genes(Data_SP, Data_P)
View(platGenes)

platGeneNameID <- platGenes$GeneNameID

# Fail: Para o mesmo GeneName, em BLASTs diferentes, temos diferentes GeneID
# Solution: Para evitar perder info, selecionamos todos por via das dúvidas
platGeneID <- unique(c(platGenes$GeneID.Ref, platGenes$GeneID.Plat))







# ==============================================================================
# TASK 2.3 - Principal Component Analysis (PCA)
# ==============================================================================
# PCA to be done:
# - with all abundance/counts (2)
# - abundance/counts subsets:
#       raw -----------------------
#       - platGeneID (2)
#       - swiss-prot fishGeneID (2)
#       - all swiss-prot (2)
#
#       deg (under- over-) --------
#       with filterByExpr ---------
#       - all zebrafish (2)
#       - all swiss-prot (2)
#
#       without filterByExpr ------
#       - all zebrafish (2)
#       - all swiss-prot (2)
#
#       dt (under- over- noSig)----
#       with filterByExpr ---------
#       - all zebrafish (2)
#       - all swiss-prot (2)
#
#       without filterByExpr ------
#       - all zebrafish (2)
#       - all swiss-prot (2)
# ------------------------------------------------------------------------------


# Data needed:
load("03-Output/01-DEG-Analysis/analysis-ready-data/abundance_tsv")
View(abundance_tsv)

# PCA data (subsets of counts \or abundance)
# tsv_datatype
# counts : 
counts <- abundance_tsv$counts
# abundance : 
abundance <- abundance_tsv$abundance





# 2.3.1 All count
PCAdata <- process_PCAdata(tsv_data = counts)
dim(PCAdata) # 266105
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-all-count.pdf")

# it's not working correctly ......... for all the cases....
#pca_res <- process_PCAloadings(pca_result, Data_SP)
#plot_pcaloadings(pca_res, save_plot = FALSE)







# 2.3.2 All abundance
PCAdata <- process_PCAdata(tsv_data = abundance)
dim(PCAdata) # 266105
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-all-abundance.pdf")


# 2.3.3 Putative plat genes | counts
PCAdata <- process_PCAdata(GeneID_vector = platGeneID, tsv_data = counts)
dim(PCAdata) # 10
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-putative-plat-genes-count.pdf")


# 2.3.4 Putative plat genes | abundance
PCAdata <- process_PCAdata(GeneID_vector = platGeneID, tsv_data = abundance)
dim(PCAdata) # 10
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-putative-plat-genes-abundance.pdf")



# 2.3.5 Fish genes | counts
PCAdata <- process_PCAdata(df1 = fish_SP, df2 = Data_Z, tsv_data = counts)
dim(PCAdata) # 2254
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-fish-genes-count.pdf")


# 2.3.6 Fish genes | abundance
PCAdata <- process_PCAdata(df1 = fish_SP, df2 = Data_Z, tsv_data = abundance)
dim(PCAdata) # 2254
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-fish-genes-count.pdf")


# to be continued .....




# ==============================================================================
# TASK 2.4 - Remove putative plat genes
# ==============================================================================

dim(Data_SP) # 15947
Data_SP <- Data_SP[!Data_SP$GeneID %in% platGeneID, ]
dim(Data_SP) # 15937
save_to_excel(Data_SP, "./03-Output/01-DEG-Analysis/analysis-ready-data/BpK-data-processed-SP[v2].xlsx")


dim(Data_Z) # 2975
Data_Z <- Data_Z[!Data_Z$GeneID %in% platGeneID, ] 
dim(Data_Z) # 2974
save_to_excel(Data_Z, "./03-Output/01-DEG-Analysis/analysis-ready-data/BpK-data-processed-Z[v2].xlsx")









# ==============================================================================
# TASK 3 - DEG analysis
# ==============================================================================
# Perform DEG analysis for swiss-prot and zebrafish datasets (with) and (without) 
# filterByExpr filter
# - explain what is done
# - explain what is done
# - explain what is done
# ------------------------------------------------------------------------------


# Data needed:

Data_SP <- read.xlsx("./03-Output/01-DEG-Analysis/analysis-ready-data/BpK-data-processed-SP[v2].xlsx")
dim(Data_SP)

Data_Z <- read.xlsx("./03-Output/01-DEG-Analysis/analysis-ready-data/BpK-data-processed-Z[v2].xlsx")
dim(Data_Z)

load("03-Output/01-DEG-Analysis/analysis-ready-data/abundance_tsv")
# counts : 
counts <- abundance_tsv$counts
# abundance : 
abundance <- abundance_tsv$abundance

# ==============================================================================
# TASK 3.1 - swiss-prot (with/without filterByExpr)
# ==============================================================================


# 3.1.1 with filterByExpr ----


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




# Get top 10 
top10_SP_with_filter <- get_top10_degs(deg_SP, use_filter = TRUE, "SP") 

# Get UniProt information (degs only)








# 3.1.1 without filterByExpr ----


deg_SP_without_filter <-run_and_save_DEG(Data_SP,
                                      "SP",
                                      model_type = "QLF",
                                      use_filter = FALSE,
                                      use_lfc = FALSE)

# Get top 10 + summary
top10_SP_without_filter <- get_top10_degs(deg_SP_without_filter$DEG_result, 
                                       use_filter = FALSE, "SP") 







# 3.1.2 continuation .... Principal Component Analysis 


# 3.1.2.1 degs | with filterByExpr | counts
deg_SP <- deg_SP_with_filter$DEG_result

PCAdata <- process_PCAdata(data = deg_SP, tsv_data = counts)
dim(PCAdata) # 1625
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-degs-with-filterByExpr-count-SP.pdf")

# 3.1.2.2 degs | with filterByExpr | abundance
PCAdata <- process_PCAdata(data = deg_SP, tsv_data = abundance)
dim(PCAdata) # 1625
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-degs-with-filterByExpr-abundance-SP.pdf")


# 3.1.2.3 degs | without filterByExpr | counts
deg_SP <- deg_SP_without_filter$DEG_result

PCAdata <- process_PCAdata(data = deg_SP, tsv_data = counts)
dim(PCAdata) # 1298
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-degs-without-filterByExpr-count-SP.pdf")

# 3.1.2.4 degs | without filterByExpr | abundance
PCAdata <- process_PCAdata(data = deg_SP, tsv_data = abundance)
dim(PCAdata) # 1298
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-degs-without-filterByExpr-abundance-SP.pdf")



# 3.1.2.5 dt | with filterByExpr | counts
dt_SP <- deg_SP_with_filter$decideTest_result

PCAdata <- process_PCAdata(data = dt_SP, tsv_data = counts)
dim(PCAdata) # 9645
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-dt-with-filterByExpr-count-SP.pdf")

# 3.1.2.6 dt | with filterByExpr | abundance
PCAdata <- process_PCAdata(data = dt_SP, tsv_data = abundance)
dim(PCAdata) # 9645
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-dt-with-filterByExpr-abundance-SP.pdf")


# 3.1.2.7 dt | without filterByExpr | counts
dt_SP <- deg_SP_without_filter$decideTest_result

PCAdata <- process_PCAdata(data = dt_SP, tsv_data = counts)
dim(PCAdata) # 15937
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-dt-without-filterByExpr-count-SP.pdf")



# 3.1.2.8 dt | without filterByExpr | abundance
PCAdata <- process_PCAdata(data = dt_SP, tsv_data = abundance)
dim(PCAdata) # 15937
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-dt-without-filterByExpr-abundance-SP.pdf")




# 3.1.3 volcano plots 

# 3.1.3.1 with filterByExpr

View(deg_SP_with_filter$DEG_result)
names(deg_SP_with_filter$decideTest_result)

dt_SP <- deg_SP_with_filter$decideTest_result

res <- process_volcano_data(dt_SP[["Sp.RAvsTR"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-RAvsTR-with-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["Sp.RAvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-RAvsRF-with-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["Sp.TRvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-TRvsRF-with-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["Su.RAvsTR"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-RAvsTR-with-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["Su.RAvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-RAvsRF-with-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["Su.TRvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-TRvsRF-with-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["RA.SpvsSu"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-RiaAveiro-SpvsSu-with-filterByExpr-SP.pdf")

#res <- process_volcano_data(dt_SP[["TR.SpvsSu"]])
#plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE)

res <- process_volcano_data(dt_SP[["RF.SpvsSu"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-RiaFormosa-SpvsSu-with-filterByExpr-SP.pdf")


# 3.1.3.2 without filterByExpr

View(deg_SP_without_filter$DEG_result)
names(deg_SP_without_filter$decideTest_result)

dt_SP <- deg_SP_without_filter$decideTest_result

res <- process_volcano_data(dt_SP[["Sp.RAvsTR"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-RAvsTR-without-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["Sp.RAvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-RAvsRF-without-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["Sp.TRvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-TRvsRF-without-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["Su.RAvsTR"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-RAvsTR-without-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["Su.RAvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-RAvsRF-without-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["Su.TRvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-TRvsRF-without-filterByExpr-SP.pdf")

res <- process_volcano_data(dt_SP[["RA.SpvsSu"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-RiaAveiro-SpvsSu-without-filterByExpr-SP.pdf")

#res <- process_volcano_data(dt_SP[["TR.SpvsSu"]])
#plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE)

res <- process_volcano_data(dt_SP[["RF.SpvsSu"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-RiaFormosa-SpvsSu-without-filterByExpr-SP.pdf")





# TO BE DONE!!!

# 3.1.4 heatmap plot

# 3.1.4.1 with filterByExpr

# 3.1.4.2 without filterByExpr




# ==============================================================================
# TASK 3.2 - Zebrafish (with/without filterByExpr)
# ==============================================================================


# 3.2.1 with filterByExpr ------------------------------------------------------



deg_Z_with_filter <-run_and_save_DEG(Data_Z,
                                      "Z",
                                      model_type = "QLF",
                                      use_filter = TRUE,
                                      use_lfc = FALSE)


# Get top 10 + summary
top10_Z_with_filter <- get_top10_degs(deg_Z_with_filter$DEG_result, 
                                       use_filter = TRUE, "Z") 






# 3.2.2 without filterByExpr ---------------------------------------------------

deg_Z_without_filter<-run_and_save_DEG(Data_Z,
                                       "Z",
                                       model_type = "QLF",
                                       use_filter = FALSE,
                                       use_lfc = FALSE)

# Get top 10 + summary
top10_Z_without_filter <- get_top10_degs(deg_Z_without_filter$DEG_result, 
                                      use_filter = FALSE, "Z") 










# 3.2.2 continuation .... Principal Component Analysis 


# 3.2.2.1 degs | with filterByExpr | counts
deg_Z <- deg_Z_with_filter$DEG_result
  
PCAdata <- process_PCAdata(data = deg_Z, tsv_data = counts)
dim(PCAdata) # 407
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-degs-with-filterByExpr-count-Z.pdf")

# 3.2.2.2 degs | with filterByExpr | abundance
PCAdata <- process_PCAdata(data = deg_Z, tsv_data = abundance)
dim(PCAdata) # 407
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-degs-with-filterByExpr-abundance-Z.pdf")


# 3.2.2.3 degs | without filterByExpr | counts
deg_Z <- deg_Z_without_filter$DEG_result

PCAdata <- process_PCAdata(data = deg_Z, tsv_data = counts)
dim(PCAdata) # 404
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-degs-without-filterByExpr-count-Z.pdf")


# 3.2.2.4 degs | without filterByExpr | abundance
PCAdata <- process_PCAdata(data = deg_Z, tsv_data = abundance)
dim(PCAdata) # 404
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-degs-without-filterByExpr-abundance-Z.pdf")




# 3.2.2.5 dt | with filterByExpr | counts
dt_Z <- deg_Z_with_filter$decideTest_result

PCAdata <- process_PCAdata(data = dt_Z, tsv_data = counts)
dim(PCAdata) # 2473
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-dt-with-filterByExpr-count-Z.pdf")


# 3.2.2.6 dt | with filterByExpr | abundance
PCAdata <- process_PCAdata(data = dt_Z, tsv_data = abundance)
dim(PCAdata) # 2473
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-dt-with-filterByExpr-abundace-Z.pdf")


# 3.2.2.7 dt | without filterByExpr | counts
dt_Z <- deg_Z_without_filter$decideTest_result

PCAdata <- process_PCAdata(data = dt_Z, tsv_data = counts)
dim(PCAdata) # 2974  
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-dt-without-filterByExpr-count-Z.pdf")


# 3.2.2.8 dt | without filterByExpr | abundance
PCAdata <- process_PCAdata(data = dt_Z, tsv_data = abundance)
dim(PCAdata) # 2974
pca_result <- runPCA(PCAdata)
plot_pca(pca_result, save_plot = TRUE, 
         output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/pcaPlot-dt-without-filterByExpr-abundance-Z.pdf")








# 3.1.3 volcano plots 

# 3.1.3.1 with filterByExpr

View(deg_Z_with_filter$DEG_result)
names(deg_Z_with_filter$decideTest_result)

dt_Z <- deg_Z_with_filter$decideTest_result

res <- process_volcano_data(dt_Z[["Sp.RAvsTR"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-RAvsTR-with-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["Sp.RAvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-RAvsRF-with-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["Sp.TRvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-TRvsRF-with-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["Su.RAvsTR"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-RAvsTR-with-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["Su.RAvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-RAvsRF-with-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["Su.TRvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-TRvsRF-with-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["RA.SpvsSu"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-RiaAveiro-SpvsSu-with-filterByExpr-Z.pdf")

#res <- process_volcano_data(dt_Z[["TR.SpvsSu"]])
#plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE)

res <- process_volcano_data(dt_Z[["RF.SpvsSu"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-RiaFormosa-SpvsSu-with-filterByExpr-Z.pdf")


# 3.1.3.2 without filterByExpr

View(deg_Z_without_filter$DEG_result)
names(deg_Z_without_filter$decideTest_result)

dt_Z <- deg_Z_without_filter$decideTest_result

res <- process_volcano_data(dt_Z[["Sp.RAvsTR"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-RAvsTR-without-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["Sp.RAvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-RAvsRF-without-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["Sp.TRvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Spring-TRvsRF-without-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["Su.RAvsTR"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-RAvsTR-without-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["Su.RAvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-RAvsRF-without-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["Su.TRvsRF"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-Summer-TRvsRF-without-filterByExpr-Z.pdf")

res <- process_volcano_data(dt_Z[["RA.SpvsSu"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-RiaAveiro-SpvsSu-without-filterByExpr-Z.pdf")

#res <- process_volcano_data(dt_Z[["TR.SpvsSu"]])
#plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE)

res <- process_volcano_data(dt_Z[["RF.SpvsSu"]])
plot_volcanoPlot(res$volcano_data, res$highlighted, save_plot = TRUE,
                 output_file = "./03-Output/01-DEG-Analysis/annotation-results/figures/volcanoPlot-RiaFormosa-SpvsSu-without-filterByExpr-Z.pdf")




# ==============================================================================
# TASK ? - Swiss-prot vs. Zebrafish analysis
# ==============================================================================

SPvsZ <- common_deg_SP_Z(deg_SP, deg_Z)

#save_to_excel(SPvsZ$common, "./03-Output/01-DEG-Analysis/DEG-results/with-filterByExpr/degs-common.xlsx")
#save_to_excel(SPvsZ$unique_SP, "./03-Output/01-DEG-Analysis/DEG-results/with-filterByExpr/degs-unique-SP.xlsx")
#save_to_excel(SPvsZ$unique_Z, "./03-Output/01-DEG-Analysis/DEG-results/with-filterByExpr/degs-unique-Z.xlsx")



SPvsZ$common$Sp.RAvsTR
SPvsZ$common$Sp.RAvsRF
SPvsZ$common$Sp.TRvsRF
SPvsZ$common$Su.RAvsTR
SPvsZ$common$Su.RAvsRF
SPvsZ$common$Su.TRvsRF
SPvsZ$common$RF.SpvsSu
SPvsZ$common$RA.SpvsSu

SPvsZ$unique_SP$Sp.RAvsTR
SPvsZ$unique_SP$Sp.RAvsRF
SPvsZ$unique_SP$Sp.TRvsRF
SPvsZ$unique_SP$Su.RAvsTR
SPvsZ$unique_SP$Su.RAvsRF
SPvsZ$unique_SP$Su.TRvsRF
SPvsZ$unique_SP$RF.SpvsSu
SPvsZ$unique_SP$RA.SpvsSu

SPvsZ$unique_Z$Sp.RAvsTR
SPvsZ$unique_Z$Sp.RAvsRF
SPvsZ$unique_Z$Sp.TRvsRF
SPvsZ$unique_Z$Su.RAvsTR
SPvsZ$unique_Z$Su.RAvsRF
SPvsZ$unique_Z$Su.TRvsRF
SPvsZ$unique_Z$RF.SpvsSu
SPvsZ$unique_Z$RA.SpvsSu









# ==============================================================================
# TASK ? - Gene set enrichment analysis (GSEA)
# ==============================================================================

# steps:

# tips: 
# Zebra fish (without filtering low gene counts)
# Important !! we want all (over-, under-, Nosig) so we need to use decideTest output



# pro tip:
# 


# plot:


# problems:


# references:
#
# ==============================================================================


## Data (with filterByExpr):
gseaobj <- read_excel("./03-Output/01-DEG-Analysis/DEG-results/with-filterByExpr/decideTest-Z.xlsx")
View(gseaobj)

# ----- to use multiGSEA package, we need a df with "Accession", "logFC", "PValue" columns -----
GSEA_data <- process_GSEAdata(gseaobj)
save_to_excel(GSEA_data, "./03-Output/02-Pathway-Enrichment/analysis-ready-data/with-filterByExpr/GSEA-data.xlsx")


pathways <- perform_GSEA(GSEA_data)
#save_to_excel(pathways$ES_processed_result, "./03-Output/02-Pathway-Enrichment/GSEA-results/without-filterByExpr/pathways.xlsx")




# 2974 each df -> output 0 enrich pathways (wtf is going on?!)
GSEA_data <- process_GSEAdata(dt_Z)
pathways <- perform_GSEA(GSEA_data)



## Data (with filterByExpr):
uniprot_taxaobj <- read_excel("./03-Output/02-Pathway-Enrichment/preprocessed-data/with-filterByExpr/uniprot-taxaobj-dtZ.xlsx")
View(uniprot_taxaobj)


# ----- to use multiGSEA package, we need a df with "Accession", "logFC", "PValue" columns -----
GSEA_data <-  process_GSEAdata(uniprot_taxaobj)
save_to_excel(GSEA_data, "./03-Output/02-Pathway-Enrichment/analysis-ready-data/with-filterByExpr/GSEA-data.xlsx")


pathways <- perform_GSEA(GSEA_data)
save_to_excel(pathways$ES_processed_result, "./03-Output/02-Pathway-Enrichment/GSEA-results/with-filterByExpr/pathways.xlsx")











# ==============================================================================
# TASK n - STRING analysis
# ==============================================================================

# String (version 12.0):
# software link:
# https://string-db.org/cgi/input?sessionId=bAHDVX8em5Is&input_page_show_search=on

# steps:
# search -> Multiple proteins -> List Of Name -> Organisms -> Advance Settings 
#                                                              |-> Required score 
#                                                                   |-> high confidence (0.700) - professor choice

# tips:
# List Of Names -> one-per-line, use Gene.Names..primary from UniProt, if you try to use Accessions, it doesn't work! 
#               -> do not use all the Gene.Names..primary! It only work up to 2.000 entry (BIG ACHO)
#               -> try for instance:
#                      - only the degs (per-contrast)
#                      - only the under-expressed (per-contrast)
#                      - only the over-expressed (per-contrast)
#                      - only possible specific genes of interests selected (per-contrast)
#               -> "per contrast" means that you have to search individually for each contrast, not all together
#               -> you cannot use it if the organisms are all different, as in BLAST against Swiss-Prot! it fails  
# Organisms -> use auto-detect, or the actual organism, in my case "Danio rerio"

# references:


# ==============================================================================


# Data needed:




# Running UniProt may take approximately 5min (degs) to 3hours (dt). 
# Alternatively, you can skip this step and read the saved file instead.
# ------------------------------------------------------------------------------
# DataZ decide test without filterByExpr
#decideTest_Z <- deg_Z_without_filter$decideTest_result
#SpvsSu <- decideTest_Z[c("RA.SpvsSu", "TR.SpvsSu", "RF.SpvsSu" )]
#Data_SpvsSu <- get_uniprot_taxainfo(SpvsSu)
#Su <- decideTest_Z[c("Su.RAvsTR", "Su.RAvsRF", "Su.TRvsRF" )]
#Data_Su <- get_uniprot_taxainfo(Su)
#Sp <- decideTest_Z[c("Sp.RAvsTR", "Sp.RAvsRF", "Sp.TRvsRF" )]
#Data_Sp <- get_uniprot_taxainfo(Sp)
#uniprot_taxaobj <- c(Data_SpvsSu,Data_Su,Data_Sp)
#View(uniprot_taxaobj)
#save_to_excel(uniprot_taxaobj, "./03-Output/02-Pathway-Enrichment/preprocessed-data/without-filterByExpr/uniprot-taxaobj-dtZ.xlsx")

# DataZ decide test with filterByExpr
decideTest_Z <- deg_Z_with_filter$decideTest_result
SpvsSu <- decideTest_Z[c("RA.SpvsSu", "TR.SpvsSu", "RF.SpvsSu" )]
Data_SpvsSu <- get_uniprot_taxainfo(SpvsSu)
Su <- decideTest_Z[c("Su.RAvsTR", "Su.RAvsRF", "Su.TRvsRF" )]
Data_Su <- get_uniprot_taxainfo(Su)
Sp <- decideTest_Z[c("Sp.RAvsTR", "Sp.RAvsRF", "Sp.TRvsRF" )]
Data_Sp <- get_uniprot_taxainfo(Sp)
uniprot_taxaobj <- c(Data_SpvsSu,Data_Su,Data_Sp)
View(uniprot_taxaobj)
save_to_excel(uniprot_taxaobj, "./03-Output/02-Pathway-Enrichment/preprocessed-data/with-filterByExpr/uniprot-taxaobj-dtZ.xlsx")



# DataZ degs without filterByExpr
#deg_Z <- deg_Z_without_filter$DEG_result
#SpvsSu <- deg_Z[c("RA.SpvsSu", "TR.SpvsSu", "RF.SpvsSu" )]
#Data_SpvsSu <- get_uniprot_taxainfo(SpvsSu)
#Su <- deg_Z[c("Su.RAvsTR", "Su.RAvsRF", "Su.TRvsRF" )]
#Data_Su <- get_uniprot_taxainfo(Su)
#Sp <- deg_Z[c("Sp.RAvsTR", "Sp.RAvsRF", "Sp.TRvsRF" )]
#Data_Sp <- get_uniprot_taxainfo(Sp)
#uniprot_taxaobj <- c(Data_SpvsSu,Data_Su,Data_Sp)
#View(uniprot_taxaobj)
#save_to_excel(uniprot_taxaobj, "./03-Output/02-Pathway-Enrichment/preprocessed-data/without-filterByExpr/uniprot-taxaobj-degZ.xlsx")

# DataZ degs with filterByExpr
#deg_Z <- deg_Z_with_filter$DEG_result
#SpvsSu <- deg_Z[c("RA.SpvsSu", "TR.SpvsSu", "RF.SpvsSu" )]
#Data_SpvsSu <- get_uniprot_taxainfo(SpvsSu)
#Su <- deg_Z[c("Su.RAvsTR", "Su.RAvsRF", "Su.TRvsRF" )]
#Data_Su <- get_uniprot_taxainfo(Su)
#Sp <- deg_Z[c("Sp.RAvsTR", "Sp.RAvsRF", "Sp.TRvsRF" )]
#Data_Sp <- get_uniprot_taxainfo(Sp)
#uniprot_taxaobj <- c(Data_SpvsSu,Data_Su,Data_Sp)
#View(uniprot_taxaobj)
#save_to_excel(uniprot_taxaobj, "./03-Output/02-Pathway-Enrichment/preprocessed-data/with-filterByExpr/uniprot-taxaobj-degZ.xlsx")
# ------------------------------------------------------------------------------

 


# prepare data to use in string software:

## Data (with filterByExpr):
uniprot_taxaobj_with <- read_excel("./03-Output/02-Pathway-Enrichment/preprocessed-data/with-filterByExpr/uniprot-taxaobj-dtZ.xlsx")

# ----- to use STRING, we need the column "Gene.Names..primary." -----
STRING_data <- process_STRINGdata(uniprot_taxaobj_with) 
save_to_excel(STRING_data, "./03-Output/02-Pathway-Enrichment/analysis-ready-data/with-filterByExpr/STRING-data.xlsx")


## Data (with filterByExpr):
uniprot_taxaobj_without <- read_excel("./03-Output/02-Pathway-Enrichment/preprocessed-data/without-filterByExpr/uniprot-taxaobj-dtZ.xlsx")

# ----- to use STRING, we need the column "Gene.Names..primary." -----
STRING_data <- STRINGdata_processing(uniprot_taxaobj_without) 
save_to_excel(STRING_data, "./03-Output/02-Pathway-Enrichment/analysis-ready-data/without-filterByExpr/STRING-data.xlsx")












