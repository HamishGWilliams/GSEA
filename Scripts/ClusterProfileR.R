## Metadata ----

## Outline of script:
# This script contains the code used to perform Gene Set Enrichment Analysis
# (GSEA), as well as over representation analysis (Functional Enrichment) for 
# the Differentially expressed gene list and the entire actinia genome
# from the short-experiment dataset (No oil vs Oil)

## Date Started: 03/05/2023

## Created By: Hamish Williams 

## For use by: Hamish Williams, Marius Wenzel, David Fisher

## !! Analysis for Chapter 1 of research !!


## 1. Building, installing and loading required packages: ----

# Build function...
check_packages_ClusterProfileR <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE)) {
      print(paste(pkg, "is not installed, installing package..."))
      install.packages(pkg)
    } else {
      suppressMessages(library(pkg, character.only = TRUE))
      print(paste(pkg, "is installed and loaded"))
    }
  }
}

## Install ClusterProfiler: !! Unhash as needed !!
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("clusterProfiler")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("AnnotationForge")
#
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#
#BiocManager::install("enrichplot", force = TRUE)

# Create list of packages needed:
pkg_list_ClusterProfileR<-c("tidyverse", 
                            "tibble", 
                            "stats",
                            "ggplot2", 
                            "dplyr", 
                            "tidyr", 
                            "ggfortify", 
                            "knitr",
                            "clusterProfiler",
                            "AnnotationForge",
                            "enrichplot",
                            "data.table")

# Run Function::
check_packages_ClusterProfileR(pkg_list_ClusterProfileR)
# Double Check everything has loaded:
check_packages_ClusterProfileR(pkg_list_ClusterProfileR)

# Open the clusterProfiler manual:
browseVignettes("clusterProfiler")
browseVignettes("AnnotationForge")

## 2. Importing Data ----
# Import GO term file:
{
DEG_GO_Terms <- read.table("Data/Attempt 1/DEG_GO_terms.txt", 
                           header = T, 
                           sep = '\t', 
                           stringsAsFactors = FALSE)

# Set the path to the zip file
zip_file <- "Data/Attempt 1/genome_GO_terms.zip"

# Specify the destination directory where you want to extract the contents
destination_dir <- "Data/Attempt 1/"

# Unzip the file
unzip(zipfile = zip_file, exdir = destination_dir)


genome_GO_Terms <- read.table("Data/Attempt 1/genome_GO_terms.txt", 
                           header = T, 
                           sep = '\t', 
                           stringsAsFactors = FALSE)

DEG_BP_MF_GO_Terms <- read.table("Data/Attempt 1/DEG_BP_MF_GOs.txt", 
                           header = T, 
                           sep = '\t', 
                           stringsAsFactors = FALSE)

genome_BP_MF_GO_Terms <- read.table("Data/Attempt 1/genome_BP_MF_merged.txt",
                              header = T, 
                              sep = '\t', 
                              stringsAsFactors = FALSE)
}

# Import blastp file:
{
DEG_matching_file <- read.table("Data/Attempt 1/DEG_matching_file.txt", 
                                header = F, 
                                sep = '\t', 
                                stringsAsFactors = FALSE)
  
# Set the path to the zip file
zip_file <- "Data/Attempt 1/genome_matching_file.zip"

# Specify the destination directory where you want to extract the contents
destination_dir <- "Data/Attempt 1/"

# Unzip the file
unzip(zipfile = zip_file, exdir = destination_dir)

genome_matching_file <- read.table("Data/Attempt 1/genome_matching_file.txt", 
                                header = F, 
                                sep = '\t', 
                                stringsAsFactors = FALSE)
}


# Change the name of the columns to be identical
{# DEGS
colnames(DEG_GO_Terms)
names(DEG_GO_Terms)[names(DEG_GO_Terms) == "From"] <- "ID"
colnames(DEG_GO_Terms)

colnames(DEG_matching_file)
names(DEG_matching_file)[names(DEG_matching_file) == "V13"] <- "ID"
colnames(DEG_matching_file)

colnames(DEG_BP_MF_GO_Terms)
names(DEG_BP_MF_GO_Terms)[names(DEG_BP_MF_GO_Terms) == "From"] <- "ID"
names(DEG_BP_MF_GO_Terms)[names(DEG_BP_MF_GO_Terms) == "Gene.Names"] <- "GENENAMES"
names(DEG_BP_MF_GO_Terms)[names(DEG_BP_MF_GO_Terms) == "Gene.Ontology..biological.process."] <- "GO:BP"
names(DEG_BP_MF_GO_Terms)[names(DEG_BP_MF_GO_Terms) == "Gene.Ontology..molecular.function."] <- "GO:MF"
colnames(DEG_BP_MF_GO_Terms)

# Genome
colnames(genome_matching_file)
names(genome_matching_file)[names(genome_matching_file) == "V13"] <- "ID"
colnames(genome_matching_file)

colnames(genome_BP_MF_GO_Terms)
names(genome_BP_MF_GO_Terms)[names(genome_BP_MF_GO_Terms) == "From"] <- "ID"
names(genome_BP_MF_GO_Terms)[names(genome_BP_MF_GO_Terms) == "Gene.Names"] <- "GENENAMES"
names(genome_BP_MF_GO_Terms)[names(genome_BP_MF_GO_Terms) == "Gene.Ontology..biological.process."] <- "GO:BP"
names(genome_BP_MF_GO_Terms)[names(genome_BP_MF_GO_Terms) == "Gene.Ontology..molecular.function."] <- "GO:MF"
colnames(genome_BP_MF_GO_Terms)
}

## 3. Merge DEG files ----

## Merge DEG Files:
merged_df <- merge(DEG_matching_file, DEG_BP_MF_GO_Terms, by = "ID")
colnames(merged_df)
merged_df <- subset(merged_df[, c("V1","ID","GO:BP","GO:MF")])
names(merged_df)[names(merged_df) == "V1"] <- "GENE"

## DEG BP only:
  # Select BP column:
    DEG_BP <- merged_df[,c("GENE","ID","GO:BP")]
    DEG_BP <- DEG_BP[!apply(DEG_BP == "", 1, any), ]
    DEG_BP <- DEG_BP[!duplicated(DEG_BP), ]
  # Extract GO terms:
    GENES_DEG_BP <- DEG_BP
    GENES_DEG_BP$GO <- ""
    bracketed_rows <- grepl("\\[.*\\]", GENES_DEG_BP$`GO:BP`)
    GENES_DEG_BP$GO[bracketed_rows] <- gsub(".*\\[(.*)\\].*", "\\1", GENES_DEG_BP$`GO:BP`[bracketed_rows])
    GENES_DEG_BP <- GENES_DEG_BP[, c("GO","GENE")]
    GENES_DEG_BP <- GENES_DEG_BP[!duplicated(GENES_DEG_BP), ]
    
    ## Extract Terms:
    TERMs_DEG_BP <- DEG_BP
    TERMs_DEG_BP$GO <- ""
    bracketed_rows <- grepl("\\[.*\\]", TERMs_DEG_BP$`GO:BP`)
    TERMs_DEG_BP$GO[bracketed_rows] <- gsub(".*\\[(.*)\\].*", "\\1", TERMs_DEG_BP$`GO:BP`[bracketed_rows])
    TERMs_DEG_BP$`GO:BP`[bracketed_rows] <- gsub("\\[.*\\]", "", TERMs_DEG_BP$`GO:BP`[bracketed_rows])
    TERMs_DEG_BP <- TERMs_DEG_BP[, c("GO","GO:BP")]
    TERMs_DEG_BP <- TERMs_DEG_BP[!duplicated(TERMs_DEG_BP), ]
    names(TERMs_DEG_BP)[names(TERMs_DEG_BP) == "GO:BP"] <- "GENENAMES_BP"

## DEG MF only:
    # Select BP column:
    DEG_MF <- merged_df[,c("GENE","ID","GO:MF")]
    DEG_MF <- DEG_MF[!apply(DEG_MF == "", 1, any), ]
    DEG_MF <- DEG_MF[!duplicated(DEG_MF), ]
    # Extract GO terms:
    GENES_DEG_MF <- DEG_MF
    GENES_DEG_MF$GO <- ""
    bracketed_rows <- grepl("\\[.*\\]", GENES_DEG_MF$`GO:MF`)
    GENES_DEG_MF$GO[bracketed_rows] <- gsub(".*\\[(.*)\\].*", "\\1", GENES_DEG_MF$`GO:MF`[bracketed_rows])
    GENES_DEG_MF <- GENES_DEG_MF[, c("GO","GENE")]
    GENES_DEG_MF <- GENES_DEG_MF[!duplicated(GENES_DEG_MF), ]
    
    ## Extract Terms:
    TERMs_DEG_MF <- DEG_MF
    TERMs_DEG_MF$GO <- ""
    bracketed_rows <- grepl("\\[.*\\]", TERMs_DEG_MF$`GO:MF`)
    TERMs_DEG_MF$GO[bracketed_rows] <- gsub(".*\\[(.*)\\].*", "\\1", TERMs_DEG_MF$`GO:MF`[bracketed_rows])
    TERMs_DEG_MF$`GO:MF`[bracketed_rows] <- gsub("\\[.*\\]", "", TERMs_DEG_MF$`GO:MF`[bracketed_rows])
    TERMs_DEG_MF <- TERMs_DEG_MF[, c("GO","GO:MF")]
    TERMs_DEG_MF <- TERMs_DEG_MF[!duplicated(TERMs_DEG_MF), ]
    names(TERMs_DEG_MF)[names(TERMs_DEG_MF) == "GO:MF"] <- "GENENAMES_MF"

## Extract DEGs list:
z <- merged_df[,1]
z <- unique(z)
DEGs <- z
        

## 4. Merge GENOME files ----

# Genome
merged_df <- merge(genome_matching_file, genome_BP_MF_GO_Terms, by = "ID")
colnames(merged_df)
names(merged_df)[names(merged_df) == "V1"] <- "GENE"
merged_df <- subset(merged_df[, c("GENE","ID","GO:BP","GO:MF")])
merged_df <- merged_df[!apply(merged_df == "", 1, any), ]
merged_df <- merged_df[!duplicated(merged_df), ]

## Genome BP terms only
  # Remove unwanted term (MF)
    aGENOME_BP <- merged_df[,c("GENE","GO:BP")]
  # Separate the Ontologies from each other onto separate rows
    aGENOME_BP <- separate_rows(aGENOME_BP, "GO:BP", sep="; ")
  # Drop any results which contain N/As:
    aGENOME_BP <- aGENOME_BP[!apply(aGENOME_BP == "N/A", 1, any), ]
  # Remove duplicate rows
    aGENOME_BP <- aGENOME_BP[!duplicated(aGENOME_BP), ]
    
  ## Create GENEs_GENOME file:
    GENEs_GENOME_BP <- aGENOME_BP
    GENEs_GENOME_BP$GO <- ""
    bracketed_rows <- grepl("\\[.*\\]", GENEs_GENOME_BP$`GO:BP`)
    GENEs_GENOME_BP$GO[bracketed_rows] <- gsub(".*\\[(.*)\\].*", "\\1", GENEs_GENOME_BP$`GO:BP`[bracketed_rows])
    GENEs_GENOME_BP <- GENEs_GENOME_BP[, c("GO","GENE")]
    
   ## Create TERMs_Genome file:
    TERMs_GENOME_BP <- aGENOME_BP
    TERMs_GENOME_BP$GO <- ""
    bracketed_rows <- grepl("\\[.*\\]", TERMs_GENOME_BP$`GO:BP`)
    TERMs_GENOME_BP$GO[bracketed_rows] <- gsub(".*\\[(.*)\\].*", "\\1", TERMs_GENOME_BP$`GO:BP`[bracketed_rows])
    TERMs_GENOME_BP$`GO:BP`[bracketed_rows] <- gsub("\\[.*\\]", "", TERMs_GENOME_BP$`GO:BP`[bracketed_rows])
    TERMs_GENOME_BP <- TERMs_GENOME_BP[, c("GO","GO:BP")]
    TERMs_GENOME_BP <- TERMs_GENOME_BP[!duplicated(TERMs_GENOME_BP), ]
    names(TERMs_GENOME_BP)[names(TERMs_GENOME_BP) == "GO:BP"] <- "GENENAMES_BP"

## Genome MF terms only:
    # Remove unwanted term (MF)
    aGENOME <- merged_df[,c("GENE","GO:MF")]
    # Seperate the Ontologies from eachother onto seperate rows
    aGENOME_MF <- separate_rows(aGENOME, "GO:MF", sep="; ")
    # Drop any results which contain N/As:
    aGENOME_MF <- aGENOME_MF[!apply(aGENOME_MF == "N/A", 1, any), ]
    # Remove duplicate rows
    aGENOME_MF <- aGENOME_MF[!duplicated(aGENOME_MF), ]
    
    ## Create GENEs_GENOME file:
    GENEs_GENOME_MF <- aGENOME_MF
    GENEs_GENOME_MF$GO <- ""
    bracketed_rows <- grepl("\\[.*\\]", GENEs_GENOME_MF$`GO:MF`)
    GENEs_GENOME_MF$GO[bracketed_rows] <- gsub(".*\\[(.*)\\].*", "\\1", GENEs_GENOME_MF$`GO:MF`[bracketed_rows])
    GENEs_GENOME_MF <- GENEs_GENOME_MF[, c("GO","GENE")]
    
    ## Create TERMs_Genome file:
    TERMs_GENOME_MF <- aGENOME_MF
    TERMs_GENOME_MF$GO <- ""
    bracketed_rows <- grepl("\\[.*\\]", TERMs_GENOME_MF$`GO:MF`)
    TERMs_GENOME_MF$GO[bracketed_rows] <- gsub(".*\\[(.*)\\].*", "\\1", TERMs_GENOME_MF$`GO:MF`[bracketed_rows])
    TERMs_GENOME_MF$`GO:MF`[bracketed_rows] <- gsub("\\[.*\\]", "", TERMs_GENOME_MF$`GO:MF`[bracketed_rows])
    TERMs_GENOME_MF <- TERMs_GENOME_MF[, c("GO","GO:MF")]
    TERMs_GENOME_MF <- TERMs_GENOME_MF[!duplicated(TERMs_GENOME_MF), ]
    names(TERMs_GENOME_MF)[names(TERMs_GENOME_MF) == "GO:MF"] <- "GENENAMES_MF"

## 4. Import DE results ----
# Import DE results for the short experiment:
DE_results <- read.table("Data/Attempt 1/DE_results_names.txt", 
                         header = T, 
                         sep = '\t', 
                         stringsAsFactors = TRUE)

# feature 1: numeric vector
geneList = DE_results[,2]

# feature 2: named vector
names(geneList) = as.character(DE_results[,1])

# feature 3: decreasing order;
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

## 5. Run GSEA ----

# Create list of packages needed:
pkg_list_ClusterProfileR<-c("stats",
                            "clusterProfiler")

# Run Function::
check_packages_ClusterProfileR(pkg_list_ClusterProfileR)


# Run GSEA of DEGs and Genome:
yBP <- enricher(DEGs, 
              TERM2GENE = GENEs_GENOME_BP,
              TERM2NAME = TERMs_GENOME_BP,
              minGSSize=1)
yBP@result
DEG_BP_results <- yBP@result
DEG_BP_results

yMF <- enricher(DEGs, 
                TERM2GENE = GENEs_GENOME_MF,
                TERM2NAME = TERMs_GENOME_MF,
                minGSSize=1,
                pvalueCutoff = 0.1)
yMF@result
DEG_MF_results <- yMF@result

xBP<- GSEA(geneList,
          TERM2GENE = GENEs_GENOME_BP,
          TERM2NAME = TERMs_GENOME_BP,
          pvalueCutoff = 0.2)
xBP@result
GENOME_BP_results <- xBP@result
GENOME_BP_results

xMF <- GSEA(geneList,
          TERM2GENE = GENEs_GENOME_MF,
          TERM2NAME = TERMs_GENOME_MF,
          pvalueCutoff = 0.2)
xMF@result
GENOME_MF_results <- xMF@result



?`simplify,enrichResult-method`
?dropGO
?GSEA
?enricher

enrichplot::dotplot(yMF) + ggtitle("Molecular Function") +
  enrichplot::dotplot(yBP) + ggtitle("Biological Processes")

enrichplot::dotplot(xMF) + ggtitle("Molecular Function") +
  enrichplot::dotplot(xBP) + ggtitle("Biological Processes")


  # works by dropping GO terms based on specificity on hierarchy

  # generate term2names file to better visualize results
      # gene names + GO term names
        # 1st: term
        # 2nd: name of gene

    # After everything, collect the BP and MF annoatations, and run the analyses
    # again just with these groups of terms to reduce the strictness of the
    # p-value adjustments

# For kegg, you would need to blast against nematostella vectensis and actinia tenebrosa, 
# Take best hit
# make vector for DEGS of KEGG equivalent IDs
# keytype = uniprot
 # organism = model organism KEGG id name (NVE, ATE)

?enrichKEGG


Vector <- as.vector(DEG_GSEA_dataset[,3])
head(Vector)

write.table(Vector, "clipboard", row.names = FALSE)


# Run another BLAST for DEGs with only SWISS-PROT data to get more precise gene names


## Troubleshooting ----
## Matching the geneList vector to the GENOME/TERM dataframes ----
missing_genes <- setdiff(rownames(geneList), GENEs_GENOME_BP$GENE)
if (length(missing_genes) > 0) {
  print(paste("Missing genes in GENEs_GENOME_BP dataframe:", paste(missing_genes, collapse = ", ")))
}

##  It seems that the genes in the geneList are covered in the GENEs_GENOME_BP dataframe,
##  Which essentially means that in this direction, there arent any GENEs mising


# Check if gene names in GENEs_GENOME_BP dataframe match geneList
missing_genes <- setdiff(GENEs_GENOME_BP$GENE, rownames(geneList))
if (length(missing_genes) > 0) {
  print(paste("Missing genes in geneList:", paste(missing_genes, collapse = ", ")))
}

##  However, when we run it in the opposite direction, we see that there is a 
##  a large number of genes which are not covered in the geneList vector
##  I suspect that this is due to the multiple GOs per Gene formatting of the
##  dataframe. I see only 2 solutions which will solve this problem:

    # 1. Reduce the GOs per gene to only the top matching GO per gene
        # Pros: Will solve this issue with the data, and potentially for the 
        #       entire function

        # Cons: Will lose some of the other information which relates to the 
        #       potential function of the genes

    # 2. Condense the GOs per gene to a single row.
        # Pros: would allow the inclusion of other potentially useful GOs
        #       per gene

        # Cons: Will be visually very difficult to understand the results
        #       , producing hard to read figures. And potentially may not work
        #       with the analysis

## Solution 1: 

## Genome BP terms only
# Remove unwanted term (MF)
aGENOME_BP <- merged_df[,c("GENE","GO:BP")]
# Separate the Ontologies from each other onto separate rows
aGENOME_BP <- separate_rows(aGENOME_BP, "GO:BP", sep="; ")
# Drop any results which contain N/As:
aGENOME_BP <- aGENOME_BP[!apply(aGENOME_BP == "N/A", 1, any), ]
# Remove duplicate rows
aGENOME_BP <- aGENOME_BP[!duplicated(aGENOME_BP), ]

# Subset dataframe to include only unique instances of gene names
aGENOME_BP <- aGENOME_BP[!duplicated(aGENOME_BP$GENE), ]

## Create GENEs_GENOME file:
GENEs_GENOME_BP <- aGENOME_BP
GENEs_GENOME_BP$GO <- ""
bracketed_rows <- grepl("\\[.*\\]", GENEs_GENOME_BP$`GO:BP`)
GENEs_GENOME_BP$GO[bracketed_rows] <- gsub(".*\\[(.*)\\].*", "\\1", GENEs_GENOME_BP$`GO:BP`[bracketed_rows])
GENEs_GENOME_BP <- GENEs_GENOME_BP[, c("GO","GENE")]

## Create TERMs_Genome file:
TERMs_GENOME_BP <- aGENOME_BP
TERMs_GENOME_BP$GO <- ""
bracketed_rows <- grepl("\\[.*\\]", TERMs_GENOME_BP$`GO:BP`)
TERMs_GENOME_BP$GO[bracketed_rows] <- gsub(".*\\[(.*)\\].*", "\\1", TERMs_GENOME_BP$`GO:BP`[bracketed_rows])
TERMs_GENOME_BP$`GO:BP`[bracketed_rows] <- gsub("\\[.*\\]", "", TERMs_GENOME_BP$`GO:BP`[bracketed_rows])
TERMs_GENOME_BP <- TERMs_GENOME_BP[, c("GO","GO:BP")]
TERMs_GENOME_BP <- TERMs_GENOME_BP[!duplicated(TERMs_GENOME_BP), ]
names(TERMs_GENOME_BP)[names(TERMs_GENOME_BP) == "GO:BP"] <- "GENENAMES_BP"

# Check if gene names in geneList match GENEs_GENOME_BP dataframe
missing_genes <- setdiff(rownames(geneList), GENEs_GENOME_BP$GENE)
if (length(missing_genes) > 0) {
  print(paste("Missing genes in GENEs_GENOME_BP dataframe:", paste(missing_genes, collapse = ", ")))
}

# Check if gene names in GENEs_GENOME_BP dataframe match geneList
missing_genes <- setdiff(GENEs_GENOME_BP$GENE, rownames(geneList))
if (length(missing_genes) > 0) {
  print(paste("Missing genes in geneList:", paste(missing_genes, collapse = ", ")))
}

    ## It's clear at this point that the likely issue is due to the unequal 
    ## number of values in the two dataframes (geneList and the dataframes)
    ## Therefore, I'm going to try to subset the geneList to only genes
    ## which match the genes found from the blastp results for the genome 
    ## dataframes

# Get the matching gene names between geneList and GENEs_GENOME_BP
rownames_genome_df <- as.data.frame(GENEs_GENOME_BP)
rownames(rownames_genome_df) <- (GENEs_GENOME_BP$GENE)
matching_genes <- intersect(row.names(geneList), row.names(rownames_genome_df))

# Convert geneList to data.frame
geneList_df <- data.frame(geneList)
geneList_df$GENE <- rownames(geneList_df)


# Merge geneList_df with GENEs_GENOME_BP by gene_name
merged_data <- merge(geneList_df, GENEs_GENOME_BP, by = "GENE", all.x = TRUE)

# Filter out non-intersecting genes
intersect_genes <- merged_data$GENE[!is.na(merged_data$GO)]

# Create a new vector with intersecting genes from geneList
intersect_vector <- geneList[intersect_genes]
head(intersect_vector)
intersect_vector <- sort(intersect_vector)

# Check if gene names in GENEs_GENOME_BP dataframe match geneList
missing_genes <- setdiff(GENEs_GENOME_BP$GENE, rownames(intersect_vector))
if (length(missing_genes) > 0) {
  print(paste("Missing genes in intersect_vector:", paste(missing_genes, collapse = ", ")))
}

## Check if Gene Ontologies are missing between TERMs and GENEs:
{
# Check if gene ontology IDs in GENEs_GENOME_BP match TERMs_GENOME_BP dataframe
missing_terms <- setdiff(GENEs_GENOME_BP$GO, TERMs_GENOME_BP$GO)
if (length(missing_terms) > 0) {
  print(paste("Missing gene ontology IDs in TERMs_GENOME_BP dataframe:", paste(missing_terms, collapse = ", ")))
}

# Check if gene ontology IDs in TERMs_GENOME_BP match GENEs_GENOME_BP dataframe
missing_terms <- setdiff(TERMs_GENOME_BP$GO, GENEs_GENOME_BP$GO)
if (length(missing_terms) > 0) {
  print(paste("Missing gene ontology IDs in GENEs_GENOME_BP dataframe:", paste(missing_terms, collapse = ", ")))
}
}
# Checked, and GOs match!


xBP<- GSEA(geneList,
           TERM2GENE = GENEs_GENOME_BP,
           TERM2NAME = TERMs_GENOME_BP,
           pvalueCutoff = 0.1)
xBP@result
GENOME_BP_results <- xBP@result
GENOME_BP_results
