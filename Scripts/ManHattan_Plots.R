### Metadata ----
# Script contains code and data for constructing Manhattan Plots for the 
# Differential Expression results. 

## Example Manhattan Plot Code: ----

# Install Package:
install.packages("qqman")
library(qqman)

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P" )
  # This is the basics of how the manhattan plot command works:
    # gwasResults = the dataframe, containing...

x <- gwasResults
head(x)

      # 4 columns:
        # SNP = single nucleotide polymorphism
        # CHR = Chromosome
        # BP = Base Pair location
        # P = P value

?manhattan

    # The command then calls for several features of the dataframe, in this case...
      # the chromosome location...
      # the base pair number location...
      # the SNP name...
      # the P value of the SNP...

# In terms of this being relevant to my work, the chromosome location can likely
# be dropped, as this is just used to denote additional information that may
# be interesting to have. 

# Base pair data forms the primary X-Axis input data. So this is a must 

# SNP can be replaced with "gene_name", as this variable is used to distinguish
# between different genes.

# The P value forms the primary Y-axis data in the example. This
# could even be replaced with LFC values to give an interesting figure 
# showing the changes in relative expression level of a particular gene across
# the entire anemone genome.

## Genes of interest: ----

# We can even specify a list of genes which form the significantly differentially
# expressed genes:

snpsOfInterest # Simply need a list of names which you which to highlight

# Highlight the listed names in the plot:
manhattan(gwasResults, highlight = snpsOfInterest)

# We can even choose to highlight points with specific Pvalues:
manhattan(gwasResults, annotatePval = 0.01)


## QQ plot: ----
# Just a simple QQ plot, but can be useful to test the results against an
# expected "chance" outcome. Which can give us an idea of how likely the 
# results are due to chance versus to response to our stimulus.

## Additional Options: ---
# The package also gives some options for conducting the work in ggplot and
# even generating an intractable plot of results. These additional result
# formats can be found at: https://r-graph-gallery.com/101_Manhattan_plot.html

### Constructing a Manhattan Plot for DE results: ----

# In order to get a data frame prepped and ready to generate a nice Manhattan
# plot, we need the following column data;

    # GENE: the gene names of all the genes aligned
    # BP: the base pair locations of each gene
    # LFC: the log fold change value for each gene
    # P: The P value for the DE of genes from the analyses

# The GENE, LFC and P variables will be easy to extract, as there is already
# data frames in other Projects which contain this data. The tricky part will
# be extracting the right BP locations for each of the genes.

# Some of the start/end locations of each gene are not consistent between 
# different samples. How to select which start and end BP locations
# should be used for each gene is not clear. Will need to decide on a method
# to tackle this issue.
    # I will likely need to simply extract the very first value from each 
    # data point, as trying to average any of the values will simply lead
    # to the possibility of overlap of gene locations, which may cause issues
    # later down the line in terms of analyses.

# The other issue with this variable is whether to take the start or end
# number as the BP variable data for the Manhattan data frame object.
# I think that either way, the result will roughly be the same, as the start 
# of the next gene is also the end of the last gene. In this case, the starting
# number for the gene BP location is likely to be the best value of the two
# to take for this variable column.

## Separating the lowest value start locations for each gene ----

# create filename
filename <- "Data/Attempt 2/A_Equina_Counts.txt"

# read in data
data <- read.table(filename, sep="\t", header = TRUE)

# Extract only the smallest value of the "Start" column
data$Start <- sapply(strsplit(as.character(data$Start), ";"), function(x) {
  x <- as.numeric(x)
  x <- unique(x)  # Remove duplicates
  min_val <- min(x)
  x <- x[x <= min_val]
  paste(x, collapse = ";")  # Combine remaining values with a semicolon delimiter
})

# Subset data to only include gene name and start columns
data_subset <- subset(data[, c(1,3)])

# Remove original Data file
rm(data)

## Currently, the Names of the gene IDs do not match the other datafile,
## so we need to make them the same:

# Remove the ".TU" component of names:
# Check if the 'Geneid' column exists in the data frame
  # if ("Geneid" %in% colnames(data_subset)) {
  #   # Replacing ".TU" with "." in the "Geneid" column
  #   data_subset$Geneid <- gsub("\\.TU", ".", data_subset$Geneid)
  #   
  #   # Replacing "path" with "mrna" in the "Geneid" column
  #   data_subset$Geneid <- gsub("path", "mrna", data_subset$Geneid)
  # }

# Rename the "data_subset" file to an appropriate name:
bp_locations <- data_subset
rm(data_subset)

## Import DESeq2 results: ----
filename <- "Data/Attempt 2/DEresultsshort.csv"
data <- read.table(filename, sep=",", header = TRUE)

# rename first row
colnames(data)
names(data)[names(data) == "X"] <- "Geneid"
names(data)[names(data) == "log2FoldChange"] <- "LFC"
colnames(data)

# subset only relevant columns:
data_subset <- subset(data[, c(1,3,6,7)])
colnames(data_subset)
rm(data)
# rename data frame:
DESeq_results <- data_subset
rm(data_subset)

# Merging datafiles:
head(DESeq_results)
head(bp_locations)

merged_files <- merge(DESeq_results, bp_locations, by = "Geneid")
colnames(merged_files)
names(merged_files)[names(merged_files) == "Start"] <- "bp"
colnames(merged_files)

# Change properties of variables:
str(merged_files)
merged_files$bp <- as.numeric(merged_files$bp)


    ## Repeat for long experiment

filename <- "Data/Attempt 2/DEresultslong.csv"
data <- read.table(filename, sep=",", header = TRUE)

# rename first row
colnames(data)
names(data)[names(data) == "X"] <- "Geneid"
names(data)[names(data) == "log2FoldChange"] <- "LFC"
colnames(data)

# subset only relevant columns:
data_subset <- subset(data[, c(1,3,6,7)])
colnames(data_subset)
rm(data)
# rename data frame:
DESeq_results_long <- data_subset
rm(data_subset)

# Merging datafiles:
head(DESeq_results_long)
head(bp_locations)

merged_files_long <- merge(DESeq_results_long, bp_locations, by = "Geneid")
colnames(merged_files_long)
names(merged_files_long)[names(merged_files_long) == "Start"] <- "bp"
colnames(merged_files_long)

# Change properties of variables:
str(merged_files_long)
merged_files_long$bp <- as.numeric(merged_files_long$bp)

## Manhattan Plot----
library(tidyverse)
library(ggplot2)
library(plotly)

# First plot:
  # LFC~BP, coloured with Padj - SHORT
ggplot(merged_files, aes(x = bp, y = LFC, color = padj)) +
  geom_point(aes(color = ifelse(padj < 0.1, "Significant", "Non-Significant")), size = 4) +
  scale_color_manual(values = c("Non-Significant" = "darkgray", "Significant" = "orange")) +
  labs(x = "Base Pair Position", y = "Log-fold Change", color = "Padj") +
  scale_x_continuous(limits = c(0,3000000), breaks = seq(0,3000000, by = 500000)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6, by = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed",
             linewidth = 1.5) +
  geom_hline(yintercept = -1, colour = "blue", linetype = "dashed",
             linewidth = 1.5)

dev.copy(png, file = "Results/Iteration 2/LFC_BP_Short_1st.png")
dev.off()
dev.copy(svg, file = "Results/Iteration 2/LFC_BP_Short_1st.svg")
dev.off()

# LFC~BP, coloured with Padj - LONG
ggplot(merged_files_long, aes(x = bp, y = LFC, color = padj)) +
  geom_point(aes(color = ifelse(padj < 0.1, "Significant", "Non-Significant")), size = 4) +
  scale_color_manual(values = c("Non-Significant" = "darkgray", "Significant" = "purple")) +
  labs(x = "Base Pair Position", y = "Log-fold Change", color = "Padj") +
  scale_x_continuous(limits = c(0,3000000), breaks = seq(0,3000000, by = 500000)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6, by = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed",
             linewidth = 1.5) +
  geom_hline(yintercept = -1, colour = "blue", linetype = "dashed",
             linewidth = 1.5)

dev.copy(png, file = "Results/Iteration 2/LFC_BP_Long_1st.png")
dev.off()
dev.copy(svg, file = "Results/Iteration 2/LFC_BP_Long_1st.svg")
dev.off()


    ## Too much noise, removing NA padj values
NA_rm_merged_files <- merged_files[complete.cases(merged_files),]
NA_rm_merged_files_long <- merged_files_long[complete.cases(merged_files_long),]

# 2nd plot
  # removed NA values - SHORT

ggplot(NA_rm_merged_files, aes(x = bp, y = LFC, color = padj)) +
  geom_point(aes(color = ifelse(padj < 0.1, "Significant", "Non-Significant")), size = 4) +
  scale_color_manual(values = c("Non-Significant" = "darkgrey", "Significant" = "orange")) +
  labs(x = "Base Pair Position", y = "Log-fold Change", color = "Padj") +
  scale_x_continuous(limits = c(0,3000000), breaks = seq(0,3000000, by = 500000)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6, by = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed",
             linewidth = 1.5) +
  geom_hline(yintercept = -1, colour = "blue", linetype = "dashed",
             linewidth = 1.5)
dev.copy(png, file = "Results/Iteration 2/LFC_BP_Short_2nd.png")
dev.off()
dev.copy(svg, file = "Results/Iteration 2/LFC_BP_Short_2nd.svg")
dev.off()

# removed NA values - LONG

ggplot(NA_rm_merged_files_long, aes(x = bp, y = LFC, color = padj)) +
  geom_point(aes(color = ifelse(padj < 0.1, "Significant", "Non-Significant")), size = 4) +
  scale_color_manual(values = c("Non-Significant" = "darkgrey", "Significant" = "purple")) +
  labs(x = "Base Pair Position", y = "Log-fold Change", color = "Padj") +
  scale_x_continuous(limits = c(0,3000000), breaks = seq(0,3000000, by = 500000)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6, by = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed",
             linewidth = 1.5) +
  geom_hline(yintercept = -1, colour = "blue", linetype = "dashed",
             linewidth = 1.5)

dev.copy(png, file = "Results/Iteration 2/LFC_BP_Long_2nd.png")
dev.off()
dev.copy(svg, file = "Results/Iteration 2/LFC_BP_Long_2nd.svg")
dev.off()

    ## A lot of non-significant values, removing them...

sig_only_merged_files <- subset(merged_files[merged_files$padj < 0.1,])
sig_only_merged_files <- sig_only_merged_files[complete.cases(sig_only_merged_files),]
count(unique(merged_files[merged_files$padj < 0.1,]))

sig_only_merged_files_long <- subset(merged_files_long[merged_files_long$padj < 0.1,])
sig_only_merged_files_long <- sig_only_merged_files_long[complete.cases(sig_only_merged_files_long),]
count(unique(merged_files_long[merged_files_long$padj < 0.1,]))

# 3rd Plot
  # Removing padj < 0.1... - SHORT

ggplot(sig_only_merged_files, aes(x = bp, y = LFC, color = padj)) +
  geom_point(aes(color = ifelse(padj < 0.1, "Significant", "Non-Significant")), size = 4) +
  scale_color_manual(values = c("Non-Significant" = "darkgrey", "Significant" = "orange")) +
  labs(x = "Base Pair Position", y = "Log-fold Change", color = "Padj") +
  scale_x_continuous(limits = c(0,3000000), breaks = seq(0,3000000, by = 500000)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6, by = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed",
             linewidth = 1.5) +
  geom_hline(yintercept = -1, colour = "blue", linetype = "dashed",
             linewidth = 1.5)

dev.copy(png, file = "Results/Iteration 2/LFC_BP_Short_3rd.png")
dev.off()
dev.copy(svg, file = "Results/Iteration 2/LFC_BP_Short_3rd.svg")
dev.off()

# Removing padj < 0.1... - LONG

ggplot(sig_only_merged_files_long, aes(x = bp, y = LFC, color = padj)) +
  geom_point(aes(color = ifelse(padj < 0.1, "Significant", "Non-Significant")), size = 4) +
  scale_color_manual(values = c("Non-Significant" = "darkgrey", "Significant" = "purple")) +
  labs(x = "Base Pair Position", y = "Log-fold Change", color = "Padj") +
  scale_x_continuous(limits = c(0,3000000), breaks = seq(0,3000000, by = 500000)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6, by = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed",
             linewidth = 1.5) +
  geom_hline(yintercept = -1, colour = "blue", linetype = "dashed",
             linewidth = 1.5)

dev.copy(png, file = "Results/Iteration 2/LFC_BP_Long_3rd.png")
dev.off()
dev.copy(svg, file = "Results/Iteration 2/LFC_BP_SLong_3rd.svg")
dev.off()

## Combined short + Long

ggplot() +
  geom_point(data = sig_only_merged_files, 
             aes(x = bp, y = LFC, color = "Short"),
             size = 4) +
  geom_point(data = sig_only_merged_files_long, 
             aes(x = bp, y = LFC, color = "Long"),
             size = 4) +
  scale_color_manual(values = c("Long" = "purple", "Short" = "orange")) +
  labs(x = "Base Pair Position", y = "Log-fold Change", color = "Experiment") +
  scale_x_continuous(limits = c(0,3000000), breaks = seq(0,3000000, by = 500000)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6, by = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed",
             linewidth = 1.5) +
  geom_hline(yintercept = -1, colour = "blue", linetype = "dashed",
             linewidth = 1.5)

dev.copy(png, file = "Results/Iteration 2/Combined_plot.png")
dev.off()
dev.copy(svg, file = "Results/Iteration 2/Combined_plot.svg")
dev.off()

      ## Interactable version

Figure_3 <- ggplot(data = sig_only_merged_files, aes(text = paste("Gene:", Geneid))) +
                geom_point(data = sig_only_merged_files, 
                           aes(x = bp, y = LFC, color = "Short"),
                           size = 4) +
                geom_point(data = sig_only_merged_files_long, 
                           aes(x = bp, y = LFC, color = "Long"),
                           size = 4) +
                scale_color_manual(values = c("Long" = "purple", "Short" = "orange")) +
                labs(x = "Base Pair Position", y = "Log-fold Change", color = "Padj") +
                scale_x_continuous(limits = c(0,3000000), breaks = seq(0,3000000, by = 500000)) +
                scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6, by = 1)) +
                theme_minimal() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed",
             linewidth = 1.5) +
  geom_hline(yintercept = -1, colour = "blue", linetype = "dashed",
             linewidth = 1.5)

Figure_3 <- ggplotly(Figure_3)
print(Figure_3)


## Adding on GOs to each data point of plot ----

# To add on the descriptions of the gene ontologies for each gene, we
# need a combined dataframe which includes a description for the GOs in a 
# separate column.

# Do we already have a data frame in our arsenal to do this? First we need to
# understand where our data has come from, which blastp analyses was performed
# on it, that way we can determine which GSEA results we should use to 
# add the descriptions of the gene ontologies to the plot.

# for the "merged_files" dataframe, this was constructed by combining the 
# DESeq_results file and the bp_locations file. The BP_locations file has not
# had any analyses performed on it, it is simply an extract of the gene positions
# The DESeq_results is the full list of results of the Differential Expression
# Analyses for the SHORT experiment performed. 

# Neither of these datasets contain Gene Ontology Data, so therefore cannot 
# already be prescribed the Descriptions of the GOs. Therefore, the first logical
# step is to determine what blastp results we should use to attribute the 
# GOs for each gene. 

# The blastp performed was that containing both SWISS-Prot and cnidarian TREMBL
# data (as of current, I am running a new analyses with only cnidarian TREMBL
# data). The top 10 hits for each blastp hit were taken, exported to UNIPROT
# to find their GO's, then imported back to UNIX to combine the files, which 
# gave a file which contained combined GO and gene data.

  ## Try to combine the GOs from the Genome Data...
merged_genome_GOs <- merge(GENEs_GENOME_BP, TERMs_GENOME_BP, by = "GO")

  ## With this now merged dataframe, try to merge GENEs from this 
  ## dataframe to those in the merged_files dataframe.
colnames(merged_genome_GOs)
names(merged_genome_GOs)[names(merged_genome_GOs) == "GENE"] <- "Geneid"
head(merged_genome_GOs)

colnames(merged_files)
head(merged_files)
data <-merged_files

if ("Geneid" %in% colnames(data)) {
  # Replacing ".TU" with "." in the "Geneid" column
  data$Geneid <- gsub("\\.TU.", ".", data$Geneid)
  
  # Replacing "path" with "mrna" in the "Geneid" column
  data$Geneid <- gsub("path", "mrna", data$Geneid)
  
  # Printing the modified data frame
  head(data)
} else {
  cat("Error: 'Geneid' column not found in the data frame.\n")
}
merged_files<-data

merge_attempt <- merge(merged_genome_GOs, merged_files, by = "Geneid") 

## I think the issue before was that the text replace code was giving the 
## merged_files dataframe two ".."'s instead of just one: "."

head(merge_attempt)

## Attempting to merge the GO and Description columns together
combined_df <- aggregate(cbind(GO, GENENAMES_BP) ~ Geneid + LFC + pvalue + padj + bp, data = merge_attempt, FUN = function(x) paste(unique(x), collapse = "; "))
str(combined_df)
combined_df$bp <- as.numeric(combined_df$bp)
sig_only_combined_df <- subset(combined_df[combined_df$padj < 0.1,])
sig_only_combined_df <- sig_only_combined_df[complete.cases(sig_only_combined_df),]
count(unique(combined_df[combined_df$padj < 0.1,]))

## Making a ggplot with this new dataframe...
ggplot(combined_df, aes(x = bp, y = LFC)) +
  geom_point(aes(color = ifelse(padj < 0.1, "Significant", "Non-Significant")), size = 2.5) +
  scale_color_manual(values = c("Non-Significant" = "purple", "Significant" = "orange")) +
  labs(x = "Base Pair Position", y = "Log-fold Change", color = "Padj") +
  scale_x_continuous(limits = c(0,3000000), breaks = seq(0,3000000, by = 500000)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6, by = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed") +
  geom_hline(yintercept = -1, colour = "blue", linetype = "dashed")

test <- ggplot(sig_only_combined_df, aes(x = bp, y = LFC, color = padj, text = paste(Geneid))) +
  geom_point(aes(color = ifelse(padj < 0.1, "Significant", "Non-Significant")), size = 2.5) +
  scale_color_manual(values = c("Non-Significant" = "purple", "Significant" = "orange")) +
  labs(x = "Base Pair Position", y = "Log-fold Change", color = "Padj") +
  scale_x_continuous(limits = c(0,3000000), breaks = seq(0,3000000, by = 500000)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6, by = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 1, colour = "red", linetype = "dashed") +
  geom_hline(yintercept = -1, colour = "blue", linetype = "dashed")

test <- ggplotly(test)
print(test)


ggplot(merged_files, aes(x = bp, y = -log10(pvalue), color = LFC)) +
  geom_point(aes(color = ifelse(LFC > 1 | LFC < -1, "Upregulated", "Down-Regulated")), size = 2) +
  scale_color_manual(values = c("Down-Regulated" = "blue", "Upregulated" = "red")) +
  labs(x = "Base Pair Position", y = "-log10(pvalue)", color = "LFC") +
  theme_minimal()

ggplot(sig_only_merged_files, aes(x = bp, y = -log10(pvalue), color = LFC)) +
  geom_point(aes(color = ifelse(LFC > 1 | LFC < -1, "Upregulated", "Down-Regulated")), size = 2) +
  scale_color_manual(values = c("Down-Regulated" = "blue", "Upregulated" = "red")) +
  labs(x = "Base Pair Position", y = "-log10(pvalue)", color = "LFC") +
  theme_minimal()

ggplot(merged_files, aes(x = bp, y = -log10(pvalue))) +
  geom_point(size = 2) +
  labs(x = "Base Pair Position", y = "-log10(pvalue)", color = "LFC") +
  theme_minimal()

ggplot(merged_files, aes(x = bp, y = -log10(padj))) +
  geom_point(size = 2) +
  labs(x = "Base Pair Position", y = "-log10(pvalue)", color = "LFC") +
  theme_minimal()


GOI_Short <- read.csv("Data/Attempt 2/GOI_Short.csv", header = FALSE)
GOI_Short <- GOI_Short[-1,]

highlited_names <- GOI_Short
names <- merged_files$Geneid
df <- data.frame(names = names)
df$highlight <- ifelse(df$names %in% highlited_names, "DEGs", "Normal")

ggplot(merged_files, aes(x = bp, y = -log10(padj))) +
  geom_point(aes(color = ifelse(Geneid %in% highlited_names, "DEGs", "Normal")), size = 2) +
  labs(x = "Base Pair Position", y = "-log10(pvalue)", color = "LFC") +
  theme_minimal()

## Merging GENEs and TERMS files for figure data point annotation ----
GENES_DEG_BP <- read.csv("Data/Attempt 2/GENES_DEG_BP.csv", sep = ",", header = TRUE)
TERMS_DEG_BP <- read.csv("Data/Attempt 2/TERMs_DEG_BP.csv", sep = ",", header = TRUE)
MERGED_DEG_BP <- merge(GENES_DEG_BP, TERMS_DEG_BP, by = "GO")
