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

# The other issue with this variable is whether to take the start or end
# number as the BP variable data for the Manhattan data frame object.
# I think that either way, the result will roughly be the same, as the start 
# of the next gene is also the end of the last gene. In this case, the starting
# number for the gene BP location is likely to be the best value of the two
# to take for this variable column.


