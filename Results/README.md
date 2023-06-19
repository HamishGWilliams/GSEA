# Overview
This folder contains the results for each of the respective iteratinos of results generated for the GSEA and OrA analyses. Each iteration of analyses will be seperated by directory, which should make identifyinng which results are relative to which iteration and therefore databases easier.

# Iteration 1
This directory contains the results of the first iteration of analyses. This analyses comprises of the first analyses completed, including the Swiss-prot database, and a portion of the Trembl database (cnidarians). This iteration resulted in non-relevant GOs being identified from the Swiss-Prot database inclusion. 

From this result, the next iteration will remove the swiss-prot database from the blastp analyses of the differential expression analyses results. 

Included in this directory is relevant figures which were outputted. The main figure contains a side-by-side plot of the biological processes and molecular function GOs identified from the iteration results, with circle data points denoting the gene ratio of the GO on the x axis, the count of genes, and their respective adjusted p values for each GO.

# Iteration 2 (*in progress*)
This is the second iteration of the analyses, directly following iteration 1. The swiss-prot database will be excluded from this analyses, as it was identifying non-relevant GOs in the final result. ALthough Swiss-prot provides an overall better dataset in terms of accuracy and reliability of the proteins identified and their relative functions, it lacks the specificity to the target species (*Actinia equina*). It was therefore decided to complete the analyses again, but with removing this dataset from the database used. In its place, several other more closely related taxanomic groups will be added to the database of organisms used for the blastp analyses.
