# ATTEMPT 1
This Data folder will contain all of the relevant data used for the first iteration of GSEA and OrA completed. This includes the combined Swiss-Prot and Tremble database from the UNIPROT databases ("https://www.uniprot.org"). 

## RESULTS
This iteration has lead to some unexpected results (*see the results section for this data*), where Gene Ontologies (GOs) such as "Face morphogenesis", very anthropomorphic traits, have been identified as the key GOs. This is clearly not suitable for the work that we are completing, as the target species is a Beadlet Anemone (*Actinia equina*). Therefore, a new iteration of the analysis will be complete. 

## Next Iteration
The SWISS-PROT database is causing some issue with the result, specifically with atributing non-relevant GOs to the Differentiall expressed genes (DEGs). TO resolve this issue, this next iteration of data and results will exclude the Swiss-Prot database. However, this does come with its own set of benefits and cons:

**CONS**:
- No longer contain highly researched and established proteins with known functions and accurate structural complexity
- Results may be less reliable as their accuracy may be reduced from less researched proteins

**PROS**:
- Will remove redundant GOs such s "Face morphogenesis", which are inappropriate for the target species of choice
- Can include a wider range of more closely related taxanomic groups, increasing likelihood of finding more relevant GOs
- Already have the Scripts generated to complete this analysis both in UNIX and on R, therefore not much time will be lost from performing this analysis.

The new datasets will be located in **Attempt 2** of the **Data** and **Results** Folders.
