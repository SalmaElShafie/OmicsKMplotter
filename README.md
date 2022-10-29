# OmicsKMplotter

OmicsKMplotter draws Kaplan Meier survival curves from all TCGA dataset types (mRNA, miRNA, rnaseq, mutations, RPPA, methylation) and all TCGA cancer types classifying patients into groups using any number of user-selected features/targets from the dataset.
Patients are divided into categories High/Low corresponding to high/low expression of each feature where the cutoff is user-selected to be either mean, median or 75% quartile; with values of target feature below this cutoff representing members of the low 'Cat', and those above represnting the high 'Cat'. 
For example if 2 genes were selected as targets, there would be 4 groups on the KM plot: 1.high expression of both genes 2. high expression of first but low expression of second 3. low expression of first gene but high expression of second 4. low expression of both genes.
So, for any n number of user-selected targets, n^2 categories and therefore n^2 KM lines are expected to be on the provided KM plot. In other words, it is possible to check survival analysis using a signature of high/low expression of features of interest, not individually but all together. For example, if a study reported a signature profile of high expression of 2 certain genes and low expression of two other genes to be associated with chemotherpy resistance in cell lines, the effect of this specific profile signature can be checked on TCGA cohorts and their survival using this package in a simple KM plot. the only dataset which is not categorized based on mean/median/quantile values, is mutations where the categories per feature represent the number of mutations (mutation load) for that feature, and again all user-selected features are assayed simultaneously in the target cancertype cohort therefore the number of KM lines represents the different mutational load combinations (patterns) observed in this cohort

## cutomizable options

1. The KM plot can be drawn by categorizing patients based on expression with the mean (default) being the cutoff between the 2 categories, or median, or 75% quantile, which is selected by the user in the argument cutoff="mean", "median" or "quantile"

2. Based on the number of targets/features and therefore categories on the KM plot, the legend could become lengthy and not readable on the graph, therefore using the argument legend="labels" (default) or "numbers", addresses this

## installing OmicsKMplotter

devtools::install_github('SalmaElShafie/OmicsKMplotter')

library(OmicsKMplotter)

## Output examples

OmicsKMplotter(cancertype= "OV", datasettype = "mRNA", UserList=c("ELMO2", "CREB3L1", "PNMA1"))

![mrna](https://user-images.githubusercontent.com/92435273/198752404-ef15e5b6-a1e3-4d06-aa04-7d15dcb5b0b7.jpeg)

OmicsKMplotter(cancertype= "OV", datasettype = "RPPA", UserList=c("ACC1", "AR", "ACVRL1"))

![rppa](https://user-images.githubusercontent.com/92435273/198752464-1cef46b8-8c56-4b39-96b7-e3c8efde27bc.jpeg)

OmicsKMplotter(cancertype= "OV", datasettype = "miRNASeq", UserList=c("hsa-let-7d","hsa-let-7e", "hsa-let-7a3"), legend="numbers")

![mirna](https://user-images.githubusercontent.com/92435273/198752585-859cd625-9248-4c6d-8bd6-06aebf01d1fd.png)
