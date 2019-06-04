# miRNA_Polarization
Code and additional data for miRNA andpolarization studies from Kupffer cells and PBMCs performed in the Nagy Lab

This repository contains all additional code and files used in the following manuscript:

miRNAs involved in M1/M2 hyperpolarization are clustered and coordinately expressed in Alcoholic Hepatitis

List of scripts with descriptions:

PBMC_miRNAseq_DE.R - R script to calculate differentially expressed miRNAs. Uses EdgeR and plots a volcano plot with ggplot2. (Figure 1)

PBMC_qPCR_Scatterplot.R - R script to draw scatterplots for pair-wise comparisons of miRNA expression from qPCR data. This will draw both the healthy control and alcoholic hepatitis data as separately colored dots with correlation lines for each group. Script has code for both comparisons in paper, but can also be modified for any comparison. (Figure 6A and Figure 6C)

PBMC_qPCR_Heatmap.R - R script to draw a heatmap of all pair-wise miRNA expression comparisons from qPCR data of healthy controls and alcoholic hepatitis patients. (Figure 6B)


List of files with descriptions:

PBMC_miRNA_qPCR.txt - All qPCR data calculated by the ddCT method (average CT of miRNA - average CT of SNORD68 control). HC - Healthy Control, AH - Alcoholic Hepatitis, phenotype labeled under columns "pheno"

PBMC_miRNA_corr.txt - Table of Pearson Correlation Coefficients for all qPCR data, organized such that all Healthy Control comparisons are above the X=Y line and all Alcoholic Hepatitis comparisons are below it.

All_Counts.txt - Raw count data for all annotated miRNAs in the rat genome.
