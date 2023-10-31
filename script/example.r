## Settings
working_dir <- "C:/Users/shihb/OneDrive - Lancaster University/work/project/2023/20231026_skinsig"
setwd(working_dir)

# input data
skinsig <- read.csv("data/skinsig.csv")
group_annotation <- read.csv("data/group.csv")
group_annotation <- group_annotation$group
in_count_mx <- read.csv("data/count_mx.csv", row.names=1)

## install required libraries
if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager") } 
if (!require("limma", quietly = TRUE)) { BiocManager::install("limma") } 
if (!require("edgeR", quietly = TRUE)) { BiocManager::install("edgeR") } 
if (!require("ggplot2", quietly = TRUE)) { install.packages("ggplot2") } 


## load required library
library(edgeR)
library(limma)
library(ggplot2)
source("script/function_skinsig_plot.r")


## Required inputs:
# - count_mx = count matrix (gene symbol in row.names and duplicates are not allowed; missing values not allowed in the count matrix, unique sample id in the column)
# - group = vector indicating test or control for the samples (must be in the same order as the count matrix); can only do 2 group comparison
# - ctrl_name = the value in group that indicates control samples (default = "ctrl")
## Optional inputs:
# - datatype = the data type (default = "RNAseq")
# - proportion_same_direction = the proportion of genes regulated in the same direction in gene set enrichment analysis (default = 0.75)
# - proportion_same_direction = FDR threshold for gene set enrichment analysis (default = 0.01)


test <- skinsig_plots(in_count_mx, group_annotation, proportion_same_direction = 0.8, fdr_threshold = 0.01)

test$enrichment
test$correlation
test$heatmap


