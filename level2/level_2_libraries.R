## ***************************
##
## Script name: level_2_libraries.R 
## Purpose of script: Libraries needed for level 2   
## Author: mg  
## Date Created: 2024-08-14
##
## ***************************
##
## Notes: These are the packages used in the execution of level 2. In 00_dependencies.R are specified other packages needed for their installation as dependencies. 
##   
## ***************************


# Plotting and visualization 
library(ggplot2)
library(ggsci)

# Data manipulation 
library(plyr)
library(dplyr)
library(patchwork)
library(tidyverse)
library(magrittr)
library(Matrix)


# Spatial transcriptomics analysis 
library(Seurat)
library(glmGamPoi)
library(presto) # this package speeds up the computational time
library(clustermole) # cell type prediction
library(iTALK)
library(CellChat)
library(liana)


# Enrichment analysis 
library(msigdbr)
library(clusterProfiler)


# rmd rendering 
library(DT)

