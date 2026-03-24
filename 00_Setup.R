### Author Qamar Feddah
### Setup

library(Seurat)
library(tidyverse)
library(DESeq2)
library(viridis)
library(EnhancedVolcano)
library(ggrepel)

set.seed(333)
setwd("~/ConM_project/Analysis/07_confounder_assessment+DESeq2_ConM/")

dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/plots", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE)
