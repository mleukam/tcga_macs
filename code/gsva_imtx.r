## -----------------------------------------------------

## R script to get GSVA scores on all TCGA RNAseq reads from solid tumors

## Expects cleaned, normalized, log-transformed, expression matrix ready for GSVA analysis
##    Originally defined using HTseq reads from GDC
##    Expression input to GSVA may be expressionset object or matrix
## Expects gene lists as a list of character vectors for each gene set
##    Gene list input must be a list
## Outputs matrix of GSVA scores for each sample (columns) and gene set (rows)

## M. Leukam
## 11/2019

## -----------------------------------------------------

# clear workspace
rm(list = ls())

# load packages
library("tidyverse")
library("GSVA")
library("parallel")

# read in data
total_reads <- read_csv("/gpfs/data/kline-lab/tcga_macs/data/imtx_matrix.csv")
gset <- readRDS("/gpfs/data/kline-lab/tcga_macs/data/gset_ids_plusnew.rds")

# set up expression matrix
rownames <- total_reads %>% pull(gene)
expr_matrix <- total_reads %>%
  dplyr::select(-gene) %>%
  as.matrix()
rownames(expr_matrix) <- rownames

# get gvsa scores
imtx_es <- gsva(expr = expr_matrix, gset.idx.list = gset, annotation = NULL, method = "gsva", verbose = TRUE)

# write out results
saveRDS(imtx_es, "/gpfs/data/kline-lab/tcga_macs/output/gsva_imtx_plusuro_results.rds")
