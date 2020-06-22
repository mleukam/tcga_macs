###############################
## Shell script for tximport ##
###############################

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t' 

## Load compiler
module load gcc/6.2.0

## Load module
module load R/3.4.1

## Navigate to directory containing R script
cd /gpfs/data/kline-lab/tcga_macs/code

## Call R script
Rscript gsva_imtx.r && exit 0

## Exit if error
exit 1
