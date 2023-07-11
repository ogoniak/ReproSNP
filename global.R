setwd("./")

#----------------------------------------------------------------------------------------------------
# PACKAGES + FUNCTIONS
#----------------------------------------------------------------------------------------------------
library(shiny)
library(DT)
library(stringr)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(dplyr)
library(formattable)
library(shinyFeedback)
library(liftOver)
library(shinyBS)
library(fresh)
library(shinyvalidate)
library(shinyalert)
library(biomaRt)
library(data.table)

source("FUNCTIONS/bgen2vcf.R")
source("FUNCTIONS/fisher.R")
source("FUNCTIONS/functions.R")
source("FUNCTIONS/liftover.R")
source("FUNCTIONS/mergeTables.R")
source("FUNCTIONS/annotation.R")

#----------------------------------------------------------------------------------------------------
# GLOBAL DATA
#----------------------------------------------------------------------------------------------------

### DIRECTORIES AND FILES
tmp_dir = "/path/to/TMP"
refseq_curated = "/path/to/ncbiRefSeqCurated_filtered.txt"
gene_symbols = "/path/to/gene_symbols.txt"
phenotable = read.table("/path/to/phenotable.tsv", header=T, sep="\t")

### REFSEQ
header = c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
RefSeq = read.table(refseq_curated, col.names=header, sep="\t")
refseq_name2 = read.table(gene_symbols, header=F)[,1]

### OTHER VARIABLES
dbsnp_url = "https://www.ncbi.nlm.nih.gov/snp/?term=(___CHROM___>%5BChromosome%5D)%20AND%20___BP___%5BBase%20Position%5D"
clinvar_url = "https://www.ncbi.nlm.nih.gov/clinvar?term=(___CHROM___%5BChromosome%5D)%20AND%20___BP___%5BBase%20Position%5D"

### COLOR THEME
mytheme <- create_theme(
  adminlte_color(
    light_blue = "#618C68",
    olive = "#8EB074",
    green = "#8EB074",
    lime = "#8A9A97"
  ),
  adminlte_sidebar(
    dark_bg = "#586F6B",
    dark_hover_bg = "#756B80",
    dark_color = "#F5F0F6",
  ),
  adminlte_global(
    content_bg = "#F9F7FA",
    box_bg = "#E4DFE9", 
    info_box_bg = "#F9F7FA"
  )
)
