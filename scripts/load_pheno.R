library(dplyr)
library(stringr)
library(Matrix)
library(fastmatch)


if (!exists('geno')) stop('Load genotype files first')
if (!exists('path.key')) stop('Specify mapping file between subject IDs in geno and pheno file as: path.key')
if (!exists('path.pheno')) stop('Specify path to phenotype file as: path.pheno')
if (!exists('data.field.pheno')) stop('Specify data field of interest in phenotype file  as: data.field.pheno')
if (!exists('code.cases.pheno')) stop('Specify encoding for cases (e.g. 2 for red hair) as: code.cases.pheno')
if (!exists('code.na.pheno')) stop('Specify encoding for NAs (e.g. c(6,-1) for hair color) as: code.na.pheno')

###############################################################################
## Load selected patient phenotypes and match to genotype data
###############################################################################

# Load patient id info from genotype files
geno.id <- rownames(geno) 

# Key for genotype/phenotype patient ids
key <- read.csv(path.key)
app.pheno <- colnames(key)[1]     #IDs for pheno file are assumed to be in first column
app.geno <- colnames(key)[2]      #IDs for geno file are assumed to be in second column

# Match phenotype ids to selected patients genotype ids
id.key <- filter(key, get(app.geno) %in% geno.id) %>%
  arrange(match(get(app.geno), geno.id))

# Load phenotypes for selected patients
pheno <- fread(path.pheno, select=c(app.pheno, data.field.pheno)) %>%
  filter(get(app.pheno) %in% id.key$app.pheno) %>%
  arrange(match(get(app.pheno), id.key$app.pheno))

pheno <- pheno$data.field.pheno

# Remove NAs from data
phen.na <- pheno %in% code.na.pheno | is.na(pheno)     
pheno <- as.numeric(pheno == code.cases.pheno)

if(length(phen.na) > 0){
  pheno <- pheno[!phen.na]
  geno <- geno[!phen.na,]
}