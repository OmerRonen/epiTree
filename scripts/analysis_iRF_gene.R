library(iRF)
library(ranger)
library(data.table)
library(dplyr)
library(stringr)
library(Matrix)
library(fastmatch)
library(glmnet)
library(RSQLite)
library("grid")
library(ape)
library(seriation)
library(gridExtra)
library(rgl)
library(pryr)
library(argparse)

# Define the command-line arguments
parser <- ArgumentParser(description = "EpiTree Analysis CLI")
parser$add_argument("--predx", type = "character", help = "Path to imputed gene expression data")
parser$add_argument("--db", type = "character", help = "Path to the database")
parser$add_argument("--pheno", type = "character", help = "Path to phenotype file")
parser$add_argument("--out", type = "character", help = "Path for output files")

# Parse the command-line arguments
args <- parser$parse_args()

# Assign arguments to variables
path.predx_0 <- args$predx
path.db <- args$db
path.pheno <- args$pheno
path.out <- args$out
name.out <- 'analysis_iRF_gene.Rdata'      #name of output files


# Construct path for lasso and ranger output files
path.lasso <- file.path(path.out, "lasso.Rdata")
path.ranger <- file.path(path.out, "ranger.Rdata")
#
# ###############################################################################
# ## Specify paths
# ###############################################################################
#
# path.predx <- "/accounts/campus/omer_ronen/projects/epiTree/results/expression/chr_6_predicted_expression.txt"      #path to imputed gene expression data
# path.db <- "/accounts/campus/omer_ronen/projects/epiTree/data/ctimp_Brain_Cortex.db"      #path to data base
path.excl <- "/accounts/campus/omer_ronen/projects/epiTree/data/geneExclusionList.txt"      #optional list of genes that should be included
#
# path.pheno <- "/accounts/campus/omer_ronen/projects/epiTree/data/Multiple_sclerosis/pheno.csv"     #path to phenotype file
# #path.key <- "/accounts/campus/omer_ronen/projects/epiTree/data/Multiple_sclerosis/mapping.csv"     #mapping file between subject IDs in pheno (1st column) and geno (2nd column) file
# # data.field.pheno <- 'n_1747_0_0'      #data field of interest in phenotype file (here for red hair)
# # code.cases.pheno <- 2    # red hair phenotype is encoded as 2
# # code.na.pheno <- c(6, -1)     #NAs for red hair phenotype are 6, -1
#
#
# path.out <- '../results/'     #path for output files
# path.lasso <- paste0("results/lasso_",name.out)      #output file for lasso results only
# path.ranger <- paste0("results/ranger_",name.out)     #output file for ranger results only

###############################################################################
## Load training  data          
###############################################################################

#here we load training data in batches to avoid memory issues
#we load a balanced sample of 13K cases and controls
#we assume that the first 100K subjects contain at least 13K controls

print("load training data")

#load batch 1
load.id <- 1:100000     #load first 100000 subjects
geno_list = list()
j = 1
for (i in 1:22){
    path.predx <- paste0(path.predx_0, "/chr_",i,"_predicted_expression.txt")
    source(paste0('scripts/load_predixcan.R'))     #load genotype files
    # remove all zero columns from geno
    # add geno to list
    geno_list[[j]] <- geno[, colSums(geno) != 0]
    j = j + 1

}
geno <- do.call(cbind, geno_list)
# remove duplicated columns values


#source(paste0('scripts/load_pheno.R'))      #load phenotype files
pheno <- read.csv(path.pheno, header=FALSE)
# set first column as row names
rownames(pheno) <- pheno[,1]
# match pheno and geno by rownames
pheno <- pheno[match(rownames(geno), rownames(pheno)),]
colnames(pheno) <- c("id", "pheno")
pheno <- pheno$pheno

numb_cases <- sum(pheno == 1)
# load.id <-  c(which(pheno == 0)[1:13000], which(pheno == 1)[1:min(13000, numb_cases)])
# geno <- geno[load.id,]
# pheno <- pheno[load.id]

#load remaining batches
# batch.id <- 1
# while(numb_cases < 13000){
#   geno_old <- geno
#   pheno_old <- pheno
#   numb_cases_old <- numb_cases
#   load.id <- (100000 * batch_id + 1):(100000 * batch_id + 100000)
#   source(paste0('../scripts/load_predixcan.R'))
#   source(paste0('../scripts/load_pheno.R'))
#
#   numb_cases <- sum(pheno == 1)
#   load.id <-  which(pheno == 1)[1:min(13000 - numb_cases_old, numb_cases)]
#   geno <- rbind(geno_old, geno[load.id,])
#   pheno <- c(pheno_old, pheno[load.id])
#   numb_cases <- numb_cases + numb_cases_old
#   batch.id <- batch.id + 1
# }


###############################################################################
## Run lasso                           
###############################################################################

if(file.exists(path.lasso)){
  load(path.lasso)
}else{
  print("Run lasso")
  lambda <- cv.glmnet(x = geno, y = pheno, family = "binomial")$lambda.min
  lasso <- glmnet(geno, pheno, family = "binomial", lambda = lambda)
  save(file = path.lasso, lasso)
}

###############################################################################
## Run ranger  
###############################################################################

set.seed(47)
if(file.exists(path.ranger)){
  load(path.ranger)
}else{
  print("Run ranger")
  frang <- ranger(data=cbind(geno, pheno), 
                  dependent.variable.name='pheno',
                  classification=TRUE, 
                  importance='impurity',
                  num.threads=1)
  
  save(file = path.ranger, frang)
}

###############################################################################
## Run iRF  on top k RF features   
###############################################################################

KTop = c(10, 50, 100, 500, 1000, 2500)      #consider top k RF features
oobError = numeric(length(KTop))      #out-of-bag error when considering top k RF features
gnames <- colnames(geno)

geno_mat <- data.matrix(geno)
for(i in 1:length(KTop)){
  kTop <- KTop[i]
  idkp <- order(frang$variable.importance, decreasing=TRUE)[1:kTop]     #select top k RF features

  # Run iRF (without interaction selection)
  fit <- iRF(x = geno_mat[,idkp],
             y = as.factor(pheno), 
             varnames.grp = gnames[idkp],
             n.iter = 3,
             iter.return = 3,
             select.iter = FALSE,
             type = 'ranger',
             n.core = 1)
  oobError[i] <- fit$rf.list$prediction.error       #out-of-bag error for top k RF features
}

kTop <- KTop[which.min(oobError)]     #select best k based on oob error

idkp <- order(frang$variable.importance, decreasing=TRUE)[1:kTop]
geno_mat <- geno_mat[,idkp]


#extract candidate interactions via iRF
rit.param <- list(ntree=5000, depth=3, nchild=5, class.id=1, min.nd=1)      #set RIT parameters
fit <- iRF(x=geno_mat,
           y=as.factor(pheno), 
           varnames.grp=gnames[idkp],
           n.iter=3,
           iter.return=3,
           select.iter=FALSE,
           n.bootstrap=50,
           rit.param=rit.param,
           int.return = 3,
           type='ranger',
           n.core=1)

rdForest <- readForest(fit$rf.list, geno_mat)     #read forest for posthoc analysis of iRF output

save(file = paste0(path.out, "/", name.out), fit, frang, lasso, rdForest)      #save results

