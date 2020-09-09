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


###############################################################################
## Specify paths       
###############################################################################

path.predx <- "../data/name_of_output_prediXcan_txt_file_predicted_expression.txt"      #path to imputed gene expression data
path.db <- "../data/gtex_v7_Skin_Sun_Exposed_Lower_leg_imputed_europeans_tw_0.5_signif.db"      #path to data base
path.excl <- "../data/geneExclusionList.txt"      #optional list of genes that should be included

path.pheno <- "path_to_phenotype_file/app15860_standard_data_2016Nov19.txt"     #path to phenotype file
path.key <- "path_to_mapping_file/ukb15860_13721_mapping.tsv"     #mapping file between subject IDs in pheno (1st column) and geno (2nd column) file
data.field.pheno <- 'n_1747_0_0'      #data field of interest in phenotype file (here for red hair)
code.cases.pheno <- 2    # red hair phenotype is encoded as 2
code.na.pheno <- c(6, -1)     #NAs for red hair phenotype are 6, -1


path.out <- '../results/'     #path for output files
name.out <- 'analysis_iRF_gene.Rdata'      #name of output files
path.lasso <- paste0("../results/lasso_",name.out)      #output file for lasso results only
path.ranger <- paste0(home,"/irfGWES/results/ranger_",name.out)     #output file for ranger results only

###############################################################################
## Load training  data          
###############################################################################

#here we load training data in batches to avoid memory issues
#we load a balanced sample of 13K cases and controls
#we assume that the first 100K subjects contain at least 13K controls

print("load training data")

#load batch 1
load.id <- 1:100000     #load first 100000 subjects

source(paste0(home, '/irfGWES/scripts/load_predixcan.R'))     #load genotype files
source(paste0(home, '/irfGWES/scripts/load_pheno.R'))      #load phenotype files

numb_cases <- sum(pheno == 1)
load.id <-  c(which(pheno == 0)[1:13000], which(pheno == 1)[1:min(13000, numb_cases)])
geno <- geno[load.id,]
pheno <- pheno[load.id]

#load remaining batches
batch.id <- 1
while(numb_cases < 13000){
  geno_old <- geno
  pheno_old <- pheno
  numb_cases_old <- numb_cases
  load.id <- (100000 * batch_id + 1):(100000 * batch_id + 100000)
  source(paste0(home, '/irfGWES/scripts/load_predixcan.R'))
  source(paste0(home, '/irfGWES/scripts/load_pheno.R'))
  
  numb_cases <- sum(pheno == 1)
  load.id <-  which(pheno == 1)[1:min(13000 - numb_cases_old, numb_cases)]
  geno <- rbind(geno_old, geno[load.id,])
  pheno <- c(pheno_old, pheno[load.id])
  numb_cases <- numb_cases + numb_cases_old
  batch.id <- batch.id + 1
}


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

for(i in 1:length(KTop)){
  kTop <- KTop[i]
  idkp <- order(frang$variable.importance, decreasing=TRUE)[1:kTop]     #select top k RF features

  # Run iRF (without interaction selection)
  fit <- iRF(x = geno[,idkp], 
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
geno <- geno[,idkp]
gnames <- gnames[idkp]

#extract candidate interactions via iRF
rit.param <- list(ntree=5000, depth=3, nchild=5, class.id=1, min.nd=1)      #set RIT parameters
fit <- iRF(x=geno, 
           y=as.factor(pheno), 
           varnames.grp=gnames,
           n.iter=3,
           iter.return=3,
           select.iter=FALSE,
           n.bootstrap=50,
           rit.param=rit.param,
           int.return = 3,
           type='ranger',
           n.core=1)

rdForest <- readForest(fit$rf.list, geno)     #read forest for posthoc analysis of iRF output

save(file = paste0(path.out, name.out), fit, frang, lasso, rdForest)      #save results

