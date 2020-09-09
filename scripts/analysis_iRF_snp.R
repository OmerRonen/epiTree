library(iRF)
library(ranger)
library(data.table)
library(dplyr)
library(stringr)
library(Matrix)
library(fastmatch)
library(glmnet)
library(PRROC)
library(RSQLite)
library("grid")
library("ggplotify")
library(ape)
library(seriation)
library(gridExtra)
library(rgl)
library(rpart)
library(pryr)
library(purrr)
library(snpStats)

source('utilities_general.R')

###############################################################################
## Specify paths       
###############################################################################
path.plink <- paste0("path_to_output_plink_file/name_of_output_plink_file_chr",1:22)      #plink files with snp data

path.pheno <- "path_to_phenotype_file/app15860_standard_data_2016Nov19.txt"     #path to phenotype file
path.key <- "path_to_mapping_file/ukb15860_13721_mapping.tsv"     #mapping file between subject IDs in pheno (1st column) and geno (2nd column) file
data.field.pheno <- 'n_1747_0_0'      #data field of interest in phenotype file (here for red hair)
code.cases.pheno <- 2    # red hair phenotype is encoded as 2
code.na.pheno <- c(6, -1)     #NAs for red hair phenotype are 6, -1

path.out <-  '../results/'     #path for output files
name.out <- 'analysis_iRF_snp.Rdata'      #name of output files
path.lasso <- paste0("../results/lasso_",name.out)      #output file for lasso results only
path.ranger <- paste0(home,"/irfGWES/results/ranger_",name.out)     #output file for ranger results only

###############################################################################
## Load training  data          
###############################################################################
# we assume to have one plink file per chromosome
# we load a balanced sample of 13K cases and controls

chr <- 1
geno <- data.frame()

while(length(geno) == 0){
  if(file.exists(paste0(path.plink[chr], '.bim'))){   #we may not have selected SNPs from every chromosome
    print(paste('load initial chr', chr))
    snp <- read.table(paste0(path.plink[chr], ".bim"))[,2]
    is.dup <- duplicated(snp)
    geno <- read.plink(path.plink[chr], select.snps = which(!is.dup))$genotypes
    
    source('load_pheno.R')
    
    load.id <-  c(which(pheno == 0)[1:13000], which(pheno == 1)[1:13000])
    geno <- geno[load.id,]
    pheno <- pheno[load.id]
    
    #rename SNPs, replace ':' by 'chr' for chromosome, '_' by 'b' for base pair and add chr number to name
    colnames(geno) <- map_chr(colnames(geno), ~{renameSNP(.,chr)}) 
    
    gc()
  }
  chr <- chr + 1
}

while(chr <= 22){
  if(file.exists(paste0(path.plink[chr], '.bim'))){
    print(paste('load chr', chr))
    snp <- read.table(paste0(path.plink[chr], ".bim"))[,2]
    is.dup <- duplicated(snp)
    
    geno_new <- read.plink(path.plink[chr], select.snps = which(!is.dup))$genotypes
    geno_new <- geno_new[!phen.na,]
    geno_new <- geno_new[load.id,]
    colnames(geno_new) <- map_chr(colnames(geno_new), ~{renameSNP(.,chr)})
    gc()
    
    geno <- cbind(geno, geno_new)
    rm(geno_new)
    gc()
    
  }
  chr <- chr + 1
}



geno <- as(geno, 'numeric')     #convert SnpMatrix object to numeric (SNPs are no encoded as 0,1,2 and NA for missing)
mode(geno) <- 'integer'
geno <- as.data.frame(geno)
geno <- geno %>% mutate_all(~ifelse(is.na(.), median(., na.rm = TRUE), .))     #impute missing SNPs by median
geno <- as.matrix(geno)

gnames <- colnames(geno)


###############################################################################
## Run lasso                           
###############################################################################

if(file.exists(path.lasso)){
  load(path.lasso)
}else{
  print("Run lasso")
  lambda <- cv.glmnet(x = geno, y = pheno, family = "binomial", maxit = 10^3)$lambda.min
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

KTop = c(10, 50, 100, 500, 1000, 2500)       #consider top k RF features
oobError = numeric(length(KTop))        #out-of-bag error when considering top k RF features

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

kTop <- KTop[which.min(oobError)]       #select best k based on oob error

idkp <- order(frang$variable.importance, decreasing=TRUE)[1:kTop]
geno <- geno[,idkp]
gnames <- gnames[idkp]
print(dim(geno))



#extract candidate interactions via iRF
rit.param <- list(ntree=5000, depth=3, nchild=5, class.id=1, min.nd=1)
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


rdForest <- readForest(fit$rf.list, geno)       #read forest for posthoc analysis of iRF output


save(file = paste0(path.out, name.out), fit, frang, lasso, rdForest)

