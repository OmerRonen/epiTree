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
library(superheat)
library("grid")
library("ggplotify")
library(ape)
library(seriation)
library(gridExtra)
library(rgl)
library(rpart)
library(pryr)

reloadData <- TRUE      #should full data set be loaded

source('utilities_general.R')
source('utilities_tests.R')

###############################################################################
## Specify paths  
###############################################################################

res <- 'analysis_iRF_gene'      #name of results
path.res <- '../results/'

path.predx <- "../data/name_of_output_prediXcan_txt_file_predicted_expression.txt"      #path to imputed gene expression data
path.db <- "../data/gtex_v7_Skin_Sun_Exposed_Lower_leg_imputed_europeans_tw_0.5_signif.db"      #path to data base
path.excl <- "../data/geneExclusionList.txt"      #optional list of genes that should be included

path.pheno <- "path_to_phenotype_file/app15860_standard_data_2016Nov19.txt"     #path to phenotype file
path.key <- "path_to_mapping_file/ukb15860_13721_mapping.tsv"     #mapping file between subject IDs in pheno (1st column) and geno (2nd column) file
data.field.pheno <- 'n_1747_0_0'      #data field of interest in phenotype file (here for red hair)
code.cases.pheno <- 2    # red hair phenotype is encoded as 2
code.na.pheno <- c(6, -1)     #NAs for red hair phenotype are 6, -1
path.pca <- "path_to_pca_file/name_of_pca_file"       #path to file with principle components of subjects 

path.out <- '../results/'     #path for output files
name.out <- paste0(res,'.Rdata')      #name of output files
path.lasso <- paste0("../results/lasso_",name.out)      #output file for lasso results only
path.ranger <- paste0("../results/ranger_",name.out)     #output file for ranger results only

###############################################################################
## Load results 
###############################################################################

load(paste0(path.res,res,".Rdata"))

#reload full data set
if(reloadData){
  
  ###############################################################################
  ## Load training and test  data          
  ###############################################################################
  
  # here we load data in batches to avoid memory issues
  # we load a balanced sample of 15K cases and controls
  # the first 13K of each class are training samples
  # we assume that the first 100K subjects contain at least 15K controls
  
  print("load training and test data")
  
  #load batch 1
  load.id <- 1:100000     #load first 100000 subjects
  
  source(paste0('../scripts/load_predixcan.R'))     #load genotype files
  source(paste0('../scripts/load_pheno.R'))      #load phenotype files
  
  numb_cases <- sum(pheno == 1)
  load.id <-  c(which(pheno == 0)[1:15000], which(pheno == 1)[1:min(15000, numb_cases)])
  geno <- geno[load.id,]
  pheno <- pheno[load.id]
  
  #load remaining batches
  batch.id <- 1
  while(numb_cases < 15000){
    geno_old <- geno
    pheno_old <- pheno
    numb_cases_old <- numb_cases
    load.id <- (100000 * batch_id + 1):(100000 * batch_id + 100000)
    source(paste0('../scripts/load_predixcan.R'))
    source(paste0('../scripts/load_pheno.R'))
    
    numb_cases <- sum(pheno == 1)
    load.id <-  which(pheno == 1)[1:min(15000 - numb_cases_old, numb_cases)]
    geno <- rbind(geno_old, geno[load.id,])
    pheno <- c(pheno_old, pheno[load.id])
    numb_cases <- numb_cases + numb_cases_old
    batch.id <- batch.id + 1
  }
  
  ind.train <- c(which(pheno == 0)[1:13000], which(pheno == 1)[1:13000])      #specify training indeces
  geno.train <- geno[ind.train,]
  pheno.train <- pheno[ind.train]
  
  geno <- geno[-ind.train,]
  pheno <- pheno[-ind.train]
  
  gc()
  
  # Load PCAs for test data
  sel.sub <- rownames(geno)
  source('load_pca.R')
  pca <- xpc[, 1:15]      #select first 15 principle components
  
  pca <- as.matrix(pca)
  geno <- as.matrix(geno)
  geno.train <- as.matrix(geno.train)
  
  
  ###############################################################################
  ## Test lasso  
  ###############################################################################

  ypred <- predict(lasso, geno, type = "response")
  ypred.lasso <- ypred
  
  pr.curve.lasso <- pr.curve(ypred[pheno == 1], ypred[pheno == 0], curve = TRUE)
  roc.curve.lasso <- roc.curve(ypred[pheno == 1], ypred[pheno == 0], curve = TRUE)
  
  ###############################################################################
  ## Test ranger  
  ###############################################################################
  
  ypred <- predict(frang, data=geno, predict.all=TRUE)
  ypred <- rowMeans(ypred$predictions)
  ypred.ranger <- ypred
  
  pr.curve.ranger <- pr.curve(ypred[pheno == 1], ypred[pheno == 0], curve = TRUE)
  roc.curve.ranger <- roc.curve(ypred[pheno == 1], ypred[pheno == 0], curve = TRUE)
  
  ###############################################################################
  ## Test iRF 
  ###############################################################################

  # select subset of features that was used for iRF (top k features from RF with k giving best oob error)
  if(length(fit$rf.list$variable.importance) != ncol(geno)){
    ind <- sapply(names(fit$rf.list$variable.importance), function(x) which(colnames(geno) == x))
    geno.all <- geno
    geno <- geno[,ind]
    
    geno.all.train <- geno.train
    geno.train <- geno.train[ ,ind]
  }
  
  ypred <- predict(fit$rf.list, data=geno, predict.all=TRUE)
  ypred <- rowMeans(ypred$predictions)
  ypred.irf <- ypred
  
  pr.curve.irf <- pr.curve(ypred[pheno == 1], ypred[pheno == 0], curve = TRUE)
  roc.curve.irf <- roc.curve(ypred[pheno == 1], ypred[pheno == 0], curve = TRUE)
  
  save(file = paste0('../results/ypred_',res,'.Rdata'), ypred.lasso, ypred.ranger, ypred.irf)
  
  ###############################################################################
  ## Save data  
  ###############################################################################
  
  save(file = paste0('../data/data_',res,'.Rdata'), geno, pheno, geno.train, pheno.train, pca)
  
}else{
  load(file = paste0('../data/data_',res,'.Rdata'))
  
}

###############################################################################
## Filter stable intra chromosome interactions 
###############################################################################

# Load data base
source('utilities_loadDB.R')

indImp <- sapply(names(fit$rf.list$variable.importance), 
                 function(x) which(x == db.gene.chr.genename$gene))     #genes considered by iRF

int <- cleanInt(fit$interaction$int)      #remove sign from interaction name
int <- sapply(int, function(x) unlist(str_split(x, '_')))     #individual genes in interaction
ind.db.int <- sapply(int, match, db.gene.chr.genename$gene)     #indexes of interaction genes in data base
ind.geno.int <- sapply(int, match, colnames(geno))      #indexes of interaction genes in geno data matrix

ind.stab <- which(apply(fit$interaction, 1, function(x) all(x[10] >= 0.5)))     #indexes of interactions with stability score >= 0.5
ind.ichr <- which(sapply(ind.db.int, function(x) length(unique(db.gene.chr.genename$chr[x]))>1))     #indexes of intra chromosome interactions  
ind.stab.intra <- intersect(ind.stab, ind.ichr)     #indexes of stable intra chromosome interactions
ind.stab.inter <- setdiff(ind.stab, ind.ichr)       #indexes of stable inter chromosome interactions

###############################################################################
## calculate p-values  
###############################################################################

pMulti.intra <- vector('list', length(ind.stab.intra))
pMulti.inter <- vector('list', length(ind.stab.inter))

pCart.intra <- vector('list', length(ind.stab.intra))
pCart.inter <- vector('list', length(ind.stab.inter))

if(file.exists(paste0('../results/pval_all_',res,'.Rdata'))){
  load(paste0('../results/pval_all_',res,'.Rdata'))
}else{
  #compute p-values for stable intra chromosome interactions
  for(i in 1:length(ind.stab.intra)){
    print(paste("interaction", i, "out of", length(ind.stab.intra)))
    ind.int <- unlist(ind.geno.int[ind.stab.intra[i]])
    int <- unlist(fit$interaction$int[ind.stab.intra[i]])
    pMulti.intra[[i]] <- glmTestPredixMulti(geno = geno, pheno = pheno, 
                                            geno.train = geno.train, pheno.train = pheno.train, 
                                            ind.int = ind.int, fam = 'binomial', pca = pca)
    pCart.intra[[i]] <- pcsTestPredixCART(geno = geno, pheno = pheno, ind.int = ind.int, 
                                                geno.train = geno.train, pheno.train = pheno.train, 
                                                single.train = "train", pv.stat = "PCS")
  }
  
  save(file = paste0('../results/pval_all_',res,'.Rdata'), 
       pMulti.inter, pMulti.intra, 
       pCart.inter, pCart.intra, 
       ind.stab, ind.ichr, ind.stab.inter, ind.stab.intra, 
       db.gene.chr.genename)  
}


###############################################################################
##  compute p-values among bootstrap replicates   
###############################################################################
if(file.exists(paste0('../results/pval_stab_',res,'.Rdata'))){
  load(paste0('../results/pval_stab_',res,'.Rdata'))
}else{

  B = 10 #number of bootstrap replicates
  
  PMulti.intra <- vector('list', B)
  PCart.intra <- vector('list', B)
  
  set.seed(31)
  for(b in 1:B){
    print(paste("bootstrap run ", b, " out of ", B))
    
    ind <- sample(1:nrow(geno), size = nrow(geno), replace = T)     #bootstrap indexes
    genoB <- geno[ind,]
    phenoB <- pheno[ind]
    bpMulti.intra <- rep(NA, length(ind.stab.intra))
    bpCart.intra <- rep(NA, length(ind.stab.intra))

    for(i in 1:length(ind.stab.intra)){
      ind.int <- unlist(ind.geno.int[ind.stab.intra[i]])
      int <- unlist(fit$interaction$int[ind.stab.intra[i]])
      bpMulti.intra[i] <- glmTestPredixMulti(geno = genoB, pheno = phenoB, 
                                             geno.train = geno.train, pheno.train = pheno.train, 
                                             ind.int = ind.int, fam = 'binomial', pca = pca)
      bpCart.intra[i] <- pcsTestPredixCART(geno = genoB, pheno = phenoB, ind.int = ind.int, 
                                                 geno.train = geno.train, pheno.train = pheno.train, 
                                                 single.train = "train", pv.stat = "PCS")
    }
    PMulti.intra[[b]] <- bpMulti.intra
    PCart.intra[[b]] <- bpCart.intra
  }
  
  tsh.intra <- 0.05/length(ind.stab.intra)
  tsh.inter <- 0.05/length(ind.stab.inter)
  
  stab.pMulti.intra <- sapply(1:length(ind.stab.intra), function(i) mean(sapply(PMulti.intra, function(x) x[i] <= tsh.intra)))
  stab.pCart.intra <- sapply(1:length(ind.stab.intra), function(i) mean(sapply(PCart.intra, function(x) x[i] <= tsh.intra)))

  save(file = paste0('../results/pval_stab_',res,'.Rdata'),
       stab.pMulti.intra, stab.pCart.intra) 
}


###############################################################################
##  save p-value results          
###############################################################################

ind <- ind.stab.intra
pValues.intra <- data.frame(genes = sapply(ind.db.int[ind], function(x) paste(db.gene.chr.genename$name[x], collapse = " , ")),
                            chr = sapply(ind.db.int[ind], function(x) paste(db.gene.chr.genename$chr[x], collapse = " , ")),
                            pM = as.numeric(pMulti.intra), 
                            mseA_pM = sapply(pMulti.intra, function(x) attr(x, "mseA")),
                            mseH_pM = sapply(pMulti.intra, function(x) attr(x, "mseH")),
                            aucA_pM = sapply(pMulti.intra, function(x) attr(x, "aucA")),
                            aucH_pM = sapply(pMulti.intra, function(x) attr(x, "aucH")),
                            mDiffA_pM = sapply(pMulti.intra, function(x) attr(x, "mDiffA")),
                            mDiffH_pM = sapply(pMulti.intra, function(x) attr(x, "mDiffH")),
                            lr_pM = sapply(pMulti.intra, function(x) attr(x, "lr")),
                            pC = as.numeric(pCart.intra), 
                            mseA_pC = sapply(pCart.intra, function(x) attr(x, "mseA")),
                            mseH_pC = sapply(pCart.intra, function(x) attr(x, "mseH")),
                            aucA_pC = sapply(pCart.intra, function(x) attr(x, "aucA")),
                            aucH_pC = sapply(pCart.intra, function(x) attr(x, "aucH")),
                            mDiffA_pC = sapply(pCart.intra, function(x) attr(x, "mDiffA")),
                            mDiffH_pC = sapply(pCart.intra, function(x) attr(x, "mDiffH")),
                            lr_pC = sapply(pCart.intra, function(x) attr(x, "lr")),
                            stab.pMulti = stab.pMulti.intra, 
                            stab.pCart = stab.pCart.intra, 
                            fit$interaction[ind, 2:10]
)


save(file = paste0('../results/pval_',res,'.Rdata'), pValues.intra, 
     ind.stab, ind.ichr, ind.stab.inter, ind.stab.intra, db.gene.chr.genename)  
