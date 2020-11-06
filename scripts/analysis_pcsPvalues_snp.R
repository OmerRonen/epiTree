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
library(purrr)
library(snpStats)


reloadData <- TRUE      #should full data set be loaded

source('utilities_general.R')
source('utilities_tests.R')

sqlite.driver <- dbDriver("SQLite")

###############################################################################
## Specify paths  
###############################################################################
path.plink <- paste0("path_to_output_plink_file/name_of_output_plink_file_chr",1:22)      #plink files with snp data

path.pheno <- "path_to_phenotype_file/app15860_standard_data_2016Nov19.txt"     #path to phenotype file
path.key <- "path_to_mapping_file/ukb15860_13721_mapping.tsv"     #mapping file between subject IDs in pheno (1st column) and geno (2nd column) file
data.field.pheno <- 'n_1747_0_0'      #data field of interest in phenotype file (here for red hair)
code.cases.pheno <- 2    # red hair phenotype is encoded as 2
code.na.pheno <- c(6, -1)     #NAs for red hair phenotype are 6, -1

res <- 'analysis_iRF_snp'      #name of results
path.res <- '../results/'
name.out <- paste0(res,'.Rdata')
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
  # we assume to have one plink file per chromosome
  # we load a balanced sample of 15K cases and controls
  # the first 13K sample of each are the training data
  
  chr <- 1
  geno <- data.frame()
  
  while(length(geno) == 0){
    if(file.exists(paste0(path.plink[chr], '.bim'))){   #we may not have selected SNPs from every chromosome
      print(paste('load initial chr', chr))
      snp <- read.table(paste0(path.plink[chr], ".bim"))[,2]
      is.dup <- duplicated(snp)
      geno <- read.plink(path.plink[chr], select.snps = which(!is.dup))$genotypes
      
      source('load_pheno.R')
      
      load.id <-  c(which(pheno == 0)[1:15000], which(pheno == 1)[1:15000])
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
  print(paste0("Dim of pca:", dim(pca)))
  
  geno <- as.matrix(geno)
  geno.train <- as.matrix(geno.train)
  
  gnames <- colnames(geno)
  
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

getChr <- function(snpName){
  # get chromosome from SNP name
  x <- str_locate(snpName, 'ch')[1]
  if(!is.element(x, c(2,3))){
    stop('invalide snp name')
  }
  chr <- as.numeric(str_sub(snpName, 1, x-1))
}

chrGeno <- sapply(colnames(geno), getChr)

int <- cleanInt(fit$interaction$int)      #remove sign from interaction name
int <- sapply(int, function(x) unlist(str_split(x, '_')))       #individual SNPs in interaction
ind.geno.int <- sapply(int, match, colnames(geno))      #indexes of interaction SNPs in geno data matrix

ind.stab <- which(apply(fit$interaction, 1, function(x) all(x[c(10)] >= 0.5)))      #indexes of interactions with stability score >= 0.5
ind.ichr <- which(sapply(ind.geno.int, function(x) length(unique(chrGeno[x]))>1))     #indexes of intra chromosome interactions
ind.stab.intra <- intersect(ind.stab, ind.ichr)       #indexes of stable intra chromosome interactions
ind.stab.inter <- setdiff(ind.stab2, ind.ichr)      #indexes of stable inter chromosome interactions


###############################################################################
## calculate p-values  
###############################################################################

pFull.intra <- vector('list', length(ind.stab.intra))
pDom.intra <- vector('list', length(ind.stab.intra))
pRec.intra <- vector('list', length(ind.stab.intra))

pFullPCS.intra <- vector('list', length(ind.stab.intra))
pDomPCS.intra <- vector('list', length(ind.stab.intra))
pRecPCS.intra <- vector('list', length(ind.stab.intra))


if(file.exists(paste0('../results/pval_all_',res,'.Rdata'))){
  load(paste0('../results/pval_all_',res,'.Rdata'))
}else{
  #compute logistic regression p-values for stable intra chromosome interactions
  for(i in 1:length(ind.stab.intra)){
    print(paste("interaction", i, "out of", length(ind.stab.intra)))
    ind.int <- unlist(ind.geno.int[ind.stab.intra[i]])
    int <- unlist(fit$interaction$int[ind.stab.intra[i]])
    pFullPCS.intra[[i]] <- pcsTestSnp(geno = geno, pheno = pheno, ind.int = ind.int, 
                                            geno.train = geno.train, pheno.train = pheno.train, model = 'full')
    pDomPCS.intra[[i]] <- pcsTestSnp(geno = geno, pheno = pheno, ind.int = ind.int, 
                                           geno.train = geno.train, pheno.train = pheno.train, model = 'dom')
    pRecPCS.intra[[i]] <- pcsTestSnp(geno = geno, pheno = pheno, ind.int = ind.int, 
                                           geno.train = geno.train, pheno.train = pheno.train, model = 'rec')
  }
 
  #compute PCS p-values for stable intra chromosome interactions
  print('intra chromosoe')
  for(i in 1:length(ind.stab.intra)){
    print(paste("interaction", i, "out of", length(ind.stab.intra)))
    ind.int <- unlist(ind.geno.int[ind.stab.intra[i]])
    int <- unlist(fit$interaction$int[ind.stab.intra[i]])
    pFull.intra[[i]] <- glmTestSnp(geno = geno, pheno = pheno, ind.int = ind.int,
                                   geno.train = geno.train, pheno.train = pheno.train,
                                   model = 'full', pca = pca)
    pDom.intra[[i]] <- glmTestSnp(geno = geno, pheno = pheno, ind.int = ind.int,
                                  geno.train = geno.train, pheno.train = pheno.train,
                                  model = 'dom', pca = pca)
    pRec.intra[[i]] <- glmTestSnp(geno = geno, pheno = pheno, ind.int = ind.int,
                                  geno.train = geno.train, pheno.train = pheno.train,
                                  model = 'rec', pca = pca)
  
  }

  save(file = paste0('../results/pval_all_',res,'.Rdata'), pFull.intra, pDom.intra, pRec.intra,
       pFullPCS.intra, pDomPCS.intra, pRecPCS.intra,
       ind.stab, ind.ichr, ind.stab.inter, ind.stab.intra)
}


###############################################################################
##  compute p-values among bootstrap replicates   
###############################################################################

if(file.exists(paste0('../results/pval_stab_',res,'.Rdata'))){
  load(paste0('../results/pval_stab_',res,'.Rdata'))
}else{

  B = 10 #number of bootstrap replicates
  
  PFull.intra <- vector('list', B)
  PDom.intra <- vector('list', B)
  PRec.intra <- vector('list', B)
  
  PFullPCS.intra <- vector('list', B)
  PDomPCS.intra <- vector('list', B)
  PRecPCS.intra <- vector('list', B)
  
  set.seed(31)
  for(b in 1:B){
    print(paste("Bootstrap run ", b, " out of ", B))
    
    ind <- sample(1:nrow(geno), size = nrow(geno), replace = T)       #bootstrap indexes
    genoB <- geno[ind,]
    phenoB <- pheno[ind]
    
    bpFull.intra <- rep(NA, length(ind.stab.intra))
    bpDom.intra <- rep(NA, length(ind.stab.intra))
    bpRec.intra <- rep(NA, length(ind.stab.intra))
    
    bpFullPCS.intra <- rep(NA, length(ind.stab.intra))
    bpDomPCS.intra <- rep(NA, length(ind.stab.intra))
    bpRecPCS.intra <- rep(NA, length(ind.stab.intra))
    
    for(i in 1:length(ind.stab.intra)){
      ind.int <- unlist(ind.geno.int[ind.stab.intra[i]])
      int <- unlist(fit$interaction$int[ind.stab.intra[i]])
      
      bpFull.intra[i] <- glmTestSnp(geno = genoB, pheno = phenoB, ind.int = ind.int,
                                    geno.train = geno.train, pheno.train = pheno.train,
                                    pca = pca)
      bpDom.intra[i] <- glmTestSnp(geno = genoB, pheno = phenoB, ind.int = ind.int,
                                   geno.train = geno.train, pheno.train = pheno.train,
                                   model = 'dom', pca = pca)
      bpRec.intra[i] <- glmTestSnp(geno = genoB, pheno = phenoB, ind.int = ind.int,
                                   geno.train = geno.train, pheno.train = pheno.train,
                                   model = 'rec', pca = pca)
      
      bpFullPCS.intra[i] <- pcsTestSnp(geno = genoB, pheno = phenoB, ind.int = ind.int, 
                                       geno.train = geno.train, pheno.train = pheno.train, model = 'full')
      bpDomPCS.intra[i] <- pcsTestSnp(geno = genoB, pheno = phenoB, ind.int = ind.int, 
                                      geno.train = geno.train, pheno.train = pheno.train, model = 'dom')
      bpRecPCS.intra[i] <- pcsTestSnp(geno = genoB, pheno = phenoB, ind.int = ind.int, 
                                      geno.train = geno.train, pheno.train = pheno.train, model = 'rec')
    }
    
    PFull.intra[[b]] <- bpFull.intra
    PDom.intra[[b]] <- bpDom.intra
    PRec.intra[[b]] <- bpRec.intra
    
    PFullPCS.intra[[b]] <- bpFullPCS.intra
    PDomPCS.intra[[b]] <- bpDomPCS.intra
    PRecPCS.intra[[b]] <- bpRecPCS.intra
  }
  
  tsh.intra <- 0.05/length(ind.stab.intra)
  tsh.inter <- 0.05/length(ind.stab.inter)

  stab.pFull.intra <- sapply(1:length(ind.stab.intra), function(i) mean(sapply(PFull.intra, function(x) x[i] <= tsh.intra)))
  stab.pFullPCS.intra <- sapply(1:length(ind.stab.intra), function(i) mean(sapply(PFullPCS.intra, function(x) x[i] <= tsh.intra)))

  stab.pDom.intra <- sapply(1:length(ind.stab.intra), function(i) mean(sapply(PDom.intra, function(x) x[i] <= tsh.intra)))
  stab.pDomPCS.intra <- sapply(1:length(ind.stab.intra), function(i) mean(sapply(PDomPCS.intra, function(x) x[i] <= tsh.intra)))

  stab.pRec.intra <- sapply(1:length(ind.stab.intra), function(i) mean(sapply(PRec.intra, function(x) x[i] <= tsh.intra)))
  stab.pRecPCS.intra <- sapply(1:length(ind.stab.intra), function(i) mean(sapply(PRecPCS.intra, function(x) x[i] <= tsh.intra)))

  save(file = paste0('../results/pval_stab_',res,'.Rdata'),
        stab.pFull.intra, stab.pFullPCS.intra, 
        stab.pDom.intra, stab.pDomPCS.intra,
        stab.pRec.intra, stab.pRecPCS.intra) 

}


###############################################################################
##  save p-value results          
###############################################################################

   
ind <- ind.stab.intra
pValues.intra <- data.frame( fit$interaction[ind, 1:10],
                            pF = as.numeric(pFull.intra),
                            aucA_pF = sapply(pFull.intra, function(x) attr(x, "aucA")),
                            aucH_pF = sapply(pFull.intra, function(x) attr(x, "aucH")),
                            mDiffA_pF = sapply(pFull.intra, function(x) attr(x, "mDiffA")),
                            mDiffH_pF = sapply(pFull.intra, function(x) attr(x, "mDiffH")),
                            lr_pF = sapply(pFull.intra, function(x) attr(x, "lr")),
                            stab.pF = stab.pFull.intra,
                            pD = as.numeric(pDom.intra),
                            aucA_pD = sapply(pDom.intra, function(x) attr(x, "aucA")),
                            aucH_pD = sapply(pDom.intra, function(x) attr(x, "aucH")),
                            mDiffA_pD = sapply(pDom.intra, function(x) attr(x, "mDiffA")),
                            mDiffH_pD = sapply(pDom.intra, function(x) attr(x, "mDiffH")),
                            lr_pD = sapply(pDom.intra, function(x) attr(x, "lr")),
                            stab.pD = stab.pDom.intra,
                            pR = as.numeric(pRec.intra),
                            aucA_pR = sapply(pRec.intra, function(x) attr(x, "aucA")),
                            aucH_pR = sapply(pRec.intra, function(x) attr(x, "aucH")),
                            mDiffA_pR = sapply(pRec.intra, function(x) attr(x, "mDiffA")),
                            mDiffH_pR = sapply(pRec.intra, function(x) attr(x, "mDiffH")),
                            lr_pR = sapply(pRec.intra, function(x) attr(x, "lr")),
                            stab.pR = stab.pRec.intra,
                            pF_PCS = as.numeric(pFullPCS.intra),
                            aucA_pF_PCS = sapply(pFullPCS.intra, function(x) attr(x, "aucA")),
                            aucH_pF_PCS = sapply(pFullPCS.intra, function(x) attr(x, "aucH")),
                            mDiffA_pF_PCS = sapply(pFullPCS.intra, function(x) attr(x, "mDiffA")),
                            mDiffH_pF_PCS = sapply(pFullPCS.intra, function(x) attr(x, "mDiffH")),
                            lr_pF_PCS = sapply(pFullPCS.intra, function(x) attr(x, "lr")),
                            stab.pF_PCS = stab.pFullPCS.intra,
                            pD_PCS = as.numeric(pDomPCS.intra),
                            aucA_pD_PCS = sapply(pDomPCS.intra, function(x) attr(x, "aucA")),
                            aucH_pD_PCS = sapply(pDomPCS.intra, function(x) attr(x, "aucH")),
                            mDiffA_pD_PCS = sapply(pDomPCS.intra, function(x) attr(x, "mDiffA")),
                            mDiffH_pD_PCS = sapply(pDomPCS.intra, function(x) attr(x, "mDiffH")),
                            lr_pD_PCS = sapply(pDomPCS.intra, function(x) attr(x, "lr")),
                            stab.pD_PCS = stab.pDomPCS.intra,
                            pR_PCS = as.numeric(pRecPCS.intra),
                            aucA_pR_PCS = sapply(pRecPCS.intra, function(x) attr(x, "aucA")),
                            aucH_pR_PCS = sapply(pRecPCS.intra, function(x) attr(x, "aucH")),
                            mDiffA_pR_PCS = sapply(pRecPCS.intra, function(x) attr(x, "mDiffA")),
                            mDiffH_pR_PCS = sapply(pRecPCS.intra, function(x) attr(x, "mDiffH")),
                            lr_pR_PCS = sapply(pRecPCS.intra, function(x) attr(x, "lr")),
                            stab.pR_PCS = stab.pRecPCS.intra
                            
)
                            
                            
save(file = paste0('../results/pval_',res,'.Rdata'), pValues.intra,
     ind.stab, ind.ichr, ind.stab.inter, ind.stab.intra) 
