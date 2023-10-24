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
# name.out <- 'analysis_iRF_gene.Rdata'      #name of output files
path.pca <- "/accounts/campus/omer_ronen/projects/epiTree/disease_phenotypes_full_PC40.csv"       #path to file with principle components of subjects


# Construct path for lasso and ranger output files
path.lasso <- file.path(path.out, "lasso.Rdata")
path.ranger <- file.path(path.out, "ranger.Rdata")
path.excl <- "/accounts/campus/omer_ronen/projects/epiTree/data/geneExclusionList.txt"      #optional list of genes that should be included
reloadData <- TRUE      #should full data set be loaded

source('/accounts/campus/omer_ronen/projects/epiTree/scripts/utilities_general.R')
source('/accounts/campus/omer_ronen/projects/epiTree/scripts/utilities_tests.R')

###############################################################################
## Specify paths  
###############################################################################

# res <- 'analysis_iRF_gene'      #name of results
# path.res <- '../results/'

# # path.key <- "path_to_mapping_file/ukb15860_13721_mapping.tsv"     #mapping file between subject IDs in pheno (1st column) and geno (2nd column) file
# data.field.pheno <- 'n_1747_0_0'      #data field of interest in phenotype file (here for red hair)
# code.cases.pheno <- 2    # red hair phenotype is encoded as 2
# code.na.pheno <- c(6, -1)     #NAs for red hair phenotype are 6, -1

# path.out <- '../results/'     #path for output files
# name.out <- paste0(res,'.Rdata')      #name of output files
# path.lasso <- paste0("../results/lasso_",name.out)      #output file for lasso results only
# path.ranger <- paste0("../results/ranger_",name.out)     #output file for ranger results only

###############################################################################
## Load results 
###############################################################################

# load(paste0(path.res,res,".Rdata"))

#reload full data set
if(reloadData){
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
    geno <- as.data.frame(as.matrix(geno))
    geno <- geno[,!duplicated(colnames(geno), fromLast = TRUE)]
    geno <- data.matrix(geno)
    ind.train <- read.csv(paste0(path.out, "/ind.train.csv"), header=TRUE)[, 1]

    # remove duplicated columns values


    #source(paste0('scripts/load_pheno.R'))      #load phenotype files
    pheno <- read.csv(path.pheno, header=FALSE)
    # set first column as row names
    rownames(pheno) <- pheno[,1]
    # match pheno and geno by rownames
    pheno <- pheno[match(rownames(geno), rownames(pheno)),]
    colnames(pheno) <- c("id", "pheno")
    id <- pheno$id
    pheno <- pheno$pheno

    numb_cases <- sum(pheno == 1)

      data_pca <- read.csv(path.pca, header=TRUE)
      # select only subjects that are in pheno
      data_pca <- data_pca[data_pca$participant.eid %in% id,]
      # select columns Genetic.PC.1 to Genetic.PC.15
      idx <- colnames(data_pca) %in% paste0("Genetic.PC", 1:15)

  # make ind.train be 80% of the data



  pca <- data_pca[, idx]
  load(paste0(path.out, "/analysis_iRF_gene.Rdata"))
  load(paste0(path.out, "/lasso.Rdata"))
  load(paste0(path.out, "/ranger.Rdata"))



  geno.train <- geno[ind.train,]
  pheno.train <- pheno[ind.train]
  
  geno.test <- geno[-ind.train,]
  pheno.test <- pheno[-ind.train]

  geno <- geno.test
  pheno <- pheno.test
  
  #gc()
  
  # Load PCAs for test data
#   sel.sub <- rownames(geno)
#   source('load_pca.R')
#   pca <- xpc[, 1:15]      #select first 15 principle components
  
  pca <- as.matrix(pca)
    pca.train <- pca[ind.train,]
    pca.test <- pca[-ind.train,]

  # Load lasso and ranger results

#   geno <- as.data.frame(as.matrix(geno))
#   geno <- geno[!duplicated(names(geno))]
#   geno.train <- as.matrix(geno.train)


  
  
  ###############################################################################
  ## Test lasso  
  ###############################################################################

  ypred <- predict(lasso, geno.test, type = "response")
  ypred.lasso <- ypred
  
  pr.curve.lasso <- pr.curve(ypred[pheno.test == 1], ypred[pheno.test == 0], curve = TRUE)
  roc.curve.lasso <- roc.curve(ypred[pheno.test == 1], ypred[pheno.test == 0], curve = TRUE)
  
  ###############################################################################
  ## Test ranger  
  ###############################################################################
  
  ypred <- predict(frang, data=geno.test, predict.all=TRUE)
  ypred <- rowMeans(ypred$predictions)
  ypred.ranger <- ypred
  
  pr.curve.ranger <- pr.curve(ypred[pheno.test == 1], ypred[pheno.test == 0], curve = TRUE)
  roc.curve.ranger <- roc.curve(ypred[pheno.test == 1], ypred[pheno.test == 0], curve = TRUE)
  
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
  
  ypred <- predict(fit$rf.list, data=geno.test, predict.all=TRUE)
  ypred <- rowMeans(ypred$predictions)
  ypred.irf <- ypred
  
  pr.curve.irf <- pr.curve(ypred[pheno.test == 1], ypred[pheno.test == 0], curve = TRUE)
  roc.curve.irf <- roc.curve(ypred[pheno.test == 1], ypred[pheno.test == 0], curve = TRUE)
  
  save(file = paste0(path.out, "/prediction.Rdata"), ypred.lasso, ypred.ranger, ypred.irf)
  
  ###############################################################################
  ## Save data  
  ###############################################################################
  
  save(file = paste0(path.out, "/pred_data.Rdata"), geno, pheno, geno.train, pheno.train, pca)
  
}else{
  load(file = paste0(path.out, "/pred_data.Rdata"))
  
}

###############################################################################
## Filter stable intra chromosome interactions 
###############################################################################

# Load data base
source('/accounts/campus/omer_ronen/projects/epiTree/scripts/utilities_loadDB.R')

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

pcs_pval_fname <- paste0(path.out, "/pcs_pval.Rdata")

if(file.exists(pcs_pval_fname)){
  load(pcs_pval_fname)
}else{
  #compute p-values for stable intra chromosome interactions
  for(i in 1:length(ind.stab.intra)){
    print(paste("interaction", i, "out of", length(ind.stab.intra)))
    ind.int <- unlist(ind.geno.int[ind.stab.intra[i]])
    int <- unlist(fit$interaction$int[ind.stab.intra[i]])
    pMulti.intra[[i]] <- glmTestPredixMulti(geno = geno, pheno = pheno, 
                                            geno.train = geno.train, pheno.train = pheno.train, 
                                            ind.int = ind.int, fam = 'binomial', pca = pca.train)
    pCart.intra[[i]] <- pcsTestPredixCART(geno = geno, pheno = pheno, ind.int = ind.int, 
                                                geno.train = geno.train, pheno.train = pheno.train, 
                                                single.train = "train", pv.stat = "PCS")
  }
  
  save(file = pcs_pval_fname,
       pMulti.inter, pMulti.intra, 
       pCart.inter, pCart.intra, 
       ind.stab, ind.ichr, ind.stab.inter, ind.stab.intra, 
       db.gene.chr.genename)  
}


###############################################################################
##  compute p-values among bootstrap replicates   
###############################################################################
pcs_pval_fname <- paste0(path.out, "/pval_bs.Rdata")

if(file.exists(pcs_pval_fname)){
  load(pcs_pval_fname)
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

  save(file = pcs_pval_fname,
       stab.pMulti.intra, stab.pCart.intra) 
}


###############################################################################
##  save p-value results          
###############################################################################
results_pval_fname <- paste0(path.out, "/pval_results.Rdata")

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


save(file = results_pval_fname, pValues.intra,
     ind.stab, ind.ichr, ind.stab.inter, ind.stab.intra, db.gene.chr.genename)  
