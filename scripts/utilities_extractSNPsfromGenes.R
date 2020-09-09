library(data.table)
library(dplyr)
library(stringr)
library(Matrix)
library(fastmatch)
library(ape)
library(seriation)
library(gridExtra)
library(rgl)
library(pryr)
library(ensembldb)
library(IRanges)
library(EnsDb.Hsapiens.v75)

#Load Plink 
system('ml biology') 
system('ml plink/1.90b5.3') 

source('utilities_general.R')

###############################################################################
## Specify paths  
###############################################################################

res <- 'analysis_iRF_gene'      #name of results
path.db <- "../data/gtex_v7_Skin_Sun_Exposed_Lower_leg_imputed_europeans_tw_0.5_signif.db"      #path to data base
path.plink <- "path_to_plink_files/name_of_plink_files_chr"     #plink files from which SNPs should be extracted

###############################################################################
## Extract coordinates for selected genes
###############################################################################

#load data base
source('utilities_loadDB.R')

#load results from iRF
load(paste0('../results/',res,'.Rdata'))
int <- cleanInt(fit$interaction$int)
int <- sapply(int, function(x) unlist(str_split(x, '_')))
ind.db.int <- sapply(int, match, db.gene.chr.genename$gene)       #get indexes in data base from iRF genes
p <- length(fit$rf.list$variable.importance)      #number of iRF features (selected via best oob error)

#select genes
genesSelect <- as.character(db.gene.chr.genename$gene[unique(unlist(ind.db.int))])      #select iRF-interaction genes
genesSelect <- unique(c(genesSelect, 
                        names(sort(fit$rf.list$variable.importance, decreasing = T)[1:min(100, p)])))     #select top 100 iRF genes
genesSelect.chr <- db.gene.chr.genename$chr[match(genesSelect, db.gene.chr.genename$gene)]      #get chromosomes for genes


# get gene list and their coordinates from EnsDb.Hsapiens.v75 package
edb <- EnsDb.Hsapiens.v75
genes <- genes(edb)
rge <- ranges(genes)

# for each chromosome and each selected gene in chromosome get gene coordinates
coGenes <- vector('list', 22)
for(i in 1:22){
  genSelect <- genesSelect[genesSelect.chr == i]
  genSelect <- sapply(genSelect, clear_ensg)
  rgeS <- rge[names(rge) %in% genSelect]
  start(rgeS) <- start(rgeS) - 1000
  end(rgeS) <- end(rgeS) + 1000
  rgeSR <- IRanges::reduce(rgeS)
  coGen <- as.matrix(rgeSR)
  coGen[,2] <- coGen[,1] + coGen[,2] 
  coGenes[[i]] <- coGen
}

#extract SNPs for imputed data for each chromosome and for each region
for(i in 1:22){
  print(i)
  coGen  <- coGenes[[i]]
  if(nrow(coGen) == 0) next
  for(j in 1:nrow(coGen)){
   
    system(paste0('plink --bfile ', path.plink, i,' --chr ', i,' --from-bp ',
                  coGen[j,1],' --to-bp ', coGen[j, 2],' --write-snplist --out ', '../results/SNPs_from_iRF_per_region_',res, '_chr', i, '_region_', j))
  }
}

#combine SNPs from individual chromosomes and regions
system(paste0('cat ../results/SNPs_from_iRF_per_region_',res, '_chr* > SNPs_from_iRF_',res,'.snplist'))
system(paste0('rem ../results/SNPs_from_iRF_per_region_',res, '_chr*'))

