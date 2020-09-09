intPairs <- function(int) {
  # Generate all pairwise interactions for LD analysis
  require(stringr)

  # Clean interaction for consistent naming with databse
  int <- cleanInt(int)
  int <- str_replace_all(int, '\\.\\.', '-')
  int <- str_split(int, '_')[[1]]

  # Generate all pairwise combinations of a given interactions
  #return(int)
  int.pairs <- t(combn(int, 2))
  return(int.pairs)
}

getLD <- function(f) { 
  # Reads plink output to get D' and R for each pair of SNPs 
  require(data.table)
  require(stringr)

  x <- fread(file=f, skip=27, nrows=1)
  if (nrow(x) == 0) return(NULL)
  pair <- tail(str_split(f, '/')[[1]], 1) %>% str_remove('\\.log')
  out <- data.table(pair=pair,  R=x$V3, D=x$V6)
  return(out)
}

bestPair <- function(x, int) {
  # Get maximum absolute correlation between any pair of SNPs in an
  # interaction
  require(stringr)
  int <- str_split(int, '_')[[1]]
  pairs <- combn(int, 2, simplify=FALSE)
  pairs <- lapply(pairs, sort)
  pairs <- sapply(pairs, paste, collapse='_')
  xf <- filter(x, pair %in% pairs)
  if (nrow(xf) == 0) return(0)
  id <- which.max(abs(xf$R))
  return(c(xf$R[id], xf$D[id]))
}

intChr <- function(ints, bim) {
  # Determine the chromosomes that interacting SNPs appear on
  # args:
  #   xint: character vector of interactions as returned by iRF
  #   bim: bim file of SNP metadata
  require(stringr)
  require(dplyr)
  library(fastmatch)
  ints <- cleanInt(ints)
  ints <- str_replace_all(ints, '\\.\\.', '-')
  ints <- str_split(ints, '_')
  int.chr <- sapply(ints, function(i) bim$V1[bim$V2 %fin% i])
  out <- sapply(int.chr, paste, collapse='_')
  return(out)
}

loadImp <- function(path) {
  # Load in variable importances from a given path
  load(path)
  out <- data.table(snp=names(importance), imp=importance)
  return(out)
}

filterInt <- function(x, tcpe=0.5, tfsd=0.5, tmip=0.5) {
  # Filters a dataframe of interactions based on stability scores of null
  # importance metrics.
  # args:
  #   x: dataframe of interactions as returned by iRF
  #   tcpe, tfsd, tmip: thresholds for corresponding null metrics
  require(dplyr)
  out <- filter(x, sta.cpe >= tcpe, sta.fsd >= tfsd, sta.mip >= tmip)
  return(out)
}

cleanInt <- function(x) {
  # Removes sign information from interactions
  require(stringr)
  return(str_remove_all(x, '[-\\+]'))
}

intTable <- function(geno, pheno, int) {
  # Generates a table of phenotype proportions and counts for all SNP 
  # combinatios of a given interaction.
  #   args:
  #     geno: genotype matrix of SNPs (in dosage format)
  #     pheno: phenoype response vector
  #     int: character vector of interacting SNPs
  require(stringr)

  int <- unlist(str_split(int, '_'))
  xint <- geno[,int]
  idrm <- rowSums(apply(xint, MAR=2, function(z) !z %in% c(0, 1, 2)))
  xint <- xint[!idrm,]
  pheno <- pheno[!idrm]
  xint <- apply(xint, MAR=1, paste, collapse='_')
  
  tint <- table(xint)
  pint <- c(by(pheno, xint, mean))

  return(list(p=pint, n=tint))
}


geneToSnp <- function(gene, path.db, gene.wise = F){
  # Returns the collection of SNPs in a a database that correspond to a given
  # set of genes.
  # args:
  #   gene: character vector of genes (in "E..." notation) for which to get SNPs
  #   path.dp: path to gene/SNP database used for PrediXcan 
  #   gene.wise: should each gene return a seperate list of SNPs
 if(!exists('sqlite.driver')){
  require(SQLite)
  sqlite.driver <- dbDriver("RSQLite")
  }
  db <- dbConnect(sqlite.driver, dbname = path.db)
  gene.snp <- dbReadTable(db, "weights")
  dbDisconnect(db)

  gene.snp <- split(gene.snp, gene.snp$gene)
  gene.names <- unname(sapply(gene.snp, function(x) x$gene[1]))
  snp.names <- lapply(gene.snp, function(x) x$rsid)
  ind <- which(str_detect(gene.names, paste(gene, collapse = '|')))
  if(gene.wise){
    return(snp.names[ind])
  }else{
    return(unique(unlist(snp.names[ind])))
  }
}

isSubsetOfAny <- function(set, list){
  #is a given set subset of any element in a list
  return(any(sapply(list, function(x) length(setdiff(x, set)) == length(x) - length(set) && length(setdiff(x, set)) > 0)))
}

impMissingSnp <- function(x){
  #impute missing SNPs by median
  x[is.na(x)] <- median(x, na.rm = T)
  return(x)
}

clear_ensg <- function(x){
  # clear ensable ID for version number
  unlist(str_split(x, "\\."))[1]
}


cleanName <- function(x) {
  #replace all '_', '+', and '-' from name
  require(stringr)
  
  x <- str_replace_all(x, fixed("-"), "..")
  x <- str_replace_all(x, fixed("_"), "..")
  x <- str_replace_all(x, fixed("+"), "..")
  
  return(x)
}

comb_inout <- function(ind_outer, ind_inner){
  #function to combine inner and outer indeces
  ind_outer[ind_outer != 0] <- ind_inner
  ans <- c()
  for(i in 1:length(ind_outer)){
    if(ind_outer[i] == 0){
      ans <- c(ans, c(0,0))
    }
    if(ind_outer[i] == 1){
      ans <- c(ans, c(1,0))
    }
    if(ind_outer[i] == 2){
      ans <- c(ans, c(0,1))
    }
  }
  return(ans)
}


geneSNP <- function(location, chr){
  # for a given location and chromosome return (of some SNP) return gene(s) that overlap with that location
  ind <- which(geneLocation$chr == chr)
  if(length(ind) > 0){
    snpDist <- sapply(ind, function(x) intervalDist(location, geneLocation$start[x], geneLocation$end[x]))
    if(min(snpDist) == 0){
      ind_gene <- which(snpDist == 0)
      gene <- paste(geneLocation$genename[ind[ind_gene]], collapse = "_")
    }else{
      gene <- "NA"
    }
    return(gene)
  }else{
    return(Inf)
  }
}

intervalDist <- function(val, start, end){
  # distance from value val to interval [start, end]
  if(val >= start & val <= end){
    return(0)
  }else{
    return(min(abs(val - start), abs(val - end)))
  }
}



# several rename functions:

intNameChr <- function(int){
  int <- cleanInt(int)
  int <- sapply(int, function(x) unlist(str_split(x, '_')))
  ind.db.int <- sapply(int, match, db.gene.chr.genename$gene) 
  name <- db.gene.chr.genename$name[ind.db.int]
  chr <- db.gene.chr.genename$chr[ind.db.int]
  return(paste(name, "(", chr, ")", collapse = " "))
}

intName <- function(int){
  int <- cleanInt(int)
  int <- sapply(int, function(x) unlist(str_split(x, '_')))
  ind.db.int <- sapply(int, match, db.gene.chr.genename$gene) 
  name <- db.gene.chr.genename$name[ind.db.int]
  return(paste(name, collapse = " + "))
}

intName2 <- function(int){
  int <- cleanInt(int)
  int <- sapply(int, function(x) unlist(str_split(x, '_')))
  ind.db.int <- sapply(int, match, db.gene.chr.genename$gene) 
  name <- db.gene.chr.genename$name[ind.db.int]
  return(paste(name, collapse = " _ "))
}

intNameSNP <- function(int){
  int <- cleanInt(int)
  int <- sapply(int, function(x) unlist(str_split(x, '_')))
  ind.db.int <- sapply(int, match, db.gene.snp.chr.genename$snp) 
  name <- db.gene.snp.chr.genename$name[ind.db.int]
  name <- levels(name)[name]
  snp <- db.gene.snp.chr.genename$snp[ind.db.int]
  snp <- levels(snp)[snp]
  chr <- db.gene.snp.chr.genename$chr[ind.db.int]
  return(paste(paste( paste(snp, collapse = " + "),"(", paste(chr, name, collapse = " "), ")")))
}


intNameSNP2 <- function(x){
  # variation of intNameSNP for database information from bim file
  
  x <- cleanInt(x)
  x <- unlist(str_split(x, '_'))
  ind.db.int <- sapply(x, match, snpL$rsid) 
  name <- snpL$rsid[ind.db.int]
  chr <- snpL$chr[ind.db.int]
  return(paste(paste( paste(name, collapse = " + "),"(", paste(chr, collapse = " "), ")")))
}


renameSNP <- function(snpName, chr){
  #rename SNPs, replace ':' by 'chr' for chromosome, '_' by 'b' for base pair and add chr number to name
  
  snpName <- str_replace_all(snpName, ":", "ch")
  snpName <- str_replace_all(snpName, "_", "b")
  if(str_sub(snpName, 1,2) == "rs"){
    snpName <- paste0(chr, 'ch', snpName)
  }
  return(snpName)
}

