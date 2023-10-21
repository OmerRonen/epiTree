library(dplyr)
library(stringr)
library(Matrix)
library(fastmatch)

if (!exists('load.id')) stop('Specifify indices to load as: "load.id"')
if (!exists('path.predx')) stop('Specifify PrediXcan path as: "path.predx"')
if (all(diff(load.id) != 1)) stop('Only consecutive indices supported')


###############################################################################
## Load selected patient genotype data from PrediXcan
###############################################################################
print(path.predx)
# Load gene identifiers
gnames <- fread(path.predx[1], nrow=1, header=F)
gnames <- as.character(gnames)[-c(1,2)]

# Load genes to be excluded
if (exists('path.excl')) {
  print("load gene-exclusion-list")
  gene.ex <- fread(path.excl, skip=1, header=F)
  gene.ex <- gene.ex[[1]]
  
  ind.ex <- fmatch(gene.ex, gnames)
  ind.ex <- ind.ex[!is.na(ind.ex)]
}

skip <- min(load.id) - 1 + 4
N <- max(load.id) - (skip - 4)

geno <- 0
active <- numeric(0)
gchrom <- list()

for(i in 1:length(path.predx)) {
  print(paste("load PrediXcan file", i))
  xin <- fread(path.predx[i], skip=skip, nrow=N)
  if(i == 1){
    print("load subject names")
    subnames <- as.character(unlist(xin[,2]))
  }
  if(!identical(subnames, as.character(unlist(xin[,2])))){
    warning("Subjects in files not identical!")
  }
  xin <- Matrix(as.matrix(xin[,-(1:2)]), sparse=TRUE)
  # Update genotype matrix with gene expression from current chromosome
  geno <- geno + xin
  gc(verbose=FALSE)
}


# Remove genes to be excluded
if (exists('path.excl')) {
  print(paste("Removing", length(ind.ex), "out of", ncol(geno), "genes"))
  if(length(ind.ex) > 0){
    print(paste("Removed genes:", gnames[ind.ex]))
    geno <- geno[, -ind.ex]
    gnames <- gnames[-ind.ex]
  }
}

colnames(geno) <- gnames
rownames(geno) <- subnames
