library(matrixStats)
library(RSQLite)
library(data.table)
library(stringr)

sqlite.driver <- dbDriver("SQLite")

if (!exists('path.db')) stop('Specifify PrediXcan database path as: "path.db"')
if (!exists('path.bim.out')) print('Specifify output bim file with list of DB SNPs as: "path.bim.out"')

source('/accounts/campus/omer_ronen/projects/epiTree/scripts/utilities_general.R')

#loading database
db <- dbConnect(sqlite.driver, dbname = path.db)

#gene - chromosome mapping
db.gene.chr <- dbReadTable(db,"weights")[,c(1,3)]
db.gene.chr[,2] <- sapply(db.gene.chr[,2], function(x) as.integer(strsplit(x,"_")[[1]][1]) ) 
db.gene.chr <- unique(db.gene.chr)

#gene - genename mapping
db.gene.genename <- dbReadTable(db,"extra")[,c(1,2)]

#gene - genename - chromosome mapping
name <- sapply(db.gene.chr[,1], function(x) db.gene.genename[which(x == db.gene.genename[,1]), 2])
db.gene.chr.genename <- data.frame(gene = db.gene.chr[,1], name = name, chr = db.gene.chr[,2])
remove(name)

dbDisconnect(db)

#gene - snp mapping
db.gene.snp <- geneToSnp(db.gene.chr.genename$gene, path.db, gene.wise = T)

#chromosome - snp mapping
db.chr.snp <- sapply(split(db.gene.chr.genename, db.gene.chr.genename$chr), function(x) geneToSnp(x$gene, path.db, gene.wise = F))

# save list of snps per chromosome in data base
if(exists('path.bim.out')){
  write.table(unlist(db.chr.snp), file = path.bim.out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE) 
}

