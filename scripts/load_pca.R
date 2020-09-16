library(stringr)
library(dplyr)
library(data.table)
library(Matrix)


# Check for valid inputs
if (!exists('sel.sub')) stop('Specify a set of selected subjects as sel.sub')
if (!exists('path.pca')) stop('Specify pca file as path.pca')
if (!exists('app.geno')) stop('Specify application number used in pca file as app.geno')


# Load in matrix of principle components
xpc <- fread(path.pca, select=c('IDs', paste0('PC', 1:40))) %>% filter(app.geno %in% sel.sub)
xpc <- as.matrix(xpc[match(sel.sub, xpc$app.geno), 2:ncol(xpc)])
xpc <- Matrix(xpc, sparse=TRUE)
rownames(xpc) <- sel.sub
colnames(xpc) <- paste0('PC', 1:40)
