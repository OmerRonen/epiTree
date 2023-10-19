library(RSQLite)

# specify path of data base file
name.db <- "ctimp_Brain_Cortex"
path.db <- paste0("data/",name.db ,".db")


sqlite.driver <- dbDriver("SQLite")

db <- dbConnect(sqlite.driver, dbname = path.db)
db_weights <- dbReadTable(db,"weights")
db_extra <- dbReadTable(db,"extra")
dbDisconnect(db)

snpList <- db_weights$rsid

# save list of SNPs in a text file
write.table(snpList, file = paste0("data/snpList_", name.db, ".snplist"), quote = F, col.names = F, row.names = F)
