library(RSQLite)

# specify path of data base file
name.db <- "gtex_v7_Skin_Sun_Exposed_Lower_leg_imputed_europeans_tw_0.5_signif"
path.db <- paste0("data/",name.db ,".db")


sqlite.driver <- dbDriver("SQLite")

db <- dbConnect(sqlite.driver, dbname = path.db)
db_weights <- dbReadTable(db,"weights")
db_extra <- dbReadTable(db,"extra")
dbDisconnect(db)

snpList <- db_weights$rsid

# save list of SNPs in a text file
write.table(snpList, file = paste0("data/snpList_", name.db, ".snplist"), quote = F, col.names = F, row.names = F)
