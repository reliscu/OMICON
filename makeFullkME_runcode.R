source("/home/rebecca/omicon/code/make_kme/makeFullkME_fxn.R")
library(tictoc) ## Calc. time for function to run

setwd("/home/rebecca/omicon/code/make_kme")
datExpr <- as.data.frame(fread("Make_full_kME_table/Sun_GSE4290_glioma_154_ComBat.csv"))
modEig <- as.data.frame(fread("Make_full_kME_table/Module_eigengenes_18-12-19.csv"))
kme <- as.data.frame(fread("Make_full_kME_table/kME_table_18-12-19.csv"))
nThreads <- 15

tic()
fullkME <- makeFullkME(datExpr, modEig, kme, nThreads)
toc()