source("/home/rebecca/omicon/code/make_kme/makeFullkME_fxn.R")
library(tictoc) ## Calc. time for function to run
library(data.table)

setwd("/home/rebecca/omicon/code/make_kme")
datExpr <- as.data.frame(fread("Sun_GSE4290_glioma_154_ComBat.csv"))
modEig <- as.data.frame(fread("Module_eigengenes_18-12-19.csv"))
kme <- as.data.frame(fread("kME_table_18-12-19.csv"))
nThreads <- 15

# tic()
# fullkME <- makeFullkME(datExpr, modEig, kme, nThreads)
# toc()

tic()
fullkME <- makeFullkMESerial(datExpr, modEig, kme)
toc()

fwrite(fullkME, file = "Full_kME_table.csv") ## Write file to current working directory