library(data.table)
library(tictoc)

source("~/make_full_kme_fxns.R")

analyte_data <- fread("data/E-MTAB-3892_Tumor_169_ComBat.csv", data.table=F)
module_eig <- fread("data/Module_eigengenes_08-39-41.csv", data.table=F)
kme_partial <- fread("data/kME_table_08-39-41.csv", data.table=F)
n_threads <- 15
parallel <- T

tic()
full_kme <- make_full_kme(analyte_data, module_eig, kme_partial, parallel, n_threads)
toc()

fwrite(full_kme, file = "data/full_kME_table.csv")