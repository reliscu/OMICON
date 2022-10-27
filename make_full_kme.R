library(data.table)
library(tictoc)

source("~/make_full_kme_fxns.R")

analyte_data <- fread("~/E-MTAB-3892_Tumor_169_ComBat.csv", data.table=F)
module_eig <- fread("~/Module_eigengenes_08-39-41.csv", data.table=F)
kme_partial <- fread("~/kME_table_08-39-41.csv", data.table=F)
n_threads <- 15
parallel <- T

full_kme <- make_full_kme(analyte_data, module_eig, kme_partial, parallel, n_threads)

fwrite(full_kme, file = "/path/to/output/directory/Full_kME_table.csv")