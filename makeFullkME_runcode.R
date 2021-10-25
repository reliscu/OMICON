library(flexiblas)
flexiblas_load_backend("OPENBLASPTHREAD")
flexiblas_switch(2)

setwd("/home/rebecca/omicon/code/make_kme")

datExpr <- as.data.frame(fread("Make_full_kME_table/Sun_GSE4290_glioma_154_ComBat.csv"))
modEig <- as.data.frame(fread("Make_full_kME_table/Module_eigengenes_18-12-19.csv"))
kME <- as.data.frame(fread("Make_full_kME_table/kME_table_18-12-19.csv"))

makeFullkME(datExpr, 
            modEig, 
            kME)