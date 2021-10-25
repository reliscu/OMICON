library(flexiblas)
flexiblas_load_backend("OPENBLASPTHREAD")
flexiblas_switch(2)

setwd("/home/rebecca/omicon/code/make_kme")

datExpr <- as.data.frame(fread("Make_full_kME_table/Sun_GSE4290_glioma_154_ComBat.csv"))
modEig <- as.data.frame(fread("Make_full_kME_table/Module_eigengenes_18-12-19.csv"))
kME <- as.data.frame(fread("Make_full_kME_table/kME_table_18-12-19.csv"))
  
corStats <- function(x, y) {
  cor.stats <- cor.test(x, y)
  return(c(cor.stats$estimate, cor.stats$p.value))
}

makeFullkME <- function(datExpr, 
                        modEig, 
                        kME) {
  
  library(data.table)
  library(WGCNA)
  
  modEig <- modEig[match(colnames(datExpr)[-c(1, 2)], modEig$Index),]
  if(!identical(as.integer(as.numeric(colnames(datExpr)[-c(1, 2)])), modEig$Index)) {
    stop("Expression samples do not match module eigengene samples")
  }
  
  datExpr <- datExpr[match(kME$Index, datExpr[, 1]),]
  if(!identical(as.integer(datExpr[, 1]), kME$Index)) {
    stop("Expression features do match kME table features")
  }
  
  expr <- apply(t(datExpr[, -c(1, 2)]), 2, as.numeric)
  modEig <- modEig[, -c(1, 2)]
  
  kme.stats <- lapply(modEig, function(ME) t(apply(expr, 2, corStats, ME)))
  kme.stats <- as.data.frame(kme.stats) 
  colnames(kme.stats) <- paste0("kME", unlist(lapply(colnames(modEig), function(x) paste0(x, c(".cor", ".pval")))))
  
  kme.full <- cbind.data.frame(kME, kme.stats)
  return(kme.full)
  
} ## makeFullkME() {