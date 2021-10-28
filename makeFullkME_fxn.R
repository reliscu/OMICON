corStats <- function(x, y) {
  cor.stats <- cor.test(x, y)
  return(c(cor.stats$estimate, cor.stats$p.value))
}

makeFullkME <- function(datExpr, 
                        modEig, 
                        kme,
                        nThreads = NULL) {
  
  library(parallel)
  
  if (is.null(nThreads)) nThreads <- detectCores()
  
  modEig <- modEig[match(colnames(datExpr)[-c(1, 2)], modEig$Index),]
  if (!identical(as.integer(as.numeric(colnames(datExpr)[-c(1, 2)])), modEig$Index)) {
    stop("Expression samples do not match module eigengene samples")
  }
  
  datExpr <- datExpr[match(kme$Index, datExpr[, 1]),]
  if (!identical(as.integer(datExpr[, 1]), kme$Index)) {
    stop("Expression features do match kME table features")
  }
  
  expr <- apply(t(datExpr[, -c(1, 2)]), 2, as.numeric)
  modEig <- modEig[, -c(1, 2)]
  
  tic()
  kmeStats <- mclapply(modEig, function(ME) t(apply(expr, 2, corStats, ME)), mc.cores = nThreads)
  toc()
  kmeStats <- as.data.frame(kmeStats) 
  colnames(kmeStats) <- paste0("kME", mapply(function(x) paste0(x, c(".cor", ".pval")), colnames(modEig)))
  
  kmeFull <- cbind.data.frame(kme, kmeStats)
  return(kmeFull)
  
}