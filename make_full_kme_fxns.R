make_full_kme <- function(analyte_data, module_eig, kme_partial, parallel=F, n_threads=NULL){

  module_eig$Index <- as.character(module_eig$Index)
  idx <- match(colnames(analyte_data)[-c(1,2)], module_eig$Index)
  module_eig <- module_eig[idx,]
  
  if(!identical(colnames(analyte_data)[-c(1,2)], module_eig$Index)){
    stop("Expression samples do not match module eigengene samples")
  }
  
  idx <- match(analyte_data$INDEX[-c(1)], kme_partial$Index)
  kme_partial <- kme_partial[idx,]
  
  if(!identical(analyte_data$INDEX[-c(1)], kme_partial$Index)){
    stop("Expression features do match kME table features")
  }
  
  analyte_data <- remove_dataset_indices(analyte_data)
  analyte_data <- analyte_data[,-c(1)]
  module_eig <- module_eig[,-c(1,2)]
  
  eig_vecs <- lapply(1:ncol(module_eig), function(i) module_eig[,i])
  
  if(parallel){
    
    kme_stats_list <- par_cor(analyte_data, eig_vecs, n_threads)
    
  } else {
    
    kme_stats_list <- lapply(eig_vecs, function(eigengene){
      t(apply(analyte_data, 1, cor_fxn, eigengene))
    })
    
  }

  kme_stats <- do.call(cbind, kme_stats_list)
  
  colnames(kme_stats) <- paste0("kME", mapply(function(x){
    paste0(x, c(".cor", ".pval"))
  }, colnames(module_eig)))
  
  return(cbind(kme_partial, kme_stats))
  
} ## make_full_kme_par <- function(

par_cor <- function(analyte_data, eig_vecs, n_threads){
  
  require("future.apply")
  
  options(future.globals.maxSize=Inf)
  
  if(is.null(n_threads)){
    n_threads <- detectCores()
  }
  
  plan(multisession, workers=n_threads)
  
  kme_stats_list <- future_lapply(eig_vecs, FUN=function(eigengene){
    t(apply(analyte_data, 1, cor_fxn, eigengene))
  })
  
  return(kme_stats_list)
  
} ## par_cor <- function(

remove_dataset_indices <- function(analyte_data){
  
  analyte_data <- analyte_data[,-c(1)]
  colnames(analyte_data) <- analyte_data[1,]
  analyte_data <- analyte_data[-c(1),]
  for(i in 2:ncol(analyte_data)){
    analyte_data[,i] <- as.numeric(analyte_data[,i])
  }
  return(analyte_data)
  
}

cor_fxn <- function(x, y){
  cor.stats <- cor.test(x, y)
  return(c(
    cor.stats$estimate, cor.stats$p.value
  ))
}
