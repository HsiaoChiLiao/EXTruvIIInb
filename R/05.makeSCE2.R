#' Create SingleCellExperiment object containing normalised data for matched sequencing count data
#'
#' This function takes extFastruvIIInb_vanilla, extFastruvIIInb_ortho, or extFastruvIIInb_cca output as input and creates a SingleCellExperiment (SCE) object containing various metrics of normalised data.  
#'
#' @param obj  object containing output of call to ruvIII.nb or fastruvIII.nb function.
#' @param cData A data frame containing cell-level metadata. This data frame will be used as 'colData' in the SingleCellExperiment object and must be specified.
#' @param batch numeric vector containing batch information for each sample. Must correspond to columns of count matrix. Only needed if batch-specific dispersion parameter is fitted.
#' @param assays type of corrected data to return. Available types are "logcounts","logPAC" and "pearson".
#' 
#' @importFrom methods as
#' @importFrom stats cor dnbinom mad median pnbinom qnbinom quantile uniroot
#' 
#' @return A SingleCellExperiment object with normalized data added to the 'assays' slot. The log percentile-adjusted count is stored in the 'logPAC' component of the 'assays' slot, the revised logPAC is in the 'logPAC_new' component, and the Pearson residuals is in the 'pearson' component.

makeSCE2 <-
function(obj,
                     cData=NULL, batch = NULL, assays = c('pearson','logPAC', "logPAC_new")){
  # cData=meta
  # assays = c('pearson','logPAC')
  
  batch<-obj$batch
  if(is.null(cData))
    stop("'cData' must be a user-specified data frame.")
  if(any(! assays %in% c('logcounts','logPAC','pearson',"logPAC_new")))
    stop("'Wrong names of corrected assays data was specified. Available corrected data assays are logcounts,logPAC, new logPAC and pearson")
  
  # sce.obj <- SingleCellExperiment::SingleCellExperiment(assays = list(counts_m1 = as(obj$counts_modality1,"sparseMatrix"),
  #                                                                     counts_m2 = as(obj$counts_modality2,"sparseMatrix")), colData = cData)
  # Error in method(object) : all assays must have the same nrow and ncol [thinking to save as Seurat object instead]
  # CHECK HOW THEY DEAL WITH CITE-SEQ in one object.
  
  #for modality 1
  obj1 <- list(counts=obj$counts_modality1,
               a=obj$a1, a_joint=obj$a3_1,
               W=obj$W1, W_joint=obj$W3,
               gmean=obj$gmean[match(rownames(obj$counts_modality1), names(obj$gmean))],
               Mb=obj$Mb[match(rownames(obj$counts_modality1), rownames(obj$Mb)),],
               psi=obj$psi[match(rownames(obj$counts_modality1), names(obj$psi))])
  sce.obj.m1 <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as(obj1$counts,"sparseMatrix")), colData = cData)
  if(any(assays=='logcounts'))
    assays(sce.obj.m1, withDimnames=FALSE)$logcounts <- tryCatch({as( get.res2(obj1,type='logcounts',batch = batch),"sparseMatrix" )}, error = function(e) {NULL})
  if(any(assays=='pearson'))
    assays(sce.obj.m1, withDimnames=FALSE)$pearson <- tryCatch({as( get.res2(obj1,type='pearson',batch = batch), "sparseMatrix")}, error = function(e) {NULL})
  if(any(assays=='logPAC'))
    assays(sce.obj.m1, withDimnames=FALSE)$logPAC <- tryCatch({as(  log(get.res2(obj1,type='quantile',batch = batch)+1), "sparseMatrix")}, error = function(e) {NULL})
  if(any(assays=='logPAC_new'))
    assays(sce.obj.m1, withDimnames=FALSE)$logPAC_new <- tryCatch({as(  log(get.res2(obj1,type='quantile_spanorm',batch = batch)+1), "sparseMatrix")}, error = function(e) {NULL})
  
  #for modality 2
  obj2 <- list(counts=obj$counts_modality2,
               a=obj$a2, a_joint=obj$a3_2,
               W=obj$W2, W_joint=obj$W3,
               gmean=obj$gmean[match(rownames(obj$counts_modality2), names(obj$gmean))],
               Mb=obj$Mb[match(rownames(obj$counts_modality2), rownames(obj$Mb)),],
               psi=obj$psi[match(rownames(obj$counts_modality2), names(obj$psi))])
  sce.obj.m2 <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as(obj2$counts,"sparseMatrix")), colData = cData)
  if(any(assays=='logcounts'))
    assays(sce.obj.m2, withDimnames=FALSE)$logcounts <- tryCatch({as( get.res2(obj2,type='logcounts',batch = batch),"sparseMatrix" )}, error = function(e) {NULL})
  if(any(assays=='pearson'))
    assays(sce.obj.m2, withDimnames=FALSE)$pearson <- tryCatch({as( get.res2(obj2,type='pearson',batch = batch), "sparseMatrix")}, error = function(e) {NULL})
  if(any(assays=='logPAC'))
    assays(sce.obj.m2, withDimnames=FALSE)$logPAC <- tryCatch({as(  log(get.res2(obj2,type='quantile',batch = batch)+1), "sparseMatrix")}, error = function(e) {NULL})
  if(any(assays=='logPAC_new'))
    assays(sce.obj.m2, withDimnames=FALSE)$logPAC_new <- tryCatch({as(  log(get.res2(obj2,type='quantile_spanorm',batch = batch)+1), "sparseMatrix")}, error = function(e) {NULL})
  
  return(list(sce.obj.m1, sce.obj.m2))
}
