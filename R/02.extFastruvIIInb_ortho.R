#' Removing unwanted variation from matched sequencing count data (with orthogonality constraint)
#'
#' This function performs the extended fast version of ruvIIInb method (with orthogonality constraint) to remove unwanted variation from matched sequencing count data that contains features from two modalities. Currently, only Negative Binomial model for UMI data is supported. It takes raw count matrix with features (e.g., genes/transcript and proteins/ADT in CITE-seq) as rows and samples (cells) as columns. 
#' Users need to specify the set of negative controls (i.e features where the between cells variation is assumed to be solely due to the unwanted variation)
#' and the (pseudo)replicate matrix that define sets of cells (e.g., group by biology of interest) that can be considered technical replicates.
#'
#' @param Y1 raw count as DelayedMatrix object with features from modality one as rows and cells as columns
#' @param Y2 raw count as DelayedMatrix object with features from modality two as rows and cells as columns
#' @param M replicate matrix with number of rows equal to the number of cells and number of columns equal to the number of distinct sub-population of cells: M(i,j) = 1 if cell i is part of sub-population j and 0 otherwise. If a row has all zero entries, it means the corresponding cell is not assumed to belong to any of the known sub-populations (a priori unannotated cells). 
#' @param ctl1  either logical vector of length n with TRUE = controls OR a character vector of names of controls (control features for inferring modality one specific unwanted factors W1)
#' @param ctl2  either logical vector of length n with TRUE = controls OR a character vector of names of controls (control features for inferring modality two specific unwanted factors W2)
#' @param ctl3  either logical vector of length n with TRUE = controls OR a character vector of names of controls (control features for inferring joint unwanted factors W3)
#' @param k1 number of unwanted factors to be estimated for W1. Default is 1.
#' @param k2 number of unwanted factors to be estimated for W2. Default is 1.
#' @param k3 number of unwanted factors to be estimated for W3. Default is 1.
#' @param robust logical value. Whether to use Huber's robust weight in the iterative reweighted least squares (IRLS) algorithm. Default is FALSE.
#' @param ortho.W logical value. Whether the unwanted factors (W) are to be orthogonalized. The default is FALSE.
#' @param lambda.a1 smoothing parameter for regularizing regression coefficients associated with W1. The default is 0.01.
#' @param lambda.a2 smoothing parameter for regularizing regression coefficients associated with W2. The default is 0.01.
#' @param lambda.a3 smoothing parameter for regularizing regression coefficients associated with W3. The default is 0.01.
#' @param lambda.b smoothing parameter for regularizing the feature-level intercepts. The default is 16.
#' @param batch numeric vector containing batch information for each cell. Must correspond to columns of count matrix. If not supplied, cells are assumed to come from one batch.
#' @param step.fac multiplicative factor to decrease IRLS step by when log-likelihood diverges. Default is 0.5
#' @param inner.maxit the maximum number of IRLS iteration for estimating mean parameters for a given dispersion parameter. Default is 50
#' @param outer.maxit the maximum number of IRLS iteration. Default is 25 
#' @param ncores The number of cores used for parallel processing. The default number of workers (cores) is 2.
#' @param use.pseudosample whether to use pseudocells to define replicates (default is FALSE). Note that the replicates defined by the pseudocells will be used in addition to any replicates defined by the M matrix above. We recommend that use.pseudosample=TRUE be used only for data with UMI.
#' @param nc.pool the number of cells per pool (used for defining pseudocells). The default is 20. Only relevant when use.pseudosample=TRUE.
#' @param batch.disp whether to fit batch-specific dispersion parameter for each feature. The default is FALSE.
#' @param pCells.touse the proportion of a priori annotated cells used to estimate alpha and dispersion parameter (Default=20%). When pseudocells are used (use.pseudosample=TRUE), ultimately only 10% of the original subset of cells will be used to estimate alpha.
#' @param block.size the maximum number of cells for block processing when estimating W matrix. Larger block size can be quicker but requires higher RAM. Default = 5000.
#' 
#' @import matrixStats
#' @import foreach
#' @import doParallel
#' @importFrom DelayedArray DelayedArray
#' @importFrom ruvIIInb estimateDisp.par
#' @importFrom utils head
#' @importFrom stats cor dnbinom mad median pnbinom qnbinom quantile uniroot
#' @importFrom irlba irlba
#' @importFrom matlib GramSchmidt
#' 
#' @return A SingleCellExperiment object containing the raw and adjusted data in the assays slot, cell-level data including the estimated unwanted factors and the M matrix in the colData slot and the feature-level data including the coefficients associated with the unwanted factors in rowData slot. 
#' @export
extFastruvIIInb_ortho <-
function(Y1,Y2,M,
                                  ctl1,ctl2,ctl3,
                                  k1=1, k2=1, k3=1,
                                  robust=FALSE, ortho.W=FALSE,
                                  lambda.a1=0.01, lambda.a2=0.01, lambda.a3=0.01, lambda.b=16,
                                  batch=NULL, step.fac=0.5,
                                  inner.maxit=50, outer.maxit=25,
                                  ncores=2,use.pseudosample=FALSE,nc.pool=20,batch.disp=FALSE,pCells.touse=0.2,block.size=5000){
  START.time <- Sys.time()
  
  if(length(unique(c(rownames(Y1),rownames(Y2)))) != length(c(rownames(Y1),rownames(Y2)))){
    stop("Please make sure the feature names in your count matrics for both modalities are unique to each other...")
  }
  
  if(k1 != k2 | k1 != k3 | k2 != k3){
    stop("The orthogonality constraints require k1=k2=k3...")
  }
  
  #register parallel backend
  #register(BPPARAM)
  
  # setup cluster for doParallel
  doParallel::registerDoParallel(ncores)
  BPPARAM=BiocParallel::DoparParam()
  parallel <- as.logical(ncores>1)
  
  # control check
  if(!is.logical(ctl1)){
    ctl1 <- rownames(Y1) %in% ctl1
    names(ctl1) <- rownames(Y1)
  } 
  idx1 <- which(ctl1)
  nctl1<- sum(ctl1)
  #
  if(!is.logical(ctl2)){
    ctl2 <- rownames(Y2) %in% ctl2
    names(ctl2) <- rownames(Y2)
  } 
  idx2 <- which(ctl2)
  nctl2<- sum(ctl2)
  #
  if(!is.logical(ctl3)){
    ctl3 <- c(rownames(Y1),rownames(Y2)) %in% ctl3
    names(ctl3) <- c(rownames(Y1),rownames(Y2))
  } 
  idx3 <- which(ctl3)
  nctl3<- sum(ctl3)
  
  
  # check if any rows of M matrix has zero sum 
  zero.Mrows <- any(rowSums(M)==0)
  if(zero.Mrows) {
    n.unannot <- sum(rowSums(M)==0)
    print(paste0(n.unannot, ' cells have all zero entries in the M matrix. They will be considered as un-annotated cells'))
    if(n.unannot < (0.01*ncol(Y1))) 
      stop('Less than 1% of cells have known annotation. Please increase the pct of cells with known annotation to 1%')
  }
  
  # check if any cols of M matrix has zero sum 
  zero.Mcols <- any(colSums(M)==0)
  if(zero.Mcols) {
    stop(paste0('Columns ', which(colSums(M)==0),' of the M matrix has zero sum. Pls remove these columns and re-run'))
  }
  
  # force M matrix into logical matrix
  mode(M) <- 'logical'
  M       <- as.matrix(M)
  
  #format and convert M matrix into list of indices to save space
  rep.ind <- apply(M,2,which,simplify=FALSE)
  # define strata as columns of M matrix for which a cell is part of
  strata <- apply(M,1,which)
  
  # force Y1 and Y2 to DelayedMatrix object
  if(!any(class(Y1)=='DelayedMatrix'))
    Y1 <- DelayedArray::DelayedArray(Y1)
  if(!any(class(Y2)=='DelayedMatrix'))
    Y2 <- DelayedArray::DelayedArray(Y2)
  
  # Winsorize gene-by-gene and calculate size-factor
  bef <- Sys.time()
  max.val  <- ceiling(DelayedMatrixStats::rowQuantiles(Y1,prob=0.995))
  winsorize<- Y1>max.val
  Y1 <- Y1*(1-winsorize) + winsorize*max.val
  aft <- Sys.time()
  print(paste0('time to Winsorize RNA count matrix: ', difftime(aft,bef,units='secs')))
  
  # Winsorize protein-by-protein and calculate size-factor
  bef <- Sys.time()
  max.val  <- ceiling(DelayedMatrixStats::rowQuantiles(Y2,prob=0.995))
  winsorize<- Y2>max.val
  Y2 <- Y2*(1-winsorize) + winsorize*max.val
  aft <- Sys.time()
  print(paste0('time to Winsorize ADT count matrix: ', difftime(aft,bef,units='secs')))
  
  
  # remove genes with all zero counts (if any)
  zero.genes <- DelayedMatrixStats::rowSums2(Y1>0)==0
  names(zero.genes) <- rownames(Y1)
  if(any(zero.genes)) {
    Y1 <- Y1[!zero.genes,]
    ctl1 <- ctl1[!zero.genes]
    print(paste0(sum(zero.genes), ' genes with all zero counts are removed following a Winsorization step.'))
  }
  # identical(rownames(Y1), names(ctl1))
  
  # remove proteins with all zero counts (if any)
  zero.proteins <- DelayedMatrixStats::rowSums2(Y2>0)==0
  names(zero.proteins) <- rownames(Y2)
  if(any(zero.proteins)) {
    Y2 <- Y2[!zero.proteins,]
    ctl2 <- ctl2[!zero.proteins]
    print(paste0(sum(zero.proteins), ' proteins with all zero counts are removed following a Winsorization step.'))
  }
  # identical(rownames(Y2), names(ctl2))
  
  # remove genes and proteins with all zero counts (apply to control set #3)
  if(sum(zero.proteins) != 0 | sum(zero.genes) != 0){
    #if there's any zero proteins or zero genes
    ctl3 <- ctl3[-match(c(names(zero.genes)[zero.genes], names(zero.proteins)[zero.proteins]), names(ctl3))]
  }else{
    ctl3 <- ctl3
  }
  
  # combine Y1 and Y2
  Y12 <- rbind(Y1,Y2) #for est. psi
  # identical(rownames(Y12), names(ctl3))
  
  # if no batch info is supplied but batch.disp=TRUE, forced batch.disp=FALSE
  if(is.null(batch) & batch.disp) {
    print(paste0('Batch variable not supplied...cannot fit batch-specific dispersion parameter.'))
    batch.disp <- FALSE
  }
  
  # if no batch info is supplied, assumes cells come from one batch (batch=NULL - why need that)
  if(is.null(batch)) {
    print(paste0('Batch variable not supplied...assuming cells come from one batch'))
    batch <- rep(1,ncol(Y12)) 
  }
  
  # if batch is not numeric vector, stop and let user knows
  if(!is.null(batch) & !is.numeric(batch)) {
    stop(paste0('Batch variable is not a numeric vector. Please supply batch variable as numeric vector with 1,2...B values'))
  }
  
  
  # get Huber's k. (robust the estimate, not robustfy.. large values: check ruviiinb supp. - robust weight)
  # leave this for now... (4 Feb 2024)
  k.huber<- ifelse(robust,1.345,100)
  
  
  # adjust pCells.touse to allow max 3000 cells (minimise for est. alpha)
  max.pCells <- 3000/sum(rowSums(M)>0)
  if(max.pCells < pCells.touse) {
    pCells.touse <- max.pCells
    print(paste0('pCells.touse parameter is too large and has now been set to ', round(pCells.touse ,6)))
  }
  # select a subset of cells representative of sub-pops defined by M matrix (stratified random sampling)
  set.seed(123) #to make sure likelihood among different sets of NC_ADT are comparable
  subsamples <- sample(rep.ind[[1]],size=ceiling(pCells.touse*length(rep.ind[[1]])))
  if(length(rep.ind)>1) {
    for(i in 2:length(rep.ind)) {
      subsamples <- c(subsamples,sample(rep.ind[[i]],size=ceiling(pCells.touse*length(rep.ind[[i]]))))
    }
  }
  subsamples <- sort(subsamples)
  #
  Y1sub  <- as.matrix(Y1[,subsamples])
  Y2sub  <- as.matrix(Y2[,subsamples])
  Y12sub  <- as.matrix(Y12[,subsamples])
  nsub <- ncol(Y12sub) #should be the same as ncol(Y1sub) and ncol(Y2sub)
  
  #
  Msub <- as.matrix(M[subsamples,])
  rep.sub <- apply(Msub,2,which,simplify=FALSE)
  sub.batch <- batch[subsamples]
  if(!is.null(strata)) 
    strata <- strata[subsamples]
  
  #borrow from below pseudosamples code
  nbatch.org <- ifelse(is.null(sub.batch),1,length(unique(sub.batch)))
  nbatch <- nbatch.org
  subsamples.org     <- subsamples
  subsubsamples.org  <- 1:nsub
  
  # adding pseudosamples (we should make modality specific i.e., add sf_rna and sf_adt for each cell, instead of having only one value; equal weights for both RNA and ADT when forming pseudocells)
  # leave this for now... (4 Feb 2024)
  {
    # sf <- colSums(Y12sub) #input: features x cells; length(sf)=cells
    # sf <- sf/mean(sf) 
    # sample.type <- rep('sc',ncol(Y12sub))
    # # number of batches before pseudocells addition
    # nbatch.org <- ifelse(is.null(sub.batch),1,length(unique(sub.batch)))
    # # if batch.disp=FALSE,force nbatch.org==1
    # nbatch.org <- ifelse(!batch.disp,1,nbatch.org)
    # 
    # # number of single cells before pseudocells addition [have not modified 31 Jan]
    # nsub <- ncol(Ysub)
    # if(use.pseudosample) {
    #   if(is.null(strata)) {
    #     batch.ext <- sub.batch
    #     for(btc in unique(sub.batch)) {
    #       nq  <- max(round(sum(batch==btc)/nc.pool),1)
    #       qsf <-  as.numeric(gtools::quantcut(sf[sub.batch==btc],q=nq))
    #       for(q in 1:max(qsf)) {
    #         pool <- Ysub[,1:nsub][,sub.batch==btc][,qsf==q]
    #         dns.pool <- rbinom(nrow(pool),size=rowSums(pool),prob=1/ncol(pool))
    #         Ysub <- cbind(Ysub,dns.pool)
    #         sample.type <- c(sample.type,'sc')
    #         batch.ext  <- c(batch.ext,btc+max(sub.batch))
    #         Msub      <- rbind(Msub,FALSE)
    #       }
    #     }
    #     ps.mcol <- rep(FALSE,nrow(Msub))
    #     ps.mcol[seq(nsub+1,ncol(Ysub),1)] <- TRUE
    #     Msub <- cbind(Msub,ps.mcol)
    #     sub.batch <- batch.ext
    #   }
    #   
    #   if(!is.null(strata)) {
    #     batch.ext <- sub.batch
    #     for(st in unique(na.omit(strata))) {
    #       for(btc in unique(sub.batch)) {
    #         nq  <- max(round(sum(sub.batch==btc & strata==st,na.rm=TRUE)/nc.pool),1)
    #         if(sum(sub.batch==btc & strata==st, na.rm=TRUE)>1) {
    #           qsf <-  na.omit(as.numeric(gtools::quantcut(sf[sub.batch==btc & strata==st],q=nq)))
    #           for(q in 1:max(qsf)) {
    #             pool <- Ysub[,1:nsub][,sub.batch==btc & strata==st & !is.na(strata)][,qsf==q]
    #             dns.pool <- rbinom(nrow(pool),size=rowSums(pool),prob=1/ncol(pool))
    #             Ysub <- cbind(Ysub,dns.pool)
    #             sample.type <- c(sample.type,'sc')
    #             batch.ext  <- c(batch.ext,btc+max(sub.batch))
    #             ps.mrow <- rep(FALSE,ncol(Msub)) 
    #             ps.mrow[st] <- TRUE
    #             Msub    <- rbind(Msub,ps.mrow)
    #           }
    #         } #if
    #       } #for batch
    #     } # for st
    #     if(!batch.disp) 
    #       batch.ext <- ifelse(batch.ext>max(sub.batch),2,1)
    #     
    #     sub.batch <- batch.ext
    #     
    #   } # if
    # }
    # nbatch <- ifelse(!use.pseudosample,nbatch.org,length(unique(sub.batch)))
    # 
    # subsamples.org     <- subsamples
    # subsubsamples.org  <- 1:nsub
    # if(use.pseudosample) {
    #   # here, we take only 10% of the original subsampled cells + all pseudocells
    #   subsamples.org     <- subsamples
    #   subsubsamples.org  <- sort(sample(nsub,size=round(0.1*nsub)))
    #   subsamples         <- c(subsubsamples.org,(nsub+1):ncol(Ysub))
    #   
    #   Ysub <- Ysub[,subsamples]
    #   Msub <- as.matrix(Msub[subsamples,])
    #   rep.sub <- apply(Msub,2,which,simplify=FALSE)
    #   sub.batch <- sub.batch[subsamples]
    #   nsub <- ncol(Ysub)
    # }
  }
  
  # select a maximum of 100 cells per batch for estimating psi [no need to modify]
  psi.idx   <- NULL
  for(B in sort(unique(sub.batch))) 
    psi.idx  <- c(psi.idx,sample(which(sub.batch==B),size=min(100,round(0.5*sum(sub.batch==B)))))
  psi.batch  <- sub.batch[psi.idx]
  nsub.psi   <- length(psi.idx)
  # Y12sub.psi   <- Y12sub[,psi.idx]
  Y1sub.psi   <- Y1sub[,psi.idx]
  Y2sub.psi   <- Y2sub[,psi.idx]
  
  # feature specific stuff... beta, psi, intercept 
  # alpha (use subset of cells - then use alpha from subset of cells to est W, Wsub for inner algorthim... only need W for all cells and alpha need to be fixed - the one from sub cells), W
  # block size = 5000... update W for 5000 cells each time (to save time... each time, don't have that huge matrix)
  # delayarray: don't eat ram...each time only bring stuff you need to ram (to speed up)
  ### batch: biological batch ###
  ### block: computational batch ###
  
  # initial estimates (starting) of alpha1, alpha2, and W1, W2
  {
    #alpha 1
    lmu.hat1<- Z1 <- log(Y1sub+1) 
    gmean1  <- rowMeans(Z1)
    proj.Z1 <- projection(rep.sub,Zmat=Z1)[,apply(Msub,1,which)]
    RM_Z1   <- Z1 - proj.Z1
    U01     <- irlba::irlba(RM_Z1,nv=k1)$v
    # alpha(n x k) matrix
    alpha1  <- as.matrix(Z1 %*% U01 %*% Matrix::chol2inv(chol(Matrix::crossprod(U01)+lambda.a1*diag(k1)))) #435   5
    
    #alpha 2
    lmu.hat2<- Z2 <- log(Y2sub+1) 
    gmean2  <- rowMeans(Z2)
    proj.Z2 <- projection(rep.sub,Zmat=Z2)[,apply(Msub,1,which)]
    RM_Z2   <- Z2 - proj.Z2
    U02     <- irlba::irlba(RM_Z2,nv=k2)$v
    # alpha(n x k) matrix
    alpha2  <- as.matrix(Z2 %*% U02 %*% Matrix::chol2inv(chol(Matrix::crossprod(U02)+lambda.a2*diag(k2)))) #97  5
    
    # reduce outliers
    #alpha1
    a1.med <- matrixStats::colMedians(alpha1)
    a1.mad <- matrixStats::colMads(alpha1)
    for(i in 1:ncol(alpha1)) {
      alpha1[which(alpha1[,i]> (a1.med[i]+4*a1.mad[i])),i] <- a1.med[i]+4*a1.mad[i]
      alpha1[which(alpha1[,i]< (a1.med[i]-4*a1.mad[i])),i] <- a1.med[i]-4*a1.mad[i]
    }
    alpha1.c <- as.matrix(alpha1[ctl1,])
    Z1.c    <- Z1[ctl1,] - gmean1[ctl1]
    
    #alpha2
    a2.med <- matrixStats::colMedians(alpha2)
    a2.mad <- matrixStats::colMads(alpha2)
    for(i in 1:ncol(alpha2)) {
      alpha2[which(alpha2[,i]> (a2.med[i]+4*a2.mad[i])),i] <- a2.med[i]+4*a2.mad[i]
      alpha2[which(alpha2[,i]< (a2.med[i]-4*a2.mad[i])),i] <- a2.med[i]-4*a2.mad[i]
    }
    alpha2.c <- as.matrix(alpha2[ctl2,])
    Z2.c    <- Z2[ctl2,] - gmean2[ctl2]
    
    # initiate W (m x k)
    #Wsub <- Matrix::crossprod(Z.c,alpha.c) %*% solve(Matrix::crossprod(alpha.c)) 
    #W1
    W1sub      <- matrix(0,nsub,k1)
    W1sub[,1]  <- Matrix::crossprod(Z1.c,as.matrix(alpha1.c[,1])) %*% solve(lambda.a1 + Matrix::crossprod(as.matrix(alpha1.c[,1]),as.matrix(alpha1.c[,1]))) 
    if(k1>1) {
      for(j in 2:k1) {
        Z1.c <- Z1.c - outer(alpha1.c[,j-1],W1sub[,j-1])
        W1sub[,j] <- Matrix::crossprod(Z1.c,as.matrix(alpha1.c[,j])) %*% solve(lambda.a1 + Matrix::crossprod(as.matrix(alpha1.c[,j]),as.matrix(alpha1.c[,j]))) 
      }
    }
    #W2
    W2sub      <- matrix(0,nsub,k2)
    W2sub[,1]  <- Matrix::crossprod(Z2.c,as.matrix(alpha2.c[,1])) %*% solve(lambda.a2 + Matrix::crossprod(as.matrix(alpha2.c[,1]),as.matrix(alpha2.c[,1]))) 
    if(k2>1) {
      for(j in 2:k2) {
        Z2.c <- Z2.c - outer(alpha2.c[,j-1],W2sub[,j-1])
        W2sub[,j] <- Matrix::crossprod(Z2.c,as.matrix(alpha2.c[,j])) %*% solve(lambda.a2 + Matrix::crossprod(as.matrix(alpha2.c[,j]),as.matrix(alpha2.c[,j]))) 
      }
    }
    
    # reduce outliers
    #W1
    w1.med <- matrixStats::colMedians(W1sub)
    w1.mad <- matrixStats::colMads(W1sub)
    for(i in 1:ncol(W1sub)) {
      W1sub[which(W1sub[,i]> (w1.med[i]+4*w1.mad[i])),i] <- w1.med[i]+4*w1.mad[i]
      W1sub[which(W1sub[,i]< (w1.med[i]-4*w1.mad[i])),i] <- w1.med[i]-4*w1.mad[i]
    }
    #W2
    w2.med <- matrixStats::colMedians(W2sub)
    w2.mad <- matrixStats::colMads(W2sub)
    for(i in 1:ncol(W2sub)) {
      W2sub[which(W2sub[,i]> (w2.med[i]+4*w2.mad[i])),i] <- w2.med[i]+4*w2.mad[i]
      W2sub[which(W2sub[,i]< (w2.med[i]-4*w2.mad[i])),i] <- w2.med[i]-4*w2.mad[i]
    }
  }
  
  # initial estimates of alpha3 and W3
  {
    # take out effect of W1a1 and W2a2 from the data...effects  that need to be taken out prior to calculating a3
    tmp1 <- log(Y1sub+1) - Matrix::tcrossprod(alpha1,W1sub)
    tmp2 <- log(Y2sub+1) - Matrix::tcrossprod(alpha2,W2sub)
    #alpha 3 
    #current: est. from concat. Y1 and Y2; 
    #Z3 is the leftover of W1a1 and W2a2 as we have subtract the effects
    lmu.hat3<- Z3 <- rbind(tmp1, tmp2)
    gmean3  <- rowMeans(Z3) #532
    proj.Z3 <- projection(rep.sub,Zmat=Z3)[,apply(Msub,1,which)]
    RM_Z3   <- Z3 - proj.Z3 #532 1966
    U03     <- irlba::irlba(RM_Z3,nv=k3)$v #1966    5 [combo of features (rows): contain shared variation across features from different modalities?]
    # alpha(n x k) matrix
    alpha3  <- as.matrix(Z3 %*% U03 %*% Matrix::chol2inv(chol(Matrix::crossprod(U03)+lambda.a3*diag(k3)))) #532   5 [Z3 rows are different features, U03 is a mix of all features]
    
    # reduce outliers
    #alpha3
    a3.med <- matrixStats::colMedians(alpha3)
    a3.mad <- matrixStats::colMads(alpha3)
    for(i in 1:ncol(alpha3)) {
      alpha3[which(alpha3[,i]> (a3.med[i]+4*a3.mad[i])),i] <- a3.med[i]+4*a3.mad[i]
      alpha3[which(alpha3[,i]< (a3.med[i]-4*a3.mad[i])),i] <- a3.med[i]-4*a3.mad[i]
    }
    
    #somehow makes sense to sep this into 435+97 (?) - [ASK]
    alpha3_1 <- as.matrix(alpha3[match(rownames(Y1), rownames(alpha3)),])
    alpha3_2 <- as.matrix(alpha3[match(rownames(Y2), rownames(alpha3)),])
    
    alpha3.c <- as.matrix(alpha3[ctl3,]) #identical(names(ctl3), rownames(alpha3)) #T
    
    Z3.c    <- Z3[ctl3,] - gmean3[ctl3]
    
    # initiate W (m x k)
    #W3
    W3sub      <- matrix(0,nsub,k3)
    W3sub[,1]  <- Matrix::crossprod(Z3.c,as.matrix(alpha3.c[,1])) %*% solve(lambda.a3 + Matrix::crossprod(as.matrix(alpha3.c[,1]),as.matrix(alpha3.c[,1]))) 
    if(k3>1) {
      for(j in 2:k3) {
        Z3.c <- Z3.c - outer(alpha3.c[,j-1],W3sub[,j-1])
        W3sub[,j] <- Matrix::crossprod(Z3.c,as.matrix(alpha3.c[,j])) %*% solve(lambda.a3 + Matrix::crossprod(as.matrix(alpha3.c[,j]),as.matrix(alpha3.c[,j]))) 
      }
    }
    
    # reduce outliers
    #W3
    w3.med <- matrixStats::colMedians(W3sub)
    w3.mad <- matrixStats::colMads(W3sub)
    for(i in 1:ncol(W3sub)) {
      W3sub[which(W3sub[,i]> (w3.med[i]+4*w3.mad[i])),i] <- w3.med[i]+4*w3.mad[i]
      W3sub[which(W3sub[,i]< (w3.med[i]-4*w3.mad[i])),i] <- w3.med[i]-4*w3.mad[i]
    }
  }
  
  
  # initiate Mb and gmean (?): by performing poisson GLM feature-by-feature
  nW1   <- ncol(W1sub)
  nW2   <- ncol(W2sub)
  nW3   <- ncol(W3sub)
  nM   <- ncol(M)
  nm1<- nrow(Y1) #m1=modality 1
  nm2<- nrow(Y2) #m2=modality 2
  ns   <- ncol(Y1) #ncol(Y1) = ncol(Y2)
  
  bef   <- Sys.time()
  # update Mb
  W1a1     <- Matrix::tcrossprod(alpha1,W1sub) #435 1966
  W2a2     <- Matrix::tcrossprod(alpha2,W2sub) #97 1966
  W3a3     <- Matrix::tcrossprod(alpha3,W3sub) #532 1966
  
  Z.res  <- Z3 - rbind(W1a1, W2a2) - W3a3  - gmean3 #dimensions: (532 1966) - (532 1966) - (532 1966) - (532)
  Mb  <- tryCatch({ projection(rep.sub, Zmat=Z.res, Wmat=NULL, lambda=lambda.b)}, error = function(e) {NA}) #532  13
  # identical(rownames(Mb), names(ctl3)) #T
  Mb[union(union(rownames(Y1)[ctl1], rownames(Y2)[ctl2]), rownames(Y12)[ctl3]),] <- 0 #be cautious when ctl3 != concat(ctl1, ctl2) [need to revise soon]
  
  # estimate initial psi
  bef <- Sys.time() 
  #offsets
  offs <- rbind(W1a1, W2a2) + W3a3 + Mb[,apply(Msub,1,which)] 
  offs.psi <- offs[,psi.idx]
  bef  <- Sys.time()
  if(parallel){
    if(nbatch==1){
      # psi <- tryCatch({ estimateDisp.par(Y12sub.psi, as.matrix(rep(1,nsub.psi)), offset=offs.psi, tagwise=TRUE, robust=TRUE, BPPARAM=BPPARAM)$tagwise.dispersion},
      #                 error = function(e) {rep(NA,(nm1+nm2))})
      
      #rna
      psi1 <- tryCatch({ estimateDisp.par(Y1sub.psi, as.matrix(rep(1,nsub.psi)), offset=offs.psi[1:nm1,], tagwise=TRUE, robust=TRUE, BPPARAM=BPPARAM)$tagwise.dispersion},
                       error = function(e) {rep(NA,(nm1))})
      
      #adt
      psi2 <- tryCatch({ estimateDisp.par(Y2sub.psi, as.matrix(rep(1,nsub.psi)), offset=offs.psi[(nm1+1):(nm1+nm2),], tagwise=TRUE, robust=TRUE, BPPARAM=BPPARAM)$tagwise.dispersion},
                       error = function(e) {rep(NA,(nm2))})
      
      psi <- c(psi1, psi2)
    }
    if(nbatch>1){
      # psi <- tryCatch({ estimateDisp.par(Y12sub.psi[,psi.batch==1],as.matrix(rep(1,sum(psi.batch==1))),
      #                                    offset=offs.psi[,psi.batch==1],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion}, 
      #                 error = function(e) {rep(NA,(nm1+nm2))})
      # 
      # for(B in 2:max(psi.batch)) {
      #   psi <- cbind(psi,tryCatch({estimateDisp.par(Y12sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
      #                                               offset=offs.psi[,psi.batch==B],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion}, 
      #                             error = function(e) {rep(NA,(nm1+nm2))}))
      # }
      
      #rna
      psi1 <- tryCatch({ estimateDisp.par(Y1sub.psi[,psi.batch==1],as.matrix(rep(1,sum(psi.batch==1))),
                                          offset=offs.psi[1:nm1, psi.batch==1],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion},
                       error = function(e) {rep(NA,(nm1))})
      #adt
      psi2 <- tryCatch({ estimateDisp.par(Y2sub.psi[,psi.batch==1],as.matrix(rep(1,sum(psi.batch==1))),
                                          offset=offs.psi[(nm1+1):(nm1+nm2), psi.batch==1],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion},
                       error = function(e) {rep(NA,(nm2))})
      psi <- c(psi1, psi2)
      
      for(B in 2:max(psi.batch)){
        #rna
        psi1 <- tryCatch({estimateDisp.par(Y1sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
                                           offset=offs.psi[1:nm1, psi.batch==B],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion},
                         error = function(e) {rep(NA,(nm1))})
        #adt
        psi2 <- tryCatch({estimateDisp.par(Y2sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
                                           offset=offs.psi[(nm1+1):(nm1+nm2), psi.batch==B],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion},
                         error = function(e) {rep(NA,(nm2))})
        psi <- cbind(psi, c(psi1, psi2))
      }
    }
  } 
  if(!parallel){
    if(nbatch==1){
      # psi <- tryCatch({edgeR::estimateDisp(Y12sub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,tagwise=TRUE,robust=TRUE)$tagwise.dispersion}, 
      #                 error = function(e) {rep(NA,(nm1+nm2))})
      
      #rna
      psi1 <- tryCatch({edgeR::estimateDisp(Y1sub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi[1:nm1,],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},
                       error = function(e) {rep(NA,(nm1))})
      #adt
      psi2 <- tryCatch({edgeR::estimateDisp(Y2sub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi[(nm1+1):(nm1+nm2),],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},
                       error = function(e) {rep(NA,(nm2))})
      psi <- c(psi1, psi2)
    }
    if(nbatch>1){
      # psi <- tryCatch({edgeR::estimateDisp(Y12sub.psi[,psi.batch==1],as.matrix(rep(1,sum(psi.batch==1))),
      #                                      offset=offs.psi[,psi.batch==1],tagwise=TRUE,robust=TRUE)$tagwise.dispersion}, 
      #                 error = function(e) {rep(NA,(nm1+nm2))})
      # for(B in 2:max(sub.batch)) 
      #   psi <- cbind(psi,tryCatch({edgeR::estimateDisp(Y12sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
      #                                                  offset=offs.psi[,psi.batch==B],tagwise=TRUE,robust=TRUE)$tagwise.dispersion}, 
      #                             error = function(e) {rep(NA,(nm1+nm2))}))
      
      #rna
      psi1 <- tryCatch({edgeR::estimateDisp(Y1sub.psi[,psi.batch==1],as.matrix(rep(1,sum(psi.batch==1))),
                                            offset=offs.psi[1:nm1, psi.batch==1],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},
                       error = function(e) {rep(NA,(nm1))})
      #adt
      psi2 <- tryCatch({edgeR::estimateDisp(Y2sub.psi[,psi.batch==1],as.matrix(rep(1,sum(psi.batch==1))),
                                            offset=offs.psi[(nm1+1):(nm1+nm2), psi.batch==1],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},
                       error = function(e) {rep(NA,(nm2))})
      
      psi <- c(psi1, psi2)
      
      for(B in 2:max(sub.batch)){
        psi1 <- tryCatch({edgeR::estimateDisp(Y1sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
                                              offset=offs.psi[1:nm1, psi.batch==B],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},
                         error = function(e) {rep(NA,(nm1))})
        psi2 <- tryCatch({edgeR::estimateDisp(Y2sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
                                              offset=offs.psi[(nm1+1):(nm1+nm2), psi.batch==B],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},
                         error = function(e) {rep(NA,(nm2))})
        
        psi <- cbind(psi, c(psi1, psi2))
        
      }
    }
  }
  aft <- Sys.time()
  #print(paste0('time to calculate psi:',  difftime(aft,bef,units='secs')))
  psi[which(is.na(psi))] <- 0
  
  # initiate
  loglikehood.ls = trace_sub.w13.ls = trace_sub.w23.ls <- list()
  tau1_sub.ls = tau2_sub.ls <- list()
  conv <- FALSE
  iter.outer <- 0
  logl.outer <- NULL
  step <- rep(1,nsub)
  print('Start...')
  while(!conv) {
    iter.outer <- iter.outer + 1
    logl.beta <- NULL
    
    # lmu.hat    <- gmean3 + Mb[,apply(Msub,1,which)] +  rbind(Matrix::tcrossprod(alpha1,W1sub), Matrix::tcrossprod(alpha2,W2sub)) + Matrix::tcrossprod(alpha3,W3sub)
    lmu.hat    <- gmean3 + Mb[,apply(Msub,1,which)] +  rbind(Matrix::tcrossprod(alpha1,W1sub), Matrix::tcrossprod(alpha2,W2sub)) + 
      rbind(Matrix::tcrossprod(alpha3_1,W3sub), Matrix::tcrossprod(alpha3_2,W3sub)) #work for both initial and updated
    
    # weight based on NB (this part is for 1/var(e)... which is the feature-cell-specific weight)
    if(nbatch>1)
      sig.inv<- 1/(exp(-lmu.hat) + sweep(psi[,sub.batch],2,rep(1,nsub),'/')) #532 1966
    if(nbatch==1)
      sig.inv<- 1/(exp(-lmu.hat) + outer(psi,rep(1,nsub),'/'))
    wt.ctl1 <- rowMeans(sig.inv[match(rownames(Y1)[ctl1], rownames(sig.inv)),]) #weights for ctl1 
    wt.ctl2 <- rowMeans(sig.inv[match(rownames(Y2)[ctl2], rownames(sig.inv)),]) #weights for ctl2 
    wt.ctl3 <- rowMeans(sig.inv[match(rownames(Y12)[ctl3], rownames(sig.inv)),]) #weights for ctl3 
    
    # forming "tmp" log-likelihood 
    if(nbatch>1) {
      temp <- foreach(i=1:nsub, .combine=c) %dopar% {
        tmp <- dnbinom(Y12sub[,i],mu=exp(lmu.hat[,i]),size=1/psi[,sub.batch[i]],log=TRUE)
        tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
        return(sum(tmp))
      }
    } 
    if(nbatch==1) {
      temp <- foreach(i=1:nsub, .combine=c) %dopar% {
        tmp <- dnbinom(Y12sub[,i],mu=exp(lmu.hat[,i]),size=1/psi,log=TRUE)
        tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
        return(sum(tmp))
      }
    } 
    loglik <- temp
    
    ############################# 
    step <- rep(1,nsub)
    conv.beta <- FALSE
    # saving temporary best estimate
    best.psi <- psi ; best.Mb <- Mb
    best.gmean3 <- gmean3 # equal to gmean in RUVIIINB
    best.a1 <- alpha1 ; best.a2 <- alpha2 ; best.a3_1 <- alpha3_1 ; best.a3_2 <- alpha3_2 #upper for modality one and lower for modality two
    best.W1 <- W1sub ; best.W2 <- W2sub ; best.W3 <- W3sub
    best.wtctl1 <- wt.ctl1 ; best.wtctl2 <- wt.ctl2 ; best.wtctl3 <- wt.ctl3
    
    logl.beta.inner <- c()
    trace_sub.w13 = trace_sub.w23 <- c()
    tau1_sub = tau2_sub <- c()
    iter <- 1
    halving <- 0
    while(!conv.beta) {
      print(paste0('Inner iter:', iter))
      # calc working vector Z for a subset of cells 
      lmu.hat <- gmean3 + Mb[,apply(Msub,1,which)] + rbind(Matrix::tcrossprod(alpha1,W1sub), Matrix::tcrossprod(alpha2,W2sub)) + Matrix::tcrossprod(rbind(alpha3_1,alpha3_2),W3sub) 
      
      # weight based on NB
      if(nbatch>1)
        sig.inv<- 1/(exp(-lmu.hat) + sweep(psi[,sub.batch],2,rep(1,nsub),'/'))
      if(nbatch==1)
        sig.inv<- 1/(exp(-lmu.hat) +outer(psi,rep(1,nsub),'/'))
      
      #  # calculate signed deviance [robust: not modified yet... 5 Feb 2024]
      #  if(nbatch>1 & robust) {
      #    signdev <- sign(Ysub-exp(lmu.hat))* sqrt(2*(dnbinom(Ysub,mu=Ysub,size=sweep(1/psi[,sub.batch],2,rep(1,nsub),'*'),log=TRUE) - 
      # 						dnbinom(Ysub,mu=exp(lmu.hat),size=sweep(1/psi[,sub.batch],2,rep(1,nsub),'*'),log=TRUE) ))
      #  }
      #  if(nbatch==1 & robust) {
      #    signdev <- sign(Ysub-exp(lmu.hat))* sqrt(2*(dnbinom(Ysub,mu=Ysub,size=outer(1/psi,rep(1,nsub),'*'),log=TRUE) - 
      # 						dnbinom(Ysub,mu=exp(lmu.hat),size=outer(1/psi,rep(1,nsub),'*'),log=TRUE) ))
      #  }  
      
      Z      <- lmu.hat + sweep( ((Y12sub+0.01)/(exp(lmu.hat)+0.01) - 1),2,step,'*')
      
      # project Z and W into M
      
      #### step1: obtain projecion of Z and W onto the range of M and get the residual ####
      # RM_Z is (n x m) matrix, m=nsample,n=ngene ((nm1; nm2 here))
      # RM_W is (m x k) matrix
      proj.Z <- projection(rep.sub,Zmat=Z)[,apply(Msub,1,which)]
      RM_Z <- Z - proj.Z
      
      ## W1
      if(nW1>1) 
        proj.tW1 <- projection(rep.sub,Zmat=t(W1sub))
      if(nW1==1) 
        proj.tW1 <- t(as.matrix(projection.K1(rep.sub,Zmat=t(W1sub))))
      
      if(nW1>1) 
        RM_W1 <- W1sub - t(proj.tW1[,apply(Msub,1,which)])
      if(nW1==1) 
        RM_W1 <- W1sub - as.matrix(proj.tW1[,apply(Msub,1,which)])
      
      ## W2
      if(nW2>1) 
        proj.tW2 <- projection(rep.sub,Zmat=t(W2sub))
      if(nW2==1) 
        proj.tW2 <- t(as.matrix(projection.K1(rep.sub,Zmat=t(W2sub))))
      
      if(nW2>1) 
        RM_W2 <- W2sub - t(proj.tW2[,apply(Msub,1,which)])
      if(nW2==1) 
        RM_W2 <- W2sub - as.matrix(proj.tW2[,apply(Msub,1,which)])
      
      ## W3
      if(nW3>1) 
        proj.tW3 <- projection(rep.sub,Zmat=t(W3sub))
      if(nW3==1) 
        proj.tW3 <- t(as.matrix(projection.K1(rep.sub,Zmat=t(W3sub))))
      
      if(nW3>1) 
        RM_W3 <- W3sub - t(proj.tW3[,apply(Msub,1,which)])
      if(nW3==1) 
        RM_W3 <- W3sub - as.matrix(proj.tW3[,apply(Msub,1,which)])
      
      
      # store current alpha est
      alpha1.old <- alpha1
      alpha2.old <- alpha2
      alpha3_1.old <- alpha3_1 #alpha3 for modality 1
      alpha3_2.old <- alpha3_2 #alpha3 for modality 2
      
      # get alpha 
      bef=Sys.time()
      
      ### for modality 1's all features ###
      extr_m1 <- match(rownames(Y1), rownames(sig.inv))
      wt.cell_m1  <- matrixStats::colMeans2(sig.inv[extr_m1,]) #cell weights that are the average across features from modality 1 
      # prevent outlier large weight
      p98      <- quantile(wt.cell_m1,probs=0.98)
      wt.cell_m1[which(wt.cell_m1>=p98)] <- p98
      #wt.cell_m1  <- rep(1,ncol(sig.inv))
      # identical(rownames(RM_Z), rownames(sig.inv)) #T
      #alpha1 (fix other paras and update)
      b11        <- (RM_Z[extr_m1,] - t(RM_W3 %*% t(alpha3_1.old))) %*% diag(wt.cell_m1) %*% RM_W1
      alpha1    <- as.matrix(b11 %*% Matrix::chol2inv(chol(Matrix::crossprod(RM_W1*wt.cell_m1,RM_W1)+lambda.a1*diag(nW1))))
      #modality specific alpha3 (fix other paras and update)
      b13        <- (RM_Z[extr_m1,] - t(RM_W1 %*% t(alpha1.old))) %*% diag(wt.cell_m1) %*% RM_W3
      alpha3_1    <- as.matrix(b13 %*% Matrix::chol2inv(chol(Matrix::crossprod(RM_W3*wt.cell_m1,RM_W3)+lambda.a3*diag(nW3))))
      
      ### for modality 2's all features ###
      extr_m2 <- match(rownames(Y2), rownames(sig.inv))
      wt.cell_m2  <- matrixStats::colMeans2(sig.inv[extr_m2,]) #cell weights that are the average across features from modality 2 
      # prevent outlier large weight
      p98      <- quantile(wt.cell_m2,probs=0.98)
      wt.cell_m2[which(wt.cell_m2>=p98)] <- p98
      #wt.cell_m2  <- rep(1,ncol(sig.inv))
      # identical(rownames(RM_Z), rownames(sig.inv)) #T
      #alpha2 (fix other paras and update)
      b21        <- (RM_Z[extr_m2,] - t(RM_W3 %*% t(alpha3_2.old))) %*% diag(wt.cell_m2) %*% RM_W2
      alpha2    <- as.matrix(b21 %*% Matrix::chol2inv(chol(Matrix::crossprod(RM_W2*wt.cell_m2,RM_W2)+lambda.a2*diag(nW2))))
      #modality specific alpha3 (fix other paras and update)
      b23        <- (RM_Z[extr_m2,] - t(RM_W2 %*% t(alpha2.old))) %*% diag(wt.cell_m2) %*% RM_W3
      alpha3_2    <- as.matrix(b23 %*% Matrix::chol2inv(chol(Matrix::crossprod(RM_W3*wt.cell_m2,RM_W3)+lambda.a3*diag(nW3))))
      
      aft <- Sys.time()
      #print(dim(alpha))
      #print(paste0('time to calculate alpha:',  difftime(aft,bef,units='secs')))
      
      # if new alpha est has missing/inf values, revert to previous estimate
      #print(paste0('Problematic alpha entries=',sum(is.na(alpha) | is.infinite(alpha))))
      if(any(is.na(alpha1) | is.infinite(alpha1) )) alpha1 <- alpha1.old
      if(any(is.na(alpha2) | is.infinite(alpha2) )) alpha2 <- alpha2.old
      if(any(is.na(alpha3_1) | is.infinite(alpha3_1) )) alpha3_1 <- alpha3_1.old #alpha 3 for modality 1
      if(any(is.na(alpha3_2) | is.infinite(alpha3_2) )) alpha3_2 <- alpha3_2.old #alpha 3 for modality 2
      
      # reduce outliers
      #alpha1
      a1.med <- matrixStats::colMedians(alpha1)
      a1.mad <- matrixStats::colMads(alpha1)
      for(i in 1:nW1) {
        alpha1[which(alpha1[,i]> (a1.med[i]+4*a1.mad[i])),i] <- a1.med[i]+4*a1.mad[i]
        alpha1[which(alpha1[,i]< (a1.med[i]-4*a1.mad[i])),i] <- a1.med[i]-4*a1.mad[i]
      } 
      #alpha2
      a2.med <- matrixStats::colMedians(alpha2)
      a2.mad <- matrixStats::colMads(alpha2)
      for(i in 1:nW2) {
        alpha2[which(alpha2[,i]> (a2.med[i]+4*a2.mad[i])),i] <- a2.med[i]+4*a2.mad[i]
        alpha2[which(alpha2[,i]< (a2.med[i]-4*a2.mad[i])),i] <- a2.med[i]-4*a2.mad[i]
      } 
      #alpha3
      a31.med <- matrixStats::colMedians(alpha3_1)
      a31.mad <- matrixStats::colMads(alpha3_1)
      for(i in 1:nW3) {
        alpha3_1[which(alpha3_1[,i]> (a31.med[i]+4*a31.mad[i])),i] <- a31.med[i]+4*a31.mad[i]
        alpha3_1[which(alpha3_1[,i]< (a31.med[i]-4*a31.mad[i])),i] <- a31.med[i]-4*a31.mad[i]
      } 
      a32.med <- matrixStats::colMedians(alpha3_2)
      a32.mad <- matrixStats::colMads(alpha3_2)
      for(i in 1:nW3) {
        alpha3_2[which(alpha3_2[,i]> (a32.med[i]+4*a32.mad[i])),i] <- a32.med[i]+4*a32.mad[i]
        alpha3_2[which(alpha3_2[,i]< (a32.med[i]-4*a32.mad[i])),i] <- a32.med[i]-4*a32.mad[i]
      } 
      
      
      #### step2: update W using WLS of Z_c (Z for control genes) on alpha.c (alpha for control genes) ####
      W1sub.old <- W1sub 
      W2sub.old <- W2sub 
      W3sub.old <- W3sub 
      
      bef <- Sys.time()
      #below w/o modification
      {
        #W   <- BiocParallel::bpmapply(FUN=getW,Z=Z.list,signdev=signdev.bycol,wt=wtlist,zeroinf=zeroinf.list,
        #			MoreArgs=list(alpha=alpha,ctl=ctl,k.huber=k.huber,k=k),BPPARAM=BPPARAM)
        
        #Wsub <- foreach(i=1:nsub, .combine=cbind, .packages="Matrix") %dopar% { 
        # out <- rep(0,k)
        # wtvec     <- sig.inv[ctl,i] 
        # wtvec     <- wtvec/mean(wtvec)
        # out[1]      <- Rfast::lmfit(y=Z[ctl,i]-gmean[ctl],x=alpha[ctl,1],w=wtvec)$be
        # if(k>1) {
        #   for(j in 2:k) {
        #     out[j] <- Rfast::lmfit(y=Z[ctl,i]-gmean[ctl]-matrixStats::rowSums2(as.matrix(alpha[ctl,1:(j-1)]) %*% as.matrix(out[1:(j-1)])),x=alpha[ctl,j],w=wtvec)$be
        #   } 
        # }
        # as.matrix(out)
        #}
      }
      
      #W1 (using nc1)
      extr_ctl1 <- match(rownames(Y1)[ctl1], rownames(Z))
      Z1.c     <- Z[extr_ctl1,] - gmean3[extr_ctl1] - Matrix::tcrossprod(alpha3_1[match(rownames(Y1)[ctl1], rownames(alpha3_1)),],W3sub.old)
      alpha1.c <- alpha1[ctl1,,drop=FALSE]
      wt.ctl1 <- rowMeans(sig.inv[match(rownames(Y1)[ctl1], rownames(sig.inv)),])
      
      W1sub[,1]  <- Matrix::crossprod(Z1.c,as.matrix(alpha1.c[,1]*wt.ctl1)) %*% solve(lambda.a1 + Matrix::crossprod(as.matrix(alpha1.c[,1]*wt.ctl1),as.matrix(alpha1.c[,1]))) 
      if(k1>1) {
        for(j in 2:k1) {
          Z1.c <- Z1.c - outer(alpha1.c[,j-1],W1sub[,j-1])
          W1sub[,j] <- Matrix::crossprod(Z1.c,as.matrix(alpha1.c[,j]*wt.ctl1)) %*% solve(lambda.a1 + Matrix::crossprod(as.matrix(alpha1.c[,j]*wt.ctl1),as.matrix(alpha1.c[,j]))) 
        }
      }
      #W2 (using nc2)
      extr_ctl2 <- match(rownames(Y2)[ctl2], rownames(Z))
      Z2.c     <- Z[extr_ctl2,] - gmean3[extr_ctl2] - Matrix::tcrossprod(alpha3_2[match(rownames(Y2)[ctl2], rownames(alpha3_2)),],W3sub.old)
      alpha2.c <- alpha2[ctl2,,drop=FALSE]
      wt.ctl2 <- rowMeans(sig.inv[match(rownames(Y2)[ctl2], rownames(sig.inv)),])
      
      W2sub[,1]  <- Matrix::crossprod(Z2.c,as.matrix(alpha2.c[,1]*wt.ctl2)) %*% solve(lambda.a2 + Matrix::crossprod(as.matrix(alpha2.c[,1]*wt.ctl2),as.matrix(alpha2.c[,1]))) 
      if(k2>1) {
        for(j in 2:k2) {
          Z2.c <- Z2.c - outer(alpha2.c[,j-1],W2sub[,j-1])
          W2sub[,j] <- Matrix::crossprod(Z2.c,as.matrix(alpha2.c[,j]*wt.ctl2)) %*% solve(lambda.a2 + Matrix::crossprod(as.matrix(alpha2.c[,j]*wt.ctl2),as.matrix(alpha2.c[,j]))) 
        }
      }
      #W3 (using nc3)
      k <- k3 #remember in this version we need k1=k2=k3
      # extr_ctl3 <- match(rownames(Y12)[ctl3], rownames(Z))
      # Z3.c     <- Z[extr_ctl3,] - gmean3[extr_ctl3] - rbind(Matrix::tcrossprod(alpha1[which(rownames(alpha1) %in% rownames(Y12)[ctl3]),], W1sub.old),
      #                                                       Matrix::tcrossprod(alpha2[which(rownames(alpha2) %in% rownames(Y12)[ctl3]),], W2sub.old))
      # alpha3.c <- rbind(alpha3_1, alpha3_2)[match(rownames(Y12)[ctl3], c(rownames(alpha3_1), rownames(alpha3_2))),,drop=FALSE]
      # wt.ctl3 <- rowMeans(sig.inv[match(rownames(Y12)[ctl3], rownames(sig.inv)),])
      # 
      # W3sub[,1]  <- Matrix::crossprod(Z3.c,as.matrix(alpha3.c[,1]*wt.ctl3)) %*% solve(lambda.a3 + Matrix::crossprod(as.matrix(alpha3.c[,1]*wt.ctl3),as.matrix(alpha3.c[,1]))) 
      # if(k3>1) {
      #   for(j in 2:k3) {
      #     Z3.c <- Z3.c - outer(alpha3.c[,j-1],W3sub[,j-1])
      #     W3sub[,j] <- Matrix::crossprod(Z3.c,as.matrix(alpha3.c[,j]*wt.ctl3)) %*% solve(lambda.a3 + Matrix::crossprod(as.matrix(alpha3.c[,j]*wt.ctl3),as.matrix(alpha3.c[,j]))) 
      #   }
      # }
      
      extr_ctl3 <- match(rownames(Y12)[ctl3], rownames(Z))
      Z3.c     <- Z[extr_ctl3,] - gmean3[extr_ctl3] - rbind(Matrix::tcrossprod(alpha1[which(rownames(alpha1) %in% rownames(Y12)[ctl3]),], W1sub), #use W1sub instead? (since it's freshly updated above)
                                                            Matrix::tcrossprod(alpha2[which(rownames(alpha2) %in% rownames(Y12)[ctl3]),], W2sub)) #use W2sub instead? (since it's freshly updated above)
      alpha3.c <- rbind(alpha3_1, alpha3_2)[match(rownames(Y12)[ctl3], c(rownames(alpha3_1), rownames(alpha3_2))),,drop=FALSE]
      wt.ctl3 <- rowMeans(sig.inv[match(rownames(Y12)[ctl3], rownames(sig.inv)),])
      
      ##update all cells, all k at once
      {
        ##terms would be repetitively used
        z.w.a <- Matrix::crossprod(Z3.c, as.matrix(alpha3.c*wt.ctl3))
        a.w.a <- Matrix::crossprod(as.matrix(alpha3.c*wt.ctl3), as.matrix(alpha3.c))
        cp.w1sub <- crossprod(W1sub, W1sub) #crossprod using W1sub freshly updated above
        cp.w2sub <- crossprod(W2sub, W2sub) #crossprod using W2sub freshly updated above
        cp.w3sub.old <- crossprod(W3sub.old, W3sub.old) #need to use the old one
        
        ##to form the scalar below ( can speed up by applying crossprod(x,y) = t(x) %*% y )
        s.op.front <- t(as.matrix(rep(1,k))) #scalar operator front
        s.op.back <- as.matrix(rep(1,k)) #scalar operator back
        S1 <- s.op.front %*% cp.w3sub.old %*% cp.w1sub %*% solve(cp.w2sub) %*% solve(cp.w3sub.old) %*% s.op.back
        S2 <- s.op.front %*% t(W3sub.old) %*% z.w.a %*% solve(cp.w2sub) %*% solve(cp.w3sub.old) %*% s.op.back
        S3 <- s.op.front %*% cp.w3sub.old %*% a.w.a %*% solve(cp.w2sub) %*% solve(cp.w3sub.old) %*% s.op.back
        
        inv.w3w2 <- solve(cp.w2sub) %*% solve(cp.w3sub.old)
        tau1 <- (2*S1 - 2*s.op.front %*% cp.w3sub.old %*% cp.w1sub %*% inv.w3w2 %*% s.op.back + 0.000001)^(-1)*(
          - (s.op.front %*% t(W3sub.old) %*% z.w.a %*% inv.w3w2 %*% s.op.back)
          + (s.op.front %*% cp.w3sub.old %*% a.w.a %*% inv.w3w2 %*% s.op.back)
          + (S2 - S3)
        )
        tau2 <- (2*k)^(-1)*(s.op.front %*% t(W3sub.old) %*% z.w.a %*% inv.w3w2 %*% s.op.back)
        - (2*k)^(-1)*(s.op.front %*% cp.w3sub.old %*% a.w.a %*% inv.w3w2 %*% s.op.back)
        - k^(-1)*tau1*(s.op.front %*% cp.w3sub.old %*% cp.w1sub %*% inv.w3w2 %*% s.op.back)
        
        S4 <- tau1
        S5 <- (1/k)*s.op.front %*% t(W3sub.old) %*% z.w.a %*% inv.w3w2 %*% s.op.back
        S6 <- (1/k)*s.op.front %*% cp.w3sub.old %*% a.w.a %*% inv.w3w2 %*% s.op.back
        S7 <- (2*tau1/k)*s.op.front %*% cp.w3sub.old %*% cp.w1sub %*% inv.w3w2 %*% s.op.back
        
        additional.terms <- (2*as.numeric(S4)*cp.w1sub + as.numeric(S5 - S6 - S7)*cp.w2sub)
      }
      W3sub  <- z.w.a %*% solve(lambda.a3 + a.w.a + additional.terms) 
      
      tau1_sub[iter] <- as.numeric(tau1)
      tau2_sub[iter] <- as.numeric(tau2)
      trace_sub.w13[iter] <- sum(diag(W3sub %*% crossprod(W1sub, W1sub) %*% t(W3sub)))
      trace_sub.w23[iter] <- sum(diag(W3sub %*% crossprod(W2sub, W2sub) %*% t(W3sub)))
      
      aft <- Sys.time()
      #print(paste0('time to calculate W:',  difftime(aft,bef,units='secs')))
      W1sub <- as.matrix(W1sub)
      if(ncol(W1sub)!=k1) 
        W1sub <- t(W1sub)
      W2sub <- as.matrix(W2sub)
      if(ncol(W2sub)!=k2) 
        W2sub <- t(W2sub)
      W3sub <- as.matrix(W3sub)
      if(ncol(W3sub)!=k3) 
        W3sub <- t(W3sub)
      
      colNorm.W1sub <- matrixStats::colSums2(W1sub^2)
      colNorm.W2sub <- matrixStats::colSums2(W2sub^2)
      colNorm.W3sub <- matrixStats::colSums2(W3sub^2)
      # normalize W to have unit length for each col (no modification)
      #Wsub  <- sweep(Wsub,2,sqrt(colNorm.Wsub),'/')
      
      # orthogonalize W (optional)
      if(nW1>1 & ortho.W)  
        W1sub <- matlib::GramSchmidt(W1sub)
      if(nW2>1 & ortho.W)  
        W2sub <- matlib::GramSchmidt(W2sub)
      if(nW3>1 & ortho.W)  
        W3sub <- matlib::GramSchmidt(W3sub)
      
      # if new W has missing/inf values, revert to previous estimate
      #print(paste0('Problematic W entries=',sum(is.na(Wsub) | is.infinite(Wsub))))
      if(any(is.na(W1sub) | is.infinite(W1sub) )) W1sub <- W1sub.old #different from RUVIIINB code line 469
      if(any(is.na(W2sub) | is.infinite(W2sub) )) W2sub <- W2sub.old
      if(any(is.na(W3sub) | is.infinite(W3sub) )) W3sub <- W3sub.old
      
      #### step3a: update gmean #### 
      Z.res  <- Z - Mb[,apply(Msub,1,which)] - rbind(Matrix::tcrossprod(alpha1,W1sub), Matrix::tcrossprod(alpha2,W2sub)) - Matrix::tcrossprod(rbind(alpha3_1,alpha3_2),W3sub) 
      bef <- Sys.time()
      gmean3<- matrixStats::rowSums2(Z.res*sig.inv)/matrixStats::rowSums2(sig.inv)
      aft <- Sys.time()
      #print(paste0('time to calculate gmean:',  difftime(aft,bef,units='secs')))
      
      #### step 3b: update Mb ####
      Z.res  <- Z - gmean3 - rbind(Matrix::tcrossprod(alpha1,W1sub), Matrix::tcrossprod(alpha2,W2sub)) - Matrix::tcrossprod(rbind(alpha3_1,alpha3_2),W3sub) 
      Mb.old <- Mb
      bef <- Sys.time()
      Mb  <- tryCatch({ projection(rep.sub,Zmat=Z.res,Wmat=sig.inv,lambda=lambda.b)}, error = function(e) {NA})
      # set Mb=0 for control genes
      Mb[union(union(rownames(Y1)[ctl1], rownames(Y2)[ctl2]), rownames(Y12)[ctl3]),] <- 0 
      aft <- Sys.time()
      #print(paste0('time to calculate Mb:',  difftime(aft,bef,units='secs')))
      
      # if new Mb has missing/inf values, revert to previous estimate
      #print(paste0('Problematic Mb entries=',sum(is.na(Mb) | is.infinite(Mb))))
      if(any(is.na(Mb))) Mb <- Mb.old
      Mb[which(is.na(Mb) | is.infinite(Mb))] <- 0
      
      # reduce outliers
      Mb.med <- matrixStats::rowMedians(Mb)
      Mb.mad <- matrixStats::rowMads(Mb)
      for(i in 1:(nm1+nm2)) {
        Mb[i,which(Mb[i,]> (Mb.med[i]+4*Mb.mad[i]))] <- Mb.med[i]+4*Mb.mad[i]
        Mb[i,which(Mb[i,]< (Mb.med[i]-4*Mb.mad[i]))] <- Mb.med[i]-4*Mb.mad[i]
      }
      
      # calculate current logl
      lmu.hat <- gmean3 + Mb[,apply(Msub,1,which)] + rbind(Matrix::tcrossprod(alpha1,W1sub), Matrix::tcrossprod(alpha2,W2sub)) + Matrix::tcrossprod(rbind(alpha3_1,alpha3_2),W3sub) 
      
      if(nbatch>1) {
        temp <- foreach(i=1:nsub, .combine=c) %dopar% {
          tmp <- dnbinom(Y12sub[,i],mu=exp(lmu.hat[,i]),size=1/psi[,sub.batch[i]],log=TRUE)
          tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
          return(sum(tmp))
        }
      } 
      if(nbatch==1) {
        temp <- foreach(i=1:nsub, .combine=c) %dopar% {
          tmp <- dnbinom(Y12sub[,i],mu=exp(lmu.hat[,i]),size=1/psi,log=TRUE)
          tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
          return(sum(tmp))
        }
      } 
      loglik.tmp <- temp 
      
      # check degenerate case
      if(iter>=1) {
        degener<- sum(loglik)>sum(loglik.tmp)
        degener[is.na(degener)] <- TRUE
        if(degener) {
          check.gene <- loglik>loglik.tmp
          step[check.gene] <- step[check.gene]*step.fac
          #
          psi <- best.psi ; Mb <- best.Mb
          gmean3 <- best.gmean3 # equal to gmean in RUVIIINB
          alpha1 <- best.a1 ; alpha2 <- best.a2 ; alpha3_1 <- best.a3_1 ; alpha3_2 <- best.a3_2 #upper for modality one and lower for modality two
          W1sub <- best.W1 ; W2sub <- best.W2 ; W3sub <- best.W3
          wt.ctl1 <- best.wtctl1 ; wt.ctl2 <- best.wtctl2 ; wt.ctl3 <- best.wtctl3
          #
          halving <- halving + 1
          if(halving>=3) {
            loglik.tmp <- loglik
            degener <- !degener
          }
        }
        if(!degener) {
          loglik     <- loglik.tmp
          logl.beta  <- c(logl.beta,sum(loglik))
          
          logl.beta.inner[iter] <- logl.beta[length(logl.beta)]
          print(paste0('Outer Iter ',iter.outer, ', Inner iter ', iter, ' logl-likelihood:', logl.beta[length(logl.beta)]))
          # saving temporary best estimate
          #
          best.psi <- psi ; best.Mb <- Mb
          best.gmean3 <- gmean3 # equal to gmean in RUVIIINB
          best.a1 <- alpha1 ; best.a2 <- alpha2 ; best.a3_1 <- alpha3_1 ; best.a3_2 <- alpha3_2 #upper for modality one and lower for modality two
          best.W1 <- W1sub ; best.W2 <- W2sub ; best.W3 <- W3sub
          best.wtctl1 <- wt.ctl1 ; best.wtctl2 <- wt.ctl2 ; best.wtctl3 <- wt.ctl3
          #
          iter <- iter + 1
          halving <- 0
        }
      }
      
      conv.logl <- FALSE
      if(iter>2) {
        conv.logl <- ifelse( (logl.beta[length(logl.beta)] - logl.beta[length(logl.beta)-1])/abs(logl.beta[length(logl.beta)-1]) < 1e-04,TRUE,FALSE)
      }
      conv.beta <- iter>=inner.maxit  | conv.logl
    } # end of IRLS inner loop
    
    trace_sub.w13.ls[[iter.outer]] <- trace_sub.w13 ##
    trace_sub.w23.ls[[iter.outer]] <- trace_sub.w23 ##
    tau1_sub.ls[[iter.outer]] <- tau1_sub
    tau2_sub.ls[[iter.outer]] <- tau2_sub
    
    loglikehood.ls[[iter.outer]] <- logl.beta.inner ##
    
    #update psi
    logl.outer.tmp <- sum(loglik,na.rm=T)
    updt.psi <- TRUE
    if(iter.outer>1)
      updt.psi <- abs(logl.outer[length(logl.outer)]-logl.outer.tmp)/abs(logl.outer[length(logl.outer)]) > 1e-08  
    
    W1a1 <- Matrix::tcrossprod(best.a1,best.W1) 
    W2a2 <- Matrix::tcrossprod(best.a2,best.W2) 
    W3a3 <- rbind(Matrix::tcrossprod(best.a3_1,best.W3), Matrix::tcrossprod(best.a3_2,best.W3))
    offs <- rbind(W1a1, W2a2) + W3a3 + best.Mb[,apply(Msub,1,which)] 
    offs.psi <- offs[,psi.idx]
    if(updt.psi){ #needing to update psi
      if(parallel){
        if(nbatch==1){
          # psi.new <- tryCatch({ estimateDisp.par(Y12sub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,
          #                                        tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion}, error = function(e) {rep(NA,(nm1+nm2))})
          
          #rna
          psi.new1 <- tryCatch({ estimateDisp.par(Y1sub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi[1:nm1,],
                                                  tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion}, error = function(e) {rep(NA,(nm1))})
          
          #adt
          psi.new2 <- tryCatch({ estimateDisp.par(Y2sub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi[(nm1+1):(nm1+nm2),],
                                                  tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion}, error = function(e) {rep(NA,(nm2))})
          
          psi.new <- c(psi.new1, psi.new2)
        }
        if(nbatch>1){
          psi.new <- psi
          for(B in 1:max(psi.batch)){
            # psi.new[,B]<-tryCatch({estimateDisp.par(Y12sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
            #                                         offset=offs.psi[,psi.batch==B],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion},error= function(e) {rep(NA,(nm1+nm2))})
            
            #rna
            psi.new.bat1 <- tryCatch({estimateDisp.par(Y1sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
                                                       offset=offs.psi[1:nm1, psi.batch==B],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion},error= function(e) {rep(NA,(nm1))})
            #adt
            psi.new.bat2 <- tryCatch({estimateDisp.par(Y2sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
                                                       offset=offs.psi[(nm1+1):(nm1+nm2), psi.batch==B],tagwise=TRUE,robust=TRUE,BPPARAM=BPPARAM)$tagwise.dispersion},error= function(e) {rep(NA,(nm2))})
            
            psi.new.bat <- c(psi.new.bat1, psi.new.bat2)
            
            psi.new[,B] <- psi.new.bat
          } 
        }
      }
      
      if(!parallel){
        if(nbatch==1){
          # psi.new <- tryCatch({ edgeR::estimateDisp(Y12sub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,tagwise=TRUE,robust=TRUE)$tagwise.dispersion},error = function(e) {rep(NA,(nm1+nm2))})
          
          #rna
          psi.new1 <- tryCatch({ edgeR::estimateDisp(Y1sub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi[1:nm1,],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},error = function(e) {rep(NA,(nm1))})
          
          #adt
          psi.new2 <- tryCatch({ edgeR::estimateDisp(Y2sub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi[(nm1+1):(nm1+nm2),],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},error = function(e) {rep(NA,(nm2))})
          
          psi.new <- c(psi.new1, psi.new2)
        }
        if(nbatch>1){
          psi.new <- psi
          for(B in 1:max(psi.batch)){
            # psi.new[,B] <- tryCatch({ edgeR::estimateDisp(Y12sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
            #                                               offset=offs.psi[,psi.batch==B],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},
            #                         error = function(e) {rep(NA,(nm1+nm2))})
            
            #rna
            psi.new.bat1 <- tryCatch({ edgeR::estimateDisp(Y1sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
                                                           offset=offs.psi[1:nm1, psi.batch==B],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},
                                     error = function(e) {rep(NA,(nm1))})
            #adt
            psi.new.bat2 <- tryCatch({ edgeR::estimateDisp(Y2sub.psi[,psi.batch==B],as.matrix(rep(1,sum(psi.batch==B))),
                                                           offset=offs.psi[(nm1+1):(nm1+nm2), psi.batch==B],tagwise=TRUE,robust=TRUE)$tagwise.dispersion},
                                     error = function(e) {rep(NA,(nm2))})
            
            psi.new.bat <- c(psi.new.bat1, psi.new.bat2)
            
            psi.new[,B] <- psi.new.bat
            
          } 
        }
      }
      # update psi
      psi[!is.na(psi.new) & !is.infinite(psi.new)] <- psi.new[!is.na(psi.new) & !is.infinite(psi.new)]
    }
    
    # record outer logl
    W1sub <- best.W1
    W2sub <- best.W2
    W3sub <- best.W3
    alpha1 <- best.a1
    alpha2 <- best.a2
    alpha3_1 <- best.a3_1
    alpha3_2 <- best.a3_2
    Mb    <- best.Mb
    gmean3 <- best.gmean3
    # new changes
    # best.psi <- psi
    
    # recalculate logl using updated psi
    lmu.hat <- gmean3 + Mb[,apply(Msub,1,which)] + rbind(Matrix::tcrossprod(alpha1,W1sub), Matrix::tcrossprod(alpha2,W2sub)) + Matrix::tcrossprod(rbind(alpha3_1,alpha3_2),W3sub) 
    if(nbatch>1) {
      temp <- foreach(i=1:nsub, .combine=c) %dopar% {
        tmp <- dnbinom(Y12sub[,i],mu=exp(lmu.hat[,i]),size=1/psi[,sub.batch[i]],log=TRUE)
        tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
        return(sum(tmp))
      }
    } 
    if(nbatch==1) {
      temp <- foreach(i=1:nsub, .combine=c) %dopar% {
        tmp <- dnbinom(Y12sub[,i],mu=exp(lmu.hat[,i]),size=1/psi,log=TRUE)
        tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
        return(sum(tmp))
      }
    } 
    #not update logl.outer here!
    # logl.outer <- c(logl.outer, sum(temp)) #-lambda.a*sum(scale(alpha,scale=FALSE)^2)-lambda.b*sum(Mb^2))
    print(paste0('If Updating psi in Outer Iter ',iter.outer, '; will obtain logl-likelihood:', sum(temp)))
    
    # check degenerate case only for psi
    if(iter.outer>=1){
      if(iter.outer == 1){
        degener<- logl.beta.inner[length(logl.beta.inner)]>sum(temp) #the former one equal to logl.beta[length(logl.beta)]
      }else{
        degener<- logl.outer[length(logl.outer)]>sum(temp) #the former one equal to logl.beta[length(logl.beta)]
      }
      degener[is.na(degener)] <- TRUE
      
      if(!degener){ #new is better
        #using recalculated logl above
        logl.outer <- c(logl.outer, sum(temp)) #-lambda.a*sum(scale(alpha,scale=FALSE)^2)-lambda.b*sum(Mb^2))
        
        best.psi <- psi
        
        print(paste0('Updating psi in Outer Iter ',iter.outer, '; logl-likelihood:', logl.outer[length(logl.outer)]))
      }
      
      
      if(degener){ #psi from the previous iteration is better
        # recalculate logl using previous psi
        psi <- best.psi #previous psi is better
        
        lmu.hat <- gmean3 + Mb[,apply(Msub,1,which)] + rbind(Matrix::tcrossprod(alpha1,W1sub), Matrix::tcrossprod(alpha2,W2sub)) + Matrix::tcrossprod(rbind(alpha3_1,alpha3_2),W3sub) 
        if(nbatch>1) {
          temp <- foreach(i=1:nsub, .combine=c) %dopar% {
            tmp <- dnbinom(Y12sub[,i],mu=exp(lmu.hat[,i]),size=1/psi[,sub.batch[i]],log=TRUE)
            tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
            return(sum(tmp))
          }
        } 
        if(nbatch==1) {
          temp <- foreach(i=1:nsub, .combine=c) %dopar% {
            tmp <- dnbinom(Y12sub[,i],mu=exp(lmu.hat[,i]),size=1/psi,log=TRUE)
            tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
            return(sum(tmp))
          }
        } 
        logl.outer <- c(logl.outer, sum(temp)) #-lambda.a*sum(scale(alpha,scale=FALSE)^2)-lambda.b*sum(Mb^2))
        print(paste0('NOT Updating psi in Outer Iter ',iter.outer, '; logl-likelihood:', logl.outer[length(logl.outer)]))
        #iter.outer does not +1 as we are "back" to the previous iteration
      }
      
      if(iter.outer == 1){ 
        #not update psi at all (default psi is better) so only compare with default psi + converged parameters from inner loop
        conv <- ifelse( (logl.outer[iter.outer] - logl.beta.inner[length(logl.beta.inner)])/abs(logl.beta.inner[length(logl.beta.inner)]) < 1e-04 | iter.outer>=outer.maxit,TRUE,FALSE)
      }else{ #if(iter.outer>1)
        conv <- ifelse( (logl.outer[iter.outer] - logl.outer[iter.outer-1])/abs(logl.outer[iter.outer-1]) < 1e-04 | iter.outer>=outer.maxit,TRUE,FALSE)
      }
    }
    
    if(conv) {
      alpha1 <- best.a1
      alpha2 <- best.a2
      alpha3_1 <- best.a3_1
      alpha3_2 <- best.a3_2
      gmean3 <- best.gmean3
      Mb  <- best.Mb
      psi   <- best.psi
      wt.ctl1 <- best.wtctl1
      wt.ctl2 <- best.wtctl2
      wt.ctl3 <- best.wtctl3
      alpha1.c <- alpha1[ctl1,,drop=FALSE]
      alpha2.c <- alpha2[ctl2,,drop=FALSE]
      alpha3.c <- rbind(alpha3_1, alpha3_2)[match(rownames(Y12)[ctl3], c(rownames(alpha3_1), rownames(alpha3_2))),,drop=FALSE]
      
      ###
      {
        lmu.hat <- gmean3 + Mb[,apply(Msub,1,which)] + rbind(Matrix::tcrossprod(alpha1,W1sub), Matrix::tcrossprod(alpha2,W2sub)) + Matrix::tcrossprod(rbind(alpha3_1,alpha3_2),W3sub) 
        if(nbatch>1) {
          temp <- foreach(i=1:nsub, .combine=c) %dopar% {
            tmp <- dnbinom(Y12sub[,i],mu=exp(lmu.hat[,i]),size=1/psi[,sub.batch[i]],log=TRUE)
            tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
            return(sum(tmp))
          }
        } 
        if(nbatch==1) {
          temp <- foreach(i=1:nsub, .combine=c) %dopar% {
            tmp <- dnbinom(Y12sub[,i],mu=exp(lmu.hat[,i]),size=1/psi,log=TRUE)
            tmp[is.infinite(tmp) | is.na(tmp)] <- log(.Machine$double.xmin)
            return(sum(tmp))
          }
        }
        print(paste0('Both loops converged; checking the logl-likelihood again: ', sum(temp))) 
      }
      ###
      
      # estimate W for all samples (initial)
      print('Estimating W for all samples...')
      bef=Sys.time()
      #W <- irlba::irlba(Y[ctl,],nv=k)$v
      #W1, W2, and W3 updating in each block [6 Feb - check untitled4 and editing]
      W1 <- matrix(0,ncol(Y1),k1)
      W2 <- matrix(0,ncol(Y2),k2)
      W3 <- matrix(0,ncol(Y12),k3)
      W1[subsamples.org[subsubsamples.org],] <- W1sub[1:length(subsubsamples.org),] #insert estimates from Wsub, remain others as 0 for following estimation
      W2[subsamples.org[subsubsamples.org],] <- W2sub[1:length(subsubsamples.org),]
      W3[subsamples.org[subsubsamples.org],] <- W3sub[1:length(subsubsamples.org),]
      # block size variable
      block.size <- min(block.size,ncol(Y12))
      nb    <- ceiling(ncol(Y12)/block.size)
      # alpha1.c <- alpha1[ctl1,,drop=FALSE]
      # alpha2.c <- alpha2[ctl2,,drop=FALSE]
      # alpha3.c <- rbind(alpha3_1, alpha3_2)[match(rownames(Y12)[ctl3], c(rownames(alpha3_1), rownames(alpha3_2))),,drop=FALSE]
      for(block in 1:nb) {
        start.idx <- (block-1)*block.size+1 ; end.idx <- min(ns,block*block.size)
        
        ## W1
        Y1sub  <- as.matrix(Y1[ctl1,start.idx:end.idx])
        #below rethink
        lmu.hat <- gmean3[match(rownames(Y1)[ctl1], names(gmean3))] + Matrix::tcrossprod(alpha1.c,W1[start.idx:end.idx,,drop=FALSE]) + Matrix::tcrossprod(alpha3_1[match(rownames(Y1)[ctl1], rownames(alpha3_1)),],W3[start.idx:end.idx,,drop=FALSE])
        Z1.c  <- lmu.hat - gmean3[match(rownames(Y1)[ctl1], names(gmean3))] - Matrix::tcrossprod(alpha3_1[match(rownames(Y1)[ctl1], rownames(alpha3_1)),],W3[start.idx:end.idx,,drop=FALSE]) + (Y1sub+0.01)/(exp(lmu.hat)+0.01) - 1
        W1[start.idx:end.idx,1]<-Matrix::crossprod(Z1.c,as.matrix(alpha1.c[,1]*wt.ctl1)) %*% solve(lambda.a1 + Matrix::crossprod(as.matrix(alpha1.c[,1]*wt.ctl1),as.matrix(alpha1.c[,1])))
        if(k1>1) {
          for(j in 2:k1) {
            Z1.c <- Z1.c - outer(alpha1.c[,j-1],W1[start.idx:end.idx,j-1])
            W1[start.idx:end.idx,j] <- Matrix::crossprod(Z1.c,as.matrix(alpha1.c[,j]*wt.ctl1)) %*% solve(lambda.a1 + Matrix::crossprod(as.matrix(alpha1.c[,j]*wt.ctl1),as.matrix(alpha1.c[,j])))
          }
        }
        
        ## W2
        Y2sub  <- as.matrix(Y2[ctl2,start.idx:end.idx])
        #below rethink
        lmu.hat <- gmean3[match(rownames(Y2)[ctl2], names(gmean3))] + Matrix::tcrossprod(alpha2.c,W2[start.idx:end.idx,,drop=FALSE]) + Matrix::tcrossprod(alpha3_2[match(rownames(Y2)[ctl2], rownames(alpha3_2)),],W3[start.idx:end.idx,,drop=FALSE])
        Z2.c  <- lmu.hat - gmean3[match(rownames(Y2)[ctl2], names(gmean3))] - Matrix::tcrossprod(alpha3_2[match(rownames(Y2)[ctl2], rownames(alpha3_2)),],W3[start.idx:end.idx,,drop=FALSE]) + (Y2sub+0.01)/(exp(lmu.hat)+0.01) - 1
        W2[start.idx:end.idx,1]<-Matrix::crossprod(Z2.c,as.matrix(alpha2.c[,1]*wt.ctl2)) %*% solve(lambda.a2 + Matrix::crossprod(as.matrix(alpha2.c[,1]*wt.ctl2),as.matrix(alpha2.c[,1])))
        if(k2>1) {
          for(j in 2:k2) {
            Z2.c <- Z2.c - outer(alpha2.c[,j-1],W2[start.idx:end.idx,j-1])
            W2[start.idx:end.idx,j] <- Matrix::crossprod(Z2.c,as.matrix(alpha2.c[,j]*wt.ctl2)) %*% solve(lambda.a2 + Matrix::crossprod(as.matrix(alpha2.c[,j]*wt.ctl2),as.matrix(alpha2.c[,j])))
          }
        }
        
        ## W3
        Y12sub  <- as.matrix(Y12[ctl3,start.idx:end.idx]) #check if ctl3 is really formed from Y12 [CHECK]
        #below rethink
        lmu.hat <- gmean3[match(rownames(Y12)[ctl3], names(gmean3))] +
          rbind(Matrix::tcrossprod(alpha1[which(rownames(alpha1) %in% rownames(Y12)[ctl3]),], W1[start.idx:end.idx,,drop=FALSE]), #alpha1 for NC3
                Matrix::tcrossprod(alpha2[which(rownames(alpha2) %in% rownames(Y12)[ctl3]),], W2[start.idx:end.idx,,drop=FALSE])) + #alpha2 for NC3
          Matrix::tcrossprod(alpha3.c,W3[start.idx:end.idx,,drop=FALSE])
        
        Z3.c  <- lmu.hat - gmean3[match(rownames(Y12)[ctl3], names(gmean3))] - rbind(Matrix::tcrossprod(alpha1[which(rownames(alpha1) %in% rownames(Y12)[ctl3]),], W1[start.idx:end.idx,,drop=FALSE]),
                                                                                     Matrix::tcrossprod(alpha2[which(rownames(alpha2) %in% rownames(Y12)[ctl3]),], W2[start.idx:end.idx,,drop=FALSE])) +
          (Y12sub+0.01)/(exp(lmu.hat)+0.01) - 1
        
        # W3[start.idx:end.idx,1] <- Matrix::crossprod(Z3.c,as.matrix(alpha3.c[,1]*wt.ctl3)) %*% solve(lambda.a3 + Matrix::crossprod(as.matrix(alpha3.c[,1]*wt.ctl3),as.matrix(alpha3.c[,1])))
        # if(k3>1) {
        #   for(j in 2:k3) {
        #     Z3.c <- Z3.c - outer(alpha3.c[,j-1],W3[start.idx:end.idx,j-1])
        #     W3[start.idx:end.idx,j] <- Matrix::crossprod(Z3.c,as.matrix(alpha3.c[,j]*wt.ctl3)) %*% solve(lambda.a3 + Matrix::crossprod(as.matrix(alpha3.c[,j]*wt.ctl3),as.matrix(alpha3.c[,j])))
        #   }
        # }
        
        ##update all cells, all k at once
        k <- k3
        {
          ##terms would be repetitively used
          z.w.a <- Matrix::crossprod(Z3.c, as.matrix(alpha3.c*wt.ctl3))
          a.w.a <- Matrix::crossprod(as.matrix(alpha3.c*wt.ctl3), as.matrix(alpha3.c))
          cp.w1 <- crossprod(W1[start.idx:end.idx,], W1[start.idx:end.idx,]) #crossprod using W1 freshly updated above
          cp.w2 <- crossprod(W2[start.idx:end.idx,], W2[start.idx:end.idx,]) #crossprod using W2 freshly updated above
          cp.w3 <- crossprod(W3[start.idx:end.idx,], W3[start.idx:end.idx,]) #need to use the old one
          
          ##to form the scalar below ( can speed up by applying crossprod(x,y) = t(x) %*% y )
          s.op.front <- t(as.matrix(rep(1,k))) #scalar operator front
          s.op.back <- as.matrix(rep(1,k)) #scalar operator back
          S1 <- s.op.front %*% cp.w3 %*% cp.w1 %*% solve(cp.w2) %*% solve(cp.w3) %*% s.op.back
          S2 <- s.op.front %*% t(W3[start.idx:end.idx,]) %*% z.w.a %*% solve(cp.w2) %*% solve(cp.w3) %*% s.op.back
          S3 <- s.op.front %*% cp.w3 %*% a.w.a %*% solve(cp.w2) %*% solve(cp.w3) %*% s.op.back
          
          inv.w3w2 <- solve(cp.w2) %*% solve(cp.w3)
          tau1 <- (2*S1 - 2*s.op.front %*% cp.w3 %*% cp.w1 %*% inv.w3w2 %*% s.op.back + 0.000001)^(-1)*(
            - (s.op.front %*% t(W3[start.idx:end.idx,]) %*% z.w.a %*% inv.w3w2 %*% s.op.back)
            + (s.op.front %*% cp.w3 %*% a.w.a %*% inv.w3w2 %*% s.op.back)
            + (S2 - S3)
          )
          tau2 <- (2*k)^(-1)*(s.op.front %*% t(W3[start.idx:end.idx,]) %*% z.w.a %*% inv.w3w2 %*% s.op.back)
          - (2*k)^(-1)*(s.op.front %*% cp.w3 %*% a.w.a %*% inv.w3w2 %*% s.op.back)
          - k^(-1)*tau1*(s.op.front %*% cp.w3 %*% cp.w1 %*% inv.w3w2 %*% s.op.back)
          
          S4 <- tau1
          S5 <- (1/k)*s.op.front %*% t(W3[start.idx:end.idx,]) %*% z.w.a %*% inv.w3w2 %*% s.op.back
          S6 <- (1/k)*s.op.front %*% cp.w3 %*% a.w.a %*% inv.w3w2 %*% s.op.back
          S7 <- (2*tau1/k)*s.op.front %*% cp.w3 %*% cp.w1 %*% inv.w3w2 %*% s.op.back
          
          additional.terms <- (2*as.numeric(S4)*cp.w1 + as.numeric(S5 - S6 - S7)*cp.w2)
        }
        W3[start.idx:end.idx,]  <- z.w.a %*% solve(lambda.a3 + a.w.a + additional.terms)
        
      }
      
      if(ncol(W1)!=k1)
        W1 <- t(W1)
      if(ncol(W2)!=k2)
        W2 <- t(W2)
      if(ncol(W3)!=k3)
        W3 <- t(W3)
      
      # regularise W
      #W1
      W1.med  <- matrixStats::colMedians(W1sub)
      W1.mad  <- matrixStats::colMads(W1sub)
      for(k1 in 1:nW1) {
        high <- W1[,k1] > (W1.med[k1]+4.*W1.mad[k1])
        low  <- W1[,k1] < (W1.med[k1]-4.*W1.mad[k1])
        W1[,k1]      <- (1-(high | low))*W1[,k1] + high*(W1.med[k1]+4*W1.mad[k1]) + low*(W1.med[k1]-4*W1.mad[k1])
      }
      #W2
      W2.med  <- matrixStats::colMedians(W2sub)
      W2.mad  <- matrixStats::colMads(W2sub)
      for(k2 in 1:nW2) {
        high <- W2[,k2] > (W2.med[k2]+4.*W2.mad[k2])
        low  <- W2[,k2] < (W2.med[k2]-4.*W2.mad[k2])
        W2[,k2]      <- (1-(high | low))*W2[,k2] + high*(W2.med[k2]+4*W2.mad[k2]) + low*(W2.med[k2]-4*W2.mad[k2])
      }
      #W3
      W3.med  <- matrixStats::colMedians(W3sub)
      W3.mad  <- matrixStats::colMads(W3sub)
      for(k3 in 1:nW3) {
        high <- W3[,k3] > (W3.med[k3]+4.*W3.mad[k3])
        low  <- W3[,k3] < (W3.med[k3]-4.*W3.mad[k3])
        W3[,k3]      <- (1-(high | low))*W3[,k3] + high*(W3.med[k3]+4*W3.mad[k3]) + low*(W3.med[k3]-4*W3.mad[k3])
      }
      
      
      # orthogonalize W (optional)
      if(nW1>1 & ortho.W)
        W1 <- matlib::GramSchmidt(W1)
      if(nW2>1 & ortho.W)
        W2 <- matlib::GramSchmidt(W2)
      if(nW3>1 & ortho.W)
        W3 <- matlib::GramSchmidt(W3)
      
      # update W until convergence [ordering: W1 -> W2 -> W3 -> W1 -> W2 -> W3 -> ...]
      conv.W1 = conv.W2 = conv.W3 <- FALSE
      iter.W1 = iter.W2 = iter.W3 <- 0
      trace_w13 = trace_w23 <- c()
      tau1.v = tau2.v <- c()
      while(!conv.W1 | !conv.W2 | !conv.W3) { #for any one of them not yet converge
        if(!conv.W1){
          iter.W1 <- iter.W1 + 1
          
          W1.old <- W1
          for(block in 1:nb) {
            start.idx <- (block-1)*block.size+1 ; end.idx <- min(ns,block*block.size)
            Y1sub  <- as.matrix(Y1[ctl1,start.idx:end.idx])
            #below rethink
            lmu.hat <- gmean3[match(rownames(Y1)[ctl1], names(gmean3))] + Matrix::tcrossprod(alpha1.c,W1[start.idx:end.idx,,drop=FALSE]) + Matrix::tcrossprod(alpha3_1[match(rownames(Y1)[ctl1], rownames(alpha3_1)),],W3[start.idx:end.idx,,drop=FALSE])
            Z1.c  <- lmu.hat - gmean3[match(rownames(Y1)[ctl1], names(gmean3))] - Matrix::tcrossprod(alpha3_1[match(rownames(Y1)[ctl1], rownames(alpha3_1)),],W3[start.idx:end.idx,,drop=FALSE]) + (Y1sub+0.01)/(exp(lmu.hat)+0.01) - 1
            W1[start.idx:end.idx,1]<-Matrix::crossprod(Z1.c,as.matrix(alpha1.c[,1]*wt.ctl1)) %*% solve(lambda.a1 + Matrix::crossprod(as.matrix(alpha1.c[,1]*wt.ctl1),as.matrix(alpha1.c[,1])))
            if(k1>1) {
              for(j in 2:k1) {
                Z1.c <- Z1.c - outer(alpha1.c[,j-1],W1[start.idx:end.idx,j-1])
                W1[start.idx:end.idx,j] <- Matrix::crossprod(Z1.c,as.matrix(alpha1.c[,j]*wt.ctl1)) %*% solve(lambda.a1 + Matrix::crossprod(as.matrix(alpha1.c[,j]*wt.ctl1),as.matrix(alpha1.c[,j])))
              }
            }
          }
          
          if(ncol(W1)!=k1)
            W1 <- t(W1)
          
          # regularise W
          #W1
          for(k1 in 1:nW1) {
            high <- W1[,k1] > (W1.med[k1]+4.*W1.mad[k1])
            low  <- W1[,k1] < (W1.med[k1]-4.*W1.mad[k1])
            W1[,k1]      <- (1-(high | low))*W1[,k1] + high*(W1.med[k1]+4*W1.mad[k1]) + low*(W1.med[k1]-4*W1.mad[k1])
          }
          
          # orthogonalize W (optional)
          if(nW1>1 & ortho.W)
            W1 <- matlib::GramSchmidt(W1)
          ##
          crit.W1 <- mean( (abs(W1-W1.old)/abs(W1.old))^2)
          #print(round(crit.W1,10))
          conv.W1 <- crit.W1< 1e-5 | iter.W1>8
        }
        if(!conv.W2){
          iter.W2 <- iter.W2 + 1
          
          W2.old <- W2
          for(block in 1:nb) {
            start.idx <- (block-1)*block.size+1 ; end.idx <- min(ns,block*block.size)
            Y2sub  <- as.matrix(Y2[ctl2,start.idx:end.idx])
            #below rethink
            lmu.hat <- gmean3[match(rownames(Y2)[ctl2], names(gmean3))] + Matrix::tcrossprod(alpha2.c,W2[start.idx:end.idx,,drop=FALSE]) + Matrix::tcrossprod(alpha3_2[match(rownames(Y2)[ctl2], rownames(alpha3_2)),],W3[start.idx:end.idx,,drop=FALSE])
            Z2.c  <- lmu.hat - gmean3[match(rownames(Y2)[ctl2], names(gmean3))] - Matrix::tcrossprod(alpha3_2[match(rownames(Y2)[ctl2], rownames(alpha3_2)),],W3[start.idx:end.idx,,drop=FALSE]) + (Y2sub+0.01)/(exp(lmu.hat)+0.01) - 1
            W2[start.idx:end.idx,1]<-Matrix::crossprod(Z2.c,as.matrix(alpha2.c[,1]*wt.ctl2)) %*% solve(lambda.a2 + Matrix::crossprod(as.matrix(alpha2.c[,1]*wt.ctl2),as.matrix(alpha2.c[,1])))
            if(k2>1) {
              for(j in 2:k2) {
                Z2.c <- Z2.c - outer(alpha2.c[,j-1],W2[start.idx:end.idx,j-1])
                W2[start.idx:end.idx,j] <- Matrix::crossprod(Z2.c,as.matrix(alpha2.c[,j]*wt.ctl2)) %*% solve(lambda.a2 + Matrix::crossprod(as.matrix(alpha2.c[,j]*wt.ctl2),as.matrix(alpha2.c[,j])))
              }
            }
          }
          
          if(ncol(W2)!=k2)
            W2 <- t(W2)
          
          # regularise W
          #W2
          for(k2 in 1:nW2) {
            high <- W2[,k2] > (W2.med[k2]+4.*W2.mad[k2])
            low  <- W2[,k2] < (W2.med[k2]-4.*W2.mad[k2])
            W2[,k2]      <- (1-(high | low))*W2[,k2] + high*(W2.med[k2]+4*W2.mad[k2]) + low*(W2.med[k2]-4*W2.mad[k2])
          }
          
          # orthogonalize W (optional)
          if(nW2>1 & ortho.W)
            W2 <- matlib::GramSchmidt(W2)
          ##
          crit.W2 <- mean( (abs(W2-W2.old)/abs(W2.old))^2)
          #print(round(crit.W2,10))
          conv.W2 <- crit.W2< 1e-5 | iter.W2>8
        }
        if(!conv.W3){
          iter.W3 <- iter.W3 + 1
          
          W3.old <- W3
          for(block in 1:nb) {
            start.idx <- (block-1)*block.size+1 ; end.idx <- min(ns,block*block.size)
            Y12sub  <- as.matrix(Y12[ctl3,start.idx:end.idx]) #check if ctl3 is really formed from Y12 [CHECK]
            #below rethink
            lmu.hat <- gmean3[match(rownames(Y12)[ctl3], names(gmean3))] +
              rbind(Matrix::tcrossprod(alpha1[which(rownames(alpha1) %in% rownames(Y12)[ctl3]),], W1[start.idx:end.idx,,drop=FALSE]), #alpha1 for NC3
                    Matrix::tcrossprod(alpha2[which(rownames(alpha2) %in% rownames(Y12)[ctl3]),], W2[start.idx:end.idx,,drop=FALSE])) + #alpha2 for NC3
              Matrix::tcrossprod(alpha3.c,W3[start.idx:end.idx,,drop=FALSE])
            
            Z3.c  <- lmu.hat - gmean3[match(rownames(Y12)[ctl3], names(gmean3))] - rbind(Matrix::tcrossprod(alpha1[which(rownames(alpha1) %in% rownames(Y12)[ctl3]),], W1[start.idx:end.idx,,drop=FALSE]),
                                                                                         Matrix::tcrossprod(alpha2[which(rownames(alpha2) %in% rownames(Y12)[ctl3]),], W2[start.idx:end.idx,,drop=FALSE])) +
              (Y12sub+0.01)/(exp(lmu.hat)+0.01) - 1
            # W3[start.idx:end.idx,1]<-Matrix::crossprod(Z3.c,as.matrix(alpha3.c[,1]*wt.ctl3)) %*% solve(lambda.a3 + Matrix::crossprod(as.matrix(alpha3.c[,1]*wt.ctl3),as.matrix(alpha3.c[,1])))
            # if(k3>1) {
            #   for(j in 2:k3) {
            #     Z3.c <- Z3.c - outer(alpha3.c[,j-1],W3[start.idx:end.idx,j-1])
            #     W3[start.idx:end.idx,j] <- Matrix::crossprod(Z3.c,as.matrix(alpha3.c[,j]*wt.ctl3)) %*% solve(lambda.a3 + Matrix::crossprod(as.matrix(alpha3.c[,j]*wt.ctl3),as.matrix(alpha3.c[,j])))
            #   }
            # }
            
            ##update all cells, all k at once
            ##update all cells, all k at once
            k <- k3
            {
              ##terms would be repetitively used
              z.w.a <- Matrix::crossprod(Z3.c, as.matrix(alpha3.c*wt.ctl3))
              a.w.a <- Matrix::crossprod(as.matrix(alpha3.c*wt.ctl3), as.matrix(alpha3.c))
              cp.w1 <- crossprod(W1[start.idx:end.idx,], W1[start.idx:end.idx,]) #crossprod using W1 freshly updated above
              cp.w2 <- crossprod(W2[start.idx:end.idx,], W2[start.idx:end.idx,]) #crossprod using W2 freshly updated above
              cp.w3 <- crossprod(W3[start.idx:end.idx,], W3[start.idx:end.idx,]) #need to use the old one
              
              ##to form the scalar below ( can speed up by applying crossprod(x,y) = t(x) %*% y )
              s.op.front <- t(as.matrix(rep(1,k))) #scalar operator front
              s.op.back <- as.matrix(rep(1,k)) #scalar operator back
              S1 <- s.op.front %*% cp.w3 %*% cp.w1 %*% solve(cp.w2) %*% solve(cp.w3) %*% s.op.back
              S2 <- s.op.front %*% t(W3[start.idx:end.idx,]) %*% z.w.a %*% solve(cp.w2) %*% solve(cp.w3) %*% s.op.back
              S3 <- s.op.front %*% cp.w3 %*% a.w.a %*% solve(cp.w2) %*% solve(cp.w3) %*% s.op.back
              
              inv.w3w2 <- solve(cp.w2) %*% solve(cp.w3)
              tau1 <- (2*S1 - 2*s.op.front %*% cp.w3 %*% cp.w1 %*% inv.w3w2 %*% s.op.back + 0.000001)^(-1)*(
                - (s.op.front %*% t(W3[start.idx:end.idx,]) %*% z.w.a %*% inv.w3w2 %*% s.op.back)
                + (s.op.front %*% cp.w3 %*% a.w.a %*% inv.w3w2 %*% s.op.back)
                + (S2 - S3)
              )
              tau2 <- (2*k)^(-1)*(s.op.front %*% t(W3[start.idx:end.idx,]) %*% z.w.a %*% inv.w3w2 %*% s.op.back)
              - (2*k)^(-1)*(s.op.front %*% cp.w3 %*% a.w.a %*% inv.w3w2 %*% s.op.back)
              - k^(-1)*tau1*(s.op.front %*% cp.w3 %*% cp.w1 %*% inv.w3w2 %*% s.op.back)
              
              S4 <- tau1
              S5 <- (1/k)*s.op.front %*% t(W3[start.idx:end.idx,]) %*% z.w.a %*% inv.w3w2 %*% s.op.back
              S6 <- (1/k)*s.op.front %*% cp.w3 %*% a.w.a %*% inv.w3w2 %*% s.op.back
              S7 <- (2*tau1/k)*s.op.front %*% cp.w3 %*% cp.w1 %*% inv.w3w2 %*% s.op.back
              
              additional.terms <- (2*as.numeric(S4)*cp.w1 + as.numeric(S5 - S6 - S7)*cp.w2)
            }
            W3[start.idx:end.idx,]  <- z.w.a %*% solve(lambda.a3 + a.w.a + additional.terms)
            
          }
          
          ##after updating all blocks
          all.cp.w1 <- crossprod(W1, W1)
          all.cp.w2 <- crossprod(W2, W2)
          
          trace_w13[iter.W3] <- sum(diag(W3 %*% crossprod(W1, W1) %*% t(W3)))
          trace_w23[iter.W3] <- sum(diag(W3 %*% crossprod(W2, W2) %*% t(W3)))
          
          tau1.v[iter.W3] <- as.numeric(tau1) 
          tau2.v[iter.W3] <- as.numeric(tau2)
          
          
          if(ncol(W3)!=k3)
            W3 <- t(W3)
          
          # regularise W
          #W3
          for(k3 in 1:nW3) {
            high <- W3[,k3] > (W3.med[k3]+4.*W3.mad[k3])
            low  <- W3[,k3] < (W3.med[k3]-4.*W3.mad[k3])
            W3[,k3]      <- (1-(high | low))*W3[,k3] + high*(W3.med[k3]+4*W3.mad[k3]) + low*(W3.med[k3]-4*W3.mad[k3])
          }
          
          # orthogonalize W (optional)
          if(nW3>1 & ortho.W)
            W3 <- matlib::GramSchmidt(W3)
          ##
          crit.W3 <- mean( (abs(W3-W3.old)/abs(W3.old))^2)
          #print(round(crit.W3,10))
          conv.W3 <- crit.W3< 1e-5 | iter.W3>8
        }
        
        
        # fixed the sign of W
        #for(i in 1:ncol(W))
        # W[,i] <- W[,i] * sign(cor(W[subsamples.org[subsubsamples.org],i],Wsub[1:length(subsubsamples.org),i]))
      }
      aft=Sys.time()
      #print(paste0('Time to estimate W for all samples:',difftime(aft,bef,units='secs')))
    }
  } # end of outer IRLS loop
  
  # now calculate Mb
  print('Estimating Mb....')
  Mb.all <- Y12
  # for cells with annotation
  idx.annot <- unlist(rep.ind)
  Mb.all[,idx.annot] <- Mb[,apply(M[idx.annot,],1,which)]
  # for other cells
  if(length(idx.annot) < ncol(Y12) ) {
    for(i in 1:3) {  
      for(block in 1:nb) {
        start.idx <- (block-1)*block.size+1 ; end.idx <- min(ns,block*block.size)
        lmu  <- gmean3 + as.matrix(Mb.all[,start.idx:end.idx]) + rbind(Matrix::tcrossprod(alpha1,W1[start.idx:end.idx,,drop=FALSE]), Matrix::tcrossprod(alpha2,W2[start.idx:end.idx,,drop=FALSE])) + Matrix::tcrossprod(rbind(alpha3_1,alpha3_2),W3[start.idx:end.idx,,drop=FALSE])
        Z    <- as.matrix(Mb.all[,start.idx:end.idx]) + ( (as.matrix(Y12[,start.idx:end.idx])+0.01)/(exp(lmu)+0.01) - 1)
        if(nbatch>1)
          wtmat<- 1/(exp(-lmu) + psi[,batch[start.idx:end.idx]])
        if(nbatch==1)
          wtmat<- 1/(exp(-lmu) + psi)
        Mb.all[,start.idx:end.idx] <- Z * (wtmat/(wtmat+lambda.b))
      }
    }
  }
  Mb.all[,idx.annot] <- Mb[,apply(M[idx.annot,],1,which)]
  Mb.all[union(union(rownames(Y1)[ctl1], rownames(Y2)[ctl2]), rownames(Y12)[ctl3]),] <- 0
  
  #
  Mb.med  <- DelayedMatrixStats::rowMedians(Mb.all[,idx.annot])
  Mb.mad  <- DelayedMatrixStats::rowMads(Mb.all[,idx.annot])
  high    <- Mb.all > (Mb.med+4*Mb.mad)
  low     <- Mb.all < (Mb.med-4*Mb.mad)
  Mb.all  <- (1-(high | low))*Mb.all + high*(Mb.med+4*Mb.mad) + low*(Mb.med-4*Mb.mad)
  
  # calibrate alpha and W
  #for(k in 1:nW) {
  # ab.W     <- Rfast::lmfit(y=Wsub[1:length(subsubsamples.org),k],x=cbind(1,W[subsamples.org[subsubsamples.org],k]))$be
  # ab.alpha <- Rfast::lmfit(y=best.a[,k],x=cbind(1,alpha[,k]))$be
  # W[,k]    <- ab.W[1] + ab.W[2]*W[,k]
  # alpha[,k]<- ab.alpha[1] + ab.alpha[2]*alpha[,k]
  #}
  #alpha <- sweep(alpha,2,scfac.alpha,'/')
  # cor.check1 <- diag(as.matrix(cor(W1sub[1:length(subsubsamples.org),],W1[subsamples.org[subsubsamples.org],])))
  # cor.check2 <- diag(as.matrix(cor(W2sub[1:length(subsubsamples.org),],W2[subsamples.org[subsubsamples.org],])))
  # cor.check3 <- diag(as.matrix(cor(W3sub[1:length(subsubsamples.org),],W3[subsamples.org[subsubsamples.org],])))
  #equal to
  cor.check1 <- diag(as.matrix(cor(W1sub,W1[subsamples.org[subsubsamples.org],])))
  cor.check2 <- diag(as.matrix(cor(W2sub,W2[subsamples.org[subsubsamples.org],])))
  cor.check3 <- diag(as.matrix(cor(W3sub,W3[subsamples.org[subsubsamples.org],])))
  
  # subsamples.org     <- subsamples
  # subsubsamples.org  <- 1:nsub
  # W1sub      <- matrix(0,nsub,k1)
  
  
  if(use.pseudosample) { 
    psi <- psi[,1:nbatch.org]
  }
  names(psi) <- names(gmean3)
  
  END.time <- Sys.time()
  
  runtime <- (END.time - START.time)
  
  #output
  return( list("counts_modality1"=Y1, "counts_modality2"=Y2,
               "gmean"=gmean3, #only this one got updated in IRLS
               "W1"=W1, "W2"=W2, "W3"=W3, 
               "W1sub"=W1sub, "W2sub"=W2sub, "W3sub"=W3sub, 
               "a1"=alpha1, "a2"=alpha2, "a3_1"=alpha3_1, "a3_2"=alpha3_2,
               'L.a1'=lambda.a1, 'L.a2'=lambda.a2, 'L.a3'=lambda.a3,
               'corW'=list(cor.check1,cor.check2,cor.check3),
               # "M"="skip", "Mb.sub"="skip", "Mb"="skip",
               "M"=M, "Mb.sub"=Mb, "Mb"=Mb.all,
               'L.b'=lambda.b,
               "ctl1"=ctl1, "ctl2"=ctl2, "ctl3"=ctl3, 
               "psi"=psi,
               "batch"=batch,"subsamples"=subsamples.org,
               "tau1_sub.ls"=tau1_sub.ls, "tau2_sub.ls"=tau2_sub.ls,
               "trace_sub.w13.ls"=trace_sub.w13.ls, "trace_sub.w23.ls"=trace_sub.w23.ls, #trace from subset of cells
               "logl.inner.ls"=loglikehood.ls,
               "logl.outer"=logl.outer, 
               "trace_w13"=trace_w13, "trace_w23"=trace_w23, #trace from all cells
               "tau1.v"=tau1.v, "tau2.v"=tau2.v,
               "iter.w1"=iter.W1, "iter.w2"=iter.W2, "iter.w3"=iter.W3, ##knowing iterations for updating W for all cells
               "runtime"=runtime
  ) )
}

#### Helpers functions ######################################################################################
rob.wt <- function(Y,mu,psi) {
  sign.dev <- sign(Y-mu)* sqrt(2*(dnbinom(Y,mu=Y,size=1/psi,log=TRUE) - dnbinom(Y,mu=mu,size=1/psi,log=TRUE) ))
  mad.est  <- mad(sign.dev)
  sign.dev <- sign.dev/mad.est
  matrix(MASS::psi.huber(sign.dev),byrow=T,nrow=nrow(sign.dev),ncol=ncol(sign.dev))
}

projection <- function(ind,Zmat,Wmat=NULL,lambda=0) {
  if(is.null(Wmat)) 
    out <- sapply(ind,rowAvg,Z=Zmat)
  if(!is.null(Wmat)) 
    out <- sapply(ind,rowWtAvg,Z=Zmat,Wt=Wmat,lambda=lambda)
  
  out
}

projection.K1 <- function(ind,Zmat,Wmat=NULL) {
  if(is.null(Wmat)) 
    out <- sapply(ind,rowAvg.K1,Z=Zmat)
  if(!is.null(Wmat)) 
    out <- sapply(ind,rowWtAvg.K1,Z=Zmat,Wt=Wmat)
  
  out
}

rowAvg <- function(ind,Z) {
  if(length(ind)>1) {
    out <- matrixStats::rowMeans2(as.matrix(Z[,ind]), na.rm=T)
  }
  if(length(ind)==1) {
    out <- Z[,ind]
  }
  out
}

rowWtAvg <- function(ind,Z,Wt,lambda=0) {
  if(length(ind)>1) {
    out <- rowSums(as.matrix(Z[,ind]*Wt[,ind]), na.rm=T)/(matrixStats::rowSums2(as.matrix(Wt[,ind]), na.rm=T) + lambda)  
  }
  if(length(ind)==1) {
    out <- c(Z[,ind])*c(Wt[,ind])/c(Wt[,ind]+lambda)
  }
  out
}

rowAvg.K1 <- function(ind,Z) {
  if(length(ind)>1) {
    out <- matrixStats::colMeans2(as.matrix(Z[,ind]), na.rm=T)
  }
  if(length(ind)==1) {
    out <- Z[,ind]
  }
  out
}

rowWtAvg.K1 <- function(ind,Z,Wt) {
  if(length(ind)>1) {
    out <- matrixStats::colSums2(as.matrix(Z[,ind]*Wt[,ind]), na.rm=T)/matrixStats::colSums2(as.matrix(Wt[,ind]), na.rm=T)   
  }
  if(length(ind)==1) {
    out <- Z[,ind]
  }
  out
}

subsetMat <- function(index,MARGIN,mat) {
  if(MARGIN==1)
    return(mat[index,])
  if(MARGIN==2)
    return(mat[,index])
}

get.coef.wt <- function(y,wt,W) {
  Rfast::lmfit(y=y,x=cbind(1,W),w=wt)$be
}

get.coef2 <- function(RM_Z,signdev,alpha.dev,wt,RM_W,k.huber=1.345) {
  sigma.est<- mad(signdev)/0.6745
  if(sigma.est!=0)
    wtvec    <- MASS::psi.huber(signdev,k=k.huber*sigma.est) * wt 
  if(sigma.est==0)
    wtvec    <- wt 
  wtvec    <- wtvec/mean(wtvec)
  A        <- Matrix::crossprod(RM_W*wtvec,RM_W)
  b        <- Matrix::crossprod(RM_W*wtvec, as.matrix(RM_Z - RM_W %*% as.matrix(alpha.dev)))
  return(A=cbind(A,b))
}

get.adev <- function(RM_Z,signdev,wt,RM_W,lambda.a,alpha.mean,k.huber=1.345) {
  sigma.est <- mad(signdev)/0.6745
  if(sigma.est==0)
    wtvec   <- wt 
  if(sigma.est!=0)
    wtvec   <- MASS::psi.huber(signdev,k=k.huber*sigma.est) * wt 
  wtvec <- wtvec/mean(wtvec)
  A <- Matrix::crossprod(RM_W*wtvec,RM_W)
  b <- Matrix::crossprod(RM_W*wtvec,as.matrix(RM_Z - RM_W %*% alpha.mean))
  return(solve(A+lambda.a*diag(ncol(A)),b))
}

getW  <- function(Z,signdev,wt,alpha,ctl,k.huber=1.345,k=2,zeroinf=FALSE) {
  W <- rep(0,k)
  sigma.est <- mad(signdev[ctl])/0.6745
  if(sigma.est!=0)
    wtvec     <- MASS::psi.huber(signdev[ctl],k=k.huber*sigma.est) * wt[ctl] 
  if(sigma.est==0)
    wtvec     <- wt[ctl] 
  wtvec     <- wtvec/mean(wtvec)
  W[1]      <- Rfast::lmfit(y=Z[ctl],x=alpha[ctl,1],w=wtvec)$be
  if(k>1) {
    for(j in 2:k) {
      W[j] <- Rfast::lmfit(y=Z[ctl]-matrixStats::rowSums2(as.matrix(alpha[ctl,1:(j-1)]) %*% as.matrix(W[1:(j-1)])),x=alpha[ctl,j],w=wtvec)$be
    } 
  }	   
  W
}

# estimate lambda.a parameter using RIDGM approach (Dempster et al 1977).
est.lambda.a <- function(alpha.dev,signdev,wt,RM_W,k.huber=1.345) {
  k <- length(alpha.dev)
  lambda.a.root <- function(lambda.a,alpha.dev,signdev,wt,RM_W,k.huber=k.huber) {
    sigma.est <- mad(signdev)/0.6745
    wtvec   <- MASS::psi.huber(signdev,k=k.huber*sigma.est) * wt 
    wtvec <- wtvec/mean(wtvec)
    TXW<- sweep(t(RM_W),2,wtvec,'*')
    A <- TXW %*% RM_W
    if(k>1)
      diff <- sum(alpha.dev^2/diag(solve(A+lambda.a*diag(k)))) - k
    if(k==1)
      diff <- sum(alpha.dev^2 * c(A+lambda.a)) - k
    
    diff
  }
  
  lambda.a.est <- tryCatch({uniroot(lambda.a.root,interval=c(0.001,10),extendInt="yes",alpha.dev=alpha.dev,signdev=signdev,wt=wt,RM_W=RM_W,k.huber=k.huber)$root},
                           error = function(e) {NA})
  lambda.a.est
}  


# ZINB logl to estimate pi0
logl.zinb <- function(p,y,mu,size) {
  -sum(ZIM::dzinb(y,lambda=mu,k=size,omega=1/(1+exp(-p)),log=TRUE))
}

