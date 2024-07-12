library(sBIC)
library(mclust)
library(splus2R)

adjust_dtrap_cols <- function(df) {
  colnames(df)[1] <- 'sample'
  df$source_name <- sub('[0-9]*$', '', rownames(df))
  rownames(df) <- df$sample
  df$sample <- NULL
  return(df)
}

# alpha = fraction of TRAP into lysate
# beta = contamination from lysate into TRAP
# gamma = amount of target tissue recovered from TRAP (default=1)
epsilon <- function(alpha, beta, gamma) gamma * alpha + beta * (1 - alpha)
Residual <- function(alpha, beta, gamma, Rt, Tr) ( gamma * Rt - epsilon(alpha, beta, gamma) * Tr ) / ((1 - alpha) * (gamma - beta))
Transcriptome <- function(alpha, beta, gamma, Rt, Tr) ( epsilon(alpha, beta, gamma) * Tr - beta * Rt ) / (alpha * (gamma - beta))
scoreFun <- function(alpha, beta, gamma, Rt, Tr, score = rep(1, length(Rt))) {
  r <- Residual(alpha, beta, gamma, Rt, Tr)
  t <- Transcriptome(alpha, beta, gamma, Rt, Tr)
  sum(score[r > 0 & t > 0])
}
# simplified system of equations
Residual2 <- function(alpha, beta, Rt, Tr) ( Rt - alpha * Tr ) / ( 1 - alpha * beta )
Transcriptome2 <- function(alpha, beta, Rt, Tr) ( Tr - beta * Rt ) / ( 1 - alpha * beta )
# system of equations that takes into account the cv of the libraries
Residual3 <- function(alpha, beta, Rt, Tr) ( Rt - alpha * Tr ) / ( 1 - alpha * beta )
Transcriptome3 <- function(alpha, beta, Rt, Tr) ( Tr - beta * Rt ) / ( 1 - alpha * beta )
#---
Residual3cv <- function(alpha, beta, P, Rt, Tr) ( Rt - alpha * Tr * P ) / ( 1 - alpha * beta )
Transcriptome3cv <- function(alpha, beta, P, Rt, Tr) ( Tr - beta * Rt * P ) / ( 1 - alpha * beta )
# system of equations that takes into account the cv of the libraries in a stepwise manner
Residual4 <- function(alpha, beta, P, Rt, Tr, p_threshold=.1) ( Rt - alpha * Tr * 1*(P < p_threshold) ) / ( 1 - alpha * beta )
Transcriptome4 <- function(alpha, beta, P, Rt, Tr, p_threshold=.1) ( Tr - beta * Rt * 1*(P < p_threshold) ) / ( 1 - alpha * beta )

log2FCDist <- function(x, y, N=101, minVal = 0, minQt=0, maxQt=1, 
                       from = NA, to = NA, delta=.5, ...) {
  stopifnot(length(minQt)==length(maxQt))
  iters <- length(minQt)
  # 0th filter (remove zeros)
  idx <- x > minVal & y > minVal
  x <- x[idx]; y <- y[idx]
  # log2 transform
  xLog <- log2(x); yLog <- log2(y)
  # calculate x-axis range 
  idx <- (xLog+yLog) >= quantile((xLog+yLog), probs = min(minQt)) & 
    (xLog+yLog) <= quantile((xLog+yLog), probs = max(maxQt))
  xRange <- xLog[idx]; yRange <- yLog[idx]
  idx <- is.finite(xRange) & is.finite(yRange)
  xRange <- xRange[idx]; yRange <- yRange[idx]
  r <- range(yRange-xRange)
  if(is.na(from)) from <- r[1]
  if(is.na(to)) to <- r[2]
  log2fc_intercept <- seq(from,to,length.out=N)
  
  log2smoothFun <- function(minQt, maxQt) {
    # 1st filter (quantiles)
    idx <- (xLog+yLog) > quantile((xLog+yLog), probs = minQt) & (xLog+yLog) < quantile((xLog+yLog), probs = maxQt)
    xLogQT <- xLog[idx]; yLogQT <- yLog[idx]
    # 2nd filter
    idx <- is.finite(xLogQT) & is.finite(yLogQT)
    xLogQT <- xLogQT[idx]; yLogQT <- yLogQT[idx]
    res <- sapply(log2fc_intercept, function(intercept) {
      length(which(
        ((yLogQT - xLogQT) < (intercept + delta)) & ((yLogQT - xLogQT) > (intercept - delta))
      ))/length(xLogQT)
    })
    res[1] <- res[1] + length(which((yLogQT - xLogQT) <= log2fc_intercept[1] - delta))
    res[N] <- res[N] + length(which((yLogQT - xLogQT) > log2fc_intercept[N] + delta))
    return(res)
  }
  
  if(iters == 1) {
    res <- log2smoothFun(minQt, maxQt)
  } else {
    res <- sapply(1:iters, function(i) log2smoothFun(minQt[i], maxQt[i]))
  }
  return(list(x=log2fc_intercept, y=res, data=yRange-xRange))
  # maxima <- 2+as.numeric(gregexpr('0-20', paste(diff(sign(diff(res))), collapse=''))[[1]])
  # abline(v=log2fc_intercept[maxima])
}

# log2FCDist <- function(x, y, N=101, minVal = 0, minQt=0, maxQt=1, 
#                        from = NA, to = NA, delta=.5, ...) {
#   stopifnot(length(minQt)==length(maxQt))
#   iters <- length(minQt)
#   # 0th filter (remove zeros)
#   idx <- x > minVal & y > minVal
#   x <- x[idx]; y <- y[idx]
#   # log2 transform
#   xLog <- log2(x); yLog <- log2(y)
#   # calculate x-axis range 
#   idx <- (xLog+yLog) >= quantile((xLog+yLog), probs = min(minQt)) & 
#     (xLog+yLog) <= quantile((xLog+yLog), probs = max(maxQt))
#   xRange <- xLog[idx]; yRange <- yLog[idx]
#   idx <- is.finite(xRange) & is.finite(yRange)
#   xRange <- xRange[idx]; yRange <- yRange[idx]
#   r <- range(yRange-xRange)
#   if(is.na(from)) from <- r[1]
#   if(is.na(to)) to <- r[2]
#   log2fc_intercept <- seq(from,to,length.out=N)
#   
#   log2smoothFun <- function(minQt, maxQt) {
#     # 1st filter (quantiles)
#     idx <- (xLog+yLog) > quantile((xLog+yLog), probs = minQt) & (xLog+yLog) < quantile((xLog+yLog), probs = maxQt)
#     xLogQT <- xLog[idx]; yLogQT <- yLog[idx]
#     # 2nd filter
#     idx <- is.finite(xLogQT) & is.finite(yLogQT)
#     xLogQT <- xLogQT[idx]; yLogQT <- yLogQT[idx]
#     res <- sapply(log2fc_intercept, function(intercept) {
#       length(which(
#         ((yLogQT - xLogQT) < (intercept + delta)) & ((yLogQT - xLogQT) > (intercept - delta))
#       ))/length(xLogQT)
#     })
#     res[1] <- res[1] + length(which((yLogQT - xLogQT) <= log2fc_intercept[1] - delta))
#     res[N] <- res[N] + length(which((yLogQT - xLogQT) > log2fc_intercept[N] + delta))
#     return(list(res = res, data=yLogQT - xLogQT))
#   }
#   
#   if(iters == 1) {
#     res <- log2smoothFun(minQt, maxQt)
#     data <- res$data
#     res <- res$res
#   } else {
#     res <- lapply(1:iters, function(i) log2smoothFun(minQt[i], maxQt[i]))
#     data <- lapply(res, '[', 2)
#     res <- sapply(res, '[', 1)
#   }
#   return(list(x=log2fc_intercept, y=res, data=data))
#   # maxima <- 2+as.numeric(gregexpr('0-20', paste(diff(sign(diff(res))), collapse=''))[[1]])
#   # abline(v=log2fc_intercept[maxima])
# }


# log2FCDist <- function(x, y, N=101, minVal = 0, minQt=0, maxQt=1, NQt=1,
#                        from = NA, to = NA, delta=.5, ...) {
#   stopifnot(minQt>=0 & minQt<=1)
#   stopifnot(maxQt>=0 & maxQt<=1)
#   qtls <- seq(minQt, maxQt, length.out=NQt+1)
# 
#   # remove low expressed
#   idx <- x > minVal & y > minVal
#   x <- x[idx]; y <- y[idx]
#   
#   # remove out of quantile genes
#   idx <- (x+y) >= quantile((x+y), probs = min(minQt)) & 
#     (x+y) <= quantile((x+y), probs = max(maxQt))
#   xTmp <- x[idx]; yTmp <- y[idx]
#   
#   # log2 transform and removo non finite
#   xTmp <- log2(xTmp); yTmp <- log2(yTmp)
#   idx <- is.finite(xTmp) & is.finite(yTmp)
#   xTmp <- xTmp[idx]; yTmp <- yTmp[idx]
#   
#   # define range of the log2 fold changes
#   r <- range(yTmp-xTmp)
#   if(is.na(from)) from <- r[1]
#   if(is.na(to)) to <- r[2]
#   
#   log2fc_intercept <- seq(from,to,length.out=N)
#   
#   tmpFun <- function(minQt, maxQt) {
#     # 1st filter (quantiles)
#     idx <- (x+y) > quantile((x+y), probs = minQt) & (x+y) < quantile((x+y), probs = maxQt)
#     x <- x[idx]; y <- y[idx]
#     # log2 transform
#     x <- log2(x); y <- log2(y)
#     # 2nd filter
#     idx <- is.finite(x) & is.finite(y)
#     x <- x[idx]; y <- y[idx]
#     res <- sapply(log2fc_intercept, function(intercept) {
#       length(which(
#         ((y - x) < (intercept + delta)) & ((y - x) > (intercept - delta))
#       ))/length(x)
#     })
#     res[1] <- res[1] + length(which((y - x) <= log2fc_intercept[1] - delta))
#     res[N] <- res[N] + length(which((y - x) > log2fc_intercept[N] + delta))
#     return(res)
#   }
#   
#   if(iters == 1) {
#     res <- tmpFun(minQt, maxQt)
#   } else {
#     res <- sapply(1:iters, function(i) tmpFun(minQt[i], maxQt[i]))
#   }
#   
#   return(list(x=log2fc_intercept, y=res, data=yTmp-xTmp))
#   # maxima <- 2+as.numeric(gregexpr('0-20', paste(diff(sign(diff(res))), collapse=''))[[1]])
#   # abline(v=log2fc_intercept[maxima])
# }

deconvolveTRAP <- function(cts_trap, cts_lysate, 
                           # alphas=NA, betas=NA,
                           method=c('dist','mean','mito'),
                           qt = .02,
                           minVal = 5, 
                           NQt=1, minQt=0.5, plot=T) {
  
  method <- method[1]
  stopifnot(method %in% c('dist','mean','mito'))
  stopifnot(ncol(cts_trap) == ncol(cts_lysate))
  n_samples <- ncol(cts_trap)
  
  alphas <- rep(NA, ncol(cts_trap))
  betas <- rep(NA, ncol(cts_trap))
  
  for(i in 1:n_samples) {
    tmp <- log2FCDist(cts_lysate[,i], cts_trap[,i], N=128, minVal = minVal,
                      minQt = seq(minQt,1,by=(1-minQt)/NQt)[-(NQt+1)],
                      maxQt = seq(minQt,1,by=(1-minQt)/NQt)[-1])
    #--------- model fc data
    minfit <- densityMclust(tmp$data, verbose=F, G=3, modelNames='V')
    if(plot) matplot(tmp$x, tmp$y, type='l', lty=1, col=1,
            lwd=seq(1/NQt,1,by=1/NQt)*2,
            main = paste(colnames(cts_lysate)[i], 'vs', colnames(cts_trap)[i]))
    mainx <- seq(min(tmp$x), max(tmp$x), length.out=128)
    cdens <- predict(minfit, mainx, what = "cdens")
    cdens <- t(apply(cdens, 1, function(d) d*minfit$parameters$pro))
    if(plot) matlines(mainx, cdens, col=c(3,1,4), lty=3:5, lwd=1)
    if(method=='dist' | method=='mito') {
      # search for intersection of the 2 sigmoids only between the means of the 2 external gaussians
      valid_search  <- which.min((mainx - minfit$parameters$mean[1])^2):which.min((mainx - minfit$parameters$mean[2])^2)
      betas[i]  <- -mainx[valid_search][which.min((cdens[valid_search,1] - cdens[valid_search,2])^2)]
      valid_search <- which.min((mainx - minfit$parameters$mean[2])^2):which.min((mainx - minfit$parameters$mean[3])^2)
      alphas[i] <- mainx[valid_search][which.min((cdens[valid_search,3] - cdens[valid_search,2])^2)]
      # betas[i]  <- - qnorm(p=qt[1], mean = minfit$parameters$mean[1], sd = sqrt(minfit$parameters$variance$sigmasq[1]))
      # alphas[i] <- qnorm(p=1-qt[2], mean = minfit$parameters$mean[3], sd = sqrt(minfit$parameters$variance$sigmasq[3]))
    }
    if(method=='mean') {
      betas[i]  <- -minfit$parameters$mean[1]
      alphas[i] <- minfit$parameters$mean[3]
    }
    #--------- display mt transcripts
    mt_idx <- grepl('^mt-', rownames(cts_lysate)) & cts_trap[,i]>minVal;
    log2mttrap   <- log2(cts_trap[mt_idx,i])
    log2mtlysate <- log2(cts_lysate[mt_idx,i])
    #
    idxmt <- is.finite(log2mttrap) & is.finite(log2mtlysate)
    log2mttrap <- log2mttrap[idxmt]
    log2mtlysate <- log2mtlysate[idxmt]
    #
    idx2mt <- (log2mttrap + log2mtlysate) >= quantile(log2mttrap + log2mtlysate, probs=0)
    log2mt <- log2mttrap[idx2mt] - log2mtlysate[idx2mt]
    if(plot) rug(log2mt, col=2, lwd=2)
    fit <- Mclust(log2mt, model='V', G=1, verbose=F)
    mtx <- mainx
    mty <- dens('V', mtx, parameters = fit$parameters)
    if(plot) lines(mtx, mty/max(mty)*max(tmp$y)/10, col=2, lwd=2)
    
    mtycumsum <- cumsum(mty)/sum(mty)
    idxtargetmt <- which.min((mtycumsum - (1-qt[1]))^2)
    if(plot) abline(v=mtx[idxtargetmt], lty=3, col=2, lwd=2)
    
    if(method=='mito') betas[i] <- -mtx[idxtargetmt]
    
    # if(is.na(alphas[i])) {
    #   betas[i] <- -min(sapply(NQt, function(k)
    #     tmp$x[min(which(cumsum(tmp$y[,k])/sum(tmp$y[,k])>qt))]))
    #   alphas[i] <- max(sapply(NQt, function(k)
    #     tmp$x[min(which(cumsum(tmp$y[,k])/sum(tmp$y[,k])>(1-qt)))]))
    # }
    if(plot) abline(v=alphas[i], lty=1, lwd=2, col=4)
    if(plot) abline(v=-betas[i], lty=1, lwd=2, col=3)
  }
  #
  # print(alphas)
  # print(betas)
  #
  lysate_deconv <- sapply(1:n_samples, function(i)
    Residual2(alpha = 1/2^alphas[i], beta = 1/2^betas[i], #gamma = 1,
              Rt = cts_lysate[,i],
              Tr = cts_trap[,i]
    ))
  lysate_deconv[lysate_deconv<0] <- 0
  lysate_deconv <- round(lysate_deconv)
  dimnames(lysate_deconv) <- dimnames(cts_lysate)
  tscptm_deconv <- sapply(1:n_samples, function(i)
    Transcriptome2(alpha = 1/2^alphas[i], beta = 1/2^betas[i], #gamma = 1,
                   Rt = cts_lysate[,i],
                   Tr = cts_trap[,i]
    ))
  tscptm_deconv[tscptm_deconv<0] <- 0
  tscptm_deconv <- round(tscptm_deconv)
  dimnames(tscptm_deconv) <- dimnames(cts_trap)
  # colnames(tscptm_deconv_wt) <- paste(paste0('EC_WT_',days), 1:3, sep='_')
  return(list(lysate=lysate_deconv, TRAP=tscptm_deconv, alphas=alphas, betas=betas))
  
}


deconvolveTRAPens <- function(cts_trap, cts_lysate,
                              sf_trap, sf_lysate,
                           # alphas=NA, betas=NA,
                           method=c('dist','mean','mito'),
                           paired=TRUE, use.cv=FALSE, 
                           qt = .02,
                           minVal = 5, 
                           NQt=1, minQt=0.5, G=NA, plot=T) {
  
  method <- method[1]
  stopifnot(method %in% c('dist','mean','mito'))
  if(paired) stopifnot(ncol(cts_trap) == ncol(cts_lysate))
  n_samples <- ncol(cts_trap)
  stopifnot(length(G) %in% c(1,n_samples))
  if(length(G)==1) G <- rep(G, n_samples)
  #------------------
  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(
    countData = cbind(cts_trap, cts_lysate), 
    colData   = data.frame(condition=c(rep('trap', ncol(cts_trap)), rep('input', ncol(cts_lysate)))), 
    design    = ~ condition))
  dds_trap  <- dds[, dds$condition == 'trap']
  dds_input <- dds[, dds$condition == 'input']
  #------------------
  norm_cts_trap  <- counts(dds_trap, normalized=TRUE)
  norm_cts_input <- counts(dds_input, normalized=TRUE)
  #------------------
  trap_cv  <- sqrt(rowVars(norm_cts_trap))/rowMeans(norm_cts_trap)
  input_cv <- sqrt(rowVars(norm_cts_input))/rowMeans(norm_cts_input)
  #------------------ impute NaNs 
  max_cv <- max(c(trap_cv, input_cv), na.rm = T)
  trap_cv[is.na(trap_cv)]   <- max_cv
  input_cv[is.na(input_cv)] <- max_cv
  #------------------
  cvLogRatio <- log(trap_cv/input_cv)
  idxCounts <- apply(cts_trap > minVal, 1, all) | apply(cts_lysate > minVal, 1, all)
  PlogRatio <- ecdf(cvLogRatio[idxCounts])(cvLogRatio)

  alphas <- rep(NA, ncol(cts_trap))
  betas <- rep(NA, ncol(cts_trap))

  for(i in 1:n_samples) {
    if(paired) {
      x <- rowMeans(norm_cts_input)
      main_str <- paste(colnames(cts_trap)[i], 'vs', colnames(cts_lysate)[i])
    } else {
      x <- norm_cts_input[,i]
      main_str <- paste(colnames(cts_trap)[i], 'vs Input')
    }
    tmp <- log2FCDist(x, norm_cts_trap[,i], N=128, minVal = minVal,
                      minQt = seq(minQt,1,by=(1-minQt)/NQt)[-(NQt+1)],
                      maxQt = seq(minQt,1,by=(1-minQt)/NQt)[-1])
    #--------- model fc data
    if(is.na(G[i])) { # estimate the best G
      minfit <- densityMclust(tmp$data, verbose=F, # G=3, 
                              modelNames='V')
      if( length(minfit$parameters$mean ) < 3) # at least 3 gaussian needed
        minfit <- densityMclust(tmp$data, verbose=F, G=3, 
                                modelNames='V')
    } else { # use the one provided by the user
      minfit <- densityMclust(tmp$data, verbose=F, G=G[i], 
                              modelNames='V')
    }
    last <- length(minfit$parameters$mean)
    #---------
    if(plot) 
      matplot(tmp$x, tmp$y, type='l', lty=1, col=1, 
              lwd=seq(1/NQt,1,by=1/NQt)*2, main = main_str)
    mainx <- seq(min(tmp$x), max(tmp$x), length.out=128)
    cdens <- predict(minfit, mainx, what = "cdens")
    cdens <- t(apply(cdens, 1, function(d) d*minfit$parameters$pro))
    if(plot) matlines(mainx, cdens[,c(1,last)], col=c(3,4), lty=3:5, lwd=1)
    if(method=='dist' ) {
      # # search for intersection of the 2 sigmoids only between the means of the 2 external gaussians
      # valid_search  <- which.min((mainx - minfit$parameters$mean[1])^2):which.min((mainx - minfit$parameters$mean[2])^2)
      # betas[i]  <- -mainx[valid_search][which.min((cdens[valid_search,1] - cdens[valid_search,2])^2)]
      # valid_search <- which.min((mainx - minfit$parameters$mean[2])^2):which.min((mainx - minfit$parameters$mean[3])^2)
      # alphas[i] <- mainx[valid_search][which.min((cdens[valid_search,3] - cdens[valid_search,2])^2)]
      betas[i]  <- - qnorm(p=qt[1], mean = minfit$parameters$mean[1], sd = sqrt(minfit$parameters$variance$sigmasq[1]))
      alphas[i] <- qnorm(p=1-qt[1], mean = minfit$parameters$mean[last], sd = sqrt(minfit$parameters$variance$sigmasq[last]))
    }
    if(method=='mean' | method=='mito') {
      betas[i]  <- -minfit$parameters$mean[1]
      alphas[i] <- minfit$parameters$mean[last]
    }
    #--------- display mt transcripts
    mt_idx <- grepl('^mt-', rownames(cts_lysate)) & cts_trap[,i]>minVal;
    log2mttrap   <- log2(cts_trap[mt_idx,i])
    log2mtlysate <- log2(cts_lysate[mt_idx,i])
    #
    idxmt <- is.finite(log2mttrap) & is.finite(log2mtlysate)
    log2mttrap <- log2mttrap[idxmt]
    log2mtlysate <- log2mtlysate[idxmt]
    #
    idx2mt <- (log2mttrap + log2mtlysate) >= quantile(log2mttrap + log2mtlysate, probs=0)
    log2mt <- log2mttrap[idx2mt] - log2mtlysate[idx2mt]
    if(plot) rug(log2mt, col=2, lwd=2)
    fit <- Mclust(log2mt, model='V', G=1, verbose=F)
    mtx <- mainx
    mty <- dens('V', mtx, parameters = fit$parameters)
    if(plot) lines(mtx, mty/max(mty)*max(tmp$y)/10, col=2, lwd=2)

    mtycumsum <- cumsum(mty)/sum(mty)
    idxtargetmt <- which.min((mtycumsum - (1-qt[1]))^2)
    if(plot) abline(v=mtx[idxtargetmt], lty=3, col=2, lwd=2)

    if(method=='mito') betas[i] <- -mtx[idxtargetmt]

    # if(is.na(alphas[i])) {
    #   betas[i] <- -min(sapply(NQt, function(k)
    #     tmp$x[min(which(cumsum(tmp$y[,k])/sum(tmp$y[,k])>qt))]))
    #   alphas[i] <- max(sapply(NQt, function(k)
    #     tmp$x[min(which(cumsum(tmp$y[,k])/sum(tmp$y[,k])>(1-qt)))]))
    # }
    if(plot) abline(v=alphas[i], lty=1, lwd=2, col=4)
    if(plot) abline(v=-betas[i], lty=1, lwd=2, col=3)
  }
  #
  # print(alphas)
  # print(betas)
  #
  lysate_deconv <- sapply(1:n_samples, function(i)
    if(use.cv)
      Residual3cv(alpha = 1/2^alphas[i], beta = 1/2^betas[i], 
                P = PlogRatio,
                Rt = norm_cts_input[,i],
                Tr = if(paired) rowMeans(norm_cts_trap) else norm_cts_trap[,i])
    else 
      Residual3(alpha = 1/2^alphas[i], beta = 1/2^betas[i], 
                Rt = norm_cts_input[,i],
                Tr = if(paired) rowMeans(norm_cts_trap) else norm_cts_trap[,i])
    )
  lysate_deconv[lysate_deconv<0] <- 0
  lysate_deconv <- round(lysate_deconv)
  dimnames(lysate_deconv) <- dimnames(cts_lysate)
  tscptm_deconv <- sapply(1:n_samples, function(i)
    if(use.cv)
      Transcriptome3cv(alpha = 1/2^alphas[i], beta = 1/2^betas[i], 
                     P = PlogRatio,
                     Rt = if(paired) rowMeans(norm_cts_input) else norm_cts_input[,i],
                     Tr = norm_cts_trap[,i])
    else
      Transcriptome3(alpha = 1/2^alphas[i], beta = 1/2^betas[i], 
                     Rt = if(paired) rowMeans(norm_cts_input) else norm_cts_input[,i],
                     Tr = norm_cts_trap[,i])
  )
  tscptm_deconv[tscptm_deconv<0] <- 0
  tscptm_deconv <- round(tscptm_deconv)
  dimnames(tscptm_deconv) <- dimnames(cts_trap)
  return(list(lysate=lysate_deconv, TRAP=tscptm_deconv, alphas=alphas, betas=betas))

}

deconvolveTRAPddsNew <- function(dds, idx_trap, idx_lysate,
                              method=c('boundary','peak'),
                              alphas=NULL, betas=NULL,
                              paired=TRUE, use.cv=FALSE, 
                              n_left_peak=1,
                              n_right_peak=1,
                              # qt = .02,
                              minVal = 5,
                              xlim = c(-10,10),
                              min_obs = 100,
                              smoothing_factor = 10,
                              # NQt=1, 
                              minQt=0.5, 
                              # G=NA, 
                              plot=T) {
  #-------
  stopifnot(all(method %in% c('boundary','peak')))
  stopifnot(length(method) %in% 1:2)
  if(length(method)==1) method <- rep(method, 2)
  #-------  
  if(paired) stopifnot(length(idx_trap) == length(idx_lysate))
  # stopifnot(length(qt) %in% c(1,2))
  # if(length(qt)==1) qt <- rep(qt, 2)
  n_samples <- length(idx_trap)
  n_input   <- length(idx_lysate)
  #----
  stopifnot(length(n_left_peak) %in% c(1,n_samples))
  if(length(n_left_peak)==1) n_left_peak <- rep(n_left_peak, n_samples)
  #----
  stopifnot(length(n_right_peak) %in% c(1,n_samples))
  if(length(n_right_peak)==1) n_right_peak <- rep(n_right_peak, n_samples)
  #----
  stopifnot(length(min_obs) %in% c(1,n_samples))
  if(length(min_obs)==1) min_obs <- rep(min_obs, n_samples)
  #----
  if(is.null(alphas)) alphas <- rep(NA, n_samples) else stopifnot(length(alphas) %in% c(1,n_samples))
  if(is.null(betas)) betas <- rep(NA, n_samples) else stopifnot(length(betas) %in% c(1,n_samples))
  if(length(alphas) == 1) alphas <- rep(alphas, n_samples)
  if(length(betas) == 1) betas <- rep(betas, n_samples)
  #------------------
  # mito <- rep(NA, n_samples)
  cts_trap   <- counts(dds, normalized=FALSE)[,idx_trap,drop=F]
  cts_input <- counts(dds, normalized=FALSE)[,idx_lysate,drop=F]
  stopifnot('sizeFactor' %in% colnames(colData(dds)))
  norm_cts_trap  <- counts(dds, normalized=TRUE)[,idx_trap,drop=F]
  norm_cts_input <- counts(dds, normalized=TRUE)[,idx_lysate,drop=F]
  #------------------ cv estimation
  if(use.cv) {
    trap_cv  <- sqrt(rowVars(norm_cts_trap))/rowMeans(norm_cts_trap)
    input_cv <- sqrt(rowVars(norm_cts_input))/rowMeans(norm_cts_input)
    max_cv <- max(c(trap_cv, input_cv), na.rm = T) # impute NaNs 
    trap_cv[is.na(trap_cv)]   <- max_cv
    input_cv[is.na(input_cv)] <- max_cv
    cvLogRatio <- log(trap_cv/input_cv)
    idxCounts <- apply(cts_trap > minVal, 1, all) | apply(cts_input > minVal, 1, all)
    PlogRatio <- ecdf(cvLogRatio[idxCounts])(cvLogRatio)
  }
  #------------------ alpha, beta estimation
  set.seed(1) # to fix Mclust results
  for(i in 1:n_samples) {
    if(paired) {
      x_norm <- log(norm_cts_input[,i])
      x_cnts <- cts_input[,i]
      main_str <- paste(colnames(norm_cts_trap)[i], 'vs', colnames(cts_input)[i])
    } else {
      x_norm <- log(rowMeans(norm_cts_input))
      x_cnts <- rowMeans(norm_cts_input)
      main_str <- paste(colnames(norm_cts_trap)[i], 'vs AVG Input')
    }
    y_norm <- log(norm_cts_trap[,i])
    y_cnts <- cts_trap[,i]
    #--------
    delta <- 20/512 * smoothing_factor
    x_rotated <- y_norm-x_norm
    y_rotated <- x_norm+y_norm
    x_steps <- seq(xlim[1],xlim[2],length.out=512)
    idx <- x_cnts>minVal & y_cnts>minVal
    x_density <- sapply(x_steps, function(x) {
      x_left  <- length(which(x_rotated[idx]>(x-delta) & x_rotated[idx]<x))
      x_right <- length(which(x_rotated[idx]<(x+delta) & x_rotated[idx]>x))
      if((x_left + x_right) > min_obs[i]) return(x_left / x_right) else return(NA)
    })
    x_density_2 <- sapply(x_steps, function(x) 
      length(which(x_rotated[idx]>(x-delta) & x_rotated[idx]<=(x+delta))))
    #--------- model fc data
    maxima_id <- which(peaks(x_density, span = 11))
    minima_id <- which(peaks(1/x_density, span = 11))
    #---
    left_boundary <- which.min(x_density) # min(minima_id)
    left_peak <- maxima_id[maxima_id>left_boundary][n_left_peak[i]]
    rigth_boundary <- which.max(x_density)
    rigth_peak <- rev(minima_id[minima_id<rigth_boundary])[n_right_peak[i]]
    #---
    left_boundary <- x_steps[left_boundary]
    left_peak <- x_steps[left_peak]
    rigth_boundary <- x_steps[rigth_boundary]
    rigth_peak <- x_steps[rigth_peak]
    #---
    if(method[1] == 'boundary')
      if(is.na(betas[i])) betas[i]  <- -left_boundary
    if(method[2] == 'boundary')
      if(is.na(alphas[i])) alphas[i] <- rigth_boundary
    if(method[1] == 'peak')
      if(is.na(betas[i])) betas[i]  <- -left_peak
    if(method[2] == 'peak')
      if(is.na(alphas[i])) alphas[i] <- rigth_peak
    #---------
    if(plot) {
      # plot(x_steps, x_density, log='', pch='.', main = main_str)
      plot(x_rotated, y_rotated, pch='.', xlim=xlim, main = main_str, 
           col = 'grey70')
      points(x_rotated[idx], y_rotated[idx], pch='.', xlim=xlim, main = main_str, 
           col = densCols(x_rotated[idx], y_rotated[idx], colramp = viridis::viridis)) 
      abline(v=c(left_boundary, left_peak, rigth_peak, rigth_boundary), col=c(3,3,4,4), lty=c(3,1,1,3), lwd=2)
      abline(v=alphas[i], lty=1, lwd=3, col=4)
      abline(v=-betas[i], lty=1, lwd=3, col=3)
      #-----
      plot(x_steps, x_density, log='y', type='l', main = main_str); abline(h=1, lty=2)
      points(x_steps[maxima_id], x_density[maxima_id], cex=3, col=2)
      points(x_steps[minima_id], x_density[minima_id], cex=3, col=4)
      abline(v=c(left_boundary, left_peak, rigth_peak, rigth_boundary), col=c(3,3,4,4), lty=c(3,1,1,3), lwd=2)
      abline(v=alphas[i], lty=1, lwd=3, col=4)
      abline(v=-betas[i], lty=1, lwd=3, col=3)
      #-----
      plot(x_steps, x_density_2, log='', type='l', main = main_str); abline(h=1, lty=2)
      abline(v=c(left_boundary, left_peak, rigth_peak, rigth_boundary), col=c(3,3,4,4), lty=c(3,1,1,3), lwd=2)
      abline(v=alphas[i], lty=1, lwd=3, col=4)
      abline(v=-betas[i], lty=1, lwd=3, col=3)
    }
    #--------- display mt transcripts
    # mt_idx <- grepl('^mt-', rownames(cts_lysate)) & cts_trap[,i]>minVal;
    # log2mttrap   <- log2(cts_trap[mt_idx,i])
    # if(paired) 
    #   log2mtlysate <- log2(cts_lysate[mt_idx,i])
    # else
    #   log2mtlysate <- log2(rowMeans(cts_lysate[mt_idx,,drop=FALSE]))
    # #
    # idxmt <- is.finite(log2mttrap) & is.finite(log2mtlysate)
    # log2mttrap <- log2mttrap[idxmt]
    # log2mtlysate <- log2mtlysate[idxmt]
    # #
    # idx2mt <- (log2mttrap + log2mtlysate) >= quantile(log2mttrap + log2mtlysate, probs=0)
    # log2mt <- log2mttrap[idx2mt] - log2mtlysate[idx2mt]
    # if(plot) rug(log2mt, col=2, lwd=2)
    # fit <- densityMclust(log2mt, model='V', verbose=F, G=1)
    # mtx <- mainx
    # mty <- dens('V', mtx, parameters = fit$parameters)
    # if(plot) lines(mtx, mty/max(mty)*max(tmp$y)/10, col=2, lwd=2)
    # #
    # mtycumsum <- cumsum(mty)/sum(mty)
    # idxtargetmt <- which.min((mtycumsum - (1-qt[1]))^2)
    # if(plot) abline(v=mtx[idxtargetmt], lty=3, col=2, lwd=2)
    #
    # if(method=='mito') betas[i] <- -mtx[idxtargetmt]
    # mito[i] <- -mtx[idxtargetmt]
    #------- display alpha and beta on the plot
    # if(plot) abline(v=alphas[i], lty=1, lwd=2, col=4)
    # if(plot) abline(v=betas[i], lty=1, lwd=2, col=3)
  }
  #
  lysate_deconv <- sapply(1:n_input, function(i)
    if(use.cv)
      Residual3cv(alpha = exp(-alphas[i]), beta = exp(-betas[i]), 
                  P = PlogRatio,
                  Rt = norm_cts_input[,i],
                  Tr = if(paired) norm_cts_trap[,i] else rowMeans(norm_cts_trap))
    else 
      Residual3(alpha = exp(-alphas[i]), beta = exp(-betas[i]), 
                Rt = norm_cts_input[,i],
                Tr = if(paired) norm_cts_trap[,i] else rowMeans(norm_cts_trap))
  )
  lysate_deconv[lysate_deconv<0] <- 0
  lysate_deconv <- round(lysate_deconv)
  dimnames(lysate_deconv) <- dimnames(cts_input)
  tscptm_deconv <- sapply(1:n_samples, function(i)
    if(use.cv)
      Transcriptome3cv(alpha = exp(-alphas[i]), beta = exp(-betas[i]), 
                       P = PlogRatio,
                       Rt = if(paired) norm_cts_input[,i] else rowMeans(norm_cts_input),
                       Tr = norm_cts_trap[,i])
    else
      Transcriptome3(alpha = exp(-alphas[i]), beta = exp(-betas[i]), 
                     Rt = if(paired) norm_cts_input[,i] else rowMeans(norm_cts_input),
                     Tr = norm_cts_trap[,i])
  )
  tscptm_deconv[tscptm_deconv<0] <- 0
  tscptm_deconv <- round(tscptm_deconv)
  dimnames(tscptm_deconv) <- dimnames(cts_trap)
  return(list(lysate=lysate_deconv, TRAP=tscptm_deconv, alphas=alphas, betas=betas))
  
}



deconvolveTRAPcv <- function(cts_trap, cts_lysate,
                              # alphas=NA, betas=NA,
                              paired=TRUE, 
                              minVal = 5, 
                              plot=T) {
  
  if(paired) stopifnot(ncol(cts_trap) == ncol(cts_lysate))
  n_samples <- ncol(cts_trap)
  #------------------
  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(
    countData = cbind(cts_trap, cts_lysate), 
    colData   = data.frame(condition=c(rep('trap', ncol(cts_trap)), rep('input', ncol(cts_lysate)))), 
    design    = ~ condition))
  dds_trap  <- dds[, dds$condition == 'trap']
  dds_input <- dds[, dds$condition == 'input']
  #------------------
  norm_cts_trap  <- counts(dds_trap, normalized=TRUE)
  norm_cts_input <- counts(dds_input, normalized=TRUE)
  #------------------
  trap_cv  <- sqrt(rowVars(norm_cts_trap))/rowMeans(norm_cts_trap)
  input_cv <- sqrt(rowVars(norm_cts_input))/rowMeans(norm_cts_input)
  #------------------ impute NaNs 
  max_cv <- max(c(trap_cv, input_cv), na.rm = T)
  trap_cv[is.na(trap_cv)]   <- max_cv
  input_cv[is.na(input_cv)] <- max_cv
  #------------------
  cvLogRatio <- log(trap_cv/input_cv)
  idxCounts <- apply(cts_trap > minVal, 1, all) | apply(cts_lysate > minVal, 1, all)
  PlogRatio <- ecdf(cvLogRatio[idxCounts])(cvLogRatio)
  #
  lysate_deconv <- sapply(1:n_samples, function(i) {
    Rt <- norm_cts_input[,i]
    Tr <- if(paired) rowMeans(norm_cts_trap) else norm_cts_trap[,i]
    return(Rt - PlogRatio * Tr)
    })
  lysate_deconv[lysate_deconv<0] <- 0
  lysate_deconv <- round(lysate_deconv)
  dimnames(lysate_deconv) <- dimnames(cts_lysate)
  tscptm_deconv <- sapply(1:n_samples, function(i) {
    Tr <- norm_cts_trap[,i]
    Rt <- if(paired) rowMeans(norm_cts_input) else norm_cts_input[,i]
    return(Tr - PlogRatio * Rt)
    })
  tscptm_deconv[tscptm_deconv<0] <- 0
  tscptm_deconv <- round(tscptm_deconv)
  dimnames(tscptm_deconv) <- dimnames(cts_trap)
  # colnames(tscptm_deconv_wt) <- paste(paste0('EC_WT_',days), 1:3, sep='_')
  return(list(lysate=lysate_deconv, TRAP=tscptm_deconv))
  
}


# function(trap_expr, lysate_expr) {
#   stopifnot(identical(names(trap_expr), names(lysate_expr)))
#   ids <- which(trap_expr > 3 | lysate_expr > 3)
#   l2fcdata <- log2(conditionmeans[ids,'Brain_ECTRAP']/conditionmeans[ids,'Brain_lysate'])
#   l2fcdata <- l2fcdata[is.finite(l2fcdata)]
#   trap_expr <- conditionmeans[names(l2fcdata),'Brain_ECTRAP']
#   
# }
# 
# 
# 
# N <- 10
# par(mfrow=c(N+1,1), mar=c(1,0.2,0.2,0.2))
# int <- quantile(exrpmeans, probs = seq(0,1,length.out=N+1))
# for(i in N:1) {
#   id <- names(exrpmeans)[which(exrpmeans >= int[i] & exrpmeans < int[i+1])]
#   hist(l2fcdata[names(l2fcdata) %in% id], breaks = seq(-12,12,by=.1), main='', 
#        xaxt = 'n', yaxt = 'n')
#   abline(v=0, lty=1, lwd=3)
# }
# axis(1, at=seq(-12,12), 
#      labels = round(1/(2^abs(seq(-12,12)))*100,1), 
#      padj = 0)
# #
# par(mfrow=c(1,1), mar=c(5,4,4,2))
# plot(x<-log2(conditionmeans[names(l2fcdata),'Brain_lysate']),
#      y<-log2(conditionmeans[names(l2fcdata),'Brain_ECTRAP']), pch='.',
#      xlab = 'Brain_lysate', ylab = 'Brain_ECTRAP')
# abline(0,1); abline(.5,1,lty=3); abline(-.5,1,lty=3)
# abline(-6.441868,1,col=2)
# abline(6.440335,1,col=2)
# #
# res <- findLog2FCIntercepts(x,y,delta = .5)