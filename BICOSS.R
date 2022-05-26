#rm(list = ls())
#Humandata <- matrix(NA, nrow = 500, ncol =  5000)
#allele_structures <- c("A","C","T","G")
#library(parallel)
#for(i in 1:ncol(Humandata)){
#  subsamples <- sample(allele_structures,2,replace = FALSE)
#  vec1 <- sample(subsamples,nrow(Humandata),replace = TRUE)
#  vec2 <- sample(subsamples,nrow(Humandata),replace = TRUE)
#  Humandata[,i] <- paste(vec1,vec2,sep = "-")
#  Humandata[sample(1:length(Humandata[,i]),2),i] <- NA
#}
internalEnumLC <- function (qrObj, ...){
  R <- qr.R(qrObj)
  numColumns <- dim(R)[2]
  rank <- qrObj$rank
  pivot <- qrObj$pivot
  if (is.null(numColumns) || rank == numColumns) {
    list()
  }
  else {
    p1 <- 1:rank
    X <- R[p1, p1]
    Y <- R[p1, -p1, drop = FALSE]
    b <- qr(X)
    b <- qr.coef(b, Y)
    b[abs(b) < 1e-06] <- 0
    lapply(1:dim(Y)[2], function(i) c(pivot[rank + i], pivot[which(b[,
                                                                     i] != 0)]))
  }
}

findLinearCombos <- function (x){
  if (!is.matrix(x))
    x <- as.matrix(x)
  lcList <- internalEnumLC(qr(x))
  initialList <- lcList
  badList <- NULL
  if (length(lcList) > 0) {
    continue <- TRUE
    while (continue) {
      tmp <- unlist(lapply(lcList, function(x) x[1]))
      tmp <- unique(tmp[!is.na(tmp)])
      badList <- unique(c(tmp, badList))
      lcList <- internalEnumLC(qr((x[, -badList, drop = FALSE])))
      continue <- (length(lcList) > 0)
    }
  }
  else badList <- NULL
  list(linearCombos = initialList, remove = badList)
}

kinship_humans <- function(SNPs,number_cores = 1,MAF){
  requireNamespace("rrBLUP")
  level_function <- function(x,snps){
    return(length(unique(snps[,x][!is.na(snps[,x])])))
  }
  if(.Platform$OS.type == "unix"){
    lengths <- mclapply(1:ncol(SNPs),level_function,mc.cores = number_cores,snps = SNPs)
  } else{
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("level_function"),envir=environment())
    lengths <- unlist(parLapply(cl,1:ncol(SNPs),level_function,snps = SNPs))
    stopCluster(cl)
  }
  std_fun <- function(x,snps,leng){
    if(leng[x] < 2){
      ret_1 <- matrix(0,nrow = nrow(snps),ncol = 1)
      ret_1[is.na(snps[,x]),] <- NA
      return(rbind(0,ret_1))
    }else if(leng[x] < 3){
      ret_1 <- integer(nrow(snps))
      ret_1[which(snps[,x] == names(which.max(table(snps[,x][!is.na(snps[,x])]))))] <- 1
      ret_1[is.na(snps[,x])] <- NA
      return(rbind(1 - sum(ret_1,na.rm = TRUE)/sum(!is.na(snps[,x])),matrix(ret_1,ncol = 1,nrow = nrow(snps))))
    }else{
      k <- unlist(strsplit(snps[,x][!is.na(snps[,x])],"-"))
      alleles <- unique(k)
      max_allele <- names(which.max(table(k)))
      dominant <- paste0(max_allele,"-",max_allele)
      minor <- paste0(alleles[!alleles%in%max_allele],"-",alleles[!alleles%in%max_allele])
      levels_o <- unique(c(dominant,unique(snps[,x])))
      q <- factor(snps[,x], levels = levels_o)
      ret_1 <- integer(nrow(snps))
      ret_1[which(snps[,x] == minor)] <- -1
      ret_1[which(snps[,x] == dominant)] <- 1
      ret_1[is.na(snps[,x])] <- NA
      return(rbind(min(table(k)/sum(table(k))),matrix(ret_1,ncol = 1,nrow = nrow(snps))))
    }
  }
  if(.Platform$OS.type == "unix"){
    z <- mclapply(1:ncol(SNPs),std_fun,mc.cores = number_cores,snps = SNPs,leng = lengths)
  } else{
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("std_fun"),envir=environment())
    z <- parLapply(cl,1:ncol(SNPs),std_fun,snps = SNPs,leng = lengths)
    stopCluster(cl)
  }
  z <- do.call(cbind, z)
  z <- z[-1,]
  return(A.mat(z,min.MAF = MAF,n.core = number_cores))
}


standardize <- function(SNPs,number_cores = 1,selfing = TRUE,MAF = 0.01){
  requireNamespace("parallel")
  
  if(!is.numeric(number_cores)){
    stop("number_cores is not numeric")
  }
  if(selfing){
    
    std_fun <- function(x,snps){
      ret_1 <- integer(nrow(snps))
      ret_1[which(snps[,x] == names(which.max(table(snps[,x][!is.na(snps[,x])]))))] <- 1
      ret_1[is.na(snps[,x])] <- NA
      return(ret_1)
    }
    if(.Platform$OS.type == "unix"){
      z <- mclapply(1:ncol(SNPs),std_fun,mc.cores = number_cores,snps = SNPs)
    } else{
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("std_fun"),envir=environment())
      z <- parLapply(cl,1:ncol(SNPs),std_fun,snps = SNPs)
      stopCluster(cl)
    }
    z <- do.call(cbind, z)
    level_function <- function(SNPs,MAF = MAF){
      if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
        stop("Not every column of SNPs is numeric")
      }
      if(sum(apply(SNPs,2,sum) >nrow(SNPs))){
        stop("Some values of the SNP matrix are not 0 or 1")
      }
      if(!is.numeric(MAF)){
        stop("MAF is not numeric")
      }
      z <- colMeans(SNPs,na.rm = TRUE)
      SNPs <- SNPs[,z <= (1 - MAF)]
      if(sum(!(z <= (1 - MAF))) == 0){
        level_dropped <- "None"
      }else{
        level_dropped <- which(!(z <= (1 - MAF)))
      }
      return(list(SNPs = SNPs,SNPs_Dropped = level_dropped))
    }
    z <- level_function(SNPs = z, MAF = MAF)
    return(z)
  }else{
    level_function <- function(x,snps){
      return(length(unique(snps[,x][!is.na(snps[,x])])))
    }
    if(.Platform$OS.type == "unix"){
      lengths <- mclapply(1:ncol(SNPs),level_function,mc.cores = number_cores,snps = SNPs)
    } else{
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("level_function"),envir=environment())
      lengths <- unlist(parLapply(cl,1:ncol(SNPs),level_function,snps = SNPs))
      stopCluster(cl)
    }
    std_fun <- function(x,snps,leng){
      if(leng[x] < 2){
        ret_1 <- matrix(0,nrow = nrow(snps),ncol = 1)
        ret_1[is.na(snps[,x])] <- NA
        return(rbind(0,ret_1))
      }else if(leng[x] < 3){
        ret_1 <- integer(nrow(snps))
        ret_1[which(snps[,x] == names(which.max(table(snps[,x][!is.na(snps[,x])]))))] <- 1
        ret_1[is.na(snps[,x])] <- NA
        return(rbind(1 - sum(ret_1[!is.na(snps[,x])])/sum(!is.na(snps[,x])),matrix(ret_1,ncol = 1,nrow = nrow(snps))))
      }else{
        k <- unlist(strsplit(snps[,x][!is.na(snps[,x])],"-"))
        alleles <- unique(k)
        max_allele <- names(which.max(table(k)))
        dominant <- paste0(max_allele,"-",max_allele)
        minor <- paste0(alleles[!alleles%in%max_allele],"-",alleles[!alleles%in%max_allele])
        levels_o <- unique(c(dominant,unique(snps[,x])))
        q <- factor(snps[,x], levels = levels_o)
        ret_1 <- integer(nrow(snps))
        ret_1[which(snps[,x] == minor)] <- -1
        ret_1[which(snps[,x] == dominant)] <- 1
        ret_1[is.na(snps[,x])] <- NA
        return(rbind(min(table(k)/sum(table(k))),matrix(ret_1,ncol = 1,nrow = nrow(snps))))
      }
    }
    if(.Platform$OS.type == "unix"){
      z <- mclapply(1:ncol(SNPs),std_fun,mc.cores = number_cores,snps = SNPs,leng = lengths)
    } else{
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("std_fun"),envir=environment())
      z <- parLapply(cl,1:ncol(SNPs),std_fun,snps = SNPs,leng = lengths)
      stopCluster(cl)
    }
    z <- do.call(cbind, z)
    maf <- z[1,]
    z <- z[-1,]
    if(!is.numeric(MAF)){
      stop("MAF is not numeric")
    }
    SNPs_return <- z[,!(maf < MAF)]
    return(list(SNPs = SNPs_return + 1,SNPs_Dropped = (1:ncol(z))[maf < MAF]))
  }
}

#myt <- standardize(SNPs = Humandata,number_cores = 3,selfing = FALSE,MAF = 0.01)
#data("RealDataSNPs_Y")
#gyg <- standardize(SNPs = RealDataSNPs_Y[,-1],number_cores = 3,selfing = TRUE,MAF = 0.01)
#kinship <- kinship_humans(SNPs = Humandata,number_cores = 3,MAF = 0.01)

SNP_data_function_pcp_selfing <- function(x,int,pcp){
  return(cbind(int,matrix(x,ncol = 1),pcp))
}

SNP_data_function_nopcp_selfing <- function(x,int){
  return(cbind(int,matrix(x,ncol = 1)))
}

SNP_data_function_nonselfing <- function(x){
  uniq_x <- unique(x)
  vec2 <- uniq_x[!(uniq_x %in% c(min(x,na.rm = TRUE),NA))]
  a <- length(vec2)
  X <- matrix(0,nrow = length(x),ncol = a)
  for(j in 1:a){
    X[x == vec2[j],j] <- 1
    X[is.na(x),j] <- .Internal(mean(X[,j]))
  }
  return(X)
}

##### RE

log_profile_likelihood_REML <- function(x,t,y,d){
  n <- nrow(x)
  p <- ncol(x)
  estar <- 1/(1 + t*d)
  xtey <- t(x)%*%(estar*y)
  xtex <- t(x)%*%(estar*x)
  beta.mle <- solve(xtex)%*%xtey
  top <- t(y - x%*%beta.mle)%*%(estar*(y - x%*%beta.mle))
  sigma.mle <- as.numeric(top/(n - p))
  as.numeric(-.5*(n - p)*log(2*pi*sigma.mle)  - .5 * sum(log(1/estar))- (1/(2*sigma.mle)) * top - .5*log(det(xtex)))
}

log_profile_likelihood_MLE <- function(x,t,y,d){
  n <- nrow(x)
  p <- ncol(x)
  estar <- 1/(1 + t*d)
  xtey <- t(x)%*%(estar*y)
  xtex <- t(x)%*%(estar*x)
  beta.mle <- solve(xtex)%*%xtey
  top <- t(y - x%*%beta.mle)%*%(estar*(y - x%*%beta.mle))
  sigma.mle <- as.numeric(top/(n))
  as.numeric(-.5*(n)*log(2*pi*sigma.mle)  - .5 * sum(log(1/estar))- (1/(2*sigma.mle)) * top)
}

RE_p <- function(x,y,d){
  n <- nrow(x)
  p <- ncol(x)
  L <- matrix(diag(p)[-1,],nrow = p - 1,ncol = p)
  z <- optimize(log_profile_likelihood_REML,interval = c(0,100),x = x, y = y,d = d,maximum = TRUE)$maximum
  estar <- 1/(1 + z*d)
  xtx <- solve(t(x)%*%(estar*x))
  beta2 <- xtx%*%t(x)%*%(estar*y)
  top <- t(y - x%*%beta2)%*%(estar*(y - x%*%beta2))
  sigma.reml <- as.numeric(top/(n - p))
  Fstat <- ((1/sigma.reml)*as.numeric(t(L%*%beta2)%*%solve(L%*%xtx%*%t(L))%*%L%*%beta2))/(nrow(L))
  return(pf(Fstat,nrow(L),n - p,lower.tail = FALSE))
}

RE_p_fixed <- function(x,y,d,P){
  n <- nrow(x)
  p <- ncol(x)
  L <- matrix(diag(p)[-c(1,(p - P + 1):p),],nrow = p - length(c(1,(p - P + 1):p)),ncol = p)
  z <- optimize(log_profile_likelihood_REML,interval = c(0,100),x = x, y = y,d = d,maximum = TRUE)$maximum
  estar <- 1/(1 + z*d)
  xtx <- solve(t(x)%*%(estar*x))
  beta2 <- xtx%*%t(x)%*%(estar*y)
  top <- t(y - x%*%beta2)%*%(estar*(y - x%*%beta2))
  sigma.reml <- as.numeric(top/(n - p))
  Fstat <- ((1/sigma.reml)*as.numeric(t(L%*%beta2)%*%solve(L%*%xtx%*%t(L))%*%L%*%beta2))/(nrow(L))
  return(p_value = pf(Fstat,nrow(L),n - p,lower.tail = FALSE))
}

RE_p_P3D <- function(x,y,d,t){
  n <- nrow(x)
  p <- ncol(x)
  L <- matrix(diag(p)[-1,],nrow = p - 1,ncol = p)
  estar <- 1/(1 + t*d)
  xtx <- solve(t(x)%*%(estar*x))
  beta2 <- xtx%*%t(x)%*%(estar*y)
  top <- t(y - x%*%beta2)%*%(estar*(y - x%*%beta2))
  sigma.reml <- as.numeric(top/(n - p))
  Fstat <- ((1/sigma.reml)*as.numeric(t(L%*%beta2)%*%solve(L%*%xtx%*%t(L))%*%L%*%beta2))/(nrow(L))
  return(p_value = pf(Fstat,nrow(L),n - p,lower.tail = FALSE))
}

RE_p_P3D_fixed <- function(x,y,d,P,t){
  n <- nrow(x)
  p <- ncol(x)
  L <- matrix(diag(p)[-c(1,(p - P + 1):p),],nrow = p - length(c(1,(p - P + 1):p)),ncol = p)
  estar <- 1/(1 + t*d)
  xtx <- solve(t(x)%*%(estar*x))
  beta2 <- xtx%*%t(x)%*%(estar*y)
  top <- t(y - x%*%beta2)%*%(estar*(y - x%*%beta2))
  sigma.reml <- as.numeric(top/(n - p))
  Fstat <- ((1/sigma.reml)*as.numeric(t(L%*%beta2)%*%solve(L%*%xtx%*%t(L))%*%L%*%beta2))/(nrow(L))
  return(p_value = pf(Fstat,nrow(L),n - p,lower.tail = FALSE))
}

RE_BIC <- function(x,y,d){
  n <- length(d)
  z <- as.numeric(optimize(log_profile_likelihood_MLE,interval = c(0,100),x = x, y = y,d = d,maximum = TRUE)$objective)
  BIC.cor <- -2 * z + ncol(x)*log(n)
  return(BIC.cor)
}

RE_BIC_P3D <- function(x,y,t,d){
  n <- length(d)
  z <- as.numeric(log_profile_likelihood_MLE(x = x,t = t, y = y,d = d))
  BIC.cor <- -2 * z + ncol(x)*log(n)
  return(BIC.cor)
}

RE_BIC_modelselection <- function(x,y,d,pstar){
  n <- length(d)
  z <- as.numeric(optimize(log_profile_likelihood_MLE,interval = c(0,100),x = x, y = y,d = d,maximum = TRUE)$objective)
  BIC.cor <- -2*z + pstar*log(n)
  return(BIC.cor)
}

SLR_p <- function(x,y){
  n <- nrow(x)
  p <- ncol(x)
  beta.mle <- solve(t(x)%*%x)%*%t(x)%*%y
  top <- sum((y - x%*%beta.mle)^2)
  sigma.mle <- top/(n-p)
  MSTr <- (sum((y - mean(y))^2) - top)/(p - 1)
  p_value <- pf((MSTr/sigma.mle),p - 1,n - p,lower.tail = FALSE)
  return(p_value)
}

SLR_p_fixed <- function(x,y,P,SS_pc){
  n <- nrow(x)
  p <- ncol(x)
  beta_all <- solve(t(x)%*%x)%*%t(x)%*%y
  SS_resids <- sum((y - x%*%beta_all)^2)
  p_value <- pf((((SS_pc - SS_resids)/(p - P - 1))/(SS_resids/(n - p))),(p - P) - 1,n - p,lower.tail = FALSE)
  return(p_value)
}

SLR_BIC <- function(x,y){
  n <- nrow(x)
  p <- ncol(x)
  xtx <- t(x)%*%x
  beta.mle <- solve(xtx)%*%t(x)%*%y
  top <- sum((y - x%*%beta.mle)^2)
  sigma.mle <- top/(n)
  BIC <- -2 *((-(n)/2)*log(2*pi*sigma.mle) - .5*(1/sigma.mle)*top) + (p)*log(n)
  return(BIC)
}

nullmodel_tau <- function(x,y,d){
  z <- optimize(log_profile_likelihood_REML,interval = c(0,100),x = x, y = y,d = d,maximum = TRUE)$maximum
  return(z)
}

Pval_function <- function(p_vals,thresh,control){
  p_vals <- as.numeric(p_vals)
  tf_mat <- matrix(NA,ncol = 2,nrow = length(p_vals))
  #Seeing if the SNP is significant given a multiple comparison error correction and the threshold
  tf_mat[,1] <- p.adjust(p_vals,control) < thresh
  tf_mat[,2] <- p_vals
  tf_mat <- as.data.frame(tf_mat)
  colnames(tf_mat) <- c("Significant","P_values")
  #Returning the p-values as well as if they are significant or not
  return(tf_mat)
}

svd1 <- function(x){
  p <- ncol(x)
  d <- svd(x,nu = p,nv = p)$d
  return(list(max(d)/min(d),p - length(d)))
}

preselection_nopc <- function(Y,X,number_cores,P3D,frequentist,controlrate,threshold,nullprob,alterprob,kinship,selfing){
  # These lines need to be run, this essentially cleans the data set
  requireNamespace("parallel")
  ####################
  
  ##################
  
  if(is.logical(kinship)){
    #Simple Linear Regression
    nX <- nrow(X)
    if(selfing){
      SNPdata_list <- split(t(X),1:ncol(X))
      rm(X)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_nopcp_selfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nopcp_selfing,int = matrix(1,ncol = 1,nrow = nX))
      stopCluster(cl)
    }else{
      SNPdata_list <- split(t(X),1:ncol(X))
      rm(X)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_nonselfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nonselfing)
      stopCluster(cl)
      SNPdata_list <- lapply(SNPdata_list,function(x){cbind(1,x)})
    }
    
    if(frequentist){
      #Frequentist SLR
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("optimize","solve"),envir=environment())
      SNPdata_list <- unlist(parLapply(cl,SNPdata_list,SLR_p, y = Y))
      stopCluster(cl)
      Pval_function(p_vals = SNPdata_list,thresh = threshold,control = controlrate)
    }else{
      #Bayesian SLR
      w <- matrix(1,nrow = nX,ncol = 1)
      BIC.null <- SLR_BIC(w,y = Y) - 2*log(nullprob)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("solve"),envir=environment())
      SNPdata_list <- unlist(parLapply(cl,SNPdata_list,SLR_BIC, y = Y))
      stopCluster(cl)
      SNPdata_list <- SNPdata_list - 2*1*log(1 - nullprob)
      #p_vec <- (alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + nullprob*exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
      p_vec <- (exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
      order_vec <- (1:length(p_vec))[order(p_vec,decreasing = TRUE)]
      p_vec <- p_vec[order(p_vec,decreasing = TRUE)]
      FDR <- cumsum(1 - p_vec)/(1:length(p_vec))
      FDR <- FDR[order(order_vec,decreasing = FALSE)]
      p_vec <- p_vec[order(order_vec,decreasing = FALSE)]
      tf_mat <- as.data.frame(cbind(FDR < threshold,p_vec))
      colnames(tf_mat) <- c("Significant","ApprPosteriorProbs")
      return(tf_mat)
    }
  }else{
    #Kinship
    #SLR with kinship component
    nX <- nrow(X)
    
    spec.decomp <- eigen(kinship,symmetric = TRUE)
    Qt <- t(spec.decomp$vectors)
    spec.decomp$values[length(spec.decomp$values)] <- 0
    D <- spec.decomp$values
    rm(spec.decomp)
    
    Y <- Qt%*%Y
    intercept <- Qt%*%matrix(1,nrow = nX,ncol = 1)
    
    if(selfing){
      X <- Qt%*%X
      SNPdata_list <- split(t(X),1:ncol(X))
      rm(X)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_nopcp_selfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nopcp_selfing,int = intercept)
      stopCluster(cl)
    }else{
      p <- ncol(X)
      SNPdata_list <- split(t(X),1:ncol(X))
      rm(X)
      level_function_new <- function(x){
        return(length(unique(x[!is.na(x)])))
      }
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("level_function_new"),envir=environment())
      lengths <- unlist(parLapply(cl,SNPdata_list,level_function_new)) - 1
      stopCluster(cl)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_nonselfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nonselfing)
      stopCluster(cl)
      
      SNPdata_list <- do.call(cbind,SNPdata_list)
      SNPdata_list <- Qt%*%SNPdata_list
      
      SNPdata_list <- split(t(SNPdata_list),rep(1:p,lengths))
      put_together <- function(i,x,leng,int){
        return(cbind(int,matrix(x[[i]],ncol = leng[[i]],byrow = TRUE)))
      }
      names(lengths) <- as.character(1:p)
      SNPdata_list <- lapply(names(SNPdata_list),put_together,x = SNPdata_list,leng = lengths,int = intercept)
    }
    
    #Frequentist RE model
    if(frequentist){
      if(P3D){
        t <- nullmodel_tau(x = intercept,y = Y,d = D)
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("diag","solve","RE_p_P3D"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,RE_p_P3D, y = Y, d = D,t = t))
        stopCluster(cl)
      }else{
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("optimize","diag","solve","log_profile_likelihood_REML","RE_p"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,RE_p, y = Y, d = D))
        stopCluster(cl)
      }
      
      Pval_function(p_vals = SNPdata_list,thresh = threshold,control = controlrate)
    }else{
      #Bayesian RE model
      BIC.null <- RE_BIC(x = intercept,y = Y,d = D) - 2*log(nullprob)
      
      if(P3D){
        t <- nullmodel_tau(x = intercept,y = Y,d = D)
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("diag","log_profile_likelihood_MLE","log_profile_likelihood_REML","RE_BIC_P3D"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,RE_BIC_P3D, y = Y, d = D,t = t))
        stopCluster(cl)
      }else{
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("solve","log_profile_likelihood_MLE","log_profile_likelihood_REML","optimize"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,RE_BIC, y = Y, d = D))
        stopCluster(cl)
      }
      SNPdata_list <- SNPdata_list - 2*1*log(1 - nullprob)
      p_vec <- (alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + nullprob*exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
      p_vec <- (exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
      order_vec <- (1:length(p_vec))[order(p_vec,decreasing = TRUE)]
      p_vec <- p_vec[order(p_vec,decreasing = TRUE)]
      FDR <- cumsum(1 - p_vec)/(1:length(p_vec))
      FDR <- FDR[order(order_vec,decreasing = FALSE)]
      p_vec <- p_vec[order(order_vec,decreasing = FALSE)]
      tf_mat <- as.data.frame(cbind(FDR < threshold,p_vec))
      colnames(tf_mat) <- c("Significant","ApprPosteriorProbs")
      return(tf_mat)
      #tf_mat <- as.data.frame(cbind(1:length(SNPdata_list),SNPdata_list))
      #colnames(tf_mat) <- c("Index","BICs")
      #return(tf_mat)
    }
  }
}


preselection_pc <- function(Y,X,number_cores,P3D,principal_components,frequentist,controlrate,threshold,nullprob,alterprob,kinship,selfing){
  # These lines need to be run, this essentially cleans the data set
  requireNamespace("parallel")
  ####################
  if(is.null(dim(principal_components))){
    principal_components[is.na(principal_components)] <- mean(principal_components,na.rm = TRUE)
  }else{
    for(i in 1:ncol(principal_components)){
      principal_components[is.na(principal_components[,i]),i] <- mean(principal_components[,i],na.rm = TRUE)
    }
  }
  ##################
  
  if(is.logical(kinship)){
    #Simple Linear Regression
    nX <- nrow(X)
    if(selfing){
      SNPdata_list <- split(t(X),1:ncol(X))
      rm(X)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_pcp_selfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_pcp_selfing,int = matrix(1,ncol = 1,nrow = nX),pcp = principal_components)
      stopCluster(cl)
    }else{
      SNPdata_list <- split(t(X),1:ncol(X))
      rm(X)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_nonselfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nonselfing)
      stopCluster(cl)
      
      SNPdata_list <- lapply(SNPdata_list,function(x){cbind(1,x,principal_components)})
    }
    
    if(frequentist){
      #Frequentist SLR
      null_dat <- cbind(1,principal_components)
      Beta_pc <- solve(t(null_dat)%*%null_dat)%*%t(null_dat)%*%Y
      SS_pc <- sum((Y - null_dat%*%Beta_pc)^2)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("optimize","solve"),envir=environment())
      SNPdata_list <- unlist(parLapply(cl,SNPdata_list,SLR_p_fixed, y = Y,P = ncol(principal_components),SS_pc = SS_pc))
      stopCluster(cl)
      Pval_function(p_vals = SNPdata_list,thresh = threshold,control = controlrate)
    }else{
      #Bayesian SLR
      BIC.null <- SLR_BIC(cbind(1,principal_components),y = Y) - 2*log(nullprob)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("solve"),envir=environment())
      SNPdata_list <- unlist(parLapply(cl,SNPdata_list,SLR_BIC, y = Y))
      stopCluster(cl)
      SNPdata_list <- SNPdata_list - 2*1*log(1 - nullprob)
      #p_vec <- (alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + nullprob*exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
      p_vec <- (exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
      order_vec <- (1:length(p_vec))[order(p_vec,decreasing = TRUE)]
      p_vec <- p_vec[order(p_vec,decreasing = TRUE)]
      FDR <- cumsum(1 - p_vec)/(1:length(p_vec))
      FDR <- FDR[order(order_vec,decreasing = FALSE)]
      p_vec <- p_vec[order(order_vec,decreasing = FALSE)]
      tf_mat <- as.data.frame(cbind(FDR < threshold,p_vec))
      colnames(tf_mat) <- c("Significant","ApprPosteriorProbs")
      return(tf_mat)
    }
  }else{
    #Kinship
    #SLR with kinship component
    nX <- nrow(X)
    
    spec.decomp <- eigen(kinship,symmetric = TRUE)
    Qt <- t(spec.decomp$vectors)
    spec.decomp$values[length(spec.decomp$values)] <- 0
    D <- spec.decomp$values
    rm(spec.decomp)
    
    Y <- Qt%*%Y
    intercept <- Qt%*%matrix(1,ncol = 1,nrow = nX)
    principal_components <- Qt%*%principal_components
    
    if(selfing){
      X <- Qt%*%X
      SNPdata_list <- split(t(X),1:ncol(X))
      rm(X)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_pcp_selfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_pcp_selfing,int = intercept,pcp = principal_components)
      stopCluster(cl)
    }else{
      p <- ncol(X)
      SNPdata_list <- split(t(X),1:ncol(X))
      rm(X)
      level_function_new <- function(x){
        return(length(unique(x[!is.na(x)])))
      }
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("level_function_new"),envir=environment())
      lengths <- unlist(parLapply(cl,SNPdata_list,level_function_new)) - 1
      stopCluster(cl)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_nonselfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nonselfing)
      stopCluster(cl)
      
      SNPdata_list <- do.call(cbind,SNPdata_list)
      SNPdata_list <- Qt%*%SNPdata_list
      SNPdata_list <- split(t(SNPdata_list),rep(1:p,lengths))
      put_together <- function(i,x,leng,int,pcp){
        return(cbind(int,matrix(x[[i]],ncol = leng[[i]],byrow = TRUE),pcp))
      }
      names(lengths) <- as.character(1:p)
      SNPdata_list <- lapply(names(SNPdata_list),put_together,x = SNPdata_list,leng = lengths,int = intercept,pcp = principal_components)
    }
    #Frequentist RE model
    if(frequentist){
      if(P3D){
        t <- nullmodel_tau(x = cbind(intercept,principal_components),y = Y,d = D)
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("diag","solve","RE_p_P3D_fixed"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,RE_p_P3D_fixed, y = Y, d = D,t = t,P = ncol(principal_components)))
        stopCluster(cl)
      }else{
        t <- nullmodel_tau(x = cbind(intercept,principal_components),y = Y,d = D)
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("optimize","diag","solve","log_profile_likelihood_REML","RE_p_fixed"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,RE_p_fixed, y = Y, d = D, P = ncol(principal_components)))
        stopCluster(cl)
      }
      Pval_function(p_vals = SNPdata_list,thresh = threshold,control = controlrate)
    }else{
      #Bayesian RE model
      BIC.null <- RE_BIC(cbind(intercept,principal_components),y = Y,d = D) - 2*log(nullprob)
      if(P3D){
        t <- nullmodel_tau(x = cbind(intercept,principal_components),y = Y,d = D)
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("diag","log_profile_likelihood_MLE","log_profile_likelihood_REML","RE_BIC_P3D"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,RE_BIC_P3D, y = Y, d = D,t = t))
        stopCluster(cl)
      }else{
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("solve","log_profile_likelihood_MLE","log_profile_likelihood_REML","optimize"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,RE_BIC, y = Y, d = D))
        stopCluster(cl)
      }
      SNPdata_list <- SNPdata_list - 2*1*log(1 - nullprob)
      #p_vec <- (alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + nullprob*exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
      p_vec <- (exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
      order_vec <- (1:length(p_vec))[order(p_vec,decreasing = TRUE)]
      p_vec <- p_vec[order(p_vec,decreasing = TRUE)]
      FDR <- cumsum(1 - p_vec)/(1:length(p_vec))
      FDR <- FDR[order(order_vec,decreasing = FALSE)]
      p_vec <- p_vec[order(order_vec,decreasing = FALSE)]
      tf_mat <- as.data.frame(cbind(FDR < threshold,p_vec))
      colnames(tf_mat) <- c("Significant","ApprPosteriorProbs")
      return(tf_mat)
    }
    
  }
  
}

SMA <- function(Y,SNPs,number_cores = 1,P3D = TRUE,fixed_effects = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,nullprob = NULL,alterprob = NULL,kinship = FALSE,info = FALSE,selfing){
  
  principal_components <- fixed_effects
  
  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(!(is.numeric(Y))){
    stop("Y is not a numeric vector")
  }
  if(!is.numeric(threshold)){
    stop("threshold is not numeric")
  }
  if(threshold > 1 | threshold < 0){
    stop("threshold needs to be between 0 and 1")
  }
  if(!is.logical(frequentist)){
    stop("frequentist needs to be a logical value")
  }
  if(!is.numeric(number_cores)){
    stop("number_cores is not numeric")
  }
  if(!(controlrate %in% p.adjust.methods)){
    stop("control rate needs to be one of p.adjust.methods")
  }
  if(!is.null(nullprob)){
    if(!is.numeric(nullprob)){
      stop("nullprob needs to be numeric")
    }
    if(!is.numeric(alterprob)){
      stop("alterprob needs to be numeric")
    }
    if(alterprob > 1 | alterprob < 0){
      stop("nullprob needs to be between 0 and 1")
    }
    if(alterprob > 1 | alterprob < 0){
      stop("alterprob needs to be between 0 and 1")
    }
    if((nullprob + alterprob) != 1){
      stop("nullprob and alterprob need to sum to 1")
    }
  }
  
  if(is.logical(info)){
    if(is.logical(principal_components)){
      preselection_nopc(Y = Y,X = SNPs,number_cores = number_cores,P3D = P3D,frequentist = frequentist,controlrate = controlrate,threshold = threshold,nullprob = nullprob,alterprob = alterprob,kinship = kinship,selfing = selfing)
    }else{
      preselection_pc(Y = Y,X = SNPs,number_cores = number_cores,P3D = P3D,principal_components = principal_components,frequentist = frequentist,controlrate = controlrate,threshold = threshold,nullprob = nullprob,alterprob = alterprob,kinship = kinship,selfing = selfing)
    }
  }else{
    if(is.logical(principal_components)){
      results <- preselection_nopc(Y = Y,X = SNPs,number_cores = number_cores,P3D = P3D,frequentist = frequentist,controlrate = controlrate,threshold = threshold,nullprob = nullprob,alterprob = alterprob,kinship = kinship,selfing = selfing)
    }else{
      results <- preselection_pc(Y = Y,X = SNPs,number_cores = number_cores,P3D = P3D,principal_components = principal_components,frequentist = frequentist,controlrate = controlrate,threshold = threshold,nullprob = nullprob,alterprob = alterprob,kinship = kinship,selfing = selfing)
    }
    results <- cbind(t(info),results)
    results[,1] <- as.numeric(as.character(results[,1]))
    results[,2] <- as.numeric(as.character(results[,2]))
    if(frequentist){
      colnames(results) <- c("Chromosomes","Positions","Significant","P_values")
    }else{
      colnames(results) <- c("Chromosomes","Positions","Significant","ApprPosteriorProbs")
    }
    return(results)
  }
}

#library(GWAS.BAYES)

#Y <- RealDataSNPs_Y$Phenotype
#pc <- cbind(rnorm(328),rnorm(328))
#Y = Y; X = gyg$SNPs;principal_components = pc;number_cores = 3;P3D = FALSE;frequentist = TRUE;controlrate = "bonferroni";threshold = 0.05;selfing = TRUE;kinship = RealDataKinship
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = pc,number_cores = 3,P3D = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = TRUE,kinship = RealDataKinship)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = pc,number_cores = 3,P3D = FALSE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = TRUE,kinship = RealDataKinship)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = pc,number_cores = 3,P3D = TRUE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = TRUE,kinship = RealDataKinship)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = pc,number_cores = 3,P3D = TRUE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = TRUE,kinship = RealDataKinship)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = pc,number_cores = 3,P3D = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = TRUE,kinship = FALSE)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = pc,number_cores = 3,P3D = FALSE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = TRUE,kinship = FALSE)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = pc,number_cores = 3,P3D = TRUE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = TRUE,kinship = FALSE)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = pc,number_cores = 3,P3D = TRUE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = TRUE,kinship = FALSE)

#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = FALSE,number_cores = 3,P3D = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = TRUE,kinship = RealDataKinship)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = FALSE,number_cores = 3,P3D = FALSE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = TRUE,kinship = RealDataKinship)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = FALSE,number_cores = 3,P3D = TRUE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = TRUE,kinship = RealDataKinship)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = FALSE,number_cores = 3,P3D = TRUE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = TRUE,kinship = RealDataKinship)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = FALSE,number_cores = 3,P3D = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = TRUE,kinship = FALSE)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = FALSE,number_cores = 3,P3D = FALSE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = TRUE,kinship = FALSE)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = FALSE,number_cores = 3,P3D = TRUE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = TRUE,kinship = FALSE)
#gyt <- preselection(Y = Y, SNPs = gyg$SNPs,principal_components = FALSE,number_cores = 3,P3D = TRUE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = TRUE,kinship = FALSE)

#library(rrBLUP)
#pc <- cbind(rnorm(500),rnorm(500))
#Y = rnorm(500); X = myt$SNPs;principal_components = pc;number_cores = 3;P3D = FALSE;frequentist = TRUE;controlrate = "bonferroni";threshold = 0.05;selfing = FALSE;kinship = kinship
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = pc,number_cores = 3,P3D = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = FALSE,kinship = kinship)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = pc,number_cores = 3,P3D = FALSE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = FALSE,kinship = kinship)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = pc,number_cores = 3,P3D = TRUE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = FALSE,kinship = kinship)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = pc,number_cores = 3,P3D = TRUE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = FALSE,kinship = kinship)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = pc,number_cores = 3,P3D = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = FALSE,kinship = FALSE)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = pc,number_cores = 3,P3D = FALSE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = FALSE,kinship = FALSE)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = pc,number_cores = 3,P3D = TRUE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = FALSE,kinship = FALSE)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = pc,number_cores = 3,P3D = TRUE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = FALSE,kinship = FALSE)

#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = FALSE,number_cores = 3,P3D = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = FALSE,kinship = kinship)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = FALSE,number_cores = 3,P3D = FALSE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = FALSE,kinship = kinship)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = FALSE,number_cores = 3,P3D = TRUE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = FALSE,kinship = kinship)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = FALSE,number_cores = 3,P3D = TRUE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = FALSE,kinship = kinship)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = FALSE,number_cores = 3,P3D = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = FALSE,kinship = FALSE)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = FALSE,number_cores = 3,P3D = FALSE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = FALSE,kinship = FALSE)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = FALSE,number_cores = 3,P3D = TRUE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,selfing = FALSE,kinship = FALSE)
#gyt <- preselection(Y = rnorm(500), SNPs = myt$SNPs,principal_components = FALSE,number_cores = 3,P3D = TRUE,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = 0.05,selfing = FALSE,kinship = FALSE)


ga_modelselection_nopc <- function(Y,X,significant,number_cores,maxiterations,runs_til_stop,kinship,selfing,pi0){
  requireNamespace("GA")
  requireNamespace("parallel")
  requireNamespace("doParallel")
  requireNamespace("Matrix")
  requireNamespace("memoise")
  
  y1 <- paste0("SNP",1:ncol(X))[significant]
  
  #SLR w Kinship
  nX <- nrow(X)
  X <- X[,significant]
  
  spec.decomp <- eigen(kinship,symmetric = TRUE)
  Q <- spec.decomp$vectors
  Qt <- t(Q)
  spec.decomp$values[length(spec.decomp$values)] <- 0
  D <- spec.decomp$values
  rm(spec.decomp)
  
  intercept <- Qt%*%matrix(1,nrow = nX,ncol = 1)
  Y <- Qt%*%Y
  rm(Q)
  
  if(selfing){
    X <- Qt%*%X
    fitness_sub <- function(string) {
      SNPdata_list_sub <- cbind(intercept,X[,string == 1])
      P <- ncol(SNPdata_list_sub)
      
      if(Matrix::rankMatrix(SNPdata_list_sub)[1] < ncol(SNPdata_list_sub)){
        dropped_cols <- findLinearCombos(SNPdata_list_sub)$remove
        SNPdata_list_sub <- SNPdata_list_sub[,-dropped_cols]
      }
      return((-1)*(RE_BIC_modelselection(SNPdata_list_sub,Y,pstar = P,d = D) - 2*(P*log(1 - pi0) + (ncol(X) - P)*log(pi0))))
    }
  }else{
    fitness_sub <- function(string) {
      SNPdata_list_sub <- matrix(1,ncol = 1,nrow = nX)
      if(sum(string) != 0){
        for(i in 1:sum(string == 1)){
          SNPdata_list_sub <- cbind(SNPdata_list_sub,SNP_data_function_nonselfing(x = X[,which(string == 1)[i]]))
        }
      }
      P <- ncol(SNPdata_list_sub)
      
      if(Matrix::rankMatrix(SNPdata_list_sub)[1] < ncol(SNPdata_list_sub)){
        dropped_cols <- findLinearCombos(SNPdata_list_sub)$remove
        SNPdata_list_sub <- SNPdata_list_sub[,-dropped_cols]
      }
      return((-1)*(RE_BIC_modelselection(Qt%*%SNPdata_list_sub,Y,pstar = P,d = D) - 2*(P*log(1 - pi0) + (ncol(X) - P)*log(pi0))))
    }
  }
  
  suggestedsol <- rbind(0,diag(length(significant)))
  if(length(significant) > 99){
    if(selfing){
      SNPdata_list <- split(t(X),1:ncol(X))
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_nopcp_selfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nopcp_selfing,int = intercept)
      stopCluster(cl)
    }else{
      p <- ncol(X)
      SNPdata_list <- split(t(X),1:ncol(X))
      level_function_new <- function(x){
        return(length(unique(x[!is.na(x)])))
      }
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("level_function_new"),envir=environment())
      lengths <- unlist(parLapply(cl,SNPdata_list,level_function_new)) - 1
      stopCluster(cl)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_nonselfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nonselfing)
      stopCluster(cl)
      
      SNPdata_list <- do.call(cbind,SNPdata_list)
      SNPdata_list <- Qt%*%SNPdata_list
      SNPdata_list <- split(t(SNPdata_list),rep(1:p,lengths))
      put_together <- function(i,x,leng,int){
        return(cbind(int,matrix(x[[i]],ncol = leng[[i]],byrow = TRUE)))
      }
      names(lengths) <- as.character(1:p)
      SNPdata_list <- lapply(names(SNPdata_list),put_together,x = SNPdata_list,leng = lengths,int = intercept)
    }
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("solve","log_profile_likelihood_MLE","log_profile_likelihood_REML","optimize"),envir=environment())
    SNPdata_list <- unlist(parLapply(cl,SNPdata_list,RE_BIC, y = Y, d = D))
    stopCluster(cl)
    suggestedsol <- as.matrix(suggestedsol[order(SNPdata_list)[1:99],])
  }
  fitness_sub <- memoise::memoise(fitness_sub)
  ans <- GA::ga("binary", fitness = fitness_sub, nBits = ncol(X),names = y1,maxiter = maxiterations,popSize = 100,elitism = min(c(10,2^ncol(X))),parallel = 1,run = runs_til_stop,suggestions = suggestedsol,monitor = FALSE)
  memoise::forget(fitness_sub)
  dat <- cbind(ans@population,(-1)*ans@fitness)
  dat <- unique(dat)
  dat <- cbind(dat,(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))/(sum(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))))
  vec1 <- matrix(0,nrow = nrow(dat),ncol = nrow(ans@solution))
  for(i in 1:nrow(dat)){
    for(j in 1:nrow(ans@solution)){
      vec1[i,j] <- sum(dat[i,1:(ncol(dat) - 2)] == ans@solution[j,])
    }
  }
  vec <- vector()
  for(i in 1:nrow(ans@solution)){
    vec[i] <- paste0("Model ",i,": y ~ ",paste(names(which(ans@solution[i,] == 1)),collapse = " + ")," Posterior Prob. of ",round(dat[which(vec1[,i] == (ncol(dat) - 2)),ncol(dat)],digits = 4))
  }
  inclusion_prob <- cbind(y1,t(dat[,1:(ncol(dat) - 2)])%*%dat[,ncol(dat)]/sum(dat[,ncol(dat)]))
  inclusion_prob <- as.data.frame(inclusion_prob)
  colnames(inclusion_prob) <- c("SNPs","InclusionProb")
  inclusion_prob$InclusionProb <- as.numeric(inclusion_prob$InclusionProb)
  rownames(inclusion_prob) <- NULL
  
  return(list(Models = vec,Solution = ans@solution,InclusionProb = inclusion_prob))
}

ga_modelselection_pc <- function(Y,X,significant,number_cores,principal_components,maxiterations,runs_til_stop,kinship = FALSE,selfing,pi0){
  
  requireNamespace("GA")
  requireNamespace("parallel")
  requireNamespace("doParallel")
  requireNamespace("Matrix")
  requireNamespace("memoise")
  
  if(is.null(dim(principal_components))){
    principal_components[is.na(principal_components)] <- mean(principal_components,na.rm = TRUE)
  }else{
    for(i in 1:ncol(principal_components)){
      principal_components[is.na(principal_components[,i]),i] <- mean(principal_components[,i],na.rm = TRUE)
    }
  }
  
  y1 <- paste0("SNP",1:length(ncol(X)))[significant]
  
  #SLR w Kinship
  nX <- nrow(X)
  X <- X[,significant]
  
  spec.decomp <- eigen(kinship,symmetric = TRUE)
  Qt <- t(spec.decomp$vectors)
  spec.decomp$values[length(spec.decomp$values)] <- 0
  D <- spec.decomp$values
  rm(spec.decomp)
  
  intercept <- Qt%*%matrix(1,nrow = nX,ncol = 1)
  principal_components <- Qt%*%principal_components
  Y <- Qt%*%Y
  rm(Q)
  
  if(selfing){
    X <- Qt%*%X
    fitness_sub <- function(string) {
      SNPdata_list_sub <- cbind(intercept,X[,string == 1],principal_components)
      P <- ncol(SNPdata_list_sub)
      
      if(Matrix::rankMatrix(SNPdata_list_sub)[1] < ncol(SNPdata_list_sub)){
        dropped_cols <- findLinearCombos(SNPdata_list_sub)$remove
        SNPdata_list_sub <- SNPdata_list_sub[,-dropped_cols]
      }
      return((-1)*(RE_BIC_modelselection(SNPdata_list_sub,Y,pstar = P,d = D) - 2*(P*log(1 - pi0) + (ncol(X) - P)*log(pi0))))
    }
  }else{
    fitness_sub <- function(string) {
      SNPdata_list_sub <- matrix(1,ncol = 1,nrow = nX)
      if(sum(string) != 0){
        for(i in 1:sum(string == 1)){
          SNPdata_list_sub <- cbind(SNPdata_list_sub,SNP_data_function_nonselfing(x = X[,which(string == 1)[i]]))
        }
      }
      SNPdata_list_sub <- cbind(Qt%*%SNPdata_list_sub,principal_components)
      P <- ncol(SNPdata_list_sub)
      
      if(Matrix::rankMatrix(SNPdata_list_sub)[1] < ncol(SNPdata_list_sub)){
        dropped_cols <- findLinearCombos(SNPdata_list_sub)$remove
        SNPdata_list_sub <- SNPdata_list_sub[,-dropped_cols]
      }
      return((-1)*(RE_BIC_modelselection(SNPdata_list_sub,Y,pstar = P,d = D) - 2*(P*log(1 - pi0) + (ncol(X) - P)*log(pi0))))
    }
  }
  
  suggestedsol <- rbind(0,diag(length(significant)))
  if(length(significant) > 99){
    if(selfing){
      SNPdata_list <- split(t(X),1:ncol(X))
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_pcp_selfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_pcp_selfing,int = intercept,pcp = principal_components)
      stopCluster(cl)
    }else{
      p <- ncol(X)
      SNPdata_list <- split(t(X),1:ncol(X))
      level_function_new <- function(x){
        return(length(unique(x[!is.na(x)])))
      }
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("level_function_new"),envir=environment())
      lengths <- unlist(parLapply(cl,SNPdata_list,level_function_new)) - 1
      stopCluster(cl)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_nonselfing"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nonselfing)
      stopCluster(cl)
      
      SNPdata_list <- do.call(cbind,SNPdata_list)
      SNPdata_list <- Qt%*%SNPdata_list
      SNPdata_list <- split(t(SNPdata_list),rep(1:p,lengths))
      put_together <- function(i,x,leng,int,pcp){
        return(cbind(int,matrix(x[[i]],ncol = leng[[i]],byrow = TRUE)))
      }
      names(lengths) <- as.character(1:p)
      SNPdata_list <- lapply(names(SNPdata_list),put_together,x = SNPdata_list,leng = lengths,int = intercept,pcp = principal_components)
    }
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("solve","log_profile_likelihood_MLE","log_profile_likelihood_REML","optimize"),envir=environment())
    SNPdata_list <- unlist(parLapply(cl,SNPdata_list,RE_BIC, y = Y, d = D))
    stopCluster(cl)
    suggestedsol <- as.matrix(suggestedsol[order(SNPdata_list)[1:99],])
  }
  fitness_sub <- memoise::memoise(fitness_sub)
  ans <- GA::ga("binary", fitness = fitness_sub, nBits = ncol(X),names = y1,maxiter = maxiterations,popSize = 100,elitism = min(c(10,2^ncol(X))),parallel = number_cores,run = runs_til_stop,suggestions = suggestedsol,monitor = FALSE)
  memoise::forget(fitness_sub)
  dat <- cbind(ans@population,(-1)*ans@fitness)
  dat <- unique(dat)
  dat <- cbind(dat,(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))/(sum(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))))
  vec1 <- matrix(0,nrow = nrow(dat),ncol = nrow(ans@solution))
  for(i in 1:nrow(dat)){
    for(j in 1:nrow(ans@solution)){
      vec1[i,j] <- sum(dat[i,1:(ncol(dat) - 2)] == ans@solution[j,])
    }
  }
  vec <- vector()
  for(i in 1:nrow(ans@solution)){
    vec[i] <- paste0("Model ",i,": y ~ ",paste(names(which(ans@solution[i,] == 1)),collapse = " + ")," Posterior Prob. of ",round(dat[which(vec1[,i] == (ncol(dat) - 2)),ncol(dat)],digits = 4))
  }
  inclusion_prob <- cbind(y1,t(dat[,1:(ncol(dat) - 2)])%*%dat[,ncol(dat)]/sum(dat[,ncol(dat)]))
  inclusion_prob <- as.data.frame(inclusion_prob)
  colnames(inclusion_prob) <- c("SNPs","InclusionProb")
  inclusion_prob$InclusionProb <- as.numeric(inclusion_prob$InclusionProb)
  rownames(inclusion_prob) <- NULL
  
  return(list(Models = vec,Solution = ans@solution,InclusionProb = inclusion_prob))
}

postGWAS <- function(Y,SNPs,significant,number_cores = 1,principal_components = FALSE,maxiterations = 100,runs_til_stop = 10,kinship,info = FALSE,selfing,pi0){
  
  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(!(is.numeric(Y))){
    stop("Y is not a numeric vector")
  }
  if(!(is.numeric(significant))){
    stop("significant is not a numeric vector")
  }
  if(!is.numeric(number_cores)){
    stop("number_cores is not numeric")
  }
  if(!is.numeric(maxiterations)){
    stop("maxiterations needs to be numeric")
  }
  if(!is.numeric(runs_til_stop)){
    stop("runs_til_stop needs to be numeric")
  }
  
  if(length(significant) > 11){
    if(is.logical(principal_components)){
      results <- ga_modelselection_nopc(Y = Y,X = SNPs,significant = significant,number_cores = number_cores,maxiterations = maxiterations,runs_til_stop = runs_til_stop,kinship = kinship,selfing = selfing,pi0 = pi0)
    }else{
      results <- ga_modelselection_pc(Y = Y,X = SNPs,significant = significant,number_cores = number_cores,principal_components = principal_components,maxiterations = maxiterations,runs_til_stop = runs_til_stop,kinship = kinship,selfing = selfing,pi0 = pi0)
    }
    if(is.logical(info)){
      return(results)
    }else{
      info <- info[,significant]
      info <- paste0(info[1,],"-",info[2,])
      colnames(results$Solution) <- info
      results$InclusionProb$SNPs <- info
      return(results)
    }
  }else{
    if(is.logical(principal_components)){
      requireNamespace("GA")
      requireNamespace("parallel")
      requireNamespace("doParallel")
      requireNamespace("Matrix")
      requireNamespace("memoise")
      
      X <- SNPs
      y1 <- paste0("SNP",1:ncol(X))[significant]
      
      #SLR w Kinship
      
      nX <- nrow(X)
      X <- matrix(X[,significant],ncol = length(significant))
      
      spec.decomp <- eigen(kinship,symmetric = TRUE)
      Qt <- t(spec.decomp$vectors)
      spec.decomp$values[length(spec.decomp$values)] <- 0
      D <- spec.decomp$values
      rm(spec.decomp)
      
      intercept <- Qt%*%matrix(1,nrow = nX,ncol = 1)
      Y <- Qt%*%Y
      n <- ncol(X)
      
      if(selfing){
        X <- Qt%*%X
        fitness_sub <- function(string) {
          SNPdata_list_sub <- cbind(intercept,X[,string == 1])
          P <- ncol(SNPdata_list_sub)
          
          if(Matrix::rankMatrix(SNPdata_list_sub)[1] < ncol(SNPdata_list_sub)){
            dropped_cols <- findLinearCombos(SNPdata_list_sub)$remove
            SNPdata_list_sub <- SNPdata_list_sub[,-dropped_cols]
          }
          return((-1)*(RE_BIC_modelselection(SNPdata_list_sub,Y,pstar = P,d = D) - 2*(P*log(1 - pi0) + (ncol(X) - P)*log(pi0))))
        }
      }else{
        fitness_sub <- function(string) {
          SNPdata_list_sub <- matrix(1,ncol = 1,nrow = nX)
          if(sum(string) != 0){
            for(i in 1:sum(string == 1)){
              SNPdata_list_sub <- cbind(SNPdata_list_sub,SNP_data_function_nonselfing(x = X[,which(string == 1)[i]]))
            }
          }
          P <- ncol(SNPdata_list_sub)
          
          if(Matrix::rankMatrix(SNPdata_list_sub)[1] < ncol(SNPdata_list_sub)){
            dropped_cols <- findLinearCombos(SNPdata_list_sub)$remove
            SNPdata_list_sub <- SNPdata_list_sub[,-dropped_cols]
          }
          return((-1)*(RE_BIC_modelselection(Qt%*%SNPdata_list_sub,Y,pstar = P,d = D) - 2*(P*log(1 - pi0) + (ncol(X) - P)*log(pi0))))
        }
      }
      
      l <- rep(list(0:1), n)
      l <- expand.grid(l)
      colnames(l) <- y1
      ol <- as.matrix(l)
      l <- split(l,1:nrow(l))
      
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("l","fitness_sub","X","nX","D","Qt","pi0","RE_BIC_modelselection","log_profile_likelihood_MLE","log_profile_likelihood_REML","SNP_data_function_nonselfing","rankMatrix","findLinearCombos","Y"),envir=environment())
      l <- unlist(parLapply(cl,l,fitness_sub))
      stopCluster(cl)
      
      dat <- cbind(ol,l)
      dat[,ncol(dat)] <- (-1)*dat[,ncol(dat)]
      dat <- cbind(dat,(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))/(sum(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))))
      
      dat <- dat[order(dat[,ncol(dat)],decreasing = FALSE),]
      anstf <- cumsum(dat[,ncol(dat)]) > .05
      ans <- matrix(dat[anstf,1:(ncol(dat) - 2)],ncol = (ncol(dat) - 2),nrow = sum(anstf))
      colnames(ans) <- y1
      
      vec1 <- matrix(0,nrow = nrow(dat),ncol = nrow(ans))
      for(i in 1:nrow(dat)){
        for(j in 1:nrow(ans)){
          vec1[i,j] <- sum(dat[i,1:(ncol(dat) - 2)] == ans[j,])
        }
      }
      vec <- vector()
      for(i in 1:nrow(ans)){
        vec[i] <- paste0("Model ",i,": y ~ ",paste(names(which(ans[i,] == 1)),collapse = " + ")," Posterior Prob. of ",round(dat[which(vec1[,i] == (ncol(dat) - 2)),ncol(dat)],digits = 4))
      }
      inclusion_prob <- cbind(y1,t(dat[,1:(ncol(dat) - 2)])%*%dat[,ncol(dat)]/sum(dat[,ncol(dat)]))
      inclusion_prob <- as.data.frame(inclusion_prob)
      colnames(inclusion_prob) <- c("SNPs","InclusionProb")
      inclusion_prob$InclusionProb <- as.numeric(inclusion_prob$InclusionProb)
      rownames(inclusion_prob) <- NULL
      
      results <- list(Models = vec,Solution = ans,InclusionProb = inclusion_prob)
      if(is.logical(info)){
        return(results)
      }else{
        info <- info[,significant]
        info <- paste0(info[1,],"-",info[2,])
        colnames(results$Solution) <- info
        results$InclusionProb$SNPs <- info
        return(results)
      }
      
    }else{
      requireNamespace("GA")
      requireNamespace("parallel")
      requireNamespace("doParallel")
      requireNamespace("Matrix")
      requireNamespace("memoise")
      
      if(is.null(dim(principal_components))){
        principal_components[is.na(principal_components)] <- mean(principal_components,na.rm = TRUE)
      }else{
        for(i in 1:ncol(principal_components)){
          principal_components[is.na(principal_components[,i]),i] <- mean(principal_components[,i],na.rm = TRUE)
        }
      }
      
      X <- SNPs
      y1 <- paste0("SNP",1:ncol(X))[significant]
      
      #SLR w Kinship
      nX <- nrow(X)
      X <- matrix(X[,significant],ncol = length(significant))
      
      spec.decomp <- eigen(kinship,symmetric = TRUE)
      Qt <- t(spec.decomp$vectors)
      spec.decomp$values[length(spec.decomp$values)] <- 0
      D <- spec.decomp$values
      rm(spec.decomp)
      
      intercept <- Qt%*%matrix(1,nrow = nX,ncol = 1)
      principal_components <- Qt%*%principal_components
      Y <- Qt%*%Y
      n <- ncol(X)
      
      if(selfing){
        X <- Qt%*%X
        fitness_sub <- function(string) {
          SNPdata_list_sub <- cbind(intercept,X[,string == 1],principal_components)
          P <- ncol(SNPdata_list_sub)
          
          if(Matrix::rankMatrix(SNPdata_list_sub)[1] < ncol(SNPdata_list_sub)){
            dropped_cols <- findLinearCombos(SNPdata_list_sub)$remove
            SNPdata_list_sub <- SNPdata_list_sub[,-dropped_cols]
          }
          return((-1)*(RE_BIC_modelselection(SNPdata_list_sub,Y,pstar = P,d = D) - 2*(P*log(1 - pi0) + (ncol(X) - P)*log(pi0))))
        }
      }else{
        fitness_sub <- function(string) {
          SNPdata_list_sub <- matrix(1,ncol = 1,nrow = nX)
          if(sum(string) != 0){
            for(i in 1:sum(string == 1)){
              SNPdata_list_sub <- cbind(SNPdata_list_sub,SNP_data_function_nonselfing(x = X[,which(string == 1)[i]]))
            }
          }
          SNPdata_list_sub <- cbind(Qt%*%SNPdata_list_sub,principal_components)
          P <- ncol(SNPdata_list_sub)
          
          if(Matrix::rankMatrix(SNPdata_list_sub)[1] < ncol(SNPdata_list_sub)){
            dropped_cols <- findLinearCombos(SNPdata_list_sub)$remove
            SNPdata_list_sub <- SNPdata_list_sub[,-dropped_cols]
          }
          return((-1)*(RE_BIC_modelselection(SNPdata_list_sub,Y,pstar = P,d = D) - 2*(P*log(1 - pi0) + (ncol(X) - P)*log(pi0))))
        }
      }
      
      l <- rep(list(0:1), n)
      l <- expand.grid(l)
      colnames(l) <- y1
      ol <- as.matrix(l)
      l <- split(l,1:nrow(l))
      
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("l","fitness_sub","X","nX","pi0","D","Qt","RE_BIC_modelselection","SNP_data_function_nonselfing","log_profile_likelihood_MLE","log_profile_likelihood_REML","rankMatrix","findLinearCombos","Y","principal_components"),envir=environment())
      l <- unlist(parLapply(cl,l,fitness_sub))
      stopCluster(cl)
      
      dat <- cbind(ol,l)
      dat[,ncol(dat)] <- (-1)*dat[,ncol(dat)]
      dat <- cbind(dat,(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))/(sum(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))))
      
      dat <- dat[order(dat[,ncol(dat)],decreasing = FALSE),]
      anstf <- cumsum(dat[,ncol(dat)]) > .05
      ans <- matrix(dat[anstf,1:(ncol(dat) - 2)],ncol = (ncol(dat) - 2),nrow = sum(anstf))
      colnames(ans) <- y1
      
      vec1 <- matrix(0,nrow = nrow(dat),ncol = nrow(ans))
      for(i in 1:nrow(dat)){
        for(j in 1:nrow(ans)){
          vec1[i,j] <- sum(dat[i,1:(ncol(dat) - 2)] == ans[j,])
        }
      }
      vec <- vector()
      for(i in 1:nrow(ans)){
        vec[i] <- paste0("Model ",i,": y ~ ",paste(names(which(ans[i,] == 1)),collapse = " + ")," Posterior Prob. of ",round(dat[which(vec1[,i] == (ncol(dat) - 2)),ncol(dat)],digits = 4))
      }
      inclusion_prob <- cbind(y1,t(dat[,1:(ncol(dat) - 2)])%*%dat[,ncol(dat)]/sum(dat[,ncol(dat)]))
      inclusion_prob <- as.data.frame(inclusion_prob)
      colnames(inclusion_prob) <- c("SNPs","InclusionProb")
      inclusion_prob$InclusionProb <- as.numeric(inclusion_prob$InclusionProb)
      rownames(inclusion_prob) <- NULL
      
      results <- list(Models = vec,Solution = ans,InclusionProb = inclusion_prob)
      if(is.logical(info)){
        return(results)
      }else{
        info <- info[,significant]
        info <- paste0(info[1,],"-",info[2,])
        colnames(results$Solution) <- info
        results$InclusionProb$SNPs <- info
        return(results)
      }
    }
  }
}

#pc <- cbind(rnorm(328),rnorm(328))

#Y <- RealDataSNPs_Y$Phenotype
#sig1 <- which(rbinom(ncol(gyg$SNPs),1,.05) == 1)
#sig1
#gyt <- postGWAS(Y = Y,SNPs = gyg$SNPs,significant = sig1,number_cores = 3,kinship = RealDataKinship,maxiterations = 100,runs_til_stop = 10,selfing = TRUE)
#gyt <- postGWAS(Y = Y,SNPs = gyg$SNPs,significant = which(rbinom(ncol(gyg$SNPs),1,.01) == 1),number_cores = 3,kinship = FALSE,maxiterations = 100,runs_til_stop = 10,selfing = TRUE)

#gyt <- postGWAS(Y = Y,SNPs = gyg$SNPs,significant = which(rbinom(ncol(gyg$SNPs),1,.01) == 1),principal_components = pc,number_cores = 3,kinship = RealDataKinship,maxiterations = 100,runs_til_stop = 10,selfing = TRUE)
#sig2 <- which(rbinom(ncol(gyg$SNPs),1,.01) == 1)
#gyt <- postGWAS(Y = Y,SNPs = gyg$SNPs,significant = sig2,principal_components = pc,number_cores = 3,kinship = FALSE,maxiterations = 100,runs_til_stop = 10,selfing = TRUE)

#pc <- cbind(rnorm(500),rnorm(500))

#sig1 <- which(rbinom(ncol(myt$SNPs),1,.001) == 1)
#sig1
#gyt <- postGWAS(Y = rnorm(500),SNPs = myt$SNPs,significant = sig1,number_cores = 3,kinship = kinship,maxiterations = 100,runs_til_stop = 10,selfing = FALSE)
#sig2 <- which(rbinom(ncol(myt$SNPs),1,.001) == 1)
#sig2
#gyt <- postGWAS(Y = rnorm(500),SNPs = myt$SNPs,significant = sig2,number_cores = 3,kinship = FALSE,maxiterations = 100,runs_til_stop = 10,selfing = FALSE)

#sig3 <- which(rbinom(ncol(myt$SNPs),1,.001) == 1)
#sig3
#gyt <- postGWAS(Y = rnorm(500),SNPs = myt$SNPs,significant = sig3,principal_components = pc,number_cores = 3,kinship = kinship,maxiterations = 100,runs_til_stop = 10,selfing = FALSE)
#sig4 <- which(rbinom(ncol(myt$SNPs),1,.001) == 1)
#sig4
#gyt <- postGWAS(Y = rnorm(500),SNPs = myt$SNPs,significant = sig4,principal_components = pc,number_cores = 3,kinship = FALSE,maxiterations = 100,runs_til_stop = 10,selfing = FALSE)


#X <- gyg$SNPs
#number_cores = 3
#fixed <- c(450,900)


preselection_nopc_2 <- function(Y,X,fixed,number_cores,frequentist,controlrate,threshold,nullprob,alterprob,kinship = FALSE,P3D,selfing){
  requireNamespace("parallel")
  ####################
  
  dropped_cols <- NULL
  
  fixed_o <- fixed
  
  if(!selfing){
    X_fixed <- matrix(1,nrow = nrow(X),ncol = 1)
    for(i in 1:length(fixed)){
      X_fixed <- cbind(X_fixed,SNP_data_function_nonselfing(x = X[,fixed[i]]))
    }
    X_fixed <- matrix(X_fixed[,-1],ncol = ncol(X_fixed) - 1)
  }else{
    X_fixed <- matrix(X[,fixed],ncol = length(fixed))
  }
  
  if(Matrix::rankMatrix(X_fixed)[1] < ncol(X_fixed)){
    dropped_cols <- findLinearCombos(X_fixed)$remove
    X_fixed <- matrix(X_fixed[,-dropped_cols],ncol = Matrix::rankMatrix(X_fixed)[1])
    fixed <- fixed[-dropped_cols]
  }
  
  ##################
  #Kinship
  #SLR with kinship component
  nX <- nrow(X)
  
  spec.decomp <- eigen(kinship,symmetric = TRUE)
  Qt <- t(spec.decomp$vectors)
  spec.decomp$values[length(spec.decomp$values)] <- 0
  D <- spec.decomp$values
  rm(spec.decomp)
  
  intercept <- Qt%*%matrix(1,nrow = nX,ncol = 1)
  Y <- Qt%*%Y
  X_fixed <- Qt%*%X_fixed
  
  if(selfing){
    X <- Qt%*%X
    SNPdata_list <- split(t(X),1:ncol(X))
    rm(X)
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("SNP_data_function_pcp_selfing"),envir=environment())
    SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_pcp_selfing,int = intercept,pcp = X_fixed)
    stopCluster(cl)
  }else{
    p <- ncol(X)
    SNPdata_list <- split(t(X),1:ncol(X))
    rm(X)
    level_function_new <- function(x){
      return(length(unique(x[!is.na(x)])))
    }
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("level_function_new"),envir=environment())
    lengths <- unlist(parLapply(cl,SNPdata_list,level_function_new)) - 1
    stopCluster(cl)
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("SNP_data_function_nonselfing"),envir=environment())
    SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nonselfing)
    stopCluster(cl)
    
    SNPdata_list <- do.call(cbind,SNPdata_list)
    SNPdata_list <- Qt%*%SNPdata_list
    
    SNPdata_list <- split(t(SNPdata_list),rep(1:p,lengths))
    put_together <- function(i,x,leng,int,X_fix){
      return(cbind(int,matrix(x[[i]],ncol = leng[[i]],byrow = TRUE),X_fix))
    }
    names(lengths) <- as.character(1:p)
    SNPdata_list <- lapply(names(SNPdata_list),put_together,x = SNPdata_list,leng = lengths,int = intercept,X_fix = X_fixed)
  }
  
  SNPdata_list[fixed] <- replicate(n = length(fixed), expr = {cbind(intercept,X_fixed)},simplify = F)
  
  l1 <- unlist(lapply(SNPdata_list, svd1))
  rank_defient <- l1[1:length(l1) %% 2 == 0]
  condition_num <- l1[1:length(l1) %% 2 == 1]
  
  if(selfing){
    rank_defient <- which(rank_defient > 0)
    condition_num <- which(condition_num > 1e6)
    SNPdata_list[rank_defient] <- replicate(n = length(rank_defient), expr = {cbind(intercept,X_fixed)},simplify = F)
    SNPdata_list[condition_num] <- replicate(n = length(condition_num), expr = {cbind(intercept,X_fixed)},simplify = F)
    exclude <- c(fixed,rank_defient,condition_num)
  }else{
    rank_defient <- which(rank_defient > 0)
    if(length(rank_defient) > 0){
      for(i in 1:length(rank_defient)){
        if(lengths[rank_defient[i]] < 2){
          SNPdata_list[[rank_defient[i]]] <- cbind(intercept,X_fixed)
        }else{
          SNPdata_list[[rank_defient[i]]] <- SNPdata_list[[rank_defient[i]]][,-3]
        }
      }
    }
    condition_num <- which(condition_num > 1e6)
    if(length(condition_num) > 0){
      for(i in 1:length(condition_num)){
        if(lengths[condition_num[i]] < 2){
          SNPdata_list[[condition_num[i]]] <- cbind(intercept,X_fixed)
        }else{
          SNPdata_list[[condition_num[i]]] <- SNPdata_list[[condition_num[i]]][,-3]
        }
      }
    }
    exclude <- c(fixed,rank_defient[lengths[rank_defient] < 2],condition_num[lengths[condition_num] < 2])
  }
  
  #Frequentist RE model
  if(frequentist){
    if(P3D){
      q <- nullmodel_tau(x = cbind(intercept,X_fixed),y = Y,d = D)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("diag","solve","RE_p_P3D_fixed"),envir=environment())
      SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]] <- parLapply(cl,SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]],RE_p_P3D_fixed, y = Y, d = D,t = q,P = ncol(X_fixed))
      stopCluster(cl)
      SNPdata_list[exclude] <- 1
      SNPdata_list <- unlist(SNPdata_list)
    }else{
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("optimize","diag","solve","log_profile_likelihood_REML","RE_p_fixed"),envir=environment())
      SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]] <- parLapply(cl,SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]],RE_p_fixed, y = Y, d = D, P = ncol(X_fixed))
      stopCluster(cl)
      SNPdata_list[exclude] <- 1
      SNPdata_list <- unlist(SNPdata_list)
    }
    
    tf_mat <- Pval_function(p_vals = SNPdata_list[(1:length(SNPdata_list))[!(1:length(SNPdata_list) %in% exclude)]],thresh = threshold,control = controlrate)
    rownames(tf_mat) <- (1:length(SNPdata_list))[!(1:length(SNPdata_list) %in% exclude)]
    Significant <- tf_mat[tf_mat[,1] == 1,]
    return(list(fixed = fixed_o,Significant = Significant,All_P_Values = tf_mat))
  }else{
    #Bayesian RE model
    BIC.null <- RE_BIC(x = cbind(intercept,X_fixed),y = Y,d = D) - 2*log(nullprob)
    
    if(P3D){
      t <- nullmodel_tau(x = cbind(intercept,X_fixed),y = Y,d = D)
      cl <- makeCluster(1)
      clusterExport(cl,c("diag","log_profile_likelihood_MLE","log_profile_likelihood_REML","RE_BIC_P3D"),envir=environment())
      SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]] <- parLapply(cl,SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]],RE_BIC_P3D, y = Y, d = D,t = t)
      stopCluster(cl)
      SNPdata_list[exclude] <- BIC.null
      SNPdata_list <- unlist(SNPdata_list)
      SNPdata_list <- SNPdata_list - 2*1*log(1 - nullprob)
    }else{
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("solve","log_profile_likelihood_MLE","log_profile_likelihood_REML","optimize"),envir=environment())
      SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]] <- parLapply(cl,SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]],RE_BIC, y = Y, d = D)
      stopCluster(cl)
      SNPdata_list[exclude] <- BIC.null
      SNPdata_list <- unlist(SNPdata_list)
      SNPdata_list <- SNPdata_list - 2*1*log(1 - nullprob)
    }
    p_vec <- (exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
    order_vec <- (1:length(p_vec))[order(p_vec,decreasing = TRUE)]
    p_vec <- p_vec[order(p_vec,decreasing = TRUE)]
    FDR <- cumsum(1 - p_vec)/(1:length(p_vec))
    FDR <- FDR[order(order_vec,decreasing = FALSE)]
    p_vec <- p_vec[order(order_vec,decreasing = FALSE)]
    tf_mat <- as.data.frame(cbind(FDR < threshold,p_vec))
    colnames(tf_mat) <- c("Significant","ApprPosteriorProbs")
    Significant <- tf_mat[tf_mat[,1] == 1,]
    return(list(fixed = fixed_o,Significant = Significant,approxposterior = tf_mat))
  }
}

#Y <- RealDataSNPs_Y$Phenotype

#gyt <- preselection_nopc_2(Y = Y, X = gyg$SNPs,fixed = c(450,900),number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = RealDataKinship,P3D = FALSE,selfing = TRUE)
#gyt <- preselection_nopc_2(Y = Y, X = gyg$SNPs,fixed = c(450,900),number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = FALSE,P3D = FALSE,selfing = TRUE)

#gyt <- preselection_nopc_2(Y = Y, X = gyg$SNPs,fixed = c(450,900),number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = RealDataKinship,P3D = TRUE,selfing = TRUE)
#gyt <- preselection_nopc_2(Y = Y, X = gyg$SNPs,fixed = c(450,900),number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = FALSE,P3D = TRUE,selfing = TRUE)

#gyt <- preselection_nopc_2(Y = Y, X = gyg$SNPs,fixed = c(450,900),number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = RealDataKinship,P3D = FALSE,selfing = TRUE)
#gyt <- preselection_nopc_2(Y = Y, X = gyg$SNPs,fixed = c(450,900),number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = FALSE,P3D = FALSE,selfing = TRUE)

#gyt <- preselection_nopc_2(Y = Y, X = gyg$SNPs,fixed = c(450,900),number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = RealDataKinship,P3D = TRUE,selfing = TRUE)
#gyt <- preselection_nopc_2(Y = Y, X = gyg$SNPs,fixed = c(450,900),number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = FALSE,P3D = TRUE,selfing = TRUE)


#gyt <- preselection_nopc_2(Y = rnorm(500), X = myt$SNPs,fixed = c(450,900),number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = kinship,P3D = FALSE,selfing = FALSE)
#gyt <- preselection_nopc_2(Y = rnorm(500), X = myt$SNPs,fixed = c(450,900),number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = FALSE,P3D = FALSE,selfing = FALSE)

#gyt <- preselection_nopc_2(Y = rnorm(500), X = myt$SNPs,fixed = c(450,900),number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = kinship,P3D = TRUE,selfing = FALSE)
#gyt <- preselection_nopc_2(Y = rnorm(500), X = myt$SNPs,fixed = c(450,900),number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = FALSE,P3D = TRUE,selfing = FALSE)

#gyt <- preselection_nopc_2(Y = rnorm(500), X = myt$SNPs,fixed = c(450,900),number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = kinship,P3D = FALSE,selfing = FALSE)
#gyt <- preselection_nopc_2(Y = rnorm(500), X = myt$SNPs,fixed = c(450,900),number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = FALSE,P3D = FALSE,selfing = FALSE)

#gyt <- preselection_nopc_2(Y = rnorm(500), X = myt$SNPs,fixed = c(450,900),number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = kinship,P3D = TRUE,selfing = FALSE)
#gyt <- preselection_nopc_2(Y = rnorm(500), X = myt$SNPs,fixed = c(450,900),number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = FALSE,P3D = TRUE,selfing = FALSE)


preselection_pc_2 <- function(Y,X,fixed,number_cores,principal_components,frequentist,controlrate,threshold,nullprob,alterprob,kinship = FALSE,P3D,selfing){
  # These lines need to be run, this essentially cleans the data set
  requireNamespace("parallel")
  ####################
  
  dropped_cols <- NULL
  
  fixed_o <- fixed
  
  if(!selfing){
    X_fixed <- matrix(1,nrow = nrow(X),ncol = 1)
    for(i in 1:length(fixed)){
      X_fixed <- cbind(X_fixed,SNP_data_function_nonselfing(x = X[,fixed[i]]))
    }
    X_fixed <- matrix(X_fixed[,-1],ncol = ncol(X_fixed) - 1)
  }else{
    X_fixed <- matrix(X[,fixed],ncol = length(fixed))
  }
  
  if(Matrix::rankMatrix(X_fixed)[1] < ncol(X_fixed)){
    dropped_cols <- findLinearCombos(X_fixed)$remove
    X_fixed <- matrix(X_fixed[,-dropped_cols],ncol = Matrix::rankMatrix(X_fixed)[1])
    fixed <- fixed[-dropped_cols]
  }
  
  if(is.null(dim(principal_components))){
    principal_components[is.na(principal_components)] <- mean(principal_components,na.rm = TRUE)
  }else{
    for(i in 1:ncol(principal_components)){
      principal_components[is.na(principal_components[,i]),i] <- mean(principal_components[,i],na.rm = TRUE)
    }
  }
  
  ##################
  #Kinship
  #SLR with kinship component
  nX <- nrow(X)
  
  spec.decomp <- eigen(kinship,symmetric = TRUE)
  Qt <- t(spec.decomp$vectors)
  spec.decomp$values[length(spec.decomp$values)] <- 0
  D <- spec.decomp$values
  rm(spec.decomp)
  
  intercept <- Qt%*%matrix(1,nrow = nX,ncol = 1)
  X_fixed <- Qt%*%X_fixed
  principal_components <- Qt%*%principal_components
  Y <- Qt%*%Y
  
  if(selfing){
    X <- Qt%*%X
    SNPdata_list <- split(t(X),1:ncol(X))
    rm(X)
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("SNP_data_function_pcp_selfing"),envir=environment())
    SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_pcp_selfing,int = intercept,pcp = cbind(X_fixed,principal_components))
    stopCluster(cl)
  }else{
    p <- ncol(X)
    SNPdata_list <- split(t(X),1:ncol(X))
    rm(X)
    level_function_new <- function(x){
      return(length(unique(x[!is.na(x)])))
    }
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("level_function_new"),envir=environment())
    lengths <- unlist(parLapply(cl,SNPdata_list,level_function_new)) - 1
    stopCluster(cl)
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("SNP_data_function_nonselfing"),envir=environment())
    SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nonselfing)
    stopCluster(cl)
    
    SNPdata_list <- do.call(cbind,SNPdata_list)
    SNPdata_list <- Qt%*%SNPdata_list
    
    SNPdata_list <- split(t(SNPdata_list),rep(1:p,lengths))
    put_together <- function(i,x,leng,int,X_fix){
      return(cbind(int,matrix(x[[i]],ncol = leng[[i]],byrow = TRUE),X_fix))
    }
    names(lengths) <- as.character(1:p)
    SNPdata_list <- lapply(names(SNPdata_list),put_together,x = SNPdata_list,leng = lengths,int = intercept,X_fix = cbind(X_fixed,principal_components))
  }
  SNPdata_list[fixed] <- replicate(n = length(fixed), expr = {cbind(intercept,X_fixed,principal_components)},simplify = F)
  
  l1 <- unlist(lapply(SNPdata_list, svd1))
  rank_defient <- l1[1:length(l1) %% 2 == 0]
  condition_num <- l1[1:length(l1) %% 2 == 1]
  
  if(selfing){
    rank_defient <- which(rank_defient > 0)
    condition_num <- which(condition_num > 1e6)
    SNPdata_list[rank_defient] <- replicate(n = length(rank_defient), expr = {cbind(intercept,X_fixed,principal_components)},simplify = F)
    SNPdata_list[condition_num] <- replicate(n = length(condition_num), expr = {cbind(intercept,X_fixed,principal_components)},simplify = F)
    exclude <- c(fixed,rank_defient,condition_num)
  }else{
    rank_defient <- which(rank_defient > 0)
    if(length(rank_defient) > 0){
      for(i in 1:length(rank_defient)){
        if(lengths[rank_defient[i]] < 2){
          SNPdata_list[[rank_defient[i]]] <- cbind(intercept,X_fixed,principal_components)
        }else{
          SNPdata_list[[rank_defient[i]]] <- SNPdata_list[[rank_defient[i]]][,-3]
        }
      }
    }
    condition_num <- which(condition_num > 1e6)
    if(length(condition_num) > 0){
      for(i in 1:length(condition_num)){
        if(lengths[condition_num[i]] < 2){
          SNPdata_list[[condition_num[i]]] <- cbind(intercept,X_fixed,principal_components)
        }else{
          SNPdata_list[[condition_num[i]]] <- SNPdata_list[[condition_num[i]]][,-3]
        }
      }
    }
    exclude <- c(fixed,rank_defient[lengths[rank_defient] < 2],condition_num[lengths[condition_num] < 2])
  }
  
  #Frequentist RE model
  if(frequentist){
    if(P3D){
      t <- nullmodel_tau(x = Qt%*%cbind(intercept,X_fixed,principal_components),y = Y,d = D)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("diag","solve","RE_p_P3D_fixed"),envir=environment())
      SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]] <- parLapply(cl,SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]],RE_p_P3D_fixed, y = Y, d = D,t = t,P = ncol(X_fixed) + ncol(principal_components))
      stopCluster(cl)
      SNPdata_list[exclude] <- 1
      SNPdata_list <- unlist(SNPdata_list)
    }else{
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("optimize","diag","solve","log_profile_likelihood_REML","RE_p_fixed"),envir=environment())
      SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]] <- parLapply(cl,SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]],RE_p_fixed, y = Y, d = D, P = ncol(X_fixed) + ncol(principal_components))
      stopCluster(cl)
      SNPdata_list[exclude] <- 1
      SNPdata_list <- unlist(SNPdata_list)
    }
    
    tf_mat <- Pval_function(p_vals = SNPdata_list[(1:length(SNPdata_list))[!(1:length(SNPdata_list) %in% exclude)]],thresh = threshold,control = controlrate)
    rownames(tf_mat) <- (1:length(SNPdata_list))[!(1:length(SNPdata_list) %in% exclude)]
    Significant <- tf_mat[tf_mat[,1] == 1,]
    return(list(fixed = fixed_o,Significant = Significant,All_P_Values = tf_mat))
  }else{
    #Bayesian RE model
    BIC.null <- RE_BIC(cbind(intercept,X_fixed,principal_components),y = Y,d = D) - 2*log(nullprob)
    if(P3D){
      t <- nullmodel_tau(x = Qt%*%cbind(intercept,X_fixed,principal_components),y = Y,d = D)
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("diag","log_profile_likelihood_MLE","log_profile_likelihood_REML","RE_BIC_P3D"),envir=environment())
      SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]] <- parLapply(cl,SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]],RE_BIC_P3D, y = Y, d = D,t = t)
      stopCluster(cl)
      SNPdata_list[exclude] <- BIC.null
      SNPdata_list <- unlist(SNPdata_list)
      SNPdata_list <- SNPdata_list - 2*1*log(1 - nullprob)
    }else{
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("solve","log_profile_likelihood_MLE","log_profile_likelihood_REML","optimize"),envir=environment())
      SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]] <- parLapply(cl,SNPdata_list[(1:length(SNPdata_list))[!((1:length(SNPdata_list)) %in% exclude)]],RE_BIC, y = Y, d = D)
      stopCluster(cl)
      SNPdata_list[exclude] <- BIC.null
      SNPdata_list <- unlist(SNPdata_list)
      SNPdata_list <- SNPdata_list - 2*1*log(1 - nullprob)
    }
    
    p_vec <- (exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
    order_vec <- (1:length(p_vec))[order(p_vec,decreasing = TRUE)]
    p_vec <- p_vec[order(p_vec,decreasing = TRUE)]
    FDR <- cumsum(1 - p_vec)/(1:length(p_vec))
    FDR <- FDR[order(order_vec,decreasing = FALSE)]
    p_vec <- p_vec[order(order_vec,decreasing = FALSE)]
    tf_mat <- as.data.frame(cbind(FDR < threshold,p_vec))
    colnames(tf_mat) <- c("Significant","ApprPosteriorProbs")
    Significant <- tf_mat[tf_mat[,1] == 1,]
    return(list(fixed = fixed_o,Significant = Significant,approxposterior = tf_mat))
  }
}

#pc <- cbind(rnorm(328),rnorm(328))

#Y <- RealDataSNPs_Y$Phenotype

#gyt <- preselection_pc_2(Y = Y, X= gyg$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = RealDataKinship,P3D = FALSE,selfing = TRUE)
#gyt <- preselection_pc_2(Y = Y, X= gyg$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = FALSE,P3D = FALSE,selfing = TRUE)

#gyt <- preselection_pc_2(Y = Y, X= gyg$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = RealDataKinship,P3D = TRUE,selfing = TRUE)
#gyt <- preselection_pc_2(Y = Y, X= gyg$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = FALSE,P3D = TRUE,selfing = TRUE)

#gyt <- preselection_pc_2(Y = Y, X= gyg$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = RealDataKinship,P3D = FALSE,selfing = TRUE)
#gyt <- preselection_pc_2(Y = Y, X= gyg$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = FALSE,P3D = FALSE,selfing = TRUE)

#gyt <- preselection_pc_2(Y = Y, X= gyg$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = RealDataKinship,P3D = TRUE,selfing = TRUE)
#gyt <- preselection_pc_2(Y = Y, X= gyg$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = FALSE,P3D = TRUE,selfing = TRUE)

#pc <- cbind(rnorm(500),rnorm(500))

#gyt <- preselection_pc_2(Y = rnorm(500), X= myt$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = kinship,P3D = FALSE,selfing = FALSE)
#gyt <- preselection_pc_2(Y = rnorm(500), X= myt$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = FALSE,P3D = FALSE,selfing = FALSE)

#gyt <- preselection_pc_2(Y = rnorm(500), X= myt$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = kinship,P3D = TRUE,selfing = FALSE)
#gyt <- preselection_pc_2(Y = rnorm(500), X= myt$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = FALSE,P3D = TRUE,selfing = FALSE)

#gyt <- preselection_pc_2(Y = rnorm(500), X= myt$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = kinship,P3D = FALSE,selfing = FALSE)
#gyt <- preselection_pc_2(Y = rnorm(500), X= myt$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = FALSE,P3D = FALSE,selfing = FALSE)

#gyt <- preselection_pc_2(Y = rnorm(500), X= myt$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = kinship,P3D = TRUE,selfing = FALSE)
#gyt <- preselection_pc_2(Y = rnorm(500), X= myt$SNPs,fixed = c(450,900),principal_components = pc,number_cores = 3,frequentist = FALSE,nullprob = .5,alterprob = .5,threshold = .05,kinship = FALSE,P3D = TRUE,selfing = FALSE)


preselection2 <- function(Y,SNPs,fixed,number_cores = 1,principal_components = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,nullprob = NULL,alterprob = NULL,kinship = FALSE,info = FALSE,P3D,selfing){
  
  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(!(is.numeric(Y))){
    stop("Y is not a numeric vector")
  }
  if(!is.logical(kinship)){
    if(!((nrow(kinship) == nrow(SNPs)) & (ncol(kinship) == nrow(SNPs)))){
      stop("kinship does not have the dimensions nrow(SNPs) x nrow(SNPs)")
    }
  }else{
    if(kinship == TRUE){
      stop("kinship can only take on the values of FALSE or a matrix nrow(SNPs) x nrow(SNPs)")
    }
  }
  if(!is.logical(principal_components)){
    if(!(nrow(principal_components) == nrow(SNPs))){
      stop("principal_components does not have the the same number of rows as SNPs")
    }
    if(sum(apply(principal_components,2,is.numeric)) != ncol(principal_components)){
      stop("Not every column of principal_components is numeric")
    }
  }else{
    if(principal_components == TRUE){
      stop("principal_components can only take on the values of FALSE or a matrix with the same number of rows as SNPs")
    }
  }
  if(!is.logical(info)){
    if(!((nrow(info) == 2) & (ncol(info) == ncol(SNPs)))){
      stop("info does not have the dimensions 2 x ncol(SNPs)")
    }
  }else{
    if(info == TRUE){
      stop("info can only take on the values of FALSE or a matrix 2 x ncol(SNPs)")
    }
  }
  if(!is.numeric(threshold)){
    stop("threshold is not numeric")
  }
  if(threshold > 1 | threshold < 0){
    stop("threshold needs to be between 0 and 1")
  }
  if(!is.logical(frequentist)){
    stop("frequentist needs to be a logical value")
  }
  if(!is.numeric(number_cores)){
    stop("number_cores is not numeric")
  }
  if(!(controlrate %in% p.adjust.methods)){
    stop("control rate needs to be one of p.adjust.methods")
  }
  if(!is.null(nullprob)){
    if(!is.numeric(nullprob)){
      stop("nullprob needs to be numeric")
    }
    if(!is.numeric(alterprob)){
      stop("alterprob needs to be numeric")
    }
    if(alterprob > 1 | alterprob < 0){
      stop("nullprob needs to be between 0 and 1")
    }
    if(alterprob > 1 | alterprob < 0){
      stop("alterprob needs to be between 0 and 1")
    }
    if((nullprob + alterprob) != 1){
      stop("nullprob and alterprob need to sum to 1")
    }
  }
  
  if(is.logical(principal_components)){
    preselection_nopc_2(Y = Y,X = SNPs,fixed = fixed,number_cores = number_cores,frequentist = frequentist,controlrate = controlrate,threshold = threshold,nullprob = nullprob,alterprob = alterprob,kinship = kinship,P3D = P3D,selfing = selfing)
  }else{
    preselection_pc_2(Y = Y,X = SNPs,fixed = fixed,number_cores = number_cores,principal_components = principal_components,frequentist = frequentist,controlrate = controlrate,threshold = threshold,nullprob = nullprob,alterprob = alterprob,kinship = kinship,P3D = P3D,selfing = selfing)
  }
}

BICOSS <- function(Y,SNPs,number_cores = 1,fixed_effects = FALSE,threshold = 0.05,kinship = FALSE,info = FALSE,maxiterations = 100,runs_til_stop = 10,P3D = TRUE,selfing = FALSE){
  principal_components <- fixed_effects
  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(!(is.numeric(Y))){
    stop("Y is not a numeric vector")
  }
  if(!is.logical(kinship)){
    if(!((nrow(kinship) == nrow(SNPs)) & (ncol(kinship) == nrow(SNPs)))){
      stop("kinship does not have the dimensions nrow(SNPs) x nrow(SNPs)")
    }
  }else{
    if(kinship == TRUE){
      stop("kinship can only take on the values of FALSE or a matrix nrow(SNPs) x nrow(SNPs)")
    }
  }
  if(!is.logical(principal_components)){
    if(!(nrow(principal_components) == nrow(SNPs))){
      stop("principal_components does not have the the same number of rows as SNPs")
    }
    if(sum(apply(principal_components,2,is.numeric)) != ncol(principal_components)){
      stop("Not every column of principal_components is numeric")
    }
  }else{
    if(principal_components == TRUE){
      stop("principal_components can only take on the values of FALSE or a matrix with the same number of rows as SNPs")
    }
  }
  if(!is.logical(info)){
    if(!((nrow(info) == 2) & (ncol(info) == ncol(SNPs)))){
      stop("info does not have the dimensions 2 x ncol(SNPs)")
    }
  }else{
    if(info == TRUE){
      stop("info can only take on the values of FALSE or a matrix 2 x ncol(SNPs)")
    }
  }
  if(!is.numeric(threshold)){
    stop("threshold is not numeric")
  }
  if(threshold > 1 | threshold < 0){
    stop("threshold needs to be between 0 and 1")
  }
  if(!is.numeric(number_cores)){
    stop("number_cores is not numeric")
  }
  if(!is.numeric(maxiterations)){
    stop("maxiterations needs to be numeric")
  }
  if(!is.numeric(runs_til_stop)){
    stop("runs_til_stop needs to be numeric")
  }
  
  prescreen1_freq <- SMA(Y,SNPs,number_cores = number_cores,fixed_effects = principal_components,frequentist = TRUE,controlrate = "bonferroni",threshold = threshold,kinship = kinship,info = FALSE,P3D = P3D,selfing = selfing)
  true_nullprop <- convest(prescreen1_freq$P_values)
  if(true_nullprop == 1){
    true_nullprop <- 1 - 100/(ncol(SNPs))
  }
  prescreen1 <- SMA(Y,SNPs,number_cores = number_cores,fixed_effects = principal_components,nullprob = true_nullprop,alterprob = 1 - true_nullprop,frequentist = FALSE,threshold = threshold,kinship = kinship,info = info,P3D = P3D,selfing = selfing)
  
  if(sum(prescreen1$Significant == 1) == 0){
    modelselection <- "No significant in prescreen1"
    prescreen_2nd <- "No significant in prescreen1"
    prescreen2_freq <- prescreen1_freq
  }else{
    significant <- which(prescreen1$Significant == 1)
    p <- length(significant)
    modelselection_tmp <- list()
    prescreen2_tmp <- list()
    prescreen2_freq <- list()
    count <- 1
    prescreen2_freq[[count]] <- prescreen1_freq
    while((p > 0) & (count < 10)){
      modelselection_tmp[[count]] <- postGWAS(Y,SNPs,significant = significant,number_cores = number_cores,principal_components,maxiterations = maxiterations,runs_til_stop = runs_til_stop,kinship = kinship,selfing = selfing,pi0 = true_nullprop)
      if(sum(dim(modelselection_tmp[[count]]$Solution) == c(1,1)) == 2){
        fixed <- as.numeric(gsub("SNP","",colnames(modelselection_tmp[[count]]$Solution)))
      }else if(nrow(modelselection_tmp[[count]]$Solution) == 1){
        fixed <- as.numeric(gsub("SNP","",colnames(modelselection_tmp[[count]]$Solution)[which(modelselection_tmp[[count]]$Solution == 1)]))
      }else{
        nums <- vector()
        for(i in 1:nrow(modelselection_tmp[[count]]$Solution)){
          nums[i] <- as.numeric(strsplit(modelselection_tmp[[count]]$Models, " ")[[i]][length(strsplit(modelselection_tmp[[count]]$Models, " ")[[i]])])
        }
        modelselection_tmp[[count]]$Solution <- cbind(modelselection_tmp[[count]]$Solution,nums)
        modelselection_tmp[[count]]$Solution <- modelselection_tmp[[count]]$Solution[which(modelselection_tmp[[count]]$Solution[,ncol(modelselection_tmp[[count]]$Solution)] == max(modelselection_tmp[[count]]$Solution[,ncol(modelselection_tmp[[count]]$Solution)],na.rm = TRUE)),]
        if(!is.matrix(modelselection_tmp[[count]]$Solution)){
          fixed <- modelselection_tmp[[count]]$Solution[which(modelselection_tmp[[count]]$Solution[-length(modelselection_tmp[[count]]$Solution)] == max(modelselection_tmp[[count]]$Solution[-length(modelselection_tmp[[count]]$Solution)],na.rm = TRUE))]
          if(is.matrix(fixed)){
            fixed <- as.numeric(gsub("SNP","",names(which(fixed[1,] == 1))))
          }else{
            fixed <- as.numeric(gsub("SNP","",names(which(fixed == 1))))
          }
        }else{
          fixed <- modelselection_tmp[[count]]$Solution[which(rowSums(modelselection_tmp[[count]]$Solution[,-ncol(modelselection_tmp[[count]]$Solution)]) == max(rowSums(modelselection_tmp[[count]]$Solution[,-ncol(modelselection_tmp[[count]]$Solution)]),na.rm = TRUE)),]
          if(is.matrix(fixed)){
            fixed <- as.numeric(gsub("SNP","",names(which(fixed[1,] == 1))))
          }else{
            fixed <- as.numeric(gsub("SNP","",names(which(fixed == 1))))
          }
        }
      }
      prescreen2_freq[[count + 1]] <- preselection2(SNPs = SNPs,Y = Y,fixed = fixed,kinship = kinship,frequentist = TRUE,controlrate = "bonferroni",threshold = threshold,number_cores = number_cores,P3D = P3D,selfing = selfing)
      true_nullprop2 <- convest(prescreen2_freq[[count + 1]]$All_P_Values$P_values)
      if(true_nullprop2 == 1){
        approx_post <- matrix(0,nrow = ncol(SNPs),ncol = 2)
        colnames(approx_post) <- c("Significant","ApprPosteriorProbs")
        approx_post <- as.data.frame(approx_post)
        sig1 <- approx_post[approx_post$Significant == 2,]
        prescreen2_tmp[[count]] <- list(fixed = fixed,Significant = as.data.frame(sig1),approxposterior = approx_post)
        prescreen2_tmp[[count]]$fixed <- fixed
        p <- 0
      }else{
        prescreen2_tmp[[count]] <- preselection2(SNPs = SNPs,Y = Y,fixed = fixed,kinship = kinship,nullprob = true_nullprop2,alterprob = 1 - true_nullprop2,frequentist = FALSE,threshold = threshold,number_cores = number_cores,P3D = P3D,selfing = selfing)
        p <- length(unique(c(prescreen2_tmp[[count]]$fixed,as.numeric(rownames(prescreen2_tmp[[count]]$Significant == 1))))) - length(prescreen2_tmp[[count]]$fixed)
        significant <- unique(c(prescreen2_tmp[[count]]$fixed,as.numeric(rownames(prescreen2_tmp[[count]]$Significant == 1))))[order(unique(c(prescreen2_tmp[[count]]$fixed,as.numeric(rownames(prescreen2_tmp[[count]]$Significant == 1)))))]
      }
      count <- count + 1
    }
    modelselection <- modelselection_tmp
    prescreen_2nd <- prescreen2_tmp
  }
  return(list(BICOSS_SNPs = fixed))
  
}

resids_diag <- function(Y,SNPs,significant,kinship = FALSE,fixed_effects = FALSE,plot_it = TRUE,selfing){
  requireNamespace("Matrix")
  
  principal_components <- fixed_effects
  
  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(sum(apply(SNPs,2,sum) > nrow(SNPs))){
    stop("Some values of the SNP matrix are not 0 or 1")
  }
  if(!(is.numeric(Y))){
    stop("Y is not a numeric vector")
  }
  if(!(is.numeric(significant))){
    stop("significant is not a numeric vector")
  }
  if(!is.logical(kinship)){
    if(!((nrow(kinship) == nrow(SNPs)) & (ncol(kinship) == nrow(SNPs)))){
      stop("kinship does not have the dimensions nrow(SNPs) x nrow(SNPs)")
    }
  }else{
    if(kinship == TRUE){
      stop("kinship can only take on the values of FALSE or a matrix nrow(SNPs) x nrow(SNPs)")
    }
  }
  if(!is.logical(principal_components)){
    if(!(nrow(principal_components) == nrow(SNPs))){
      stop("principal_components does not have the the same number of rows as SNPs")
    }
    if(sum(apply(principal_components,2,is.numeric)) != ncol(principal_components)){
      stop("Not every column of principal_components is numeric")
    }
  }else{
    if(principal_components == TRUE){
      stop("principal_components can only take on the values of FALSE or a matrix with the same number of rows as SNPs")
    }
  }
  if(!is.logical(plot_it)){
    stop("plot_it is not a logical value")
  }
  
  X <- SNPs
  
  RE_BETA <- function(x,y,d){
    n <- nrow(x)
    p <- ncol(x)
    z <- optimize(log_profile_likelihood_REML,interval = c(0,100),x = x, y = y,d = d,maximum = TRUE)$maximum
    estar <- 1/(1 + z*d)
    xtx <- solve(t(x)%*%(estar*x))
    beta2 <- xtx%*%t(x)%*%(estar*y)
    return(list(beta2,z))
  }
  
  SLR_BETA <- function(x,y){
    beta.mle <- solve(t(x)%*%x)%*%t(x)%*%y
    return(beta.mle)
  }
  
  if(!selfing){
    X_fixed <- matrix(1,nrow = nrow(X),ncol = 1)
    for(i in 1:length(significant)){
      X_fixed <- cbind(X_fixed,SNP_data_function_nonselfing(x = X[,significant[i]]))
    }
    X_fixed <- matrix(X_fixed[,-1],ncol = ncol(X_fixed) - 1)
  }else{
    X_fixed <- matrix(X[,significant],ncol = length(significant))
  }
  
  if(Matrix::rankMatrix(X_fixed)[1] < ncol(X_fixed)){
    dropped_cols <- findLinearCombos(X_fixed)$remove
    X_fixed <- matrix(X_fixed[,-dropped_cols],ncol = Matrix::rankMatrix(X_fixed)[1])
    significant <- significant[-dropped_cols]
  }
  
  Good <- X_fixed
  oG <- Good
  
  nX <- nrow(Good)
  
  oY <- Y
  
  if(is.logical(kinship)){
    if(is.logical(principal_components)){
      #No principal or kinship
      intercept <- matrix(1,nrow = nX,ncol = 1)
      Good <- cbind(intercept,Good)
      
      if(Matrix::rankMatrix(Good)[1] < ncol(Good)){
        dropped_cols <- findLinearCombos(Good)$remove
        Good <- Good[,-dropped_cols]
      }
      beta.good <- SLR_BETA(x = Good, y = Y)
      yhat <- Good %*% beta.good
      resids <- Y - yhat
      
      if(plot_it == TRUE){
        hist(resids,main = "Histogram of Residuals",xlab = "Residuals",breaks = 20)
        return(shapiro.test(resids))
      }else{
        return(shapiro.test(resids))
      }
    }else{
      #Principal and no kinship
      intercept <- matrix(1,nrow = nX,ncol = 1)
      Good <- cbind(intercept,Good,principal_components)
      
      if(Matrix::rankMatrix(Good)[1] < ncol(Good)){
        dropped_cols <- findLinearCombos(Good)$remove
        Good <- Good[,-dropped_cols]
      }
      beta.good <- SLR_BETA(x = Good, y = Y)
      yhat <- Good %*% beta.good
      resids <- Y - yhat
      
      if(plot_it == TRUE){
        hist(resids,main = "Histogram of Residuals",xlab = "Residuals",breaks = 20)
        return(shapiro.test(resids))
      }else{
        return(shapiro.test(resids))
      }
    }
  }else{
    if(is.logical(principal_components)){
      #Kinship no principal
      spec.decomp <- eigen(kinship,symmetric = TRUE)
      spec.decomp$values[length(spec.decomp$values)] <- 0
      Q <- spec.decomp$vectors
      Qt <- t(Q)
      D <- spec.decomp$values
      rm(spec.decomp)
      
      intercept <- matrix(1,nrow = nX,ncol = 1)
      intercept <- Qt%*%intercept
      Y <- Qt%*%Y; Good <- Qt%*%Good
      
      Good <- cbind(intercept,Good)
      oG <- cbind(1,oG)
      beta.good <- RE_BETA(x = Good, y = Y, d = D)
      
      Q <- Q[,-ncol(Q)]
      Qt <- t(Q)
      D <- D[-length(D)]
      uhat <- Q%*%((1/(1 + beta.good[[2]]*D))*Qt)%*%(oY - oG %*% beta.good[[1]])
      resids <- oY - oG %*% beta.good[[1]] - uhat
      
      if(plot_it == TRUE){
        hist(resids,main = "Histogram of Residuals",xlab = "Residuals",breaks = 20)
        return(shapiro.test(resids))
      }else{
        return(shapiro.test(resids))
      }
    }else{
      #Kinship and principal
      ogprinc <- principal_components
      
      spec.decomp <- eigen(kinship,symmetric = TRUE)
      spec.decomp$values[length(spec.decomp$values)] <- 0
      Q <- spec.decomp$vectors
      Qt <- t(Q)
      D <- diag(spec.decomp$values,nrow = nX,ncol = nX)
      rm(spec.decomp)
      
      intercept <- matrix(1,nrow = nX,ncol = 1)
      intercept <- Qt%*%intercept
      Y <- Qt%*%Y; Good <- Qt%*%Good
      principal_components <- Qt%*%principal_components
      
      Good <- cbind(intercept,Good,principal_components)
      oG <- cbind(1,oG,ogprinc)
      
      beta.good <- RE_BETA(x = Good, y = Y, d = D)
      
      Q <- Q[,-ncol(Q)]
      Qt <- t(Q)
      D <- D[-ncol(D),-ncol(D)]
      uhat <- Q%*%((1/(1 + z*d))*Qt)%*%(oY - oG %*% beta.good[[1]])
      resids <- oY - oG %*% beta.good[[1]] - uhat
      
      if(plot_it == TRUE){
        hist(resids,main = "Histogram of Residuals",xlab = "Residuals",breaks = 20)
        return(shapiro.test(resids))
      }else{
        return(shapiro.test(resids))
      }
    }
  }
}
