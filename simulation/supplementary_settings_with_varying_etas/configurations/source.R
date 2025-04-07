
#' * Helpers *

node_preprocess <- function(Y, X, directed, nodes, subtract=NULL)
{
  
  #### Preprep
  directed <- as.logical(directed)
  
  Y <- as.vector(as.numeric(Y))
  if(is.null(dim(X))){  # if X is a vector
    X <- matrix(X, ncol=1)
  } else {
    X <- as.matrix(X)  # matrix
  }
  
  
  d <- length(Y)
  if(nrow(X) != d){stop("X and Y must have the same number of observations")}
  
  remove <- which(is.na(Y))
  if(length(remove) > 0){
    Y <- Y[-remove]
    X <- X[-remove,]
    if(!is.null(nodes)){
      nodes <- nodes[-remove,]
    } else {
      nodes <- node.gen(n, directed)[-remove,]
    }
  }

 
  missing <- FALSE
  cc <- 4 + 4*(directed==F)
  n <- (1+sqrt(1+cc*d))/2
  if(n != round(n)){stop("Size of Y and X must be valid for a network of some size; assuming complete network at this point")}
  
  node_list <- node.set(n, directed) 
  row_list <- lapply(node_list, function(z) cbind( dyad(z[,1], z[,2], n, directed), dyad(z[,3], z[,4], n, directed)) )
  dyads <- 1:d
    
  
  return(list(row_list=row_list))
}

node.set <- function(n.tot, directed=T)
{  

  if(directed){
    nodes.1 <- cbind(rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)],
                     rep(1:n.tot, each=n.tot)[-seq(1,n.tot^2,n.tot+1)], rep(1:n.tot, n.tot)[-seq(1,n.tot^2,n.tot+1)])
    
    nodes.2 <- nodes.3 <- nodes.4 <- nodes.5 <-  c()
    
    for (i in 1:n.tot){ 
      if (i<n.tot){
        c1 <- rep(i,(n.tot-i))
        c2 <- ((i+1):n.tot)
        nodes.2 <- rbind(nodes.2,cbind(c1,c2,c2,c1))
      }
      
      c1 <- rep(i,(n.tot - 1)*(n.tot-2)/2)
      c2 <- rep( (1:(n.tot-1))[-i], times=(n.tot-2):(1 - 1*(i==n.tot) ) )
      c3 <- c1
      c4.mat <- outer(1:n.tot,rep(1,n.tot))[-i,-i]
      c4 <- c4.mat[lower.tri(diag(n.tot-1))]
      
      nodes.3 <- rbind(nodes.3, cbind(c1,c2,c3,c4))
      nodes.4 <- rbind(nodes.4, cbind(c2,c1,c4,c3))
      
      nodes.5 <- rbind(nodes.5, cbind(c1,c2,c4,c3), cbind(c2,c1,c3,c4))   
    }
    return(list(n1=nodes.1,n2=nodes.2,n3=nodes.3,n4=nodes.4,n5=nodes.5))
    
  } else { # undirected
    
    node.list <- node.set(n.tot, directed=T)
    node.list <- lapply(node.list, function(z) z[z[,2] < z[,1] & z[,4] < z[,3], ])   # keep lower tri only
    node.list[[2]] <- unique( rbind(node.list[[2]], node.list[[3]], node.list[[4]], node.list[[5]]) )
    
    node.list <- node.list[1:2]
    return(node.list)
  }
  
}

dyad <- function(i.in, j.in, n.tot, directed=T) 
{
  if (directed==T){ 
    dyad.result <- ((i.in-1) + (j.in-1)*(n.tot-1) + (j.in > i.in)) * (i.in != j.in)
  } else if (directed==F) {
    ij <- cbind(i.in, j.in)
    ij <- t(apply(ij, 1, sort, decreasing=T))
    j.1 <- ij[,2]    
    i.1 <- ij[,1]   
    dyad.result <- ( i.1 - 1 -.5*j.1^2 + (n.tot - 1/2)*j.1 - n.tot + 1  ) * (i.1 != j.1)
  } else {stop('Logical T/F only for directed input')}
  return(dyad.result)
}



#' * Error Term: Use Gamma dist to generate, mean = 1  *

rmvgamma <- function(n,p,shape,rate = 1/scale,rho = 0,scale = 1) { 
  
  ###' * the two gamma component have the same shape and rate*
  ###' * a,b have the same marginal distribution and cov being compound symmetric *
  
  stopifnot(length(rho) == 1,length(shape) == 1,length(rate) == 1, 
            length(scale) == 1,length(n) == 1,length(p) == 1, 
            rate > 0, shape > 0, rho >= 0, rho <= 1, scale > 0, n >= 0, p > 0 
  ) 
  n = round(n); p = ceiling(p) 
  theta = rate * rho / (1-rho) 
  k = rnbinom(n,shape,rate / (rate + theta)) 
  matrix(rgamma(p * n, shape + k, rate + theta), n) 
}



#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Exchangeable_Error_Gamma_mixing
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
error_term_gene_gamma <- function(n,a1,b1,a2,b2,a0,b0,pho_ab){
  
  ###'@param n: number of cities.
  ###' Gamma: @param a1: shape, @param b1: scale, 
  ###'        @mean: a1b1, @variance: a1 * b1^2,
  ###' Epsilon: @param a2: shape, @param b2: scale,
  ###'           @mean: a2b2, @variance: a2 * b2^2,
  ###' ab: @param a0: marginal shape, @param b0: marginal rate = 1/scale,
  ###'     @param pho_ab: correlation coefficient of ab,
  ###'     * a,b have the same marginal distribution and cov being compound symmetric *
  
  ###
  ### Variables and Priors
  ###
  e_ij.save <- matrix(0,n,n)
  
  ###
  ###'@Generate_Epsilon_ij__Paras:a2,b2
  ###
  epsilon_ij.save <- matrix(rgamma(n^2, shape = a2, scale = b2), n, n)  
  diag(epsilon_ij.save) <- 0
  # epsilon.mean <- a2 * b2
  
  ###
  ###'@Generate_Gamma_ij__Paras:a1,b1
  ###
  gamma_ij.save <- matrix(rgamma(n^2, shape = a1, scale = b1), n, n)
  gamma_ij.save[upper.tri(gamma_ij.save)] <- t(gamma_ij.save)[upper.tri(gamma_ij.save)]
  diag(gamma_ij.save) <- 0
  # gamma.mean <- a1 * b1
  
  ###
  ###'@Generate_ab_ij__Paras:a0,b0,pho_ab
  ###
  ab.temp <- rmvgamma(n,2,a0,b0,pho_ab)
  a_i.save <- ab.temp[,1]
  b_j.save <- ab.temp[,2]
  
  ###
  ### Normalization Constant
  ###
  C <- 1 / (2 * a0 / b0 + a1 * b1 + a2 * b2)
  
  ###
  ### Generate Error Terms
  ###
  e_ij.save <- (matrix(rep(a_i.save, n), n, n) + matrix(rep(b_j.save, n), n, n, byrow = TRUE) + 
                  gamma_ij.save + epsilon_ij.save) * C
  diag(e_ij.save) <- 0
  
  ###
  ### Write Output 
  ###  
  list(e_ij.save = e_ij.save, gamma_ij.save = gamma_ij.save, epsilon_ij.save = epsilon_ij.save,
        a_i.save = a_i.save, b_j.save = b_j.save, C = C)
  
}



#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Error_Term_Genrating_Function_IID_gamma
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
error_term_gene_iidgamma <- function(n){
  
  ###
  ### Variables 
  ###
  e_ij.save <- matrix(rgamma(n^2, shape = 1/2, scale = 2),n,n)  # mean=1, var=2
  diag(e_ij.save) <- 0
  
  ###
  ### Write Output 
  ###  
  list(e_ij.save = e_ij.save)
  
}



#' * Error Term: Use TruncNormal dist to generate, mean = 1 *

#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Exchangeable_Error_TrunckNormal_mixing
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
error_term_gene_truNormal <- function(n,a1,b1,mu1,sigma2_1,a2,b2,mu2,sigma2_2,a0_vec,b0_vec,
                                      mu0_vec,pho_ab,sigma_a,sigma_b){ # 注：需取 b1,b2=Inf 便于下面算6para的真值
  
  ###'@param n: number of cities.
  ###' Gamma: @param a1: lower bound of truncNormal, @param b1=Inf: upper bound of truncNormal, 
  ###'        @param mu1: original mean, @param sigma2_1: original variance
  ###' Epsilon: @param a2: lower bound of truncNormal, @param b2=Inf: upper bound of truncNormal, 
  ###'          @param mu2: original mean, @param sigma2_2: original variance
  ###' ab: @param a0_vec: lower bound of a, @param b0_vec: lower bound of b, 
  ###'     @param mu0_vec: original mean vector of ab, @param pho_ab: correlation coefficient of ab,
  ###'     @param sigma_a: original standard deviation of a, @param sigma_b: original standard deviation of b
  
  ###
  ### Subroutines
  ###
  # library(tmvtnorm)
  library(truncnorm)
  library(TruncatedNormal)
  
  
  ###
  ### Variables and Priors
  ###
  e_ij.save <- matrix(0,n,n) 

  
  ###
  ###'@Generate_Epsilon_ij__Paras:a2,b2,mu2,sigma2_2
  ###
  epsilon_ij.save <- matrix(rtruncnorm(n^2, a = a2, b = b2, mean = mu2, sd = sqrt(sigma2_2)), n, n) 
  diag(epsilon_ij.save) <- 0
  
  epsilon.mean <- mu2 - sqrt(sigma2_2) * (dnorm((b2-mu2)/sqrt(sigma2_2)) - dnorm((a2-mu2)/sqrt(sigma2_2))) /
    (pnorm((b2-mu2)/sqrt(sigma2_2)) - pnorm((a2-mu2)/sqrt(sigma2_2)))

  
  ###
  ###'@Generate_Gamma_ij__Paras:a1,b1,mu1,sigma2_1
  ###
  gamma_ij.save <- matrix(rtruncnorm(n^2, a = a1, b = b1, mean = mu1, sd = sqrt(sigma2_1)), n, n)  
  gamma_ij.save[upper.tri(gamma_ij.save)] <- t(gamma_ij.save)[upper.tri(gamma_ij.save)]
  diag(gamma_ij.save) <- 0
  
  gamma.mean <- mu1 - sqrt(sigma2_1) * (dnorm((b1-mu1)/sqrt(sigma2_1)) - dnorm((a1-mu1)/sqrt(sigma2_1))) /
    (pnorm((b1-mu1)/sqrt(sigma2_1)) - pnorm((a1-mu1)/sqrt(sigma2_1)))

  ###
  ###'@Generate_ab_ij__Paras:a0_vec,b0_vec,mu0_vec,pho_ab,sigma_a,sigma_b
  ###
  mu <- mu0_vec
  sigma <- matrix(c((sigma_a)^2, pho_ab*sigma_a*sigma_b, pho_ab*sigma_a*sigma_b, (sigma_b)^2), nrow = 2, ncol = 2, byrow = TRUE)
  lb <- c(a0_vec,b0_vec)
  ab.temp <- TruncatedNormal::rtmvnorm(n = n, mu = mu, sigma = sigma, lb = lb)
  a_i.save <- ab.temp[,1]
  b_j.save <- ab.temp[,2]
  
  # ab.moments <- tmvtnorm::mtmvnorm(mean = mu, sigma = sigma, lower = lb)
  
  ###
  ### Normalization Constant
  ###
  # C <- 1 / (ab.moments$tmean[1] + ab.moments$tmean[2] + epsilon.mean + gamma.mean)
  C <- 1 / (0.5298593 + 1.8137138 + epsilon.mean + gamma.mean)
  
  ###
  ### Generate Error Terms
  ###
  e_ij.save <- (matrix(rep(a_i.save, n), n, n) + matrix(rep(b_j.save, n), n, n, byrow = TRUE) + 
                  gamma_ij.save + epsilon_ij.save) * C
  diag(e_ij.save) <- 0
  
  ###
  ### Write Output 
  ###  
  list(e_ij.save = e_ij.save, gamma_ij.save = gamma_ij.save, epsilon_ij.save = epsilon_ij.save,
        a_i.save = a_i.save, b_j.save = b_j.save, C = C)

}



#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Error_Term_Genrating_Function_IID_TruncNormal
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
error_term_gene_iidtruNormal <- function(n){
  
  library(truncnorm)
  
  ###
  ### Variables 
  ###
  e_ij.save <- matrix(rtruncnorm(n^2, a = 0, b = Inf, mean = 1, sd = 1),n,n) / 1.2876 # mean=1, var=0.3798064
  diag(e_ij.save) <- 0
  
  ###
  ### Write Output 
  ###  
  list(e_ij.save = e_ij.save)
  
}



#' * lambda_ij *

#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Lambda_ij_link_function_exp(.)
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
lamb_gene_exp <- function(n, beta, X_full, Error_term){
  
  Xbeta <- X_full$X1 * beta[1] + X_full$X2 * beta[2] + X_full$X3 * beta[3] + X_full$X4 * beta[4]
  lamb.matr <- exp(Xbeta) * Error_term$e_ij.save  #'@e_ij.save_in_matrix_form
  
  diag(lamb.matr) <- NA
  elements <- as.vector(t(lamb.matr))
  lamb.vec <- na.omit(elements)
  diag(lamb.matr) <- 0
  
  list(Xbeta = Xbeta, lamb.matr = lamb.matr, lamb.vec = lamb.vec)
  
}



#' * Y *

#### >>>>>>>>>>>>>>>
####'@Y_ij_Genrating
#### >>>>>>>>>>>>>>>
Y_gene <- function(n, Lambda_matr_and_vec){

  Y.matr <- matrix(rpois(n * n, Lambda_matr_and_vec$lamb.matr), ncol = n, nrow = n)
  
  diag(Y.matr) <- NA
  elements <- as.vector(t(Y.matr))
  Y.vec <- na.omit(elements)
  
  diag(Y.matr) <- 0

  list(Y.matr = Y.matr, Y.vec = Y.vec)
  
}



#' * \wh \beta from PMLE *

#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Pseudo-likelihood_for_Pois_Family:g=exp(.)
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pseudo_estimators_exp <- function(beta_true, X_full, Y_matr_and_vec){
  
  ###
  ### Variables to save
  ###
  Beta_pmle <- rep(0,4)
  
  ###
  ### Data Manipulation
  ###
  X1 <- X_full$x1
  X2 <- X_full$x2
  X3 <- X_full$x3
  X4 <- X_full$x4
  Y <- Y_matr_and_vec$Y.vec
  
  ####
  #### Inverse Log-likelihood Function
  ####
  inv_lglk_degen_exp <- function(p) 
    sum(exp(X1 * p[1] + X2 * p[2] + X3 * p[3] + X4 * p[4])) - 
    sum(Y *( X1 * p[1] + X2 * p[2] + X3 * p[3] + X4 * p[4]))
  
  ####
  #### PMLE of beta
  ####
  optim.temp <- optim(par = c(0,0,0,0), fn = inv_lglk_degen_exp, method = "Nelder-Mead") 
  for(i in c(1:4)){
    Beta_pmle[i] <- optim.temp$`par`[i]
  }

  ###
  ### Write Output 
  ###  
  list(Beta_pmle = Beta_pmle)  
}



#' * moment estimation of eta2, eta3, eta4, eta5 *

#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Estimate_6_para_of_Covariance_eij_for_Pois_Family:g=exp(.)
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
para_6_estimate_exp <- function(n, beta_est, X_full, Y_matr_and_vec){
  
  xi.vec = Y_matr_and_vec$Y.vec / exp(X_full$x1 * beta_est[1] + X_full$x2 * beta_est[2] + 
                                        X_full$x3 * beta_est[3] + X_full$x4 * beta_est[4])
  
  Matrix.temp <- matrix(xi.vec, (n-1), n, byrow = TRUE)
  elementsWithNA <- cbind(rep(NA, (n-1)), Matrix.temp)
  xi.matr <- matrix(c(as.vector(t(elementsWithNA)), NA), n, n, byrow = TRUE)
  
  lam_est.save <- rep(0,6)
  
  ###
  ###'@Var(e_ij),keep_it_0,no_need_to_get_mmt_estimation
  ###
  Xbeta.ha <- X_full$x1 * beta_est[1] + X_full$x2 * beta_est[2] + X_full$x3 * beta_est[3] + X_full$x4 * beta_est[4]
  c.vec <- mean(xi.vec)^2 + 1 / exp(Xbeta.ha)
  temp.vec <- xi.vec^2 - c.vec #' * <- empirical estimation, an alternative way shown below *
                               #' * used for later shorth estimation *
  c.vec_thoretical <- 1 + 1 / exp(Xbeta.ha)
  temp.vec_thoretical <- xi.vec^2 - c.vec_thoretical

  ###
  ###'@Cov(e_ij,e_ji)=eta_2_with_count=(n^2-n)
  ###
  a_2 <- xi.vec
  b_2 <- na.omit(as.vector(xi.matr))
  lam_est.save[2] <- cov(a_2, b_2)  #' * <- empirical estimation, an alternative way -> mean(a_2b_2) - 1 *
                                    #' * same approach for eta_3,4,5 *
  lam_est.save[8] <- mean(a_2 * b_2) - 1  #' * <- an alternative way -> mean(a_2b_2) - 1 *
  
  ###
  ###'@Cov(e_ij,e_il)=eta_3_with_count=(n^2-n)(n-2)
  ###
  a_3 <- rep(xi.vec, each = n-2)
  B.temp <- matrix(xi.vec, (n-1), n)
  combine.temp_3 <- B.temp[-1,]
  for(i in 2:(n-1)){
    combine.temp_3 <- rbind(combine.temp_3, B.temp[-i,])
  }
  b_3 <- as.vector(combine.temp_3)
  lam_est.save[3] <- cov(a_3, b_3)
  lam_est.save[9] <- mean(a_3 * b_3) - 1
  
  ###
  ###'@Cov(e_ij,e_kj)=eta_4_with_count=(n^2-n)(n-2)
  ###
  xi.vec.transp <- na.omit(as.vector(xi.matr))
  a_4 <- rep(xi.vec.transp, each = n-2)
  D.temp <- matrix(xi.vec.transp, (n-1), n)
  combine.temp_4 <- D.temp[-1,]
  for(i in 2:(n-1)){
    combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
  }
  b_4 <- as.vector(combine.temp_4)
  lam_est.save[4] <- cov(a_4, b_4)
  lam_est.save[10] <- mean(a_4 * b_4) - 1
  
  ###
  ###'@Cov(e_ij,e_ki)((e_ij,e_jk))
  ###
  lam_est.save[5] <- cov(a_3, b_4)
  lam_est.save[11] <- mean(a_3 * b_4) - 1
  
  ###
  ### Write Output 
  ###  
  list(lam_6_est = lam_est.save, eta1mmt = temp.vec, eta1mmt_thoretical = temp.vec_thoretical) 
  
}



#' * hybrid Shorth estimation of eta1 (strictly > 0) *

#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Set_group_for_cross_validation
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CVgroup <- function(k, datasize){  
  
  cvlist <- list()
  n <- rep(1:k, ceiling(datasize / k))[1: datasize]    
  temp <- sample(n, datasize)   
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp == x])  
  return(cvlist)
  
}


#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Positive_shorth_estimation_based_on_mmt_estimation_of_eta_1
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
shorthPosi <- function(l, k){
  
  ###'@param l: list of numbers
  ###'@param k: number of points contained in target intervals
  
  l <- sort(l)
  N <- length(l)
  needle <- 1
  
  if(needle + k > N){
    stop("Try a smaller value of tunning parameter", call. = FALSE)
  }
  
  while(-l[needle] > l[needle + k -1]){
    needle <- 1 + needle
    
    if(needle + k > N){
      stop("Try a smaller value of tunning parameter", call. = FALSE)
    }
  }
  
  l.omit.nega <- needle
  
  len <- l[l.omit.nega + k] - l[l.omit.nega + 1] 
  position <- l.omit.nega + 1
  
  for(i in (l.omit.nega + 2):(N - k)){
    if(l[i + k - 1] - l[i] <= len){
      len <- l[i + k - 1] - l[i]
      position <- i
    }
  }
  
  return(c((l[position] / 2 + l[position + k - 1] / 2), len))
  
}



#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Adding_cross_validation_Positive_shorth_estimation_based_on_mmt_estimation_of_eta_1
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cvShorth <- function(v, Ctune = c(0)){
  
  Ctune <- sort(Ctune)
  
  fold <- 10
  datasize <- length(v)
  cvlist <- CVgroup(k = fold, datasize = datasize)  
  
  data <- v
  Cmin <- 2 / log(datasize)
  Cmax <- sum(v > - max(v)) / log(datasize)
  
  if(Ctune[1] == 0){
    C_list <-  seq(Cmin, Cmax / 2 , length.out = 50) 
  } else{
    C_list <- Ctune
  }
  
  cv.res <- data.frame()   
  
  for(C_tune in C_list){  
    
    fold.mse.list <- rep(0, fold)
    k.init <- 0
    
    for (i in 1:fold){
      train <- data[-cvlist[[i]]] 
      test <- data[cvlist[[i]]]
      
      k.train <- ceiling(C_tune * log(length(train)))
      k.init <- k.train
      
      k.shorth <- shorthPosi(train, k = k.train)
      fold.mse.list[i] <- mean((test - k.shorth[1])^2)
      
    }
    
    temp <- data.frame(cbind("C_tune" = C_tune, "k" = k.init, k.mse = mean(fold.mse.list)))
    cv.res <- rbind(cv.res, temp)   
    
  }
  
  k_final <- ceiling(cv.res$C_tune[which.min(cv.res$k.mse)] * log(length(v)))
  
  fit <- shorthPosi(v, k = k_final)
  
  return(fit[1]) 
  
}



#### >>>>>>>>>>>>>>>>>>>>>>>>>
####'@Hybrid_shorth_estimation
#### >>>>>>>>>>>>>>>>>>>>>>>>>
posiEstimator_hyCV <- function(v, Ctune = c(0)){
  
  v_mmt <- mean(v)
  
  if (v_mmt <= 0){
    v_hyCvShorth <- cvShorth(v)
  } else {
    v_hyCvShorth <- v_mmt
  }
  
  list(hyCvShorth = v_hyCvShorth)
  
}



#' * Estimation of covariance matrix of \wh \beta *

#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@One_step_estimation_of_JIJ_for_single_simulation-exp()_link
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

AinvBAinv_exp_one_step <- function(n, beta_est, lambda, X_full){
  
  library(netregR)
  
  X.inner = matrix(c(X_full$x1, X_full$x2, X_full$x3, X_full$x4), 4, (n^2 - n), byrow = TRUE)
  
  A.temp = X.inner %*% diag(exp(X_full$x1 * beta_est[1] + X_full$x2 * beta_est[2] + 
                                  X_full$x3 * beta_est[3] + X_full$x4 * beta_est[4])) %*% t(X.inner)
  
  X.Sigma = t(X.inner %*% diag(exp(X_full$x1 * beta_est[1] + X_full$x2 * beta_est[2] +
                                  X_full$x3 * beta_est[3] + X_full$x4 * beta_est[4])))
  
  temp <- node_preprocess(rep(1, n^2 - n), t(X.inner), TRUE, NULL)
  row.list <- temp$row_list

  meat.1 <- lambda[1] * crossprod(X.Sigma)

  lambda_new <- lambda
  lambda_new[3] <- lambda[4]
  lambda_new[4] <- lambda[3]
  meat.2 <- matrix(0,4,4)
  for (i in 2:5){
    r <- matrix(row.list[[i]], ncol = 2)
    if(nrow(r) > 0){
      meat.2 <- meat.2 + crossprod(X.Sigma[r[,1], ,drop=F], X.Sigma[r[,2], ,drop=F]) * lambda_new[i]  +
                         crossprod(X.Sigma[r[,2], ,drop=F], X.Sigma[r[,1], ,drop=F]) * lambda_new[i]
    }
  }
  meat.out <- meat.2 + meat.1
  
  Cov.save_simple  = chol2inv(chol(A.temp)) %*% (meat.out + A.temp) %*% chol2inv(chol(A.temp))
  
  return(Cov.save_simple)
}












