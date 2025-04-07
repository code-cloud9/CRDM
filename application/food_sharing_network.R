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
  # Generate node sets of various overlapping dyad pairs
  # Return list of each set of nodes, with null for first set (diagonal)
  
  # if(!directed){stop("node.set() not yet coded for undirected")}
  
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


#' * \wh \beta from PMLE *

#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@Pseudo-likelihood_for_Pois_Family:g=exp(.)
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pseudo_estimators_exp_with_off_set <- function(X, y, offset){
  
  ###
  ### Variables to save
  ###
  if (is.null(ncol(X))){
    p <- 1
  } else {
    p <- ncol(X)
  }
  Beta_pmle <- rep(0,p)
  
  ###
  ### Data Manipulation
  ###
  X <- data.matrix(X)
  
  
  ####
  #### Inverse Log-likelihood Function
  ####
  if (p == 1){
    inv_lglk_degen_exp <- function(para,X,y){
      return(sum(exp(X * para + offset)) - sum(y * (X * para + offset)))
    }
  } else {
    inv_lglk_degen_exp <- function(para,X,y){
      return(sum(exp(X %*% para + offset)) - sum(y * (X %*% para + offset)))
    }
  }
  
  
  ####
  #### PMLE of beta
  ####
  if (p == 1){
    optim.temp <- optimize(f = inv_lglk_degen_exp, interval = c(0, 10), y = y, X = X) 
  } else {
    optim.temp <- optim(par = rep(0,p), fn = inv_lglk_degen_exp, method = "Nelder-Mead", y = y, X = X) 
  }
  
  
  if (p == 1){
    Beta_pmle[1] <- optim.temp$minimum
  } else {
    for(i in c(1:p)){
      Beta_pmle[i] <- optim.temp$`par`[i]
    }
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
para_6_estimate_exp_with_off_set <- function(n, beta_est, X, y, offset){
  
  ###
  ### Data Manipulation
  ###
  if (is.null(ncol(X))){
    p <- 1
  } else {
    p <- ncol(X)
  }
  X <- data.matrix(X)
  
  xi.vec <- y / exp(X %*% beta_est + offset)
  
  Matrix.temp <- matrix(xi.vec, (n-1), n, byrow = TRUE)
  elementsWithNA <- cbind(rep(NA, (n-1)), Matrix.temp)
  xi.matr <- matrix(c(as.vector(t(elementsWithNA)), NA), n, n, byrow = TRUE)
  
  lam_est.save <- rep(0,12)
  
  ###
  ###'@Var(e_ij),keep_it_0,no_need_to_get_mmt_estimation
  ###
  Xbeta.ha <- X %*% beta_est + offset
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
AinvBAinv_exp_one_step_with_off_set <- function(n, beta_est, lambda, X, offset){
  
  ###
  ### Data Manipulation
  ###
  if (is.null(ncol(X))){
    p <- 1
  } else {
    p <- ncol(X)
  }
  X <- data.matrix(X)
  
  library(netregR)
  
  X.inner = t(X) 
  
  A.temp = X.inner %*% diag(exp(as.vector(unlist(X %*% beta_est)) + offset)) %*% t(X.inner)
  
  X.Sigma = t(X.inner %*% diag(exp(as.vector(unlist(X %*% beta_est)) + offset)))
  
  
  temp <- node_preprocess(rep(1, n^2 - n), t(X.inner), TRUE, NULL)
  row.list <- temp$row_list
  
  meat.1 <- lambda[1] * crossprod(X.Sigma)
  
  lambda_new <- lambda
  lambda_new[3] <- lambda[4]
  lambda_new[4] <- lambda[3]
  meat.2 <- matrix(0,p,p)
  for (i in 2:5){
    r <- matrix(row.list[[i]], ncol = 2)
    if(nrow(r) > 0){
      meat.2 <- meat.2 + crossprod(X.Sigma[r[,1], ,drop=F], X.Sigma[r[,2], ,drop=F]) * lambda_new[i]  +
        crossprod(X.Sigma[r[,2], ,drop=F], X.Sigma[r[,1], ,drop=F]) * lambda_new[i]
    }
  }
  meat.out <- meat.2 + meat.1
  
  Cov.save_simple <- chol2inv(chol(A.temp)) %*% (meat.out + A.temp) %*% chol2inv(chol(A.temp))  
  
  return(Cov.save_simple)
}


#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
####'@CI_based_on_Our_method_empirical
#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CI_with_offset <- function(df){
  
  ###'@param data: column 1: y
  ###'@param data: column 2~p+1: X
  # second column in df is offset
  
  dimen <- ncol(df)
  
  PMLE_summar.example <- pseudo_estimators_exp_with_off_set(X = df[,3:dimen], y = as.vector(unlist(df[,1])), offset = as.vector(unlist(df[,2])))

  eta_6.moment <- para_6_estimate_exp_with_off_set(n = city_size, beta_est = PMLE_summar.example$Beta_pmle, X = df[,3:dimen], y = as.vector(unlist(df[,1])), offset = as.vector(unlist(df[,2])))
  
  #' * empiri_CVshorth_eigenCorrection *
  
  #'@_eta1_empiri_CVshorth
  empi_beta_hat <- eta_6.moment$eta1mmt
  
  #'@_To_avoid_posiEstimator_hyCV_having_no_solution_caused_by_10fold_CV
  sss <- 0
  while(sss == 0){
    try(fit.empi_beta_hat <- posiEstimator_hyCV(v = empi_beta_hat, Ctune = c(0)), silent = TRUE)
    try(sss <- fit.empi_beta_hat$hyCvShorth, silent = TRUE)
  }
  
  #'@_calculate_smallest_eigen-value
  a <- eta_6.moment$lam_6_est[2] / fit.empi_beta_hat$hyCvShorth
  b <- eta_6.moment$lam_6_est[3] / fit.empi_beta_hat$hyCvShorth
  c <- eta_6.moment$lam_6_est[4] / fit.empi_beta_hat$hyCvShorth
  d <- eta_6.moment$lam_6_est[5] / fit.empi_beta_hat$hyCvShorth
  
  al <- (city_size - 1)^2 * (b^2 + c^2) + 4 * d^2 * (city_size - 3)^2 + 2 * b * c * (1 - city_size^2 + 2 * city_size)
  bta <- a * d * (8 * city_size - 24) + (b + c) * d * (12 - 4 * city_size) + 4 * a * (a - b - c)
  
  #'@_hybrid_Shorth_estimation_result
  CV_s_fix <- fit.empi_beta_hat$hyCvShorth +
    max(fit.empi_beta_hat$hyCvShorth * c(0, -(1 + a + (city_size - 2) * (b + c) + 2 * (city_size - 2) * d),
                                         -(1 + a - b - c - 2 * d), -(1 - a - b - c + 2 * d),
                                         -(((city_size - 3) * (b + c) - 2 * d + 2) / 2 - sqrt(al + bta) / 2)))
  
  #'@_Var(beta)_mine_empiri
  eta_list_1 <- c(CV_s_fix, eta_6.moment$lam_6_est[2:6])
  AinvBAinv_1 <- AinvBAinv_exp_one_step_with_off_set(n = city_size, beta_est = PMLE_summar.example$Beta_pmle, lambda = eta_list_1, X = df[,3:dimen], offset = as.vector(unlist(df[,2])))
  diag(AinvBAinv_1)
  
  #'@_Construct_CIs
  return(list(vct = c( PMLE_summar.example$Beta_pmle - qnorm(0.975) * sqrt(diag(AinvBAinv_1)),
            PMLE_summar.example$Beta_pmle,
            PMLE_summar.example$Beta_pmle + qnorm(0.975) * sqrt(diag(AinvBAinv_1)),
            eta_list_1), Covbeta = AinvBAinv_1))
  
}




#'@<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#'@_X_elementwise

library(readr)
households <- read_csv("households.csv")

# X1: Giver – Game
x1_vec <- rep(households$hgame, each = 24)

# X2: Giver – Fish
x2_vec <- rep(households$hfish, each = 24)

# X3: Giver – Pigs
x3_vec <- rep(households$hpigs, each = 24)

# X4: Giver – Wealth
x4_vec <- rep(households$hwealth, each = 24)

# X5: Receiver – Game
A5 <- matrix(rep(households$hgame,25), 25, 25, byrow = TRUE)  
diag(A5) <- NA
elements5 <- as.vector(t(A5))
x5_vec <- na.omit(elements5)

# X6: Receiver – Fish
A6 <- matrix(rep(households$hfish,25), 25, 25, byrow = TRUE)  
diag(A6) <- NA
elements6 <- as.vector(t(A6))
x6_vec <- na.omit(elements6)

# X7: Receiver – Pigs
A7 <- matrix(rep(households$hpigs,25), 25, 25, byrow = TRUE)  
diag(A7) <- NA
elements7 <- as.vector(t(A7))
x7_vec <- na.omit(elements7)

# X8: Receiver – Wealth
A8 <- matrix(rep(households$hwealth,25), 25, 25, byrow = TRUE)  
diag(A8) <- NA
elements8 <- as.vector(t(A8))
x8_vec <- na.omit(elements8)

# X9: Receiver – Pastors
A9 <- matrix(rep(households$hpastor,25), 25, 25, byrow = TRUE)  
diag(A9) <- NA
elements9 <- as.vector(t(A9))
x9_vec <- na.omit(elements9)


#'@<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#'@_Y_plots

# **********
# Figure S.4
# **********
dyads <- read_csv("dyads.csv")
Y_matr <- matrix(NA,25,25) 
for (i in 1:300) {
    Y_matr[dyads$hidA[i],dyads$hidB[i]] <- dyads$giftsAB[i]
    Y_matr[dyads$hidB[i],dyads$hidA[i]] <- dyads$giftsBA[i]
}
elementsy <- as.vector(t(Y_matr))
y_vec <- na.omit(elementsy)
hist(y_vec, breaks = 30, main = "histgram of network")

Y_matr_plot <- Y_matr
diag(Y_matr_plot) <- 0
library('plot.matrix')
plot(Y_matr_plot, breaks = c(0,5,20,50,200), 
     col = c(rgb(210/252,237/252,252/252), 
             rgb(135/252,206/252,250/252), 
             rgb(22/252,160/252,245/252),
             rgb(6/252,91/252,144/252)),
     key = NULL, main = "heat map of network", xlab = '', ylab = '', cex.main = 1,
     axis.col = list(cex.axis = 0.5),
     axis.row = list(cex.axis = 0.5))

dyads$drel5 <- 1- (dyads$drel1 + dyads$drel2 + dyads$drel3 + dyads$drel4)


#'@<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#'@_X_pairwise

# X10: Relationship - Relatedness 1
X10_matr <- matrix(NA,25,25) 
for (i in 1:300) {
  X10_matr[dyads$hidA[i],dyads$hidB[i]] <- dyads$drel1[i]
  X10_matr[dyads$hidB[i],dyads$hidA[i]] <- dyads$drel1[i]
}

elements10 <- as.vector(t(X10_matr))
x10_vec <- na.omit(elements10)

# X11: Relationship - Relatedness 2
X11_matr <- matrix(NA,25,25) 
for (i in 1:300) {
  X11_matr[dyads$hidA[i],dyads$hidB[i]] <- dyads$drel2[i]
  X11_matr[dyads$hidB[i],dyads$hidA[i]] <- dyads$drel2[i]
}

elements11 <- as.vector(t(X11_matr))
x11_vec <- na.omit(elements11)

# X12: Relationship - Relatedness 3
X12_matr <- matrix(NA,25,25) 
for (i in 1:300) {
  X12_matr[dyads$hidA[i],dyads$hidB[i]] <- dyads$drel3[i]
  X12_matr[dyads$hidB[i],dyads$hidA[i]] <- dyads$drel3[i]
}

elements12 <- as.vector(t(X12_matr))
x12_vec <- na.omit(elements12)

# X13: Relationship - Relatedness 4
X13_matr <- matrix(NA,25,25) 
for (i in 1:300) {
  X13_matr[dyads$hidA[i],dyads$hidB[i]] <- dyads$drel4[i]
  X13_matr[dyads$hidB[i],dyads$hidA[i]] <- dyads$drel4[i]
}

elements13 <- as.vector(t(X13_matr))
x13_vec <- na.omit(elements13)


# X1333: Relationship - Relatedness 5
X1333_matr <- matrix(NA,25,25) 
for (i in 1:300) {
  X1333_matr[dyads$hidA[i],dyads$hidB[i]] <- dyads$drel5[i]
  X1333_matr[dyads$hidB[i],dyads$hidA[i]] <- dyads$drel5[i]
}

elements1333 <- as.vector(t(X1333_matr))
x1333_vec <- na.omit(elements1333)


# X14: Relationship - Distance (log transformed)
X14_matr <- matrix(NA,25,25) 
for (i in 1:300) {
  X14_matr[dyads$hidA[i],dyads$hidB[i]] <- dyads$dlndist[i]
  X14_matr[dyads$hidB[i],dyads$hidA[i]] <- dyads$dlndist[i]
}

elements14 <- as.vector(t(X14_matr))
x14_vec <- na.omit(elements14)

# X15: Relationship - Association index
X15_matr <- matrix(NA,25,25) 
for (i in 1:300) {
  X15_matr[dyads$hidA[i],dyads$hidB[i]] <- dyads$dass[i]
  X15_matr[dyads$hidB[i],dyads$hidA[i]] <- dyads$dass[i]
}

elements15 <- as.vector(t(X15_matr))
x15_vec <- na.omit(elements15)

# X16: Relationship - Giver 1 & Receiver 25
X16_matr <- matrix(NA,25,25) 
for (i in 1:300) {
  X16_matr[dyads$hidA[i],dyads$hidB[i]] <- dyads$d0125[i]
  X16_matr[dyads$hidB[i],dyads$hidA[i]] <- dyads$d0125[i]
}

elements16 <- as.vector(t(X16_matr))
x16_vec <- na.omit(elements16)


# X_offset: offset
Xoffset_matr <- matrix(NA,25,25) 
for (i in 1:300) {
  Xoffset_matr[dyads$hidA[i],dyads$hidB[i]] <- dyads$offset[i]
  Xoffset_matr[dyads$hidB[i],dyads$hidA[i]] <- dyads$offset[i]
}

elementsoffset <- as.vector(t(Xoffset_matr))
xoffset_vec <- na.omit(elementsoffset)


# check co-linearity
library(corrgram)

data_check <- data.frame("x1" = x1_vec, "x2" = x2_vec,
                         "x3" = x3_vec, "x4" = x4_vec, "x5" = x5_vec, "x6" = x6_vec,
                         "x7" = x7_vec, "x8" = x8_vec, "x9" = x9_vec, "x10" = x10_vec,
                         "x11" = x11_vec, "x12" = x12_vec, "x13" = x13_vec,"x1333" = x1333_vec, "x14" = x14_vec, 
                         "x15" = x15_vec) 

corrgram(data_check, order = FALSE, upper.panel = panel.cor, lower.panel = panel.pie, 
         main = "Correlogram intercorrelations")


data_total_with_offset <- data.frame("y" = y_vec, "offset" = xoffset_vec, "x0" = rep(1,600), "x1" = x1_vec, "x2" = x2_vec,
                       "x3" = x3_vec, "x4" = x4_vec, "x5" = x5_vec, "x6" = x6_vec,
                       "x7" = x7_vec, "x8" = x8_vec, "x9" = x9_vec, "x10" = x10_vec,
                       "x11" = x11_vec, "x12" = x12_vec, "x13" = x13_vec, "x14" = x14_vec, 
                       "x15" = x15_vec)


city_size <- 25
result <- CI_with_offset(data_total_with_offset)
vct <- result$vct
Covbeta <- result$Covbeta
Corrbeta <- cov2cor(Covbeta)

library(pheatmap)
pheatmap(Corrbeta, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)

# **********
# Figure S.4
# **********
library(corrgram)
corrgram(abs(Corrbeta), order = FALSE, type = "cor", labels = c(expression(hat(beta[0])), expression(hat(beta[1])),
                                                           expression(hat(beta[2])), expression(hat(beta[3]))
                                                           , expression(hat(beta[4])), expression(hat(beta[5]))
                                                           , expression(hat(beta[6])), expression(hat(beta[7]))
                                                           , expression(hat(beta[8])), expression(hat(beta[9]))
                                                           , expression(hat(beta[10])), expression(hat(beta[11]))
                                                           , expression(hat(beta[12])), expression(hat(beta[13]))
                                                           , expression(hat(beta[14])), expression(hat(beta[15]))), 
         upper.panel = panel.cor, lower.panel = panel.pie, 
         main = NULL, gap = 0.1)





#'@_Data_analysis_on_selected_attributes<><><><><><><><><><><><><><><><><><>

data_reduced_with_offset <- data.frame("y" = y_vec, "offset" = xoffset_vec, "x0" = rep(1,600), 
                                       "x1" = x1_vec, "x2" = x2_vec,
                                       "x3" = x3_vec, 
                                       "x5" = x5_vec, "x6" = x6_vec,
                                       "x7" = x7_vec, 
                                       "x1333" = x1333_vec+x13_vec,
                                       "x14" = x14_vec, "x15" = x15_vec) 

city_size <- 25
result <- CI_with_offset(data_reduced_with_offset)
vct <- result$vct

CI_idx <- matrix(vct[c(1:30)], 10, 3)
idx_eta_linePlot <- vct[31:35]
CI_idx <- as.data.frame(CI_idx)
CI_idx$idx <- c(0:9)  
CI_idx



#'@_beta_Plots 
# ********
# Figure 2
# ********

library(ggplot2)
# Intercept
scaleFUN <- function(x) sprintf("%.1f", x) 

df <- data.frame(
  year = c("Intercept"),
  x = factor(c(1)),
  est = c(CI_idx$V2[1]),
  upper = c(CI_idx$V3[1]),
  lower = c(CI_idx$V1[1])
)

p0=ggplot(df, aes(y = est, x = x, color = "red")) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  scale_y_continuous(labels = scaleFUN) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  facet_wrap(~ year) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.text=element_text(size=9),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# Relatedness
df <- data.frame(
  year = c("Relatedness"),
  x = factor(c(1)),
  est = c(CI_idx$V2[8]),
  upper = c(CI_idx$V3[8]),
  lower = c(CI_idx$V1[8])
)

p7=ggplot(df, aes(y = est, x = x, color = "red")) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  scale_y_continuous(labels = scaleFUN) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  facet_wrap(~ year) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.text=element_text(size=9),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# Distance
df <- data.frame(
  year = c("Distance"),
  x = factor(c(1)),
  est = c(CI_idx$V2[9]),
  upper = c(CI_idx$V3[9]),
  lower = c(CI_idx$V1[9])
)
  
p8=ggplot(df, aes(y = est, x = x, color = "red")) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  scale_y_continuous(labels = scaleFUN) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  facet_wrap(~ year) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.text=element_text(size=9),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# Association
df <- data.frame(
  year = c("Association"),
  x = factor(c(1)),
  est = c(CI_idx$V2[10]),
  upper = c(CI_idx$V3[10]),
  lower = c(CI_idx$V1[10])
)

p9=ggplot(df, aes(y = est, x = x, color = "red")) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  scale_y_continuous(labels = scaleFUN) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  facet_wrap(~ year) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.text=element_text(size=9),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# Game Fish Pig
df <- data.frame(
  year = rep("Household-level Variables", 6),
  x = factor(c(1,2,3,1,2,3)),
  group = c("Giver","Giver","Giver","Receiver","Receiver","Receiver"),
  est = c(CI_idx$V2[2], CI_idx$V2[3], CI_idx$V2[4], CI_idx$V2[5], CI_idx$V2[6], CI_idx$V2[7]),
  upper = c(CI_idx$V3[2], CI_idx$V3[3], CI_idx$V3[4], CI_idx$V3[5], CI_idx$V3[6], CI_idx$V3[7]),
  lower = c(CI_idx$V1[2], CI_idx$V1[3], CI_idx$V1[4], CI_idx$V1[5], CI_idx$V1[6], CI_idx$V1[7])
)

p16 = ggplot(df, aes(y = est, x = x, color = group)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(0.3), size = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.3), width = 0.1) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_discrete(
    name="",
    labels = c("Game", "Fish", "Pigs"),
    expand=c(0.1,0.1),
    position='bottom'
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  facet_wrap(~ year) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=9),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank())


library(patchwork)
layout3 <- "
AABBBB
CCDDEE
"
p0 + p16 + p7 + p8 + p9 + plot_layout(design = layout3) 



#'@_eta_Plots 
# ********
# Figure 2
# ********
idx_eta_linePlot <- as.data.frame(idx_eta_linePlot)
idx_eta_linePlot$idx <- c(1:5)
ggplot(data=idx_eta_linePlot, aes(x=idx, y=idx_eta_linePlot)) +
  geom_line(aes(color="red"))+
  scale_x_discrete(
    name="",
    limits=factor(c(1:5)),
    breaks=factor(c(1:5)),
    labels = c(expression(eta[1]), expression(eta[2]), expression(eta[3]), 
               expression(eta[4]), expression(eta[5])),
    expand=c(0.1,0.1),
    position='bottom'
  ) +
  ylim(-0.1,1) +
  ylab("Estimate") +
  theme_light() +
  geom_point() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", axis.text=element_text(size=12))



