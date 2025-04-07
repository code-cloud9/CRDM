source("source.R")

library(netregR)
library(sandwich)
library(lmtest)

###' * \beta setting *
beta = c(-0.5, 0.5, 0.5, 1)

###' * network size *
city_sizes <- c(50)

###' * total iteration number of e_ij for each sample size *
iter_num <- c(1:15000)

###' * total iteration number of a same design matrix X *
X_num <- 1000

design_supp <- function(n, bin_p2, mu3, sd3, mu4, sd4){
  
  ###
  ### Generate X1 Matrix and x1 vector
  ###
  X1.save <-  matrix(rnorm(n*n, -4, 1), n, n)
  diag(X1.save) <- NA
  elements1 <- as.vector(t(X1.save))
  x1.save <- na.omit(elements1)
  
  ###
  ### Generate X2 Matrix and x2 vector
  ###
  bino.temp <- as.matrix(rbinom(n, 1, bin_p2)) 
  X2.save <- bino.temp %*% t(bino.temp)
  diag(X2.save) <- NA
  elements2 <- as.vector(t(X2.save))
  x2.save <- na.omit(elements2)
  
  ###
  ### Generate X3 Matrix and x3 vector
  ###
  z3.temp <- as.matrix(rnorm(n, mu3, sd3)) 
  ones <- as.matrix(rep(1, n)) 
  X3.temp <- z3.temp %*% t(ones)
  X3.save <- abs(X3.temp - t(X3.temp))
  diag(X3.save) <- NA
  elements3 <- as.vector(t(X3.save))
  x3.save <- na.omit(elements3)
  
  ###
  ### Generate X4 Matrix and x4 vector
  ###
  X4.save <-  matrix(rnorm(n*n, mu4, sd4), n, n)
  diag(X4.save) <- NA
  elements4 <- as.vector(t(X4.save))
  x4.save <- na.omit(elements4)
  
  ###
  ### Write Output 
  ###  
  list(X1 = X1.save, X2 = X2.save, X3 = X3.save, X4 = X4.save,
       x1 = x1.save, x2 = x2.save, x3 = x3.save, x4 = x4.save)
  
}

set.seed(221)
for(iter in iter_num){  
  for(city_size in city_sizes){
    
    print(iter)
    
    #'@_e_ij
    set.seed(iter)
    Error_term.example <- error_term_gene_gamma(n = city_size, a1 = 0.2, b1 = 0.3,
                                                a2 = 0.2, b2 = 0.3, a0 = 0.2,
                                                b0 = 0.4, pho_ab = 0.9)
    
    #'@_X
    set.seed((iter-1) %/% X_num)
    X.example <- design_supp(n = city_size, 0.6, 1, 1/2, 1, 1)
    
    #'@_Lambda_ij
    Lambda_matr_and_vec.example <- lamb_gene_exp(n = city_size, beta, X.example, Error_term.example)
    
    #'@_Y
    Y_matr_and_vec.example <- Y_gene(n = city_size, Lambda_matr_and_vec.example)
    
    
    #' * Wenqin Exp link *
    
    #'@_PMLE
    PMLE_summar.example <- pseudo_estimators_exp(beta, X.example,
                                                 Y_matr_and_vec.example)
    
    #'@_eta_moment
    eta_6.moment <- para_6_estimate_exp(n = city_size,
                                        PMLE_summar.example$Beta_pmle, X.example,
                                        Y_matr_and_vec.example)
    
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
    eta_list_1 <- c(CV_s_fix, eta_6.moment$lam_6_est[2:5],0)
    AinvBAinv_1 <- AinvBAinv_exp_one_step(n = city_size, beta_est = PMLE_summar.example$Beta_pmle, lambda = eta_list_1, X_full = X.example)
    
    
    #'@_Var(beta)_naive
    eta_list_2 <- c(fit.empi_beta_hat$hyCvShorth, 0, 0, 0, 0, 0)
    AinvBAinv_2 <- AinvBAinv_exp_one_step(n = city_size, beta_est = PMLE_summar.example$Beta_pmle, lambda = eta_list_2, X_full = X.example)
    
    
    #'@_Var(beta)_oracle
    eta_list_3 <- c(2.0216837, 1.2101403, 0.9964923, 0.9964923, 0.5978954, 0)
    AinvBAinv_3 <- AinvBAinv_exp_one_step(n = city_size, beta_est = PMLE_summar.example$Beta_pmle, lambda = eta_list_3, X_full = X.example)
    
    
    
    #'@_Set_up_table
    write.out <- data.frame("beta0_mle" = PMLE_summar.example$Beta_pmle[1],
                            "beta1_mle" = PMLE_summar.example$Beta_pmle[2],
                            "beta2_mle" = PMLE_summar.example$Beta_pmle[3],
                            "beta3_mle" = PMLE_summar.example$Beta_pmle[4],
                        
                            "var0_mine_emp" = AinvBAinv_1[1,1],
                            "var1_mine_emp" = AinvBAinv_1[2,2],
                            "var2_mine_emp" = AinvBAinv_1[3,3],
                            "var3_mine_emp" = AinvBAinv_1[4,4],
                            
                            "var0_naive" = AinvBAinv_2[1,1],
                            "var1_naive" = AinvBAinv_2[2,2],
                            "var2_naive" = AinvBAinv_2[3,3],
                            "var3_naive" = AinvBAinv_2[4,4],
                            
                            "var0_oracle" = AinvBAinv_3[1,1],
                            "var1_oracle" = AinvBAinv_3[2,2],
                            "var2_oracle" = AinvBAinv_3[3,3],
                            "var3_oracle" = AinvBAinv_3[4,4],
                            
                            "X_class" = (iter-1) %/% X_num,
                            "Setting" = 4,
                            "iter_num" = iter)
    
    write.table(write.out, file = paste("settings.csv",
                                        sep = ","), append = T, sep = ',', row.names = F, col.names = F)
  }
}

