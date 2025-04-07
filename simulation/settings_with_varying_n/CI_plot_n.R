plotCP <- function(X, beta, data, ylims,yicept){
  
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  

  
  colnames(data) <- c("beta0_mle",
                      "beta1_mle",
                      "beta2_mle",
                      "beta3_mle",
                      
                      "var0_mine_thm",
                      "var1_mine_thm",
                      "var2_mine_thm",
                      "var3_mine_thm",
                      
                      "var0_mine_emp",
                      "var1_mine_emp",
                      "var2_mine_emp",
                      "var3_mine_emp",
                      
                      "var0_naive",
                      "var1_naive",
                      "var2_naive",
                      "var3_naive",
                      
                      "var0_oracle",
                      "var1_oracle",
                      "var2_oracle",
                      "var3_oracle",
                      
                      "var0_oracle_oracle",
                      "var1_oracle_oracle",
                      "var2_oracle_oracle",
                      "var3_oracle_oracle",
                      
                      "X_class",
                      "city_size",
                      "iter_num")
  
  smry <- data %>%
    group_by(city_size, X_class) %>%
    summarise(cp_mine_thm_0 = mean(ifelse(beta0_mle >= beta[1] - qnorm(0.975) * sqrt(var0_mine_thm) & beta0_mle <= beta[1] + qnorm(0.975) * sqrt(var0_mine_thm), 1, 0)),
              cp_mine_thm_1 = mean(ifelse(beta1_mle >= beta[2] - qnorm(0.975) * sqrt(var1_mine_thm) & beta1_mle <= beta[2] + qnorm(0.975) * sqrt(var1_mine_thm), 1, 0)),
              cp_mine_thm_2 = mean(ifelse(beta2_mle >= beta[3] - qnorm(0.975) * sqrt(var2_mine_thm) & beta2_mle <= beta[3] + qnorm(0.975) * sqrt(var2_mine_thm), 1, 0)),
              cp_mine_thm_3 = mean(ifelse(beta3_mle >= beta[4] - qnorm(0.975) * sqrt(var3_mine_thm) & beta3_mle <= beta[4] + qnorm(0.975) * sqrt(var3_mine_thm), 1, 0)),
              
              cp_mine_emp_0 = mean(ifelse(beta0_mle >= beta[1] - qnorm(0.975) * sqrt(var0_mine_emp) & beta0_mle <= beta[1] + qnorm(0.975) * sqrt(var0_mine_emp), 1, 0)),
              cp_mine_emp_1 = mean(ifelse(beta1_mle >= beta[2] - qnorm(0.975) * sqrt(var1_mine_emp) & beta1_mle <= beta[2] + qnorm(0.975) * sqrt(var1_mine_emp), 1, 0)),
              cp_mine_emp_2 = mean(ifelse(beta2_mle >= beta[3] - qnorm(0.975) * sqrt(var2_mine_emp) & beta2_mle <= beta[3] + qnorm(0.975) * sqrt(var2_mine_emp), 1, 0)),
              cp_mine_emp_3 = mean(ifelse(beta3_mle >= beta[4] - qnorm(0.975) * sqrt(var3_mine_emp) & beta3_mle <= beta[4] + qnorm(0.975) * sqrt(var3_mine_emp), 1, 0)),
              
              cp_naive_0 = mean(ifelse(beta0_mle >= beta[1] - qnorm(0.975) * sqrt(var0_naive) & beta0_mle <= beta[1] + qnorm(0.975) * sqrt(var0_naive), 1, 0)),
              cp_naive_1 = mean(ifelse(beta1_mle >= beta[2] - qnorm(0.975) * sqrt(var1_naive) & beta1_mle <= beta[2] + qnorm(0.975) * sqrt(var1_naive), 1, 0)),
              cp_naive_2 = mean(ifelse(beta2_mle >= beta[3] - qnorm(0.975) * sqrt(var2_naive) & beta2_mle <= beta[3] + qnorm(0.975) * sqrt(var2_naive), 1, 0)),
              cp_naive_3 = mean(ifelse(beta3_mle >= beta[4] - qnorm(0.975) * sqrt(var3_naive) & beta3_mle <= beta[4] + qnorm(0.975) * sqrt(var3_naive), 1, 0)),
              
              cp_oracle_0 = mean(ifelse(beta0_mle >= beta[1] - qnorm(0.975) * sqrt(var0_oracle) & beta0_mle <= beta[1] + qnorm(0.975) * sqrt(var0_oracle), 1, 0)),
              cp_oracle_1 = mean(ifelse(beta1_mle >= beta[2] - qnorm(0.975) * sqrt(var1_oracle) & beta1_mle <= beta[2] + qnorm(0.975) * sqrt(var1_oracle), 1, 0)),
              cp_oracle_2 = mean(ifelse(beta2_mle >= beta[3] - qnorm(0.975) * sqrt(var2_oracle) & beta2_mle <= beta[3] + qnorm(0.975) * sqrt(var2_oracle), 1, 0)),
              cp_oracle_3 = mean(ifelse(beta3_mle >= beta[4] - qnorm(0.975) * sqrt(var3_oracle) & beta3_mle <= beta[4] + qnorm(0.975) * sqrt(var3_oracle), 1, 0)),
              
              cp_oracle_oracle_0 = mean(ifelse(beta0_mle >= beta[1] - qnorm(0.975) * sqrt(var0_oracle_oracle) & beta0_mle <= beta[1] + qnorm(0.975) * sqrt(var0_oracle_oracle), 1, 0)),
              cp_oracle_oracle_1 = mean(ifelse(beta1_mle >= beta[2] - qnorm(0.975) * sqrt(var1_oracle_oracle) & beta1_mle <= beta[2] + qnorm(0.975) * sqrt(var1_oracle_oracle), 1, 0)),
              cp_oracle_oracle_2 = mean(ifelse(beta2_mle >= beta[3] - qnorm(0.975) * sqrt(var2_oracle_oracle) & beta2_mle <= beta[3] + qnorm(0.975) * sqrt(var2_oracle_oracle), 1, 0)),
              cp_oracle_oracle_3 = mean(ifelse(beta3_mle >= beta[4] - qnorm(0.975) * sqrt(var3_oracle_oracle) & beta3_mle <= beta[4] + qnorm(0.975) * sqrt(var3_oracle_oracle), 1, 0))
              
    )
  
  
  
  
  smryMelt0 <- melt(smry, id = c("city_size"), measure.vars = c("cp_mine_emp_0", "cp_naive_0", "cp_oracle_0"))
  
  
  p0 <- ggplot(smryMelt0, aes(x = factor(city_size), y = value, color = variable)) +
    geom_boxplot(outlier.size = 0.3) +
    geom_hline(yintercept = yicept[1], linetype = "dashed",
               color = "black", size = 0.5) +
    scale_color_manual(values = c("#FFCC99", "#99CCFF", "#999999"),
                       name = "",
                       labels = c("Our model", "Naive", "Oracle"),
                       drop = FALSE) +
    ylim(ylims[1], ylims[2]) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 11)) +
    xlab("Number of nodes") +
    ylab("Coverage probability")
  
  if(X == 1) {
    p0 <- p0 + ggtitle(~ paste(x['1ij']))
  } else{
    p0 <- p0 + ggtitle("1")
  }
  
  
  
  smryMelt1 <- melt(smry, id = c("city_size"), measure.vars = c("cp_mine_emp_1", "cp_naive_1", "cp_oracle_1"))
  
  
  p1 <- ggplot(smryMelt1, aes(x = factor(city_size), y = value, color = variable)) +
    geom_boxplot(outlier.size = 0.3) +
    geom_hline(yintercept = yicept[2], linetype = "dashed",
               color = "black", size = 0.5) +
    scale_color_manual(values = c("#FFCC99", "#99CCFF", "#999999"),
                       labels = c("Our model", "Our model (naive)", "Our model (oracle)"),
                       name = "",
                       drop = FALSE) +
    ylim(ylims[1], ylims[2]) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 11)) +
    xlab("Number of nodes") +
    ylab("Coverage probability")
  
  if(X == 3) {
    p1 <- p1 + ggtitle(~ paste(x['2ij']))
  } else{
    p1 <- p1 + ggtitle(~ paste(x['2i'],x['2j'],)) 
  }
  
  smryMelt2 <- melt(smry, id = c("city_size"), measure.vars = c("cp_mine_emp_2", "cp_naive_2", "cp_oracle_2"))
  
  p2 <- ggplot(smryMelt2, aes(x = factor(city_size), y = value, color = variable)) +
    geom_boxplot(outlier.size = 0.3) +
    geom_hline(yintercept = yicept[3], linetype = "dashed",
               color = "black", size = 0.5) +
    scale_color_manual(values = c("#FFCC99", "#99CCFF", "#999999"),
                       name = "",
                       labels = c("Our model", "Our model (naive)", "Our model (oracle)"),
                       drop = FALSE) +
    ylim(ylims[1], ylims[2]) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none",
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 11)) +
    xlab("Number of nodes") +
    ylab("Coverage probability")
  
  if(X == 3) {
    p2 <- p2 + ggtitle(~ paste(x['3ij']))
  } else{
    p2 <- p2 + ggtitle(~ paste("|", x['3i'] - x['3j'], "|"))
  }
  
  smryMelt3 <- melt(smry, id = c("city_size"), measure.vars = c("cp_mine_emp_3", "cp_naive_3", "cp_oracle_oracle_3"))
  
  p3 <- ggplot(smryMelt3, aes(x = factor(city_size), y = value, color = variable)) +
    geom_boxplot(outlier.size = 0.3) +
    geom_hline(yintercept = yicept[4], linetype = "dashed",
               color = "black", size = 0.5) +
    scale_color_manual(values = c("#FFCC99", "#99CCFF", "#999999"),
                       name = "",
                       labels = c("Our model", "Our model (naive)", "Our model (oracle)"),
                       drop = FALSE) +
    ylim(ylims[1], ylims[2]) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 11)) +
    xlab("Number of nodes") +
    ylab("Coverage probability")
  
  
  if(X == 1) {
    p3 <- p3 + ggtitle(~ paste(x['4ij']))
  } else{
    p3 <- p3 + ggtitle(~ paste("|", x['4ij'], "|"))
  }
  

  
  ggarrange(p0, p1, p2, p3, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")

  
}

# ********
# Figure 2
# ********
# Design-main-truncN-iid
plotCP(X = 1, beta = c(1, -0.5, -0.5, -1), data = read.csv("setting1.csv", header = FALSE), ylims = c(0.55, 1), yicept = c(0.95,0.95,0.95,0.95))

# Design-main-truncN-exch
plotCP(X = 1, beta = c(1, -0.5, -0.5, -1), data = read.csv("setting2.csv", header = FALSE), ylims = c(0.2, 1),  yicept = c(0.95,0.95,0.95,0.95))

# **********
# Figure S.1
# **********
# Design-main-gamma-iid
plotCP(X = 1, beta = c(1, -0.5, -0.5, -1), data = read.csv("setting3.csv", header = FALSE), ylims = c(0.55, 1), yicept = c(0.95,0.95,0.95,0.95))

# Design-main-gamma-exch
plotCP(X = 1, beta = c(1, -0.5, -0.5, -1), data = read.csv("setting4.csv", header = FALSE), ylims = c(0.2, 1),  yicept = c(0.95,0.95,0.95,0.95))

# **********
# Figure S.2
# **********
# Design-supp-truncN-iid
plotCP(X = 1, beta = c(-0.5, 0.5, 0.5, 1), data = read.csv("setting5.csv", header = FALSE), ylims = c(0.55, 1), yicept = c(0.95,0.95,0.95,0.95))

# Design-supp-truncN-exch
plotCP(X = 1, beta = c(-0.5, 0.5, 0.5, 1), data = read.csv("setting6.csv", header = FALSE), ylims = c(0.2, 1), yicept = c(0.95,0.95,0.95,0.95))

# **********
# Figure S.3
# **********
# Design-supp-Gamma-iid
plotCP(X = 1, beta = c(-0.5, 0.5, 0.5, 1), data = read.csv("setting7.csv", header = FALSE), ylims = c(0.55, 1), yicept = c(0.95,0.95,0.95,0.95))

# Design-supp-Gamma-exch
plotCP(X = 1, beta = c(-0.5, 0.5, 0.5, 1), data = read.csv("setting8.csv", header = FALSE), ylims = c(0.2, 1), yicept = c(0.95,0.95,0.95,0.95))

