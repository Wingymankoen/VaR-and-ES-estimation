# ES shortfall estimation, 3 approaches
ES_1 <- function(BETAS, VaR) {
  ES <- (1 + exp(BETAS[4])) * VaR
  return(ES)
}

ES_2 <- function(BETAS,
                 ret,
                 EMPES,
                 RowsR,
                 VaR,
                 EMPQNT) {
  ES <- c(EMPES)
  j_t <- EMPQNT - EMPES
  
  for (i in 1:(RowsR-1)) {
    if (ret[i] < VaR[i]) {
      j_t = BETAS[4] + BETAS[5] * (VaR[i] - ret[i]) + BETAS[6] * j_t
    } else{
      j_t = j_t
    }
    ES[i + 1] <- VaR[i + 1] - j_t
  }
  return(ES)
}

ES_3 <- function(BETAS,
                 ret,
                 EMPES,
                 RowsR,
                 VaR,
                 EMPQNT,
                 mu_3,
                 pval) {
  
  ES <- c(EMPES)
  j_t <- EMPQNT - EMPES
  for (i in 1:(RowsR-1)) {
    if(is.na(VaR[i])){
      VaR[i]<- VaR[i-1] 
    }
  }
  for (i in 1:(RowsR-1)) {
    if (ret[i] < VaR[i]) {
      .e2 <- j_t + mu_3[i] - VaR[i]
      nom<-  pval * (pval - ret[i] < VaR[i]) * (ret[i] - VaR[i])/(pval * .e2)^2 - 1/.e2
      denom<- 1/.e2^2 - 2 * (pval^3 * (pval - ret[i] < VaR[i]) * .e2 * (ret[i] - VaR[i])/(pval * 
                                                                      .e2)^4)
      j_t = BETAS[4] + BETAS[5] * (VaR[i] - ret[i]) + BETAS[6] * j_t + BETAS[7]* nom/-denom
    } else{
      j_t = j_t
    }
    ES[i + 1] <- VaR[i + 1] - j_t
  }
  return(ES)
  
}
