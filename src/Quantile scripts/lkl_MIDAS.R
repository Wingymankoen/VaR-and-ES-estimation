# Likelihood

lklfct_MIDAS <-
  function(BETAS,
           rt,
           RowsR,
           pval,
           EMPQNT,
           EMPES,
           penalty,
           es_approach,
           mu,
           rt_h ,
           d ) {
    dbeta_filter <- c()
    for (i in 1:d) {
      dbeta_filter[i] <- dbeta(i / d, shape1 = 1, shape2 = BETAS[3])
    }
    dbeta_filter_f <- dbeta_filter / sum(dbeta_filter)
    VaR <- VaR_midas_fct(BETAS, rt, d, dbeta_filter_f)
    
    if (! sum(is.na(VaR)) > 0) {
      if (es_approach == 1) {
        ES <- ES_1(BETAS, VaR)
      } else if (es_approach == 2) {
        ES <- ES_2(BETAS, rt_h, EMPES, RowsR, VaR, EMPQNT)
      } else if (es_approach == 3) {
        ES <- ES_3(BETAS, rt_h, EMPES, RowsR, VaR, EMPQNT, mu, pval)
      }
    }
    hfz <- rep(1, (RowsR - d))
    if ( sum(is.na(ES))==0 |  sum(is.na(VaR))==0 ) {
      score <-
        ((1-pval) / (mu[1:(RowsR-d)] - ES[1:(RowsR-d)])) * exp( - (rt_h[1:(RowsR-d)] - VaR[1:(RowsR-d)]) *(pval - (rt_h[1:(RowsR-d)] < VaR[1:(RowsR-d)]))  / (pval * ( mu[1:(RowsR-d)]-ES[1:(RowsR-d)] ) ))
    } else{
      score = NA
    }
    
    if ( sum(is.na(score)) > 0 |
        sum(score) == Inf | sum(score) == -Inf) {
      summed_hfz = penalty
    } else{
      summed_hfz = -sum(log(score))
    }
    
    # for (t in d:RowsR) {
    #   score_2 <-
    #       exp( - (rt_h[t] - VaR[t]) *(pval - (rt_h[t] < VaR[t]))  / (pval * ( mu[t]-ES[t]) ))
    #   score <-  ((1-pval) / (mu[t] - ES[t])) * score_2
    #   hfz[t] = log(score)
    # }
         # for (t in d:RowsR) {
    #   if (!is.nan(ES[t])  |  !is.nan(VaR[t])) {
    #     score_2 <-
    #       exp( - (rt_h[t] - VaR[t]) *(pval - (rt_h[t] < VaR[t]))  / (pval * ( mu[t]-ES[t]) ))
    #   } else{
    #     score_2 = 0
    #   }
    #   
    #   if (is.nan(ES[t]) | is.nan(score_2)) {
    #     score = 0.00001
    #   } else if (ES[t] == 0 | ES[t] > 0 | score_2 == 0) {
    #     score = 0.00001
    #   }  else{
    #     score <-  ((1-pval) / (mu[t] - ES[t])) * score_2
    #   }
    #   hfz[t] = log(score)
    # }
    # summed_hfz = -sum(hfz)
    # if (is.nan(summed_hfz) |
    #     summed_hfz == Inf | summed_hfz == -Inf) {
    #   summed_hfz = 1000000
    # }
    if(is.na(summed_hfz) | summed_hfz == Inf | summed_hfz == -Inf){
      summed_hfz<- penalty
    }
    return(summed_hfz)
  }
