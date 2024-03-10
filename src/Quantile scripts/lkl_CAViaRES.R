
lklfct_CAViaRES <-
  function(BETAS,
           rt,
           RowsR,
           pval,
           EMPQNT,
           EMPES,
           penalty,
           es_approach,
           mu) {
    VaR <- caviar_SAV(BETAS, rt, EMPQNT, RowsR)
    
    if (es_approach == 1) {
      ES <- ES_1(BETAS, VaR)
    } else if (es_approach == 2) {
      ES <- ES_2(BETAS, ret = rt, EMPES, RowsR, VaR, EMPQNT)
    } else if (es_approach == 3) {
      ES <- ES_3(BETAS, ret = rt, EMPES, RowsR, VaR, EMPQNT, mu, pval)
    }
    if ( sum(is.na(ES))==0 |  sum(is.na(VaR))==0 ) {
      score <-
        ((1-pval) / (mu[1:RowsR] - ES[1:(RowsR-1)])) * exp( - (rt[2:RowsR] - VaR[1:(RowsR-1)]) *(pval - (rt[2:RowsR] < VaR[1:(RowsR-1)]))  / (pval * ( mu[1:RowsR]-ES[1:(RowsR-1)]) ))
    } else{
      score = NA
    }

    if ( sum(is.na(score)) > 0 |
        sum(score) == Inf | sum(score) == -Inf) {
      summed_hfz = penalty
    } else{
      summed_hfz = -sum(log(score))
    }
    
    if(is.na(summed_hfz) | summed_hfz == Inf | summed_hfz == -Inf){
      summed_hfz<- penalty
    }
    return(summed_hfz)
  }