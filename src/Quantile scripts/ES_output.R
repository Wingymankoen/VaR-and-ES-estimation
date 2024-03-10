ES_output <-
  function(setup,
           BETAS,
           VaR,
           data,
           EMPES,
           EMPQNT,
           RowsR,
           mu,
           pval) {
    if (setup == 1) {
      ES <- ES_1(BETAS, VaR)
    } else if (setup == 2) {
      ES <-
        ES_2(
          BETAS,
          ret = data,
          EMPES = EMPES,
          RowsR = RowsR,
          VaR,
          EMPQNT = EMPQNT
        )
    } else{
      ES <- ES_3(
        BETAS,
        ret = data,
        EMPES = EMPES,
        RowsR = RowsR,
        VaR,
        EMPQNT = EMPQNT,
        mu = mu,
        pval = pval
      )
    }

    return(ES)
  }