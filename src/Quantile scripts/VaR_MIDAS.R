VaR_midas_fct <- function(BETAS, rt, d, dbeta_filter_f) {
  VaR <- rep(0, d)
  for (t in d:length(rt)) {
    VaR[t] <-
      BETAS[1] + BETAS[2] * sum(dbeta_filter_f * abs(rt[t:(t - d + 1)]))
  }
  return(VaR)
}