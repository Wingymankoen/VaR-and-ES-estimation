# VaR estimation CAViaR-ES
caviar_SAV <- function(BETAS, rt, EMPQNT, RowsR) {
  VaR <- c(EMPQNT)
  for (i in 1:RowsR) {
    VaR[i + 1] <-
      BETAS[1] +  BETAS[2] * VaR[i] + BETAS[3] * abs(rt[i])
  }
  return(VaR)
}