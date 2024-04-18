cond_mean <- function(returns, days_ahead, returns_arima_order){
  returns_h <- c()
  # Calculate h-day forward returns
  for (j in 1:length(returns)) {
    returns_h[j] <- sum(returns[(j + 1):(j + days_ahead)])
  }
  ar1 <- arima(returns_h, order = returns_arima_order)
  con_mean_h <-
    ar1$coef[2] + ar1$coef[1] * ar1$residuals[length(ar1$residuals) - days_ahead]
  return(con_mean_h)
  
  } 