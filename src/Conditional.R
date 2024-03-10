

# Filter ------------------------------------------------------------------
filterfct <-
  function(prob,
           no_sim,
           i = i,
           look_ahead,
           std_res,
           variables,
           resids,
           est_sig,
           returns,
           model_type,
           variables_ar,
           evt = FALSE,
           u = NULL,
           variables_gpd = NULL) {
    sampled_values_list <-
      lapply(1:no_sim, function(x)
        sample(
          std_res,
          size = look_ahead,
          replace = TRUE,
          prob = NULL
        ))
    temp_resid <-
      matrix(unlist(sampled_values_list),
             ncol = look_ahead,
             byrow = TRUE)
    init_sig <- rep(est_sig[i] ^ 2, no_sim)
    init_res <- rep(resids[i], no_sim)
    init_return <- rep(returns[i], no_sim)
    current_factors <-
      matrix(c(init_sig, init_res, init_return), ncol = 3)
    daily_returns <- matrix(NA, ncol = look_ahead, nrow = no_sim)
    if (evt) {
      temp_resid[temp_resid < -u] <-
        -rgpd(sum(temp_resid < -u),
              xi = variables_gpd[1],
              mu = u,
              beta = variables_gpd[2])
    }
    
    if (model_type == 'garch') {
      for (j in 1:look_ahead) {
        current_factors[, 1] <-
          variables[2] * current_factors[, 2] ^ 2 + variables[3] * current_factors[, 1]
        current_factors[, 2] <-
          temp_resid[, j] * sqrt(current_factors[, 1])
        current_factors[, 3] <-
          variables_ar[2] +   current_factors[, 3] * variables_ar[1] + current_factors[, 2]
        daily_returns[, j] <- current_factors[, 3]
      }
      
    } else if (model_type == 'gjr-garch') {
      for (j in 1:look_ahead) {
        current_factors[, 1] <-
          variables[2] * current_factors[, 2] ^ 2 + variables[3] * current_factors[, 1] + variables[4] * (current_factors[, 2] <
                                                                                                            0) ^ 2
        current_factors[, 2] <-
          temp_resid[, j] * sqrt(current_factors[, 1])
        current_factors[, 3] <-
          variables_ar[2] +   current_factors[, 3] * variables_ar[1]  + current_factors[, 2]
        daily_returns[, j] <- current_factors[, 3]
      }
    }
    results <- rowSums(daily_returns)
    
    VaRgarch <- quantile(results, prob, na.rm = TRUE)
    ESgarch <- mean(results[results < VaRgarch])
    return(cbind(VaRgarch, ESgarch))
  }


# Conditional -------------------------------------------------------------

conditional <- function(rt, days_ahead, starts, prob, len_r) {
  # initial parameters
  K <- 5000 # Number of simulations
  nvol <- 1
  decay <- 0.95
  evt_boundary <- 0.15
  rol_window <- 500
  
  ### Specifications
  # GARCH
  spec <- ugarchspec(
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution = "norm",
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)) ,
    fixed.pars = list('omega' = 0)
  )
  # t-GARCH
  t_spec <- ugarchspec(
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution = "sstd",
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)) ,
    fixed.pars = list('omega' = 0)
  )
  # GJR-GARCH
  gjr_spec = ugarchspec(
    distribution.model = "sstd",
    variance.model = list(model = "gjrGARCH",
                          garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    fixed.pars = list('omega' = 0)
  )
  
  # MS-GARCH
  spec_msgarch <- CreateSpec(
    variance.spec = list(model = c("sGARCH", "sGARCH")),
    distribution.spec = list(distribution = c("sstd", "sstd")),
    switch.spec = list(do.mix = FALSE, K = NULL),
    constraint.spec = list(
      fixed = list(alpha0_1 = 1e-07, alpha0_2 = 1e-07),
      regime.const = NULL
    ),
    prior = list(mean = list(), sd = list())
  )
  
  # LOOP
  # Set up initials
  garch_df <- data.frame()
  tgarch_df <- data.frame()
  gjrgarch_df <- data.frame()
  msgarch_df <- data.frame()
  evtgarch_df <- data.frame()
  counter <- 0
  vol_g <- c()
  vol_tg <- c()
  vol_gjr <- c()
  vol_msg <- c()
  vol_evtg <- c()
  coef_fit <- c()
  for (i in starts:len_r) {
    counter <- counter + 1
    print(paste0("Con %: ", round(100 * counter / len_r, 3)))
    #Rolling window
    if (i > rol_window) {
      tempdat <- rt[(i - rol_window):i,]
    } else{
      tempdat <- rt[1:i,]
    }
    temp_returns_ar1 <- arima(tempdat, order = c(1, 0, 0))
    coef_ar1 <- temp_returns_ar1$coef
    tempresid <- temp_returns_ar1$residuals
    returns_h <- c()
    # Calculate h-day forward returns
    for (j in 1:length(tempdat)) {
      returns_h[j] <- sum(tempdat[(j + 1):(j + days_ahead)])
    }
    temp_h_ar1 <- arima(returns_h, order = c(1, 0, 0))
    con_mean_h <-
      temp_h_ar1$coef[2] + temp_h_ar1$coef[1] * temp_h_ar1$residuals[length(temp_h_ar1$residuals) -
                                                                       days_ahead]
    
    # Regular GARCH
    def.fit1 <-
      ugarchfit(
        spec = spec,
        data = tempresid,
        solver = 'hybrid',
        fit.control = list(rec.init = decay)
      )
    garch_stdres <- tempresid / def.fit1@fit$sigma
    
    garch_df <- rbind(
      garch_df,
      filterfct(
        prob = prob,
        no_sim = K,
        i = length(tempresid),
        look_ahead = days_ahead,
        std_res = garch_stdres,
        variables = def.fit1@fit$coef ,
        resids = tempresid,
        est_sig =  def.fit1@fit$sigma,
        returns = tempresid,
        model_type = 'garch',
        variables_ar = coef_ar1
      )
    )
    vol_g[i - starts + 1] <- def.fit1@fit$sigma[length(tempresid)]
    
    # t-GARCH
    def.fit2 = ugarchfit(
      spec = t_spec,
      data = tempresid,
      solver = 'hybrid',
      fit.control = list(rec.init = decay)
    )
    tgarch_stdres <- tempresid / sqrt(def.fit2@fit$var)
    tgarch_df <- rbind(
      tgarch_df,
      filterfct(
        prob = prob,
        no_sim = K,
        i = length(tempresid),
        look_ahead = days_ahead,
        std_res = tgarch_stdres,
        variables = def.fit2@fit$coef ,
        resids = tempresid,
        est_sig = def.fit2@fit$sigma  ,
        returns = tempresid,
        model_type = 'garch',
        variables_ar = coef_ar1
      )
    )
    vol_tg[i - starts + 1] <- def.fit2@fit$sigma[length(tempresid)]
    
    coef_fit[i - starts + i] <- def.fit2@fit$coef[2]
    # GJR-garch
    gjr_fit = ugarchfit(
      spec = gjr_spec,
      data = tempresid,
      solver = 'hybrid',
      fit.control = list(rec.init = decay)
    )
    gjrgarch_stdres <- tempresid / sqrt(gjr_fit@fit$var)
    gjrgarch_df <- rbind(
      gjrgarch_df,
      filterfct(
        prob = prob,
        no_sim = K,
        i = length(tempresid),
        look_ahead = days_ahead,
        std_res = gjrgarch_stdres,
        variables = gjr_fit@fit$coef ,
        resids = tempresid,
        est_sig = gjr_fit@fit$sigma  ,
        returns = tempresid,
        model_type = 'gjr-garch',
        variables_ar = coef_ar1
      )
    )
    vol_gjr[i - starts + 1] <- gjr_fit@fit$sigma[length(tempresid)]
    
    #MS-GARCH
    msgarch_fit <-
      try(FitML(spec_msgarch, data = tempresid))
    
    if ("try-error" %in% class(msgarch_fit)) {
      msgarch_fit <-
        FitML(spec_msgarch, data = tempdat[1:(length(tempresid) - 1)])
    }
    risk <-
      Risk(
        object = msgarch_fit,
        nahead = 1L ,
        do.cumulative = TRUE,
        alpha = c(prob)
      )
    msg_update <- c(risk$VaR[days_ahead, 1], risk$ES[days_ahead, 1])
    msg_update <- con_mean_h + msg_update
    msgarch_df <- rbind(msgarch_df, msg_update)
    
    vol_msg[i - starts + 1] <-
      Volatility(object = msgarch_fit)[length(tempresid)]
    
    # EVT
    # Using t-GARCH residuals
    thres <- -quantile(tgarch_stdres, evt_boundary)
    gpdresid <- gpd(-tgarch_stdres, threshold = thres)
    evtgarch_df <- rbind(
      evtgarch_df,
      filterfct(
        prob = prob,
        no_sim = K,
        i = length(tempresid),
        look_ahead = days_ahead,
        std_res = tgarch_stdres,
        variables = def.fit2@fit$coef ,
        resids = tempresid,
        est_sig = def.fit2@fit$sigma  ,
        returns = tempresid,
        model_type = 'garch',
        variables_ar = coef_ar1,
        u = thres,
        variables_gpd = gpdresid$par.ests,
        evt = TRUE
      )
    )
  }
  
  outputdf <-
    cbind(garch_df, tgarch_df, gjrgarch_df, msgarch_df, evtgarch_df)
  names(outputdf) <-
    c(
      rep('GARCH', 2),
      rep('t-GARCH', 2),
      rep('GJR-GARCH', 2),
      rep('MS-GARCH', 2),
      rep('EVT-GARCH', 2)
    )
  
  vol_df <- data.frame(vol_g,
                       vol_tg,
                       vol_gjr,
                       vol_msg)
  return(list(outputdf, vol_df))
}
