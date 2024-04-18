# Conditional -------------------------------------------------------------

conditional <- function(rt, days_ahead, starts, prob, len_r, wd) {
  # Source files
  setwd(paste0(wd, '/src/Conditional scripts/'))
  files.sources = list.files()
  sapply(files.sources, source)
  setwd(wd)
  
  # initial parameters
  msgarch_daysahead <- as.integer(days_ahead)
  K <- 5000 # Number of simulations
  decay <- 0.95
  evt_boundary <- 0.15
  rol_window <- 500
  no_of_risk_measures <- 2
  # Model-specific
  arma_order <- c(0, 0)
  garch_order <- c(1, 1)
  includemean <- FALSE
  returns_arima_order <- c(1, 0, 0)
  
  ### Specifications
  # GARCH
  spec <- ugarchspec(
    mean.model = list(armaOrder = arma_order, include.mean = includemean),
    distribution = "norm",
    variance.model = list(model = "sGARCH", garch_order) 
  )
  # t-GARCH
  t_spec <- ugarchspec(
    mean.model = list(armaOrder = arma_order, include.mean = includemean),
    distribution = "sstd",
    variance.model = list(model = "sGARCH", garch_order)
  )
  # GJR-GARCH
  gjr_spec = ugarchspec(
    distribution.model = "sstd",
    variance.model = list(model = "gjrGARCH",
                          garch_order),
    mean.model = list(armaOrder = arma_order, include.mean = includemean)
  )
  
  # MS-GARCH
  spec_msgarch <- CreateSpec(
    variance.spec = list(model = c("sGARCH", "sGARCH")),
    distribution.spec = list(distribution = c("sstd", "sstd"))
  )
  
  # LOOP
  # Set up initials
  counter <- 0
  garch_df <- data.frame()
  tgarch_df <- data.frame()
  gjrgarch_df <- data.frame()
  msgarch_df <- data.frame()
  evtgarch_df <- data.frame()
  vol_g <- c()
  vol_tg <- c()
  vol_gjr <- c()
  vol_msg <- c()
  vol_evtg <- c()
  
  for (i in starts:len_r) {
    counter <- counter + 1
    print(paste0("Con %: ", round(100 * counter / len_r, 3)))
    
    #Rolling window
    if (i > rol_window) {
      tempreturns <- rt[(i - rol_window):i,]
    } else{
      tempreturns <- rt[1:i,]
    }
    
    # Demean the returns
    temp_returns_ar1 <-
      arima(tempreturns, order = returns_arima_order)
    coef_ar1 <- temp_returns_ar1$coef
    demeaned_returns <- temp_returns_ar1$residuals
    # Find the h-day ahead conditional mean
    con_mean_h <-
      cond_mean(tempreturns, days_ahead, returns_arima_order)
    
    # Regular GARCH
    def.fit1 <-
      ugarchfit(
        spec = spec,
        data = demeaned_returns,
        solver = 'hybrid',
        fit.control = list(rec.init = decay)
      )
    garch_stdres <- demeaned_returns / def.fit1@fit$sigma
    
    garch_df <- rbind(
      garch_df,
      filterfct(
        prob = prob,
        no_sim = K,
        i = length(demeaned_returns),
        look_ahead = days_ahead,
        std_res = garch_stdres,
        variables = def.fit1@fit$coef,
        est_sig =  def.fit1@fit$sigma,
        returns = demeaned_returns,
        model_type = 'garch',
        variables_ar = coef_ar1
      )
    )
    vol_g[i - starts + 1] <-
      def.fit1@fit$sigma[length(demeaned_returns)]
    
    # t-GARCH
    def.fit2 = ugarchfit(
      spec = t_spec,
      data = demeaned_returns,
      solver = 'hybrid',
      fit.control = list(rec.init = decay)
    )
    tgarch_stdres <- demeaned_returns / sqrt(def.fit2@fit$var)
    tgarch_df <- rbind(
      tgarch_df,
      filterfct(
        prob = prob,
        no_sim = K,
        i = length(demeaned_returns),
        look_ahead = days_ahead,
        std_res = tgarch_stdres,
        variables = def.fit2@fit$coef,
        est_sig = def.fit2@fit$sigma,
        returns = demeaned_returns,
        model_type = 'garch',
        variables_ar = coef_ar1
      )
    )
    vol_tg[i - starts + 1] <-
      def.fit2@fit$sigma[length(demeaned_returns)]
    
    # GJR-garch
    gjr_fit = ugarchfit(
      spec = gjr_spec,
      data = demeaned_returns,
      solver = 'hybrid',
      fit.control = list(rec.init = decay)
    )
    gjrgarch_stdres <- demeaned_returns / sqrt(gjr_fit@fit$var)
    gjrgarch_df <- rbind(
      gjrgarch_df,
      filterfct(
        prob = prob,
        no_sim = K,
        i = length(demeaned_returns),
        look_ahead = days_ahead,
        std_res = gjrgarch_stdres,
        variables = gjr_fit@fit$coef,
        est_sig = gjr_fit@fit$sigma,
        returns = demeaned_returns,
        model_type = 'gjr-garch',
        variables_ar = coef_ar1
      )
    )
    vol_gjr[i - starts + 1] <-
      gjr_fit@fit$sigma[length(demeaned_returns)]
    
    #MS-GARCH
    msgarch_fit <-
      try(FitML(spec_msgarch, data = demeaned_returns))
    
    if ("try-error" %in% class(msgarch_fit)) {
      msgarch_fit <-
        FitML(spec_msgarch, data = tempreturns[1:(length(demeaned_returns) - 1)])
    }
    risk <-
      Risk(
        object = msgarch_fit,
        nahead = msgarch_daysahead ,
        do.cumulative = TRUE,
        alpha = c(prob)
      )
    msg_update <- c(risk$VaR[days_ahead, 1], risk$ES[days_ahead, 1])
    msg_update <- con_mean_h + msg_update
    msgarch_df <- rbind(msgarch_df, msg_update)
    
    vol_msg[i - starts + 1] <-
      Volatility(object = msgarch_fit)[length(demeaned_returns)]
    
    # EVT
    # Using t-GARCH residuals
    thres <- -quantile(tgarch_stdres, evt_boundary)
    gpdresid <- gpd(-tgarch_stdres, threshold = thres)
    evtgarch_df <- rbind(
      evtgarch_df,
      filterfct(
        prob = prob,
        no_sim = K,
        i = length(demeaned_returns),
        look_ahead = days_ahead,
        std_res = tgarch_stdres,
        variables = def.fit2@fit$coef,
        est_sig = def.fit2@fit$sigma,
        returns = demeaned_returns,
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
      rep('GARCH', no_of_risk_measures),
      rep('t-GARCH', no_of_risk_measures),
      rep('GJR-GARCH', no_of_risk_measures),
      rep('MS-GARCH', no_of_risk_measures),
      rep('EVT-GARCH', no_of_risk_measures)
    )
  
  vol_df <- data.frame(vol_g,
                       vol_tg,
                       vol_gjr,
                       vol_msg)
  return(list(outputdf, vol_df))
}
