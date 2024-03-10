


quantile_fct <-
  function(orig_returns,
           wd,
           len_r,
           prob,
           days_ahead,
           starts) {
    # Parameters
    ## General
    d <- 50 # number of lagged daily returns to take into account
    setup <- 1# Setup can be 1,2,3
    penalty <- 1e6
    re_estimate <- 1 # Changed
    size <- starts
    total_trys <- 1000
    nInitialCond <- 10
    met <- "L-BFGS-B"
    control <- list(maxit = 100)
    
    #load quantile formulas
    setwd(paste0(wd, '/src/Quantile scripts/'))
    files.sources = list.files()
    sapply(files.sources, source)
    
    # Set returns CAViar-ES
    returns_ar1 <-
      arima(orig_returns[1:len_r, ], order = c(1, 0, 0))
    returns_1day <- returns_ar1$residuals
    con_mean_1day <-
      returns_ar1$coef[2] + returns_ar1$coef[1] * orig_returns[1:len_r, ]
    
    # Calculate h-day forward returns MIDAS
    returns_h <- c()
    for (i in 1:(len_r - days_ahead)) {
      returns_h[i] <- sum(orig_returns[(i + 1):(i + days_ahead), ])
    }
    returns_h_ar1 <- arima(returns_h, order = (c(1, 0, 0)))
    returns_h <- returns_h_ar1$residuals
    
    temp_h <- c()
    for (i in 1:(len_r - days_ahead)) {
      temp_h[i] <- sum(orig_returns[(i + 1 - 1):(i + days_ahead - 1), ])
    }
    temp_h_ar1 <- arima(temp_h, order = c(1, 0, 0))
    con_mean_h <- temp_h_ar1$coef[2] + temp_h_ar1$coef[1] * temp_h
    
    # set a counter
    counter <- 0
    
    # Loop --------------------------------------------------------------------
    
    for (i in (len_r - size)) {
      counter <- counter + 1
      print(paste0("Quantile %: ", round(100 * counter / len_r, 3)))
      output <- matrix(NA,
                       nrow = (len_r - size - days_ahead),
                       ncol = 4)
      # Initialise data CES
      tempdat_1day <- returns_1day[1:(i + size)]
      l_tempdat <- length(tempdat_1day)
      emp_qnt <-
        quantile(tempdat_1day[(l_tempdat - size):l_tempdat] , prob)
      emp_es <-  mean(tempdat_1day[tempdat_1day < emp_qnt])
      
      #### CAVIAR-ES ####
      if (i %% re_estimate == 0) {
        ## Caviar-ES parameters
        init_ces <- sv_ces(setup, total_trys)
        lowb_ces <- lowerb_ces(setup)
        upb_ces <- upperb_ces(setup)
        
        ####### CAViaR-ES ######
        RQfval_ces <-
          apply(
            init_ces,
            1,
            lklfct_CAViaRES,
            rt = tempdat_1day,
            RowsR = l_tempdat,
            pval = prob,
            EMPQNT = emp_qnt,
            EMPES = emp_es,
            mu = con_mean_1day,
            penalty = penalty,
            es_approach = setup
          )
        BestInitialCond_ces <-
          init_ces[order(RQfval_ces), ][1:nInitialCond, ]
        RQoptim_ces <- cbind(NA, BestInitialCond_ces, NA)
        
        for (j in 1:nrow(BestInitialCond_ces)) {
          # Changed
          # initial optimization procedure
          diff_ces <- 1
          par_to_upd <- BestInitialCond_ces[j,]
          rounds <- 0
          while (diff_ces > 0.1 & rounds < 10) {
            rounds <- rounds + 1
            old_par <- par_to_upd
            vOptim_ces <- optim(
              par = par_to_upd,
              # Changed [j,]
              lklfct_CAViaRES,
              rt = tempdat_1day,
              RowsR = l_tempdat,
              pval = prob,
              EMPQNT = emp_qnt,
              EMPES = emp_es,
              penalty = penalty,
              es_approach = setup,
              mu = con_mean_1day,
              method = "Nelder-Mead",
              control = control
            )
            par_to_upd <- vOptim_ces$par
            # 2nd part
            vOptim_ces <- optim(
              par = par_to_upd,
              # Changed [j,]
              lklfct_CAViaRES,
              rt = tempdat_1day,
              RowsR = l_tempdat,
              pval = prob,
              EMPQNT = emp_qnt,
              EMPES = emp_es,
              penalty = penalty,
              es_approach = setup,
              mu = con_mean_1day,
              method = met,
              lower = lowb_ces,
              upper = upb_ces,
              control = control
            )
            par_to_upd <- vOptim_ces$par
            diff_ces <- sum(abs(par_to_upd - old_par))
          }
          RQoptim_ces[j, 1] <- vOptim_ces$value
          RQoptim_ces[j, 2:(ncol(init_ces) + 1)] <- par_to_upd
          
        }
      }
      
      params_caviares <-
        RQoptim_ces[order(RQoptim_ces[, 1]), ][1, 2:(ncol(init_ces) + 1)]
      VaR_caves <-
        caviar_SAV(params_caviares,
                   rt = tempdat_1day,
                   EMPQNT = emp_qnt,
                   RowsR = l_tempdat) * sqrt(days_ahead)
      
      ES_caves <-
        ES_output(
          setup = setup,
          BETAS = params_caviares,
          VaR = VaR_caves,
          data = tempdat_1day,
          EMPES = emp_es,
          EMPQNT = emp_qnt ,
          RowsR = l_tempdat,
          mu = con_mean_1day,
          pval = prob
        )
      output[, 1] <- c(VaR_caves[size:(len_r - days_ahead - 1)])
      output[, 2] <-  ES_caves[size:(len_r - days_ahead - 1)]
      
      ####### MIDAS #######
      rt_midas <-
        tempdat_1day[1:(l_tempdat - days_ahead)]
      rt_h_midas <-
        returns_h[1:(l_tempdat - days_ahead)]
      l_tempdat_midas <- length(rt_midas)
      
      if (i %% re_estimate == 0) {
        ## MIDAS parameters & returns
        init_midas <- sv_midas(setup, total_trys)
        lowb_midas <- lowerb_midas(setup)
        upb_midas <- upperb_midas(setup)
        
        RQfval_midas <-
          apply(
            init_midas,
            1,
            lklfct_MIDAS,
            rt = rt_midas,
            RowsR = l_tempdat_midas,
            pval = prob,
            d = d,
            EMPQNT = emp_qnt,
            EMPES = emp_es,
            penalty = penalty,
            es_approach = setup,
            mu = con_mean_h,
            rt_h = rt_h_midas
          )
        BestInitialCond_midas <-
          init_midas[order(RQfval_midas), ][1:nInitialCond, ]
        RQoptim_midas <- cbind(NA, BestInitialCond_midas, NA)
        
        # initial optimization procedure
        for (j in 1:nrow(BestInitialCond_midas)) {
          diff_midas <- 1
          par_to_upd <- BestInitialCond_midas[j,]
          rounds <- 0
          while (diff_midas > 0.1 & rounds < 10) {
            old_par <- par_to_upd
            rounds <- rounds + 1
            # Nelder Mead
            vOptim_midas <- optim(
              par = par_to_upd ,
              #BestInitialCond_midas[j,],
              lklfct_MIDAS,
              rt = rt_midas,
              RowsR = l_tempdat_midas,
              pval = prob,
              d = d,
              EMPQNT = emp_qnt,
              EMPES = emp_es,
              penalty = penalty,
              es_approach = setup,
              mu = con_mean_h,
              rt_h = rt_h_midas,
              method = 'Nelder-Mead',
              control = control
            )
            par_to_upd <- vOptim_midas$par
            
            # L-BFGS
            vOptim_midas <- optim(
              par = par_to_upd ,
              #BestInitialCond_midas[j,],
              lklfct_MIDAS,
              rt = rt_midas,
              RowsR = l_tempdat_midas,
              pval = prob,
              d = d,
              EMPQNT = emp_qnt,
              EMPES = emp_es,
              penalty = penalty,
              es_approach = setup,
              mu = con_mean_h,
              rt_h = rt_h_midas,
              method = met,
              lower = lowb_midas,
              upper = upb_midas,
              control = control
            )
            par_to_upd <- vOptim_midas$par
            diff_midas <- sum(abs(par_to_upd - old_par))
            
          }
          RQoptim_midas[j, 1] <- vOptim_midas$value
          RQoptim_midas[j, 2:(ncol(init_midas) + 1)] <-
            vOptim_midas$par
        }
      }
      
      params_midas <-
        RQoptim_midas[order(RQoptim_midas[, 1]), ][1, 2:(ncol(init_midas) + 1)]
      d_betas <- c()
      for (j in 1:d) {
        d_betas[j] <- dbeta(j / d, shape1 = 1, shape2 = params_midas[3])
      }
      d_betas <- d_betas / sum(d_betas)
      
      # Get VaR and ES one step further by using the parameters
      VaR_MIDAS <-
        VaR_midas_fct(params_midas, rt = rt_midas[1:(l_tempdat)], d, d_betas)
      ES_MIDAS <-
        ES_output(
          setup = setup,
          BETAS = params_midas,
          VaR = VaR_MIDAS,
          data = rt_h_midas,
          EMPES = emp_es,
          EMPQNT = emp_qnt ,
          RowsR = l_tempdat_midas,
          mu = con_mean_h,
          pval = prob
        )
      output[, 3] <- c(VaR_MIDAS[size:(l_tempdat_midas - 1)])
      output[, 4] <- ES_MIDAS[size:(l_tempdat_midas - 1)]
      
      outputdf <-
        data.frame(output)
      names(outputdf) <-
        c('VaRcaviares', 'EScaviares', 'VaRMIDAS', 'ESMIDAS')
    }
    
    return(outputdf)
  }
