

# Functions ---------------------------------------------------------------

fc_opt <- function(param, ret, var_f, es_f, pval, no_m) {
  weights_fc <- param
  weights_q <- c(weights_fc[1:(no_m)])
  weights_es <- c(weights_fc[(no_m + 1):(length(weights_fc))])
  
  var_ct <- as.matrix(var_f) %*% as.vector(weights_q)
  es_ct <-
    var_ct + (as.matrix(es_f) - as.matrix(var_f)) %*% as.vector(weights_es)
  
  
  score <-
    -1 / (es_ct) * (es_ct - var_ct + (ret < var_ct) * (var_ct - ret) / pval) +   log(-es_ct) + (1 - log(1 - pval))
  return(sum(score))
}

equal_fct <- function(param, ret, var_f, es_f, pval, no_m) {
  weights_fc <- param
  weights_q <- c(weights_fc[1:(no_m)])
  weights_es <- c(weights_fc[(no_m + 1):(length(weights_fc))])
  
  return(c(sum(weights_q), sum(weights_es)))
}

full_procedure <- function(var_f, es_f, returns, pval) {
  no_models <- dim(var_f)[2]
  weights_q <- rep(1 / no_models, no_models)
  weights_es <- rep(1 / no_models, no_models)
  params <- c(weights_q, weights_es)
  lowb <- rep(0, length(params))
  upb <- rep(1, length(params))
  
  opt <- solnp(
    params,
    fc_opt,
    ret = returns,
    var_f = var_f,
    es_f = es_f,
    no_m = no_models,
    pval = pval,
    eqfun = equal_fct,
    eqB = c(1, 1),
    LB = lowb,
    UB = upb,
    control = list(trace = 0)
  )
  
  return(opt$par)
}

calc_var <- function(weights, qfc) {
  qct <-
    as.matrix(qfc, ncol = dim(qfc)[2]) * as.matrix(weights, ncol = dim(weights)[2])
  qct <- rowSums(qct)
  return(qct)
}

calc_es <- function(weights, qfc, esfc, qct) {
  qfc <- as.matrix(qfc, ncol = dim(qfc)[2])
  esfc <- as.matrix(esfc, ncol = dim(esfc)[2])
  esct <-
    qct + rowSums((esfc - qfc) * as.matrix(weights, ncol = dim(weights)[2]))
  return(esct)
}

test_na <- function(x) {
  if (sum(is.na(x)) > 0) {
    for (i in 1:dim(x)[2]) {
      for (j in 1:dim(x)[1]) {
        if (is.na(x[j, i])) {
          x[j, i] = x[j - 1, i]
        }
      }
    }
  }
  return(x)
}

# Model -------------------------------------------------------------------


FC_fct <- function(rt, wd, len_r, prob, days_ahead, starts) {
  returns_h <- c()
  # Calculate h-day forward returns
  for (i in starts:(len_r - days_ahead)) {
    returns_h[i - starts + 1] <- sum(rt[(i + 1):(i + days_ahead), ])
  }
  mean_rt <- mean(returns_h)
  
  # Unconditional -----------------------------------------------------------
  
  lloc_fc <- paste(wd, '/Output/' , sep = '')
  lloc_fc <- paste(lloc_fc, as.character(days_ahead), sep = '')
  lloc_fc <- paste(lloc_fc, '/' , sep = "")
  lloc_fc <- paste(lloc_fc, bank , sep = "")
  lloc_fc_var <- paste(lloc_fc, '/Unconditional_VaR.csv' , sep = "")
  lloc_fc_es <- paste(lloc_fc, '/Unconditional_ES.csv' , sep = "")
  VaR_models <- read.csv(lloc_fc_var)
  ES_models <- read.csv(lloc_fc_es)
  VaR_models <- test_na(VaR_models)
  ES_models <- test_na(ES_models)
  total_models <- dim(VaR_models)[2]
  
  
  output <- list()
  counter <- 0
  loop_t <- min(dim(VaR_models)[1], length(returns_h))
  for (t in 1:loop_t) {
    counter <- counter + 1
    print(paste0("FC Unc %: ", round(100 * counter / len_r, 3)))
    temp_returns <- returns_h[1:t] - mean_rt
    temp_fc_var <- VaR_models[1:t, ] - mean_rt
    temp_fc_es <- ES_models[1:t, ] - mean_rt
    output[[t]] <-
      full_procedure(var_f = temp_fc_var, es_f = temp_fc_es, temp_returns, prob)
  }
  
  weights_unc <- do.call(rbind.data.frame, output)
  VaR_unc <-
    calc_var(weights_unc[, 1:total_models] , VaR_models[(1 + days_ahead):dim(VaR_models)[1], ])
  ES_unc <-
    calc_es(weights_unc[, (total_models + 1):dim(weights_unc)[2]] , VaR_models[(1 +
                                                                                  days_ahead):dim(VaR_models)[1], ], ES_models[(1 + days_ahead):dim(VaR_models)[1], ], qct = VaR_unc)
  
  # Conditional -------------------------------------------------------------
  lloc_fc_var <- paste(lloc_fc, '/Conditional_VaR.csv' , sep = "")
  lloc_fc_es <- paste(lloc_fc, '/Conditional_ES.csv' , sep = "")
  VaR_models <- read.csv(lloc_fc_var)
  ES_models <- read.csv(lloc_fc_es)
  VaR_models <- test_na(VaR_models)
  ES_models <- test_na(ES_models)
  total_models <- dim(VaR_models)[2]
  
  output <- list()
  counter <- 0
  loop_t <- length(returns_h)
  for (t in 1:loop_t) {
    counter <- counter + 1
    print(paste0("FC Con %: ", round(100 * counter / len_r, 3)))
    temp_returns <-  returns_h[1:t] - mean_rt
    temp_fc_var <- VaR_models[1:t, ] - mean_rt
    temp_fc_es <- ES_models[1:t, ] - mean_rt
    output[[t]] <-
      full_procedure(var_f = temp_fc_var, es_f = temp_fc_es, temp_returns, prob)
  }
  
  weights_con <- do.call(rbind.data.frame, output)
  VaR_con <-
    calc_var(weights_con[, 1:total_models] , VaR_models[(1 + days_ahead):dim(VaR_models)[1], ])
  ES_con <-
    calc_es(weights_con[, (total_models + 1):dim(weights_con)[2]] , VaR_models[(1 +
                                                                                  days_ahead):dim(VaR_models)[1], ], ES_models[(1 + days_ahead):dim(VaR_models)[1], ], qct = VaR_con)
  
  # Quantile ----------------------------------------------------------------
  lloc_fc_var <- paste(lloc_fc, '/Quantile_VaR.csv' , sep = "")
  lloc_fc_es <- paste(lloc_fc, '/Quantile_ES.csv' , sep = "")
  VaR_models <- read.csv(lloc_fc_var)
  ES_models <- read.csv(lloc_fc_es)
  VaR_models <- test_na(VaR_models)
  ES_models <- test_na(ES_models)
  total_models <- dim(VaR_models)[2]
  
  output <- list()
  counter <- 0
  loop_t <- min(dim(VaR_models)[1], length(returns_h))
  for (t in 1:loop_t) {
    counter <- counter + 1
    print(paste0("FC Quan %: ", round(100 * counter / len_r, 3)))
    temp_returns <-  returns_h[1:t] - mean_rt
    temp_fc_var <- VaR_models[1:t, ] - mean_rt
    temp_fc_es <- ES_models[1:t, ] - mean_rt
    output[[t]] <-
      full_procedure(var_f = temp_fc_var, es_f = temp_fc_es, temp_returns, prob)
  }
  
  weights_q <- do.call(rbind.data.frame, output)
  VaR_q <-
    calc_var(weights_q[-(1:days_ahead), 1:total_models] , VaR_models[(1 +
                                                                        days_ahead):dim(VaR_models)[1], ])
  ES_q <-
    calc_es(weights_q[-(1:days_ahead), (total_models + 1):dim(weights_q)[2]] , VaR_models[(1 +
                                                                                             days_ahead):dim(VaR_models)[1], ], ES_models[(1 + days_ahead):dim(VaR_models)[1], ], qct = VaR_q)
  
  # All ---------------------------------------------------------------------
  
  lloc_fc <- paste(wd, '/Output/' , sep = '')
  lloc_fc <- paste(lloc_fc, as.character(days_ahead), sep = '')
  lloc_fc <- paste(lloc_fc, '/' , sep = "")
  lloc_fc <- paste(lloc_fc, bank , sep = "")
  
  lloc_fc_var <- paste(lloc_fc, '/Unconditional_VaR.csv' , sep = "")
  lloc_fc_es <- paste(lloc_fc, '/Unconditional_ES.csv' , sep = "")
  VaR_models_unc <- read.csv(lloc_fc_var)
  ES_models_unc <- read.csv(lloc_fc_es)
  
  lloc_fc_var <- paste(lloc_fc, '/Conditional_VaR.csv' , sep = "")
  lloc_fc_es <- paste(lloc_fc, '/Conditional_ES.csv' , sep = "")
  VaR_models_con <- read.csv(lloc_fc_var)
  ES_models_con <- read.csv(lloc_fc_es)
  
  lloc_fc_var <- paste(lloc_fc, '/Quantile_VaR.csv' , sep = "")
  lloc_fc_es <- paste(lloc_fc, '/Quantile_ES.csv' , sep = "")
  VaR_models_q <- read.csv(lloc_fc_var)
  ES_models_q <- read.csv(lloc_fc_es)
  
  
  shortest_l <- dim(VaR_models_q)[1]
  VaR_models_all <-
    cbind(VaR_models_unc[1:shortest_l, ], VaR_models_con[1:shortest_l, ], VaR_models_q)
  ES_models_all <-
    cbind(ES_models_unc[1:shortest_l, ], ES_models_con[1:shortest_l, ], ES_models_q)
  
  VaR_models_all <- test_na(VaR_models_all)
  ES_models_all <- test_na(ES_models_all)
  total_models_all <- dim(VaR_models_all)[2]
  
  output <- list()
  counter <- 0
  loop_t <- min(dim(VaR_models)[1], length(returns_h))
  
  for (t in 1:loop_t) {
    counter <- counter + 1
    print(paste0("FC All %: ", round(100 * counter / loop_t, 3)))
    temp_returns <-  returns_h[1:t] - mean_rt
    temp_fc_var <- VaR_models_all[1:t, ] - mean_rt
    temp_fc_es <- ES_models_all[1:t, ] - mean_rt
    output[[t]] <-
      full_procedure(var_f = temp_fc_var, es_f = temp_fc_es, temp_returns, prob)
  }
  
  weights_all <- do.call(rbind.data.frame, output)
  VaR_all <-
    calc_var(weights_all[-(1:days_ahead), 1:total_models_all] , VaR_models_all[(1 +
                                                                                  days_ahead):dim(VaR_models_all)[1], ])
  ES_all <-
    calc_es(weights_all[-(1:days_ahead), (total_models_all + 1):dim(weights_all)[2]] , VaR_models_all[(1 +
                                                                                                         days_ahead):dim(VaR_models_all)[1], ], ES_models_all[(1 + days_ahead):dim(VaR_models_all)[1], ], qct = VaR_all)
  
  
  # Output ------------------------------------------------------------------
  shortest_l <- length(VaR_all)
  outputdf_var <-
    cbind(VaR_unc[1:shortest_l], VaR_con[1:shortest_l], VaR_q[1:shortest_l],  VaR_all)
  outputdf_es <-
    cbind(ES_unc[1:shortest_l], ES_con[1:shortest_l], ES_q[1:shortest_l], ES_all)
  names(outputdf_var) <-
    c('Unconditional', 'Conditional', 'Quantile', 'All')
  names(outputdf_es) <-
    c('Unconditional', 'Conditional', 'Quantile', 'All')
  return(list(outputdf_var, outputdf_es))
}
