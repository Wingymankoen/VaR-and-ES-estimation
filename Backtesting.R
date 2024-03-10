# First install rstudioapi to find the main.R file location
initial_package <- 'rstudioapi'
if (!requireNamespace(initial_package, quietly = TRUE)) {
  install.packages(initial_package, dependencies = TRUE)
  library(initial_package, character.only = TRUE)
} else {
  library(initial_package, character.only = TRUE)
}
main_wd <- (dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages & functions ----------------------------------------------------

library('Matrix')
library('readxl')
library('GAS')
library('rugarch')
library('esback')
library('esreg')
library('Formula')
library('SparseM')
library('quantreg')
library('MCS')
library('knitr')

# Load data ---------------------------------------------------------------

# Parameters
starts <- 504
days_ahead <- 1
prob <- 0.05
bank <- 'MS' # MS, HSBC, S&P500
category <- 'Unconditional'

# Load returns
fileinput <- paste(main_wd, '/Data', sep = '')
fileinput <- paste(fileinput, bank, sep = '/')
fileinput <- paste(fileinput, '.csv', sep = '')
data <- read.csv(fileinput, sep = ";")

returns_h <- c()
# Calculate h-day forward returns
for (i in starts:(dim(data)[1] - days_ahead)) {
  returns_h[i - starts + 1] <- sum(data[(i + 1):(i + days_ahead), ])
}
returns_h <- returns_h[1:(length(returns_h) - 1)]
#returns_h<- returns_h[(days_ahead+1):length(returns_h)] # FOR FC

# Load VaR & ES
direct <-
  paste(main_wd, '/Output', sep = '')
direct <-
  paste(direct,
        as.character(days_ahead),
        sep = '/')
direct <- paste(direct, bank, sep = '/')
direct <- paste(direct, category, sep = '/')

var_est <-
  read.csv(paste(direct, '_VaR.csv', sep = ''))
es_est <-
  read.csv(paste(direct, '_ES.csv', sep = ''))

if (sum(is.na(es_est)) > 0) {
  for (i in 1:dim(es_est)[2]) {
    for (j in 1:dim(es_est)[1]) {
      if (is.na(es_est[j, i])) {
        es_est[j, i] = es_est[j - 1, i]
      }
    }
  }
}
if (sum(is.na(var_est)) > 0) {
  for (i in 1:dim(var_est)[2]) {
    for (j in 1:dim(var_est)[1]) {
      if (is.na(var_est[j, i])) {
        var_est[j, i] = var_est[j - 1, i]
      }
    }
  }
}
# VaR & ES Functions ------------------------------------------------------

# VaR Testing
var_test <- function(ret, val_at_risk, pval) {
  no_models <- 1
  test_output <- matrix(NA, nrow = 6, ncol = no_models)
  hit_rate <- c()
  
  for (j in 1:no_models) {
    tests <- BacktestVaR(data = ret,
                         VaR = val_at_risk[1:length(ret), j],
                         alpha = pval)
    test_output[1:6, j] <-
      c(
        tests$LRuc[1],
        tests$LRuc[2],
        tests$LRcc[1],
        tests$LRcc[2],
        tests$DQ$stat,
        tests$DQ$pvalue
      )
    hit_rate[j] <-
      sum(ret < val_at_risk[1:length(ret), j]) / length(ret)
  }
  output <- rbind(test_output, hit_rate)
  colnames(output) <- names(val_at_risk)
  output <- round(output, 4)
  return(output)
}

# ES Testing
es_test <- function(ret, val_at_risk, exp_shf, pval) {
  no_models <- 1#dim(val_at_risk)[2]
  test_output <- matrix(NA, nrow = 5, ncol = no_models)
  hit_rate <- c()
  for (j in 1:no_models) {
    test <-  ESTest(
      alpha = pval,
      actual = ret,
      ES = exp_shf[1:length(ret), j],
      VaR = val_at_risk[1:length(ret), j],
    )
    test_output[1:3, j] <-
      c(test$expected.exceed, test$actual.exceed, test$p.value)
    
    test_output[4:5, j] <-
      c(
        esr_backtest(
          r = ret,
          q = val_at_risk[1:length(ret), j],
          e = exp_shf[1:length(ret), j],
          alpha = pval,
          version = 2
        )$pvalue_twosided_asymptotic,
        esr_backtest(
          r = ret,
          q = val_at_risk[1:length(ret), j],
          e = exp_shf[1:length(ret), j],
          alpha = pval,
          version = 1
        )$pvalue_twosided_asymptotic
      )
    
  }
  
  return(test_output)
}

# Model Confidence Set ----------------------------------------------------
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

# Gather all models
lloc_fc <- paste(main_wd, '/Output/', sep = '')
lloc_fc <- paste(lloc_fc, as.character(days_ahead), sep = '')
lloc_fc <- paste(lloc_fc, bank , sep = "/")

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

lloc_fc_var <- paste(lloc_fc, '/FC_VaR.csv' , sep = "")
lloc_fc_es <- paste(lloc_fc, '/FC_ES.csv' , sep = "")
VaR_models_fc <- read.csv(lloc_fc_var)
ES_models_fc <- read.csv(lloc_fc_es)

shortest_l <- dim(VaR_models_fc)[1]
extra <- 1 + days_ahead
VaR_models_all <-
  cbind(VaR_models_unc[extra:(shortest_l + extra - 1), ],
        VaR_models_con[extra:(shortest_l + extra - 1), ],
        VaR_models_q[extra:(shortest_l + extra - 1), ],
        VaR_models_fc)
ES_models_all <-
  cbind(ES_models_unc[extra:(shortest_l + extra - 1), ], ES_models_con[extra:(shortest_l +
                                                                                extra - 1), ], ES_models_q[extra:(shortest_l + extra - 1), ], ES_models_fc)

VaR_models_all <- test_na(VaR_models_all)
ES_models_all <- test_na(ES_models_all)
total_models_all <- dim(VaR_models_all)[2]


n_obs <-  dim(VaR_models_all)[1]
n_models <- total_models_all

returns_h <- c()
# Calculate h-day forward returns
for (i in starts:(length(data) - days_ahead)) {
  returns_h[i - starts + 1] <- sum(data[(i + 1):(i + days_ahead)])
}
returns_h <- returns_h[2:(length(returns_h))]
m_rth <- mean(returns_h)
returns_h <- returns_h - m_rth
VaR_models_all <- VaR_models_all - m_rth
ES_models_all <- ES_models_all - m_rth

loss_fct <- function(vlar, exs, rets, pval, days_a) {
  W <- 4
  if (days_a == 1) {
    outp <- matrix(NA, nrow = n_obs, ncol = n_models)
    n_models <- dim(vlar)[2]
    n_obs <- dim(vlar)[1]
    for (j in 1:n_models) {
      outp[, j] <-
        pval * (exs[, j] ^ 2 / 2 + W * vlar[, j] ^ 2 / 2 - vlar[, j] * exs[, j]) + (rets < vlar[, j]) *
        (-exs[, j] * (rets - vlar[, j]) + W * (rets ^ 2 - vlar[, j] ^ 2) / 2)
    }
  }
  if (days_a == 5) {
    outp <-
      pval * (exs ^ 2 / 2 + W * vlar ^ 2 / 2 - vlar * exs) + (rets < vlar) *
      (-exs * (rets - vlar) + W * (rets ^ 2 - vlar ^ 2) / 2)
    
  }
  
  # Use a different loss function
  #(rets < vlar[,j] - pval)*vlar[,j] - (rets < vlar[,j])*rets + (exs[,j])/(1+exp(exs[,j]))* (exs[,j]-vlar[,j]+(rets < vlar[,j])*(vlar[,j]-rets)/pval) + log(2/(1+exp(exs[,j])))
  #-1 / (exs[,j]) * (exs[,j] - vlar[,j] + (rets < vlar[,j]) * (vlar[,j] - rets) / pval) +   log(-exs[,j]) + (1 - log(1 - pval))
  # pval * (exs[,j]^2 /2+ 4*vlar[,j]^2/2 - vlar[,j]*exs[,j]) + (rets < vlar[,j])*(-exs[,j]*(rets-vlar[,j]) +4*(rets^2-vlar[,j]^2)/2  )
  
  # }
  
  return(outp)
}

if (days_ahead == 1) {
  loss_final <-
    loss_fct(
      vlar = VaR_models_all,
      exs = ES_models_all,
      rets = returns_h,
      pval = prob,
      days_a = days_ahead
    )
}
if (days_ahead == 5) {
  loss_final <-
    matrix(NA,
           nrow = length(returns_h[seq(1, length(returns_h), days_ahead)]) - 1,
           ncol = dim(VaR_models_all)[2])
  for (m in 1:dim(VaR_models_all)[2]) {
    loss <-
      matrix(NA, nrow = length(returns_h[seq(1, length(returns_h), days_ahead)]) -
               1, ncol = days_ahead)
    for (j in 1:days_ahead) {
      returns_h_temp = returns_h[seq(j, length(returns_h), days_ahead)]
      var_est_temp = VaR_models_all[seq(j, nrow(VaR_models_all), days_ahead), m]
      es_est_temp = ES_models_all[seq(j, nrow(ES_models_all), days_ahead), m]
      loss[, j] <-
        loss_fct(
          vlar = var_est_temp[1:dim(loss)[1]],
          exs = es_est_temp[1:dim(loss)[1]],
          rets = returns_h_temp[1:dim(loss)[1]],
          pval = prob,
          days_a = days_ahead
        )
    }
    loss <- data.frame(loss)
    loss_final[, m] <- apply(loss, 1, median, na.rm = T)
  }
  
}

loss_final <- na.omit(loss_final)
boots <- 1000
res1 <- MCSprocedure(loss_final, alpha = 0.2, B = boots)
