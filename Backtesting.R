# First install rstudioapi to find the main.R file location
# List of packages to install
packages_to_install <- c(
  "rstudioapi",
  "Matrix",
  "readxl",
  'GAS',
  'rugarch',
  'esback',
  'esreg',
  'dplyr',
  'Formula',
  'SparseM' ,
  'quantreg',
  'MCS' ,
  'knitr'
)

# Install and load packages without messages
invisible(lapply(packages_to_install, function(package) {
  {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package, dependencies = TRUE)
    }
    library(package, character.only = TRUE)
  }
}))

# Get wd
main_wd <- (dirname(rstudioapi::getActiveDocumentContext()$path))

# Functions ---------------------------------------------------------------

test_na<- function(x){
  if(sum(is.na(x)) > 0 ){
    for(i in 1: dim(x)[2]){
      for(j in 1:dim(x)[1]){
        if( is.na(x[j,i] ) ){
          x[j,i]= x[j-1,i]
        }
      }
    } 
  }
  return(x)
}

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
rets <- read.csv(fileinput)

# Model Confidence Set ----------------------------------------------------

lloc_fc <- paste(main_wd, '/Output/', sep = '')
lloc_fc <- paste(lloc_fc, as.character(days_ahead), sep = '')
lloc_fc <- paste(lloc_fc, bank , sep = "/")

# Get a list of all CSV files in the folder
csv_files <-
  list.files(path = lloc_fc,
             pattern = "\\.csv$",
             full.names = TRUE)

# Read all CSV files into a list
invisible(lapply(csv_files, function(file) {
  # Extract the file name without extension
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the CSV file
  data <- read.csv(file)
  
  # Assign the data frame to an object with the file name
  assign(file_name, data, envir = .GlobalEnv)
}))

# Find the shortest length of the FC, as it combines all other datasets
shortest_l <- dim(FC_VaR)[1]
extra <- 1 + days_ahead
VaR_models_all <-
  cbind(Unconditional_VaR[extra:(shortest_l + extra - 1),],
        Conditional_VaR[extra:(shortest_l + extra - 1),],
        Quantile_VaR[extra:(shortest_l + extra - 1),],
        FC_VaR)

ES_models_all <-
  cbind(Unconditional_ES[extra:(shortest_l + extra - 1),], Conditional_ES[extra:(shortest_l +
                                                                               extra - 1),], Quantile_ES[extra:(shortest_l + extra - 1),], FC_ES)
VaR_models_all <- test_na(VaR_models_all)
ES_models_all <- test_na(ES_models_all)
total_models_all <- dim(VaR_models_all)[2]
n_obs <-  dim(VaR_models_all)[1]
n_models <- total_models_all

# Calculate h-day forward returns
returns_h <- c()
for (i in starts:(dim(rets)[1] - days_ahead)) {
  returns_h[i - starts + 1] <- sum(rets[(i + 1):(i + days_ahead),])
}
returns_h <- returns_h[2:(length(returns_h))]

# Lower the returns, VaR and ES by the mean of the returns
m_rth <- mean(returns_h)
returns_h <- returns_h - m_rth
VaR_models_all <- VaR_models_all - m_rth
ES_models_all <- ES_models_all - m_rth

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
# For the 5 days ahead, we need to find the median of the MCS for improved
# Stability of the results
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
