set.seed(1234)

# Initial -----------------------------------------------------------------

# Set initial parameters, general
# Each type of models have their own parameters.
# These can be found in their category files, e.g. Unconditional.R
bank <- 'MS' # Set which bank data is used
prob <-
  0.05 # Set the boundary value to look at. Usual values are 0.01 and 0.05
days_ahead <-
  1 # Set days to look ahead for. 1 and 5 have been looked into
start_count <- 252 # Minimum number of data points to use

# First install rstudioapi to find the main.R file location
initial_package <- 'rstudioapi'
if (!requireNamespace(initial_package, quietly = TRUE)) {
  install.packages(initial_package, dependencies = TRUE)
  library(initial_package, character.only = TRUE)
} else {
  library(initial_package, character.only = TRUE)
}
main_wd <- (dirname(rstudioapi::getActiveDocumentContext()$path))

# Source files/functions
r_files <-
  list.files(paste0(main_wd, '/src/'),
             pattern = "\\.R$",
             full.names = TRUE)
lapply(r_files, source)

# Install packages
read_packages(paste0(main_wd, '/requirements.txt'))

# Read data
df <- read_data(bank, main_wd)
no_of_data <- dim(df)[1]
setwd(main_wd)

# Modelling ---------------------------------------------------------------

# Unconditional models
unc_output <- unconditional(
  data = df,
  starts = start_count,
  prob = prob,
  days_ahead = days_ahead,
  len_r = no_of_data
)

# Conditional models
con_output <- conditional(
  rt = df,
  starts = start_count,
  prob = prob,
  days_ahead = days_ahead,
  len_r = no_of_data,
  wd = main_wd
)

# Quantile models
quan_output <-
  quantile_fct(
    orig_returns = df,
    wd = main_wd,
    len_r = no_of_data,
    prob = prob,
    days_ahead = days_ahead,
    starts = start_count
  )

# Output all models
export_data(
  data = unc_output,
  wd = main_wd,
  bank = bank,
  days_ahead = days_ahead,
  category = 'Unconditional'
)


# Export data -------------------------------------------------------------
export_data(
  data = con_output[[1]],
  wd = main_wd,
  bank = bank,
  days_ahead = days_ahead,
  category = 'Conditional'
)

export_data(
  data = con_output[[2]],
  wd = main_wd,
  bank = bank,
  days_ahead = days_ahead,
  category = 'Conditional',
  time_var = TRUE
)

export_data(
  data = quan_output,
  wd = main_wd,
  bank = bank,
  days_ahead = days_ahead,
  category = 'Quantile'
)

# FC
fc_output <-
  FC_fct(
    rt = df ,
    wd = main_wd,
    len_r = no_of_data,
    prob = prob,
    days_ahead = days_ahead
    ,
    starts = start_count
  )

export_data(
  data = cbind(fc_output[[1]], fc_output[[2]]),
  wd = main_wd,
  bank = bank,
  days_ahead = days_ahead,
  category = 'FC'
)
