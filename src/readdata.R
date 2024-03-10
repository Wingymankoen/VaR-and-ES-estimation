read_data <- function(bank, wd) {
  fileinput <- paste0(paste0(wd, '/Data/'), bank)
  fileinput <- paste0(fileinput, '.csv')
  data <- read.csv(fileinput, sep = ";")
  return(data)
}