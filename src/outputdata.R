export_data <-
  function(data,
           wd,
           bank,
           days_ahead,
           category,
           time_var = FALSE,
           weights = FALSE) {
    lloc <- paste0(wd, '/Output/')
    if (time_var) {
      lloc <- paste0(lloc, as.character(days_ahead))
      lloc <- paste0(lloc, '/')
      lloc <- paste0(lloc, 'Time-variation/')
      lloc <- paste0(lloc, bank)
      lloc <- paste0(lloc, '/Volatilities.csv')
      write.csv(data,
                lloc,
                row.names = FALSE)
    } else{
      lloc <- paste0(lloc, as.character(days_ahead))
      lloc <- paste0(lloc, '/')
      lloc <- paste0(lloc, bank)
      lloc <- paste0(lloc, '/')
      lloc <- paste0(lloc, category)
      lloc_var <- paste0(lloc, "_VaR.csv")
      lloc_es <- paste0(lloc, "_ES.csv")
      if (category == 'FC') {
        write.csv(data[, ncol(data) / 2],
                  lloc_var,
                  row.names = FALSE)
        write.csv(data[, (ncol(data) / 2) + 1],
                  lloc_es,
                  row.names = FALSE)
      } else{
        write.csv(data[, seq(1, ncol(data), by = 2)],
                  lloc_var,
                  row.names = FALSE)
        write.csv(data[, seq(2, ncol(data), by = 2)],
                  lloc_es,
                  row.names = FALSE)
      }
    }
    
  }