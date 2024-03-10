

# Functions ---------------------------------------------------------------

fit_tdist <- function(input) {
  out <- fitdistrplus::fitdist(input, "t.scaled",
                               start = list(
                                 df = 3,
                                 mean = mean(input),
                                 sd = sd(input)
                               ))
  return(out)
}

unconditional <- function(data, starts, prob, days_ahead, len_r) {
  # Parameters
  size <- 252
  evt_boundary <- 0.15
  
  # Initialize
  outp_list <- list()
  counter <- 0
  
  for (i in starts:(len_r)) {
    counter <- counter + 1
    print(paste0("Unc %: ", round(100 * counter / len_r, 3)))
    
    output <- c()
    tempdat <- data[(i - size + 1):i,]
    thres <- -quantile(tempdat, evt_boundary)
    
    #normal
    norm <-  fitdistrplus::fitdist(tempdat, "norm")
    output[1:2] <-
      c(
        cvar::VaR(
          qnorm,
          prob,
          mean = norm$estimate[1],
          sd = norm$estimate[2]
        ),
        cvar::ES(
          qnorm,
          prob,
          mean = norm$estimate[1],
          sd = norm$estimate[2]
        )
      )
    #t-dist
    try_t <- try(fit_tdist(tempdat))
    if ("try-error" %in% class(try_t)) {
      tdist <- tdist
    } else{
      tdist <- fit_tdist(tempdat)
    }
    output[3:4] <-
      c(
        cvar::VaR(
          qlst,
          prob,
          df = tdist$estimate[1],
          mu = tdist$estimate[2],
          sigma = tdist$estimate[3]
        ),
        cvar::ES(
          qlst,
          prob,
          df = tdist$estimate[1],
          m = tdist$estimate[2],
          s = tdist$estimate[3]
        )
      )
    
    # EVT
    try_gpd <- try(gpd(-tempdat, threshold = thres))
    if ("try-error" %in% class(try_gpd)) {
      gpddist <- gpddist
    } else{
      gpddist <- gpd(-tempdat, threshold = thres)
    }
    vargpd <-
      thres +  (gpddist$par.ests[2] / gpddist$par.ests[1]) * ((length(tempdat) /
                                                                 gpddist$n.exceed * (prob)) ^ (-gpddist$par.ests[1]) - 1)
    esgpd <-
      vargpd / (1 - gpddist$par.ests[1])  + (gpddist$par.ests[2] - gpddist$par.ests[1] *
                                               thres) / (1 - gpddist$par.ests[1])
    output[5:6] <- c(vargpd, esgpd)
    
    #nonpar
    nonpar <- density(as.numeric(tempdat), n = 2 ^ 15)
    g <- function(x, z, bw, p)
      sum(pnorm(x - z, sd = bw)) / length(z) - p
    varnp <-
      uniroot(
        g,
        range(nonpar$x) + c(-1, 1),
        z = tempdat,
        bw = nonpar$bw,
        p = 1 - prob
      )$root
    ESnp <-
      sum(nonpar$x[nonpar$x < -varnp] * nonpar$y[nonpar$x < -varnp]) / sum(nonpar$y[nonpar$x < -varnp])
    output[7:8] <- c(varnp,-ESnp)
    
    # Output
    outp_list[[i - size + 1]] <- output
  }
  outputdf <- do.call(rbind.data.frame, outp_list)
  outputdf <- -outputdf * sqrt(days_ahead)
  names(outputdf) <-
    c(
      'Normal',
      'Normal',
      "Student's t",
      "Student's t",
      'EVT',
      'EVT',
      'Non-parametric',
      'Non-parametric'
    )
  
  return(outputdf)
}