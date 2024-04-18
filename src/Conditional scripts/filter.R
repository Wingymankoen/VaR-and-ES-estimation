# Filter ------------------------------------------------------------------
filterfct <-
  function(prob,
           no_sim,
           i = i,
           look_ahead,
           std_res,
           variables,
           est_sig,
           returns,
           model_type,
           variables_ar,
           evt = FALSE,
           u = NULL,
           variables_gpd = NULL) {
    # Sample the standardized residuals for the simulation
    sampled_values_list <-
      lapply(1:no_sim, function(x)
        sample(
          std_res,
          size = look_ahead,
          replace = TRUE,
          prob = NULL
        ))
    sampled_residuals <-
      matrix(unlist(sampled_values_list),
             ncol = look_ahead,
             byrow = TRUE)
    # Set starting values for the simulation
    init_sig <- rep(est_sig[i] ^ 2, no_sim)
    init_res <- rep(returns[i], no_sim)
    init_return <- rep(returns[i], no_sim)
    
    # The factors which will be updated each iteration
    current_factors <-
      matrix(c(init_sig, init_res, init_return), ncol = 3)
    # Create an empty matrix to fill with daily returns
    daily_returns <- matrix(NA, ncol = look_ahead, nrow = no_sim)
    # For evt approach, resample the residuals to fit the evt distribution
    if (evt) {
      sampled_residuals[sampled_residuals < -u] <-
        -rgpd(sum(sampled_residuals < -u),
              xi = variables_gpd[1],
              mu = u,
              beta = variables_gpd[2])
    }
    
    # Simulate for garch methods
    if (model_type == 'garch') {
      for (j in 1:look_ahead) {
        current_factors[, 1] <-
          variables[2] * current_factors[, 2] ^ 2 + variables[3] * current_factors[, 1]
        current_factors[, 2] <-
          sampled_residuals[, j] * sqrt(current_factors[, 1])
        current_factors[, 3] <-
          variables_ar[2] +   current_factors[, 3] * variables_ar[1] + current_factors[, 2]
        daily_returns[, j] <- current_factors[, 3]
      }
      
      # Simulate for gjr-garch methods
    } else if (model_type == 'gjr-garch') {
      for (j in 1:look_ahead) {
        current_factors[, 1] <-
          variables[2] * current_factors[, 2] ^ 2 + variables[3] * current_factors[, 1] + variables[4] * (current_factors[, 2] <
                                                                                                            0) ^ 2
        current_factors[, 2] <-
          sampled_residuals[, j] * sqrt(current_factors[, 1])
        current_factors[, 3] <-
          variables_ar[2] +   current_factors[, 3] * variables_ar[1]  + current_factors[, 2]
        daily_returns[, j] <- current_factors[, 3]
      }
    }
    
    # Take row sums, to get the h-day ahead returns
    results <- rowSums(daily_returns)
    # Estimate VaR and ES
    VaRgarch <- quantile(results, prob, na.rm = TRUE)
    ESgarch <- mean(results[results < VaRgarch])
    return(cbind(VaRgarch, ESgarch))
  }
