
# CES ---------------------------------------------------------------------

sv_ces<- function(setup, total_trys ){
  if (setup == 1) {
    init_CES <- matrix(runif(total_trys * 4), ncol = 4)
    init_CES[, 4] <-
      runif(total_trys, min = -6, max = 1)
    
  } else if (setup == 2) {
    init_CES <- matrix(runif(total_trys * 6), ncol = 6)
    
  } else if (setup==3) {
    init_CES <- matrix(runif(total_trys * 7), ncol = 7)
  }
  
  init_CES[, 1] <-
    runif(total_trys,  min = -100, max = -1)
  init_CES[, 3] <-
    runif(total_trys,  min = -1, max = 1)
  
  return(init_CES)
}

lowerb_ces<- function(setup){
  if (setup == 1) {
  lowb <- c(rep(-1,3), rep(-6, 1))
  } else if(setup==2){
    lowb <- c(rep(-1,3), rep(0.0001, 3)) 
  } else if (setup==3) {
    lowb <- c(rep(-1,3), rep(0.0001, 4)) 
  }
  lowb[1]<- -Inf
  return(lowb)
}

upperb_ces<- function(setup){
  if (setup == 1) {
    upb <- c(rep(1, 3), Inf)
  } else if(setup==2){
    upb <- c(rep(1, 6))
    upb[4]<- Inf
  } else if (setup==3) {
    upb <- c(rep(1, 7))
    upb[4]<- Inf
  }
  return(upb)
}

# MIDAS -------------------------------------------------------------------

sv_midas <- function(setup, total_trys){
  if (setup == 1) {
    init_m <- matrix(runif(total_trys * 4), ncol = 4)
    init_m[, 4] <-
      runif(total_trys, min = -6, max = -1)
    
  } else if (setup == 2) {
    init_m <- matrix(runif(total_trys * 6), ncol = 6)
    
  } else if (setup==3) {
    init_m <- matrix(runif(total_trys * 7), ncol = 7)
    
  }
  init_m[, 1] <-
    runif(total_trys, min = -100, max = -1)
  init_m[, 2] <-
    runif(total_trys, min = -10, max = 1)
  init_m[, 3] <-
    runif(total_trys, min = 1, max = 10)
  
  return(init_m)
}

lowerb_midas<- function(setup){
  if (setup == 1) {
    lowb <- c(rep(-1, 3),-6)
  } else if(setup==2){
    lowb <- c(rep(-1, 3), rep(0.0001, 3)) 
  } else if (setup==3) {
    lowb <- c(rep(-1, 3), rep(0.0001, 4)) 
  }
  lowb[3]<- 1
  lowb[2]<- -Inf
  lowb[1]<- -Inf
  return(lowb)
}

upperb_midas<- function(setup){
  if (setup == 1) {
    upb <- c(rep(1, 4))
    upb[4]<- Inf
  } else if(setup==2){
    upb <- c(rep(1, 6))
    upb[4]<- Inf
    upb[5]<- Inf
  } else if (setup==3) {
    upb <- c(rep(1, 7))
    upb[4]<- Inf
    upb[5]<- Inf
  }
  upb[1]<- Inf
  upb[3]<- Inf
  upb[2] <- Inf
  return(upb)
}