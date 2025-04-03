solve_G05 <- function(cobraSeed) {
  ## creating an instance for G05 
  G05<-COP$new("G05")
  
  ## Initializing SACOBRA  
  cobra <- cobraInit(xStart=G05$xStart, fn=G05$fn, fName=G05$name, 
                     lower=G05$lower, upper=G05$upper, feval=70,
                     solu=as.vector(G05$solu),   
                     initDesign = "LHS", cobraSeed=cobraSeed)
  cobra$squares = T
  # browser()

  ## Run SACOBRA optimizer 
  cobra <- cobraPhaseII(cobra)
  
  ## The true solution is at G05$solu = c(679.9453175, 1026.0671351, 0.1188764, -0.3962336)
  ## with objective 5126.498 
  ## The solution found by SACOBRA:
  print(getXbest(cobra))          # 679.9380064 1026.0749485 0.1188816 -0.3962311
  print(getFbest(cobra))          # 5126.498 
  cobra$errGCOP = abs(cobra$df$Best-G05$fn(G05$solu)[1])
  plot(cobra$errGCOP,log="y",type="l",
       ylab="error",xlab="iteration",main=G05$name)
  return (cobra)
}

solve_G06 <- function(cobraSeed) {
  ## creating an instance for G06 
  G06<-COP$new("G06")
  
  ## Initializing SACOBRA  
  cobra <- cobraInit(xStart=G06$xStart, fn=G06$fn, fName=G06$name, 
                     lower=G06$lower, upper=G06$upper, feval=40,
                     solu=as.vector(G06$solu),   
                     initDesign = "LHS", cobraSeed=cobraSeed)
  cobra$squares = F
  # browser()
  cobra$sac$RS = T    # temporarily, to check equivalence to Python side
  cobra$sac$RS_rep = T
  
  ## Run SACOBRA optimizer 
  cobra <- cobraPhaseII(cobra)
  
  ## The true solution is at G06$solu = c(14.0950, 0.84296)
  ## with objective -6961.814 
  ## The solution found by SACOBRA:
  print(getXbest(cobra))          # 14.0950000  0.8429608
  print(getFbest(cobra))          # -6961.814 
  cobra$errGCOP = abs(cobra$df$Best-G06$fn(G06$solu)[1])
  plot(cobra$errGCOP,log="y",type="l",
       ylab="error",xlab="iteration",main=G06$name)
  return (cobra)
}

multi_G06 <- function(runs, cobraSeed) {
  # consecutive runs to collect some statistics:
  ptm <- proc.time()
  fin_err_list = c()
  for (run in 1:runs) {
    cobra = solve_G06(cobraSeed+run)
    # print(cobra$df)
    fin_err = tail(cobra$errGCOP,1)
    fin_err_list = c(fin_err_list, fin_err)
    print(fin_err)
  }
  print(sprintf("%10.5f %s",(proc.time()-ptm)[3]/runs,"sec per run"))
  print(c("min", min(fin_err_list), "max", max(fin_err_list)))
  print(c("median", median(fin_err_list)))
  return (cobra)
}

solve_G24 <- function() {
  ##creating an instance for G24 problem 
  G24<-COP$new("G24")
  
  ##initializing SACOBRA 
  cobra <- cobraInit(xStart=G24$lower, fName=G24$name, 
                     solu=as.vector(G24$solu),   
                     fn=G24$fn, lower=G24$lower, upper=G24$upper, feval=25)

  ## Run SACOBRA optimizer 
  cobra <- cobraPhaseII(cobra)
  
  ## The true solution is at G24$solu = c(2.329520, 3.178493) 
  ## with objective -5.508013
  ## The solution found by SACOBRA:
  print(getXbest(cobra))          # 2.329725 3.176549
  print(getFbest(cobra))          # -5.506275
  cobra$errGCOP = abs(cobra$df$Best-G24$fn(G24$solu)[1])
  plot(cobra$errGCOP,log="y",type="l",
       ylab="error",xlab="iteration",main=G24$name)
  return (cobra)
}

cobra = solve_G05(42)
# cobra = solve_G06(42)
# cobra = multi_G06(10, 390)
# cobra = solve_G24()

