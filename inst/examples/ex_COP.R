solve_G24 <- function() {
  ##creating an instance for G24 problem 
  G24<-COP$new("G24")
  
  ##initializing SACOBRA 
  cobra <- cobraInit(xStart=G24$lower, fName=G24$name, 
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

solve_G06 <- function(cobraSeed) {
  ## creating an instance for G06 
  G06<-COP$new("G06")
  
  ## Initializing SACOBRA  
  cobra <- cobraInit(xStart=G06$xStart, fn=G06$fn, fName=G06$name, 
                     lower=G06$lower, upper=G06$upper, feval=40,
                     initDesign = "RAND_REP", cobraSeed=cobraSeed)
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

# cobra = solve_G24()
cobra = solve_G06(42)
# cobra = multi_G06(10, 390)

