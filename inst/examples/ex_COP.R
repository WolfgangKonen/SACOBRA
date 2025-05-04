solve_G03 <- function(cobraSeed, dimension=20) {
  G03<-COP$new("G03", dimension)
  
  cobra <- cobraInit(xStart=G03$x0, fn=G03$fn, fName=G03$name, 
                     lower=G03$lower, upper=G03$upper, feval=250,
                     solu=as.vector(G03$solu),   
                     initDesign = "LHS", cobraSeed=cobraSeed, conTol=0)
  cobra$squares = T

  cobra <- cobraPhaseII(cobra)
  
  ## The true solution is at G03$solu = vector (z,z,...) along main diagonal with z=1/sqrt(d)
  ## with objective -1.0000 
  ## The solution found by SACOBRA: similar
  ## with objective -1.0000      (error 1e-3 or smaller)
  cobra <- print_and_plot(cobra,c(1e-8,1e-1))
  return (cobra)
}

solve_G05 <- function(cobraSeed) {
  G05<-COP$new("G05")
  
  cobra <- cobraInit(xStart=G05$x0, fn=G05$fn, fName=G05$name, 
                     lower=G05$lower, upper=G05$upper, feval=70,
                     finalEpsXiZero = TRUE,
                     solu=as.vector(G05$solu),   
                     initDesign = "RAND_REP", initDesPoints = 15,  cobraSeed=cobraSeed, conTol=0)
  cobra$squares = T
  browser()
  
  cobra <- cobraPhaseII(cobra)
  
  ## The true solution is at G05$solu = c(679.9453175, 1026.0671351, 0.1188764, -0.3962336)
  ## with objective 5126.498 
  ## The solution found by SACOBRA: 679.9380064 1026.0749485 0.1188816 -0.3962311
  ## with objective 5126.498      (error 5e-6 or smaller)
  cobra <- print_and_plot(cobra,c(1e-8,1e-1))
  return (cobra)
}

solve_G06_R <- function(cobraSeed) {
  G06<-COP$new("G06")
  idp = 5
  
  cobra <- cobraInit(xStart=G06$x0, fn=G06$fn, fName=G06$name, squares=F,
                     lower=G06$lower, upper=G06$upper, feval=40,
                     #epsilonInit = 0.0, epsilonMax = 0.0,
                     finalEpsXiZero = TRUE,
                     conTol = 1e-7,
                     solu=as.vector(G06$solu),   
                     initDesign = "RAND_REP", initDesPoints=idp, cobraSeed=cobraSeed)
  cobra$squares = T
  # cobra$sac$RS = T    # temporarily, to check equivalence to Python side
  # cobra$sac$RS_rep = T
  
  cobra <- cobraPhaseII(cobra)
  
  ## The true solution is at G06$solu = c(14.0950, 0.84296)
  ## with objective -6961.814 
  ## The solution found by SACOBRA: 14.0950000  0.8429608
  ## with objective -6961.814     (error 2e-4 or smaller)
  cobra <- print_and_plot(cobra,c(1e-7,1e-0))
  return (cobra)
}

solve_G06 <- function(cobraSeed) {
  G06<-COP$new("G06")
  
  cobra <- cobraInit(xStart=G06$x0, fn=G06$fn, fName=G06$name, 
                     lower=G06$lower, upper=G06$upper, feval=80,
                     #epsilonInit = 0.0, epsilonMax = 0.0,
                     finalEpsXiZero = FALSE,
                     conTol = 1e-7,
                     solu=as.vector(G06$solu),   
                     initDesign = "LHS", cobraSeed=cobraSeed)
  cobra$squares = F
  cobra$sac$RS = T    # temporarily, to check equivalence to Python side
  cobra$sac$RS_rep = T
  
  cobra <- cobraPhaseII(cobra)
  
  ## The true solution is at G06$solu = c(14.0950, 0.84296)
  ## with objective -6961.814 
  ## The solution found by SACOBRA: 14.0950000  0.8429608
  ## with objective -6961.814     (error 2e-4 or smaller)
  cobra <- print_and_plot(cobra,c(1e-7,1e-0))
  return (cobra)
}

solve_G11 <- function(cobraSeed) {
  G11<-COP$new("G11")
  
  cobra <- cobraInit(xStart=G11$x0, fn=G11$fn, fName=G11$name, 
                     lower=G11$lower, upper=G11$upper, feval=70,
                     solu=as.vector(G11$solu),   
                     initDesign = "LHS", cobraSeed=cobraSeed, conTol=0)
  cobra$squares = T
  cobra$equHandle$equEpsFinal = 1e-10

  cobra <- cobraPhaseII(cobra)
  
  ## The true solution is at G11$solu = c(0.7071068, 0.5) = c(1/sqrt(2), 1/2)
  ## with objective 0.75 
  ## The solution found by SACOBRA: 0.7071070 0.5000003
  ## with objective 0.75      (error 1e-10 or smaller)
  cobra <- print_and_plot(cobra,c(1e-11,1e-2))
  return (cobra)
}

solve_G13 <- function(cobraSeed) {
  G13<-COP$new("G13")
  
  cobra <- cobraInit(#xStart=G13$solu[1,]+0.01, fn=G13$fn, fName=G13$name, 
                     xStart=G13$x0, fn=G13$fn, fName=G13$name, 
                     lower=G13$lower, upper=G13$upper, feval=500,
                     epsilonInit = 0.0, epsilonMax = 0.0,
                     RBFrho = 2.5, RBFrhoDec = 2.0,
                     solu=as.matrix(G13$solu),  conTol=1e-6,  # conTol=0, #
                     initDesign = "LHS", initDesPoints = 100, cobraSeed=cobraSeed)
  cobra$squares = T
  cobra$equHandle$equEpsFinal = 1e-7 # 2e-5 # 1e-7 # 
  cobra$muGrow = 100
  
  cobra <- cobraPhaseII(cobra)

  ## The true solution is at G13$solu = c(-1.7171435947203,1.5957097321519,1.8272456947885,-0.7636422812896,-0.7636439027742) and permutations
  ## with objective 0.05394984 
  ## The solution found by SACOBRA: -1.7046923  1.5812652  1.8500992 -0.7768188 -0.7531995
  ## with objective 0.05404521      (error < 1e-3, but only with p=50%. The rest gets trapped in local minimum with error = 0.385)
  cobra <- print_and_plot(cobra,c(1e-10,1e-0))
  return (cobra)
}

solve_G24 <- function(cobraSeed) {
  G24<-COP$new("G24")
  
  cobra <- cobraInit(xStart=G24$lower, fName=G24$name, 
                     fn=G24$fn, lower=G24$lower, upper=G24$upper, 
                     solu=as.vector(G24$solu),   
                     cobraSeed=cobraSeed, feval=45)

  cobra <- cobraPhaseII(cobra)
  
  ## The true solution is at G24$solu = c(2.329520, 3.178493) 
  ## with objective -5.50801327
  ## The solution found by SACOBRA: 2.329725 3.176549
  ## with objective -5.50801317   (error 1e-7 or smaller)
  cobra <- print_and_plot(cobra,c(1e-6,1e-0))
  return (cobra)
}

print_and_plot <- function(cobra,ylim) {
  testit::assert("[print_and_plot] cobra$solu is NULL", !is.null(cobra$solu))
  cat("xbest = \n")
  print(getXbest(cobra))          
  cat(sprintf("fbest = %12.6e\n",getFbest(cobra)))
  # if there are multiple solutions, they must have all the same objective value obj which we infer from first solution:
  firstSolu <- ifelse (is.matrix(solu), solu[1,], solu)
  obj <- cobra$fn(firstSolu)[1];   
  fn_solu = ifelse(is.matrix(cobra$solu),cobra$fn(cobra$solu[1,])[1],cobra$fn(cobra$solu)[1])
  cobra$errGCOP = cobra$df$Best-fn_solu
  plot(abs(cobra$errGCOP),log="y",type="l", ylab="error",xlab="iteration",
       ylim=ylim,main=sprintf("%s, best: %.3f", cobra$fName,getFbest(cobra)))
  return (cobra)
}

multi_gfnc <- function(gfnc, runs, cobraSeed) {
  # consecutive runs to collect some statistics:
  ptm <- proc.time()
  finErrList = c()
  for (run in 1:runs) {
    cobra = gfnc(cobraSeed+run)
    finErr = tail(cobra$errGCOP,1)
    finErrList = c(finErrList, finErr)
    cat(sprintf("final error: %e\n", finErr))
  }
  print(sprintf("%10.5f %s",(proc.time()-ptm)[3]/runs,"sec per run"))
  print(summary(abs(finErrList)))
  cobra$finErrList = finErrList
  return (cobra)
}

#cobra = solve_G03(42, 5)
cobra = solve_G05(42)
#cobra = solve_G06_R(42)
#cobra = solve_G06(42)
#cobra = solve_G11(42)
#cobra = solve_G13(385)
#cobra = solve_G24(42)
#cobra = multi_gfnc(solve_G05, 6, 390)

