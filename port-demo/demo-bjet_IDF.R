#library(R.matlab) 
## creating an MATLAB instance for MOPTA:
Matlab$startServer()	# start MATLAB server
matlab <- Matlab()	  # create MATLAB client object
isOpen <- open(matlab)

# compile a function from Don Jones automotive benchmark (see example.m):
evaluate(matlab,"cd('C:/Users/wolfgang/Documents/GitHub/2025JonesBenchmarks/Benchmarks/businessjet')")
evaluate(matlab,"[bjet,frel,grel,hrel,xl,xu,xopt] = businessjet_benchmark(true,true);")
data <- getVariable(matlab, c("xopt", "xl", "xu"))
solu = data$xopt
lower = data$xl
upper = data$xu
evaluate(matlab,"load('x0.mat','x0')");
data <- getVariable(matlab, c("x0"))
xStart = data$x0
# xStart = solu + 0.001
dim = length(xStart)
# evaluate(matlab,"[f,g,h] = automotive(x0)")

# define the problem function for SACOBRA, which calls the MATLAB client inside
# and which ensures that argument x is a matrix.
# It calls 'bjet' with x on the MATLAB side and returns the results
fn <- function(x) {
  testit::assert("length(x) is not 41", length(x)==41)
  if (is.vector(x)) x = matrix(x, nrow=1, ncol=41)
  # x has to be for MATLAB 'bjet' an (1 x 41)-matrix, NOT a vector (!)
  setVariable(matlab, x = x)      # transport x to MATLAB
  evaluate(matlab,"[f,g,h] = bjet(x)")
  data <- getVariable(matlab, c("f", "g", "h"))
  xEval = c(data$f, data$g)
  testit::assert("length(data$h) is not 10", length(data$h)==10)
  for (i in 1:length(data$h)) {
    xEval = c(xEval, equ=data$h[i])
  }
  return (xEval)
}

solve_bjet <- function(cobraSeed) {
  
  ## Initializing SACOBRA  
  equHandle <- defaultEquMu()
  equHandle$initType = "useGrange" # default "TAV" produces factor 15 too large init value for mu
  cobra <- cobraInit(xStart=xStart, fn=fn, fName="bjet_IDF", 
                     lower=lower, upper=upper, 
                     #feval=800, initDesPoints=83,     # for a quick run
                     feval=1000, initDesPoints=183,     # for a quick run
                     # feval=500, # 1860, #initDesPoints = 249, # the default 2*d+1 
                     initDesign = "LHS",  # "OPTIMIZED" does not work with MATLAB
                     skipPhaseI = T,
                     # initDesFactory2 = T, 
                     cobraSeed=cobraSeed, verboseIter=2,
                     epsilonInit = 0.0, epsilonMax = 0.0,
                     equHandle = equHandle, 
                     muGrow=100, mu4inequality = TRUE,
                     seqFeval=5000, seqOptimizer="COBYLA", conTol=2e-3, # 1e-10,
                     XI=c(0.3,0.01, 0.001, 0.0), # 0.0005), # distance requirement parameter
                     # sac=list(aDRC=F),   # no automatic DRC
                     solu=as.vector(solu),   
                     repairInfea=TRUE,
                     ri=list(repairMargin=0.2 ,OLD=FALSE))
  cobra$squares = T
  # browser()
  cobra$sac$RS = T
  cobra$sac$RS_rep = F
  cobra$sac$RSmin = 0.15
  cobra$equHandle$dec=1.2471           # mu should decay in 53 steps to 1e-5
  #cobra$equHandle$dec = 1.5 # 1.073534 # mu should decay in about 200 time steps to 1e-5
  cobra$equHandle$equEpsFinal = 1e-5   # set it larger than the default 1e-7
  cobra$equHandle$refineMaxit = 2000
  # cobra$trueFuncForSurrogates = T    # does not work

  ## Run SACOBRA optimizer 
  cobra <- startCobra(cobra)   # cobraPhaseII(cobra)

  ## The true solution is at solu = [0.048743 60000.000000 ...]
  ## with objective 32.63534 
  ## xStart (load('x0.mat','x0')) has objective 50.38838 and is infeasible
  ## The solution found by SACOBRA, if we start with solu+0.001:
  print(getXbest(cobra))          # 0.0484 59980.7728 ....
  print(getFbest(cobra))          # 32.79498 
  ## and if we start with solu+0.01:
  ##print(getXbest(cobra))        # 0.050232 59452.458810 ....
  ##print(getFbest(cobra))        # 33.15005
  cobra$errGCOP = abs(cobra$df$Best-fn(solu)[1])
  plot(cobra$errGCOP,log="y",type="l", ylim=c(0.5,12),
       ylab="error",xlab="iteration",main=sprintf("%s, best: %.3f", cobra$fName,getFbest(cobra)) )
  return (cobra)
}


cobra = solve_bjet(42+4)
# good solutions for seed=42,43
# only bad solutions (no feasible point found) for seed=44,45


solu2 = forwardRescale(solu,cobra);
dlA = distLine(as.vector(solu2),cobra$A)
mydf2 = data.frame(nViol=cobra$df$nViolations, maxViol=cobra$df$maxViolation, XI=cobra$df$XI, distSolu=dlA, y=cobra$df$y)

close(matlab)   # close the MATLAB client (will also shutdown MATLAB server)
print(matlab)  

