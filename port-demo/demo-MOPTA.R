library(R.matlab) 
## creating an MATLAB instance for MOPTA:
Matlab$startServer()	# start MATLAB server
matlab <- Matlab()	  # create MATLAB client object
isOpen <- open(matlab)

# compile a function from Don Jones automotive benchmark (see example.m):
evaluate(matlab,"cd('C:/Users/wolfgang/Documents/GitHub/2025JonesBenchmarks/Benchmarks/automotive')")
evaluate(matlab,"[automotive,frel,grel,hrel,xl,xu,xopt,x0] = automotive_benchmark();")
data <- getVariable(matlab, c("xopt", "x0"))
solu = data$xopt
xStart = data$x0
# xStart = solu + 0.001
dim = length(xStart)
lower = rep(0,dim)
upper = rep(1,dim)
# evaluate(matlab,"[f,g,h] = automotive(x0)")

# define the problem function for SACOBRA, which calls the MATLAB client inside
# and which ensures that argument x is a matrix.
# It calls 'automotive' with x on the MATLAB side and returns the results
fn <- function(x) {
  testit::assert("length(x) is not 124", length(x)==124)
  if (is.vector(x)) x = matrix(x, nrow=1, ncol=124)
  # x has to be for MATLAB 'automotive' an (1 x 124)-matrix, NOT a vector (!)
  setVariable(matlab, x = x)      # transport x to MATLAB
  evaluate(matlab,"[f,g,h] = automotive(x)")
  data <- getVariable(matlab, c("f", "g", "h"))
  xEval = c(data$f, data$g)
  return (xEval)
}

solve_mopta <- function(cobraSeed) {
  
  ## Initializing SACOBRA  
  browser()
  cobra <- cobraInit(xStart=xStart, fn=fn, fName="MOPTA", 
                     lower=lower, upper=upper, 
                     feval=90, initDesPoints=50,     # for a quick run
                     # feval=500, # 1860, #initDesPoints = 249, # the default 2*d+1 
                     initDesign = "LHS",  # "OPTIMIZED" does not work with MATLAB
                     # initDesFactory2 = T, 
                     cobraSeed=cobraSeed, verboseIter=2,
                     seqFeval=3000,
                     XI=c(0.3,0.01, 0.001, 0.0005), # distance requirement parameter
                     # sac=list(aDRC=F),   # no automatic DRC
                     repairInfea=TRUE ,
                     ri=list(repairMargin=0.2 ,OLD=FALSE))
  cobra$squares = T
  # browser()
  cobra$sac$RS = F    
  cobra$sac$RS_rep = T

  ## Run SACOBRA optimizer 
  cobra <- cobraPhaseII(cobra)

  ## The true solution is at solu ("solu.csv")
  ## with objective 222.6124 
  ## xStart ("xStart.csv") has objective 251.07
  ## The solution found by SACOBRA:
  print(getXbest(cobra))          # 14.0950000  0.8429608
  print(getFbest(cobra))          # -6961.814 
  cobra$errGCOP = abs(cobra$df$Best-fn(solu)[1])
  plot(cobra$errGCOP,log="y",type="l",
       ylab="error",xlab="iteration",main=G06$name)
  return (cobra)
}

cobra = solve_mopta(42)
close(matlab)   # close the MATLAB client (will also shutdown MATLAB server)
print(matlab)  

