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
  equHandle$initType = "useGrange"
  cobra <- cobraInit(xStart=xStart, fn=fn, fName="bjet_IDF", 
                     lower=lower, upper=upper, 
                     feval=800, initDesPoints=83,     # for a quick run
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
  cobra$sac$RS = F
  # cobra$sac$RS_rep = F
  cobra$equHandle$initType="useGrange" # default "TAV" produces factor 15 too large init value for mu
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

check_bjet <- function(cobra, ind4 = 344) {
  #
  # Step 1: Check whether we can recompute numViol from Gres:
  #         numViol should be equal to rs1 (or rs2) and it is.
  #
  temp=cobra$Gres
  npts=nrow(temp)
  equInd = cobra$equIndex 
  temp[,equInd] = abs(temp[,equInd])
  # rs1 holds for each row 1:npts the number of constraints that are violated 
  # (in their *artificial* feasibility)
  rs1 = sapply(1:npts , function(i){ 
    temp[i,equInd] <<- temp[i,equInd] - cobra$currentEps[i] 
    # IMPORTANT: use "<<-" to change temp on the check_bjet level (!)
    return(sum(temp[i,]  > cobra$conTol))   
  })
  cond1 = temp > cobra$conTol
  cond2=ifelse(temp > cobra$conTol,1,0)
  rs2 = rowSums(cond2)
  testit::assert(length(rs1)==length(cobra$df$nViolations))
  testit::assert(all(rs1==cobra$df$nViolations))
  testit::assert(all(rs2==cobra$df$nViolations))
  # data frame nv holds: rs1=row sums of cond1, rs2=..of cond2, nv=number of (constr.) violations, currentEps=mu for equ constr.
  #                      Condition rs1 == rs2 == nv should be true
  nv=data.frame(rs1=rs1, rs2=rs2, nv=cobra$df$nViolations, currentEps=cobra$currentEps)
  # View(nv)
  
  # Step 2: Matrix cond1 (or cond2) calculated in step 1 contains the info, which point (row)
  #         has which constraint (column) violated?
  #         Take in cs1 the colSums of cond1 for the last half of points: Which constraints
  #         are often and which constraints are never violated?
  #
  cs1 = colSums(cond1[(npts/2):npts,])
  cs2 = colSums(cond2[(npts/2):npts,])
  testit::assert(all(cs1==cs2))
  print(cs1)
  print(which(cs1>2))
  quant=sapply(1:88,function(i) {quantile(cobra$Gres[400:800,i],0.75)})
  plot(cs1,quant,col="blue")
  points(cs1[79:88],quant[79:88],col="red")
  # Inequality constraints that are often violated have 3rd quantile nv$quant near 0.0 (probably active).
  # Inequality constraints that are seldom violated have in most cases nv$quant near -0.75 (probabliy inactive).
  nv=data.frame(quant=quant,cs1=cs1, c_solu=cobra$fn(cobra$solu)[-1])
  viol=which(nv$cs1[1:78]>20)
  print(nv$c_solu[viol])
  # the fact that the printed numbers are all (with one exception) very small in magnitude 
  # shows that the inequality constraints often violated in our run belong to those 
  # which are *active* at the true solution. 
  # View(nv)

  # 
  # Step 3: Is the pattern different if we concentrate on the infill points with less then 6 violations?
  #         Answer: No.
  #
  cond1_1 = cond1[(npts/2):npts,]
  rs1_1 = rs1[(npts/2):npts]
  ind6 = which(rs1_1 < 6)
  cs1_1 = colSums(cond1_1[ind6,])
  #print(cs1_1)
  #print(which(cs1_1>2))
  # 
  # Step 4: Check the trueA calculation: Take the given index ind4 for step 4 and check whether the 
  #         number of constraints where trueA > conTol coincides with rs1[ind4]
  #
  trueA=cobra$fn(cobra$A[ind4,])[-1]
  trueA[equInd] = abs(trueA[equInd]) - cobra$currentEps[ind4]
  pos_trueA = which(trueA > cobra$conTol)
  testit::assert(length(pos_trueA)==rs1[ind4])
  browser()
  
  # return (cobra)
}

cobra = solve_bjet(42+4)
# good solutions for seed=42,43
# only bad solutions (no feasible point found) for seed=44,45

chk = check_bjet(cobra)

solu2 = forwardRescale(solu,cobra);
dlA = distLine(as.vector(solu2),cobra$A)
mydf2 = data.frame(nViol=cobra$df$nViolations, maxViol=cobra$df$maxViolation, XI=cobra$df$XI, distSolu=dlA, y=cobra$df$y)

close(matlab)   # close the MATLAB client (will also shutdown MATLAB server)
print(matlab)  

