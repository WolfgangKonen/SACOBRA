## solve G01 problem

## creating an instance for G01 
G01<-COP$new("G01")

set.seed(1)
xStart<-runif(d,min=G01$lower,max=G01$upper)

## Initializing cobra
cobra <- cobraInit(xStart=xStart, fn=G01$fn, fName=G01$name, lower=G01$lower, upper=G01$upper,
                   feval=70, initDesPoints=3*d, DOSAC=1, cobraSeed=1)

cobra$skipPhaseI = FALSE # TRUE #
cobra$seqOptimizer = "COBYLA" # "HJKB" 
## From the possible settings for cobra$seqOptimizer:
## COBYLA (default),DEOPTIM work best (error = 1e-6), but DEOPTIM terribly slow
## HJKB,NMKB work fine (error = 3e-4)
## ISRESCOBY works medium (error = 2e-5...2, depending on seed)
## ISRESN, ISRES do not work (error = 8...14)

## Running cobra optimizer
cobra <- startCobra(cobra)
## The true solution is at G01$solu = c(rep(1,9),rep(3,3),1)
## where the true optimum is G01$fn(G01$solu)[1] = optim = -15
optim = G01$fn(G01$solu)[1]
## The solutions from SACOBRA is close to this:
print(getXbest(cobra))
print(getFbest(cobra))

## Plot the resulting error (best-so-far feasible optimizer result - true optimum)
## on a logarithmic scale:
fprintf("%s %10.4e\n","error: ",getFbest(cobra)-optim)
# ---
#fb=cobra$df$Best
#fb[1:(which(cobra$df$feasible)[1]-1)] <- NA  # invalidate iterates before the 1st feasible point
# --- this is now done at the end of cobraPhaseII

plot(cobra$df$Best-optim,log="y",type="l",ylab="error",xlab="iteration",main=fName)
#print(cobra$df)

