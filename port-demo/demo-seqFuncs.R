# Demo for inner funcs subProb2 and gCOBRA (in innerFuncs.R)
#
# The printed results are copied into test_inner_funcs in testPhaseII.py
# and compared with the computation on the Python side.

require(pracma)   # needed for fprintf

nobs =  10       # number of observations for trainCubicRBF

demo_inner_funcs <- function(nobs) {
  cat("*** Starting demo_inner_funcs ... ***\n")

  d = 2
  xStart = c(2.5, 2.4)
  xflat = matrix(c(  1.3,  4.1,
                    -4.5, -2.3,
                     2.4, -2.1,
                    -3.9, -1.6), byrow=T,nrow=4,ncol=2)
  for (r in 1:nrow(xflat)) {
    xflat[r,]<-sapply(1:dimension, 
                   function(i){scales::rescale(xflat[r,i],to=c(-1,1),from=c(-5,5))})
  }
  # print(xflat)
  # browser()
  for (ro in c(1,0.5)) { 
    cat("ro = ",ro,"\n")
    fn <- function(x) {
      c(3 * sum(x ** 2),  sum(x)-1)
    }
    cobra <- cobraInit(xStart=xStart, fName="sphere", initDesign = "RAND_R",
                       fn=fn, lower=c(-5,-5), upper=c(5,5), feval=25)
    cobra$sac$TFRange = 500
    cobra$printP = TRUE
    cobra$squares = FALSE
    
    #globalOptCounter<-1
    #gama<-cobra$XI[((globalOptCounter) %% length(cobra$XI))+1]  
    cobra$ro <- ro # 1.0 # gama*cobra$l      
    cobra$EPS <- EPS <- cobra$epsilonInit 

    cobra = trainSurrogates(cobra)
    
    print(max(cobra$SurrogateInput) - min(cobra$SurrogateInput))
    print(cobra$SurrogateInput)
    yfit = predict.RBFinter(cobra$fitnessSurrogate, xflat)
    ycon = predict.RBFinter(cobra$constraintSurrogates, xflat)
    print(yfit)
    print(ycon)
    
    ### construct **now** cobra$gCOBRA_c from gCOBRA
    cobra$gCOBRA_c <- function(x){  # inequality constraints for nloptr::cobyla
      return(-gCOBRA(x,cobra))      # new convention h_i <= 0.
    } # Note that gCOBRA (in innerFuncs.R) contains the distance requirement calculation.
    
    yfit2 = subProb2(xflat,cobra)
    ycon2 = matrix(nrow=4,ncol=2)
    for (i in 1:nrow(xflat)) {
      ycon2[i,] = -gCOBRA(xflat[i,],cobra)
    }
    print(yfit2)
    print(ycon2)
  }
}

demo_inner_funcs(nobs)
G01 = COP$new("G01")
print(G01$solu+1)
print(G01$fn(G01$solu+1))
