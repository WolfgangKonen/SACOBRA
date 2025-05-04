# Demo for train_surrogates and calcPEffect
#
# The printed results are copied into test_train_surr in testPhaseII.py
# and compared with the computation on the Python side.

require(pracma)   # needed for fprintf

nobs =  10       # number of observations for trainCubicRBF

demo_train_surr <- function(nobs) {
  cat("*** Starting demo_train_surr ... ***\n")

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
    #browser()

  }
}

demo_RS_EPS <- function() {
  cat("*** Starting demo_RS_EPS ... ***\n")
  
  d = 2
  xStart = c(2.5, 2.4)
  fn <- function(x) {
    c(3 * sum(x ** 2),  -(sum(x)-1))
  }
  for (rstype in c("SIGMOID", "CONSTANT")) {
    myseed = ifelse(rstype=="SIGMOID", 42, 52)
    cobra <- cobraInit(xStart=xStart, fName="sphere", initDesign = "RAND_R",
                       fn=fn, lower=c(-5,-5), upper=c(5,5), feval=20, 
                       cobraSeed = myseed)
    cobra$sac$TFRange = 500
    cobra$sac$RS_rep = TRUE
    cobra$sac$RStype = rstype
    cobra$printP = TRUE
    cobra$squares = FALSE
    cobra$trueFuncForSurrogates = TRUE
    
    ## Run SACOBRA optimizer 
    cobra <- cobraPhaseII(cobra)
    
    print(cobra$df$RS)
    print(cobra$df2$EPS)
    print(cobra$A[6:8,])
    print(cobra$A[20])
    print(getXbest(cobra))
  }
  
  return (cobra)
  
}

demo_phase2 <- function() {
  cat("*** Starting demo_phase2 ... ***\n")
  
  d = 2
  xStart = c(2.5, 2.4)
  fn <- function(x) {
    c(3 * sum(x ** 2),  -(sum(x)-1))
  }
  cobra <- cobraInit(xStart=xStart, fName="sphere", initDesign = "RAND_R",
                     fn=fn, lower=c(-5,-5), upper=c(5,5), feval=20)
  cobra$sac$TFRange = 500
  cobra$printP = TRUE
  cobra$squares = FALSE
  
  ## Run SACOBRA optimizer 
  cobra <- cobraPhaseII(cobra)
  
}

# demo_train_surr(nobs)
cobra = demo_RS_EPS()
#demo_phase2()
