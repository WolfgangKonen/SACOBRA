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
  for (fnfac in c(1,10,100)) {
    fn <- function(x) {
      c(3 * fnfac * sum(x ** 2),  sum(x)-1)
    }
    cobra <- cobraInit(xStart=xStart, fName="sphere", initDesign = "RAND_R",
                       fn=fn, lower=c(-5,-5), upper=c(5,5), feval=25)
    cobra$sac$TFRange = 500
    cobra$printP = TRUE
    cobra$squares = FALSE
    print(max(cobra$Fres) - min(cobra$Fres))
    
    cobra = trainSurrogates(cobra)
    
    print(max(cobra$SurrogateInput) - min(cobra$SurrogateInput))
    print(cobra$SurrogateInput)
    yfit = predict.RBFinter(cobra$fitnessSurrogate, xflat)
    ycon = predict.RBFinter(cobra$constraintSurrogates, xflat)
    print(yfit)
    print(ycon)

    for (i in 1:nrow(xflat)) {
      xNew = xflat[i,]
      xNewEval = fn(xNew)
      # print(c(i, getPredY1(xNew,cobra$fitnessSurrogate1,cobra), getPredY1(xNew,cobra$fitnessSurrogate2,cobra)))
      cobra = calcPEffect(cobra,xNew,xNewEval)  # appends to err1 and err2
    }
    print(cobra$err1/cobra$err2)
    print(c(cobra$PLOG, cobra$pEffect))
  }
}

demo_train_surr(nobs)
