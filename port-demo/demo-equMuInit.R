# Unit test test_init_design_R from testCobraInit.py ported to R
# ["C:\Users\wolfgang\Documents\GitHub\SACOBRA_Py\test\testCobraInit.py"]
#
require(pracma)   # needed for fprintf

single_equMuInit <- function(hnfac) {
  cat("\n*** Starting demo_equMuInit ... ***\n")
  d = 2
  equ = defaultEquMu()
  xStart = c(2.5, 2.4)
  initTVec = c("useGrange", "TAV", "TMV", "EMV")
  muInitVec = c()
  
  # define a problem that has one equality constraint: 
  fn <- function(x){c(obj=3*sum(x^2), equ=sum(x*hnfac+1)-2, equ=sum(x*hnfac)-1)}
  
  for (initT in initTVec) {
    equ$initType = initT
    cobra <- cobraInit(xStart=xStart, fName="sphere", fn=fn, 
                       lower=c(-5,-5), upper=c(5,5), feval=25,
                       equHandle=equ, initDesign = "RAND_R"
                       )
    #print(cobra$Gres)
    #fprintf("%12.8f\n",cobra$Gres) # just to get more digits after decimal point
    muInit = tail(cobra$currentEps,1)
    cat(initT, muInit)
    muInitVec = c(muInitVec, initT=muInit)
  }
  
  cat("\n")
  muInitVec=t(as.matrix(muInitVec))
  colnames(muInitVec) <- initTVec
  return(muInitVec)

  # These numbers are nearly exactly (np.allclose) reproduced by test_adCon_R on  
  # the Python side 
}

demo_equMuInit <- function() {
  muInitVec1 = single_equMuInit(1)
  #print(muInitVec1)
  muInitVec2 = single_equMuInit(10)
  fprintf("%12.8f,",muInitVec1) # just to get more digits after decimal point
  cat("\n")
  fprintf("%12.8f,",muInitVec2) 
  cat("\n")
  muInitVec = rbind(muInitVec1, muInitVec2)
  print(muInitVec)
  return (muInitVec)
}

demo_currentEps <- function() {
  cat("\n*** Starting demo_currentEps ... ***\n")
  d = 2
  idp = 6
  equ = defaultEquMu()
  xStart = c(2.5, 2.4)
  epsTVec = c("expFunc", "SAexpFunc", "funcDim", "funcSDim", "Zhang", "CONS")
  muMat = c()
  
  # define a problem that has one equality constraint: 
  fn <- function(x){c(obj=3*sum(x^2), equ=sum(x+1)-2, equ=sum(x)-1)}
  
  for (epsT in epsTVec) {
    equ$epsType = epsT
    cobra <- cobraInit(xStart=xStart, fName="sphere", fn=fn, 
                       lower=c(-5,-5), upper=c(5,5), feval=15,
                       equHandle=equ, initDesign = "RAND_R", initDesPoints = idp
    )
    #print(cobra$Gres)
    #fprintf("%12.8f\n",cobra$Gres) # just to get more digits after decimal point
    cobra <- cobraPhaseII(cobra)
    
    muVec = cobra$currentEps[(idp+1):(idp+5)]
    cat(epsT, muVec,"\n")
    muMat = rbind(muMat,muVec)
  }
  
  cat("\n")
  rownames(muMat) <- epsTVec
  colnames(muMat) <- (idp+1):(idp+5)
  print(muMat)
  return(muMat)
}

demo_refine1 <- function() {
  cat("\n*** Starting demo_refine1 ... ***\n")
  d = 2
  idp = 6
  equ = defaultEquMu()
  xStart = c(2.5, 2.4)
  epsTVec = c("expFunc")

  
  for (hnfac in c(1)) {   # c(1,10,100)
    # define a problem that has one equality constraint: 
    fn <- function(x){c(obj=3*sum(x^2), equ=sum(x*hnfac)-2)}

    equ$refine = F
    equ$refinePrint = TRUE
    cobra <- cobraInit(xStart=xStart, fName="sphere", fn=fn, squares=F,
                       lower=c(-5,-5), upper=c(5,5), feval=40,
                       equHandle=equ, initDesign = "RAND_R", initDesPoints = idp
    )
    cobra$finalEpsXiZero = TRUE
    #print(cobra$Gres)
    #fprintf("%12.8f\n",cobra$Gres) # just to get more digits after decimal point
    cobra <- cobraPhaseII(cobra)
    
    print(getXbest(cobra))
    print(getFbest(cobra))
    
  }
  
  return(cobra)
}

demo_refine2 <- function() {
  cat("\n*** Starting demo_refine2 ... ***\n")
  d = 2
  idp = 6
  equ = defaultEquMu()
  xStart = c(2.5, 2.4)

  
  for (hnfac in c(1)) {   # c(1,10,100)
    # define a problem that has two equality constraints: 
    fn <- function(x){c(obj=3*sum(x^2), equ=x[1]-x[2]-1, equ=sum(x*hnfac)-2)}

    equ$refine = F
    equ$refinePrint = TRUE
    feval = ifelse(hnfac==1, 30, 40)
    cobra <- cobraInit(xStart=xStart, fName="sphere", fn=fn, squares=F,
                       lower=c(-5,-5), upper=c(5,5), feval=feval,
                       equHandle=equ, initDesign = "RAND_R", initDesPoints = idp
    )
    cobra$finalEpsXiZero = TRUE
    #print(cobra$Gres)
    #fprintf("%12.8f\n",cobra$Gres) # just to get more digits after decimal point
    cobra <- cobraPhaseII(cobra)
    
    print(getXbest(cobra))
    print(getFbest(cobra))
  }
  
  return(cobra)
}



#muInitVec = demo_equMuInit()
#demo_currentEps()
cobra = demo_refine2()
