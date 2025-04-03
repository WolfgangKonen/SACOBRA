## creating an instance for MOPTA
path = "C:/Users/wolfgang/Documents/GitHub/2025JonesBenchmarks/Benchmarks/automotive/"
solu =  as.numeric(read.table(file=paste0(path,"solu.csv"), dec=".",sep=","))
#xStart = as.numeric(read.table(file=paste0(path,"xStart.csv"), dec=".",sep=","))
xStart = solu + 0.001
dim = length(xStart)
lower = rep(0,dim)
upper = rep(1,dim)
fn <- function(x) {
  write.table(x,paste0(path,"infill.csv"), row.names=F, col.names=F, dec=".", sep=",")
  fileEval = paste0(path,"infillEval.csv")
  xEval = c()
  while (length(xEval)==0) {      
    if (file.exists(fileEval)) {
      xEval <- as.numeric(read.table(fileEval, dec=".",sep=","))
      testit::assert("length(xEval) is not 69", length(xEval)==69)
      file.rename(fileEval, sub(".csv","",fileEval))
      # file.remove(fileEval)
    } else {
      Sys.sleep(0.2)     # wait for MATLAB process to create "infillEval.csv"
    }
  }
  return (xEval)
}

solve_mopta <- function(cobraSeed) {
  
  ## Initializing SACOBRA  
  cobra <- cobraInit(xStart=xStart, fn=fn, fName="MOPTA", 
                     lower=lower, upper=upper, feval=30,
                     initDesign = "RAND_R", initDesPoints = 10, 
                     initDesFactory2 = T, cobraSeed=cobraSeed)
  cobra$squares = T
  # browser()
  cobra$sac$RS = F    
  cobra$sac$RS_rep = T

  ## Run SACOBRA optimizer 
  cobra <- cobraPhaseII(cobra)
  browser()
  
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

