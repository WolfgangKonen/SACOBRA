######################################################################################
#startCobra
#
#' Start COBRA (constraint-based optimization) phase I and/or phase II
#' 
#' Start COBRA (constraint-based optimization) phase I and/or phase II for object \code{cobra} 
#' 
#' @param cobra initialized COBRA object, i.e. the return value from \code{\link{cobraInit}}
#'
#' @return \code{cobra},  an object of class COBRA
#'  
#' @examples
#' 
#' ## solve G01 problem
#' 
#' ## creating an instance for G01 
#' G01<-COP$new("G01")
#' 
#' set.seed(1)
#' xStart<-runif(d,min=G01$lower,max=G01$upper)
#' 
#' ## Initializing cobra
#' cobra <- cobraInit(xStart=xStart, fn=G01$fn, fName=G01$name, lower=G01$lower, upper=G01$upper,
#'                     feval=70, initDesPoints=3*d, DOSAC=1, cobraSeed=1)
#' 
#' ## Running cobra optimizer                    
#' cobra <- startCobra(cobra)
#' ## The true solution is at G01$solu = c(rep(1,9),rep(3,3),1)
#' ## where the true optimum is G01$fn(G01$solu)[1] = optim = -15
#' optim = G01$fn(G01$solu)[1]
#' ## The solutions from SACOBRA is close to this:
#' print(getXbest(cobra))
#' print(getFbest(cobra))
#' 
#' ## Plot the resulting error (best-so-far feasible optimizer result - true optimum)
#' ## on a logarithmic scale:
#' optim = -15
#' plot(cobra$df$Best-optim,log="y",type="l",ylab="error",xlab="iteration",main=fName)
#'                             
#' @seealso   \code{\link{cobraInit}}, \code{\link{cobraPhaseI}}, \code{\link{cobraPhaseII}}
#' @export
######################################################################################
startCobra <- function(cobra){
  feasibleSolutionExists<-(0 %in% cobra$numViol | cobra$skipPhaseI)
  if(feasibleSolutionExists){    # If there is any point with no constraint violation
    cobra <- cobraPhaseII(cobra)
    return(cobra)
  }else{
    
    verboseprint(verbose=2,important=FALSE,"Starting COBRA PHASE I ")
    res1 <- cobraPhaseI(cobra)
    cobra  <- cobraPhaseII(res1)
    return(cobra)
  }
  
}
