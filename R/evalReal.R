# evalReal.R
#
#' Evaluate new iterate on real function(s)
#'
#' Helper for \code{\link{cobraPhaseII}}: The new iterate \code{xNew}, which was found by optimization
#' on the surrogate models, is evaluated on the real function \code{cobra$fn}. In the case of equality
#' constraints, \code{evalReal} does the additional refine step (see Details).
#'
#' If \code{cobra$equHandle$active==TRUE}, then \code{xNew} is first \strong{refined}: The artificially feasible
#' solution \code{xNew} is replaced by a refined solution \code{ev1$xNew}. 
#' \code{ev1$xNew} is created by using \code{optim} to minimize the function 
#'    \deqn{ \sum_i{\max(0,g_i(x))} + \sum_j{h_j^2(x) } } 
#' Ideally, the refined solution \code{ev1$xNew} should be \strong{on} the equality constraints 
#' (within machine accuracy), but there is no guarantee that \code{optim} reaches this desired result.
#'
#' @param cobra an object of class COBRA, this is a (long) list containing all settings
#'        from \code{\link{cobraPhaseII}}
#' @param ev1   a list, initially empty, gradually filled by calls to \code{evalReal}  
#' @param xNew   the new point, see \code{\link{cobraPhaseII}}
#' @param fValue   fitness value estimated for  \code{xNew} 
#' @param feval    function evaluations on surrogates needed by COBRA optimizer
#' @param optimConv   see \code{\link{cobraPhaseII}} 
#' @param optimTime   see \code{\link{cobraPhaseII}} 
#' @param currentEps   artificial current margin for the equality constraints: A point is said to be 
#'          \strong{artificially feasible}, if \eqn{h_j(x) - currentEps \le 0,    -h_j(x) - currentEps \le 0,} 
#'          for all equality constraints and if it is feasible in the inequality constraints.  
#' @param fitnessSurrogate  [\code{cobra$fitnessSurrogate}]   see \code{\link{cobraPhaseII}} 
#'
#' @return \code{ev1}, a list with the following \code{n}-dim vectors ( \code{n} = number of
#'      iterations, the last element is from the new iterate / point \code{xNew} ):
#'      \item{\code{predY}}{  prediction of  \code{fitnessSurrogate} at \code{xNew}}
#'      \item{\code{predVal}}{ \code{fvalue} (fitness + penalty in case of NMKB et al.)}
#'      \item{\code{feval}}{  function evaluations on surrogates needed by COBRA optimizer }
#'      \item{\code{optimizerConvergence}}{ see \code{\link{cobraPhaseII}} }
#'      \item{\code{optimizationTime}}{ see \code{\link{cobraPhaseII}} }
#'      \item{\code{predC}}{ prediction of  \code{cobra$constraintSurrogates} at \code{xNew} }
#'      \item{\code{feas}}{  TRUE, if \code{xNew} is feasible for the current constraints}
#'      \item{\code{feasPred}}{ TRUE, if \code{xNew} is feasible for \code{cobra$constraintSurrogates} }
#'      In addition, \code{ev1} has these elements:
#'      \item{\code{xNew}}{  \code{d}-dim vector, the new point, refined in the case of equality handling}
#       \item{\code{x_1}}{  -- in the end the same as xNew, only internal -- }
#'      \item{\code{xNewEval}}{  \code{cobra$fn(xNew)}, an \code{(1+nConstraints)}-dim vector 
#'                  (objective,constraints) }
#'      \item{\code{newNumViol}}{  scalar, the number of constraint violations (above 
#'                  \code{cobra$conTol}) on true constraints from \code{xNewEval} }
#'      \item{\code{newNumPred}}{  scalar, the number of constraint violations (above 
#'                  \code{cobra$conTol}) on constraint surrogates for \code{xNew} }
#'      \item{\code{newMaxViol}}{  scalar, the maximum constraint violation (with  
#'                  \code{currentEps} subtracted) on true constraints from \code{xNewEval} }
#'      \item{\code{trueMaxViol}}{  scalar, the maximum constraint violation (w/o  
#'                  \code{currentEps} subtracted) on true constraints from \code{xNewEval} }
#'                  
#'      If \code{cobra$equHandle$active==TRUE}, then the last four values are for \code{xNew} after 
#'      the refine step. In this case, the first three elements \code{newNumViol}, \code{newNumPred},
#'      and \code{newMaxViol} refer to the artificially enlarged equality constraints, i.e.
#'          \deqn{ h_j(x) - currentEps \le 0,    -h_j(x) - currentEps \le 0, }
#'      and the true inequality constraints \eqn{max(0,g_i(x))}. The last element \code{trueMaxViol}  
#'      measures the maximum violation among the true equality constraints \eqn{|h_j(x)|} and the  
#'      true inequality constraints \eqn{max(0,g_i(x))}. 
#'      
# @keywords internal   
#' @seealso  \code{\link{cobraPhaseII}}
#'      
evalReal <- function(cobra,ev1,xNew,fValue,feval,optimConv,optimTime,currentEps
                     ,fitnessSurrogate=cobra$fitnessSurrogate) {
  fn=cobra$fn
  ev1$xNew <- pmax(xNew,cobra$lower)  
  ev1$xNew <- pmin(xNew,cobra$upper)  
  if(cobra$equHandle$active){
    if(cobra$equHandle$refine & 
       attr(ev1,"state")=="optimized")    # do refine step only after surrogate optimizer (not after repair or TR)
    {
      # if (cobra$trueMaxViol[cobra$ibest] > cobra$equHandle$equEpsFinal) {
      if (cobra$trueMaxViol[cobra$ibest] >= 0) {    # testwise, enforce refine in every step
        #cat("[evalReal] starting refine ...\n")
        #
        # refine step, special for equality constraints: 
        # Due to currentEps, xNew will not fulfill the equality constraints exactly. 
        # Search with optim (i.e. CG, BFGS or similar) a point near to xNew which fulfills 
        # the equality constraints: Minimize the square sum of s_i(x) where s_i is the 
        # surrogate model for the ith constraint. The solution cg$par from optim replaces xNew.
        # 
        # The old version (deprecated): minimize equality violations only
        #myf <- function(x)sum(interpRBF(x,cobra$constraintSurrogates)[cobra$equIndex]^2); 
        #
        # The new version: minimize *all* constraint violations (inequalities + equalities):
        #
        #        sum( max(0,g_i(x))^2 )  + sum ( h_j(x)^2 )
        #
        myf <- function(x){
          conR=interpRBF(x,cobra$constraintSurrogates);
          sum(c(max(0,conR[-cobra$equIndex])^2,conR[cobra$equIndex]^2));
        }
        
        if (cobra$trueFuncForSurrogates)
          myf <- function(x){
            conR=cobra$fn(x)[-1];
            sum(c(max(0,conR[-cobra$equIndex])^2,conR[cobra$equIndex]^2)); 
          }
        
        # /WK/ the name "cg"'" is just a reminiscent that we initially used optim with method="CG" (conjugate gradient)
        cg = optim(ev1$xNew,myf,lower=cobra$lower,upper=cobra$upper,method="L-BFGS-B",
                   control=list(maxit=cobra$equHandle$refineMaxit))  
        # /WK/ bug fix: lower and upper added (otherwise cg$par might get out of bounds) 
        #      --> this requires method="L-BFGS-B"
        if (cg$convergence==1) { # iteration limit maxit has been reached
          warning(sprintf("Refine step: optim did not converge within maxit=%d iterations. %s",
                          cobra$equHandle$refineMaxit,
                          "Consider to increase cobra$equHandle$refineMaxit"))
        }
        if(! (cg$convergence %in% c(0,1))){print(cg$message)}
        if(any(is.na(cg$par))){browser()}
        
        #### NOTE: sometimes we see convergence code 52 with message from optim
        #### "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH" but the result looks very o.k. 
        #### Therefore we comment the following warning out
        #if (cg$convergence>0) {
        #  warning(sprintf("Refine step: optim did not converge, code=%d, message=%s",cg$convergence,cg$message))
        #}
        ev1$x_1 = cg$par
        
        # ------ The following is only debug/diagnostics: ----------------
        # /WK/ The debug-printout of the three lines below shows that optim succeeds in all cases to 
        #      turn the equality mismatch in cgbefore, which may be large (>1) in initial iterates,
        #      down to cg$value = 1e-8 or better. And cgtrue, the evaluation of cg$par on the true
        #      equality constraints (not the RBF models) gives usually 1e-8 or better.
        cgbefore = myf(ev1$xNew)
        conTrue = cobra$fn(ev1$x_1)[-1]                     # true constraints after optim
        cgtrue =  sum(c(max(0,conTrue[-cobra$equIndex])^2,conTrue[cobra$equIndex]^2))
        #### printout, only for debugging
        if(is.na(cgtrue))print("The solution returned by refine mechanism is not defined for the objective function, please check if the given lower and upper limits are correct")
        cat(sprintf("cg-values (before,after,true) = (%g, %g, %g)\n",cgbefore,cg$value,cgtrue))
        #
        # this is just more detailed constraint information, which may be inspected in browser:
        equInd = cobra$equIndex
        if (!cobra$trueFuncForSurrogates) {
          conB=interpRBF(ev1$xNew,cobra$constraintSurrogates) # constraint surrogates before optim
          conA=interpRBF(ev1$x_1,cobra$constraintSurrogates)  # constraint surrogates after optim 
        } else {
          conB=cobra$fn(ev1$xNew)[-1]   # true constraints before optim
          conA=cobra$fn(ev1$x_1)[-1]    # true constraints after optim 
          
        }
        trueB=cobra$fn(ev1$xNew)[-1]
        trueA=cobra$fn(ev1$x_1)[-1]
        conB[equInd] = abs(conB[equInd]) - currentEps
        conA[equInd] = abs(conA[equInd]) - currentEps
        trueB[equInd] = abs(trueB[equInd]) - currentEps
        trueA[equInd] = abs(trueA[equInd]) - currentEps
        pos_conB = which(conB > cobra$conTol)
        pos_conA = which(conA > cobra$conTol)
        pos_trueB = which(trueB > cobra$conTol)
        pos_trueA = which(trueA > cobra$conTol)
        ev1$nv_conB = length(pos_conB)     # will be added to df2 in updateSaveCobra
        ev1$nv_conA = length(pos_conA)
        ev1$nv_trueB = length(pos_trueB)
        ev1$nv_trueA = length(pos_trueA)
        check_df = data.frame(conB=conB, conA=conA, trueB=trueB, trueA=trueA,
                              diffB=pmax(0,trueB)-pmax(0,conB),
                              diffA=pmax(0,trueA)-pmax(0,conA))
        
        ### --- only for debugging the cgbefore peculiarity:
        #if (cgbefore < 1e-15) {
        #  browser()
        #}
        # ----------------------------------------------------------------
        
        ev1$xNew = ev1$x_1            # set ev1$xNew to be the infill point AFTER refine
        attr(ev1,"state") <- "refined"
        
      } # if (currentEps < ...)
    } # if (cobra$equHandle$refine)
  } # if(cobra$equHandle$active)
  
  newPredY <- getPredY0(ev1$xNew,fitnessSurrogate,cobra)
  if (cobra$trueFuncForSurrogates) newPredY<-fn(ev1$xNew)[1]   
  ev1$predY <- c(ev1$predY,newPredY)      # bug fix: now predY is the fitness surrogate value /WK/
  ev1$predVal <- c(ev1$predVal,fValue)   # fitness + penalty (in case of NMKB et al.) /WK/
  ev1$feval <- c(ev1$feval,feval)
  ev1$optimizerConvergence <- c(ev1$optimizerConvergence,optimConv)
  ev1$optimizationTime <- c(ev1$optimizationTime, optimTime )
  newPredC<-c()
  if(cobra$CONSTRAINED){
    newPredC <- interpRBF(ev1$xNew,cobra$constraintSurrogates)
    ev1$predC <- rbind(ev1$predC,newPredC)
  }
  
  ### Old version was
  ###   ev1$predC <- rbind(ev1$predC,constraintPrediction)
  ### but this would not be correct if we do the refine step and furthermore constraintPrediction
  ### is only known inside cobraPhaseII() /WK/04/2016
  
  #xNewTemp <- ev1$xNew
  ev1$xNewEval<-fn(ev1$xNew)
  if(cobra$CA$active && cobra$TFlag){
     xNewT<-(ev1$xNew-cobra$tCenter)%*%cobra$TM
     xNewT<-xNewT+cobra$tCenter
     ev1$xNewEvalT<-fn(xNewT)
  }
  
  if(!cobra$CONSTRAINED){
    ev1$newNumViol<-0
    ev1$newMaxViol<-0
    ev1$trueMaxViol<-0
    ev1$feas=T
  }
  if(cobra$CONSTRAINED)
    if (cobra$equHandle$active){
      temp<-ev1$xNewEval[-1]
      # /WK/the new version: we check whether 
      #
      #          g_i(x) <= 0,  h_j(x) - currentEps <= 0,    -h_j(x) - currentEps <= 0
      #
      # for the approximation newPredC with cobra$constraintSurrogates and set ev1$newNumViol to the
      # number of violated constraints.
      # NOTE that temp is also used for ev1$newMaxViol below.
      temp<-c(temp,-temp[cobra$equIndex])
      equ2Index <- c(cobra$equIndex,cobra$nConstraints+(1: length(cobra$equIndex)))
      temp[equ2Index] <- temp[equ2Index] - currentEps
      ev1$newNumViol<-length(which(temp > cobra$conTol)) # number of constraint violations for new point # 0 changed to conTol
      
      # just a debug check:
      if (attr(ev1,"state") %in% c("optimized","refined")) {
        if (ev1$newNumViol != ev1$nv_trueA) {
          browser()
        }
        # cat(attr(ev1,"state"), ev1$newNumViol, ev1$nv_trueA)
        testit::assert(ev1$newNumViol == ev1$nv_trueA)
        # (If state is not "optimized" at start of evalReal, then the branch that 
        #  computes nv_trueA will not been executed and the assertion would fail. Otherwise it should hold.)
      }

      ev1$feas = c(ev1$feas, ev1$newNumViol == 0 )
      #  # /WK/ this OLD currentFeas-calculation has a bug for mixed equality-inequality constraints
      #  #      (the inequality constraints are compared to currentEps, but should be compared to 0)
      #  temp[cobra$equIndex]<-abs(temp[cobra$equIndex])
      #  ev1$newNumViol<-length(which(temp-currentEps > cobra$conTol)) # number of constraint violations for new point #0 change to conTol
      
      
      # /WK/ brought here currentEps and equ2Index into play as well
      ptemp<- c(newPredC,-newPredC[cobra$equIndex])
      equ2Index <- c(cobra$equIndex,cobra$nConstraints+(1: length(cobra$equIndex)))
      ptemp[equ2Index] <- ptemp[equ2Index] - currentEps
      ev1$newNumPred<-length(which(ptemp > cobra$conTol)) # the same on constraint surrogates
      ev1$feasPred = c(ev1$feasPred, ev1$newNumPred < 1 ) 
      
      #WK: changed ev1$newMaxViol back to hold the artificial constraint max violation (currentEps-
      #    margin for equality constraints). This is one condition for entering repair (cobraPhaseII)
      M <-  max(0,max(temp))  # maximum violation
      try(if(M <= cobra$conTol){M<=0 })
      if(class(.Last.value)[1]=="try-error"){
        print("an error occurred in evalReal.R")
        browser()}
      else{
        # print("no exception")  
      }
      ev1$newMaxViol <- M  
      #SB: added the following lines, because it is also interesting to observe or save the information 
      #about the real maximum violation instead of the maximum distance to the artificial constraints 
      temp<-ev1$xNewEval[-1]
      temp[cobra$equIndex]<-abs(temp[cobra$equIndex])
      ev1$trueNumViol<-length(which(temp > cobra$conTol))
      #browser()
      M <-  max(0,max(temp*cobra$GRfact))  # maximum violation
      #  M <-  max(0,max(temp))  # maximum violation
      if(M <= cobra$conTol) M=0 
      ev1$trueMaxViol <- M  
      
      #--only debug
      #if (ev1$newMaxViol==0) browser()  
    }else{  # i.e. if (!cobra$equHandle$active)
      ev1$newNumViol<-length(which(ev1$xNewEval[-1] > cobra$conTol)) # number of constraint violations for new point #0 change to conTol
      ev1$feas = c(ev1$feas, ev1$newNumViol == 0 )
      ev1$newNumPred<-length(which(newPredC > cobra$conTol)) # the same on constraint surrogates
      ev1$feasPred = c(ev1$feasPred, ev1$newNumPred == 0 ) 
      
      if((max(0,max((ev1$xNewEval[-1])) )) > cobra$conTol){  # maximum violation
        ev1$newMaxViol<-max(0,max((ev1$xNewEval[-1])) )  
      }else{
        ev1$newMaxViol<-0
      }
      ev1$trueMaxViol = ev1$newMaxViol   # w/o artificial band around eq constraints, trueMaxViol and newMaxViol are the same
    } # (cobra$equHandle$active)
  
  return (ev1)
}
