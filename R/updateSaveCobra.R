# updateSaveCobra.R
#
#' Update and save cobra
#'
#' Helper for \code{\link{cobraPhaseII}}: make some assertion tests, 
#' update elements in object \code{cobra}, including data frames \code{df} and \code{df2}, 
#' and - if \code{cobra$saveIntermediate==TRUE} - save cobra in subdir \code{results/}.
#' Most importantly \code{cobra$xbest, cobra$fbest, cobra$ibest} are updated. They characterize the 
#' best feasible point (least violating point if no feasible point was found so far) and
#' influence the next starting point.
#'
#' Note: the elements \code{A, Fres, Gres} of \code{cobra} are set in \code{updateInfoAndCounters},
#' an internal function of \code{\link{cobraPhaseII}}.
#'
#' @param cobra an object of class COBRA, this is a (long) list containing all settings
#'        from \code{\link{cobraPhaseII}}
#' @param ev1   a list filled by calls to \code{\link{evalReal}}. We need here the elements xNew,
#'        feas, feasPred, feval, optimizerConvergence, optimizationTime, predY, predVal  
#' @param subMin   see \code{\link{cobraPhaseII}} 
#' @param sigmaD   see \code{\link{cobraPhaseII}} 
#' @param penaF   see \code{\link{cobraPhaseII}}
#' @param gama   see \code{\link{cobraPhaseII}} 
#' @param EPS   see \code{\link{cobraPhaseII}}
#' @param fitFuncPenalRBF   helper function from \code{\link{cobraPhaseII}, only needed for predSoluPenal}
#' @param distRequirement   helper function from \code{\link{cobraPhaseII}}
#' @param fitnessSurrogate  the model used for \code{predSoluFunc}
#'
#' @return \code{cobra}, an object of class COBRA, 
#'    enhanced here by the following elements (among others):
#'      \item{\code{df}}{  data frame with summary of the optimization run (see below)}
#'      \item{\code{df2}}{  data frame with additional summary information (see below)}
#'      \item{\code{dftr}}{  data frame with additional summary information for TR (see \code{\link{cobraPhaseII}})}
#'      \item{\code{fbest}}{ the best feasible objective value found }
#'      \item{\code{xbest}}{ the point in input space yielding the best feasible objective value }
#'      \item{\code{ibest}}{ the corresponding iteration number (row of \code{cobra$df}, of \code{cobra$A}) }
#'      \item{\code{fbestArray}}{ vector of all \code{fbest} }
#'      \item{\code{xbestArray}}{ vector of all \code{xbest} }
#'      
#'   The data frame \strong{cobra$df} contains one row per iteration with columns 
#'   \describe{
#'      \item{iter}{ iteration index }
#'      \item{y}{   true objective value \code{Fres} }
#'      \item{predY}{  surrogate objective value. Note: The surrogate may be trained on  
#'            plog-transformed training data, but \code{predY} is transformed back to the original 
#'            objective range. NA for the initial design points.}
#'      \item{predSolu}{  surrogate objective value at best-known solution \code{cobra$solu}, if given. 
#'            If \code{cobra$solu} is NULL, take the current point instead. Note: The surrogate may be trained on  
#'            plog-transformed training data, but \code{predSolu} is transformed back to the original 
#'            objective range. NA for the initial design points.}
#'      \item{feasible}{ boolean indicating the feasibiltiy of infill point }
#'      \item{feasPred}{ boolean indicating if each infill point is feasible for \code{cobra$constraintSurrogates} }
#'      \item{nViolations}{ number of violated constraints }
#'      \item{maxViolation}{ maximum constraint violation. }
#'      \item{FEval}{  number of function evaluations in sequential optimizer. NA if it was a repair step }
#'      \item{Best}{  ever-best feasible objective value \code{fbest}. As long as there is 
#'            no feasible point, take among those with minimum number of violated constraints the
#'            one with minimum Fres. }
#'      \item{optimizer}{ e.g. "COBYLA"  }
#'      \item{optimizationTime}{  in sec}
#'      \item{conv}{ optimizer convergence code }
#'      \item{distSolu}{ distance of the current point (row of \code{cobra$A}) to the true solution 
#'            \code{cobra$solu} in rescaled space. If there is more than one solution, take the one
#'            which has the minimum distance element (since this is the solution to which the 
#'            current run converges). If \code{cobra$solu} is NULL, use instead the optimal 
#'            solution found so far (\code{submin$par})}
#'      \item{distOrig}{ same as \code{distSolu}, but in original space  }
#'      \item{distX0}{ same as \code{distSolu}, but with distance to the initial starting point \code{cobra$x0}}
#'      \item{XI}{  the DRC element used in the current iteration }
#'      \item{seed}{ the seed used in each run }
#'   }
#'   The columns \code{feasPred, nViolations, maxViolation} are only present for constrained
#'   problems (i.e. if \code{cobra$CONSTRAINED=T}).
#'
#'   The data frame \strong{cobra$df2} contains one row per phase-II-iteration with columns 
#'   \describe{
#'      \item{iter}{ iteration index}
#'      \item{predY}{  surrogate objective value. Note: The surrogate may be trained on  
#'            plog-transformed training data, but \code{predY} is transformed back to the original 
#'            objective range. NA for the initial design points.}
#'      \item{predVal}{   surrogate objective value + penalty }
#'      \item{predSolu}{   surrogate objective value at true solution (see \code{cobra$df$predSolu}) }
#'      \item{predSoluPenal}{   surrogate objective value + penalty at true solution (only diagnostics)}
#'      \item{sigmaD}{ the sigmaD element used in the current iteration (see \code{\link{cobraInit}})  }
#'      \item{penaF}{ penalty factor used in the current iteration (see \code{\link{cobraInit}}) }
#'      \item{XI}{  the DRC element used in the current iteration }
#'      \item{EPS}{ the current used margin for constraint function modeling (see \code{epsilonInit} in \code{\link{cobraInit}} ) }
#'   }
#'   
#' @seealso   \code{\link{cobraPhaseII}}, \code{\link{cobraInit}}
#' @export
#'      
updateSaveCobra <- function(cobra,ev1,subMin,sigmaD,penaF,gama,EPS,
                            fitFuncPenalRBF,distRequirement,
                            fitnessSurrogate=cobra$fitnessSurrogate)
{
  xNew=ev1$xNew;
  feas=ev1$feas;
  feasPred=ev1$feasPred;
  feval=ev1$feval;
  optimizerConvergence=ev1$optimizerConvergence;
  predY=ev1$predY;
  predVal=ev1$predVal;
  df_RS  <-  cobra$df$RS
  if (is.null(df_RS)) df_RS = c(rep(FALSE,cobra$initDesPoints))
  diff<-nrow(cobra$A)-length(predY)
  #deprecated
  # if(cobra$LocalExpensiveSearch && nrow(cobra$A)>length(predY)){
  # 
  #   ev1$feas<-feas<-c(feas,rep(NA,diff))
  #   ev1$feasPred<-feasPred<-c(feasPred,rep(NA,diff))
  #   ev1$feval<-feval<-c(feval,rep(NA,diff))
  #   ev1$optimizerConvergence<-optimizerConvergence<-c(optimizerConvergence,rep(NA,diff))
  #   ev1$predY<-predY<-c(predY,rep(NA,diff))
  #   ev1$predVal<-predVal<-c(predVal,rep(NA,diff)) 
  #   ev1$optimizationTime<-c(ev1$optimizationTime,rep(NA,diff))
  #   df_RS  <-  c(cobra$df$RS,rep(NA,diff))
  # }
  is_RS = (cobra$RSDONE[length(cobra$RSDONE)]=="RS")
  df_RS  <-  c(df_RS,is_RS)
  if (cobra$DEBUG_XI) {
    df_fxStart <- c(cobra$df$fxStart,cobra$fn(cobra$xStart)[1])
    df_fxbest <-  c(cobra$df$fxbest,cobra$fn(cobra$xbest)[1])
    df_RS2 <-  c(cobra$df$RS2,cobra$DEBUG_RS)    
  }
  
  if (cobra$WRITE_XI) {
    ro<-gama*cobra$l
    sumViol <- distRequirement(xNew,cobra$fitnessSurrogate,ro)$sumViol
    if (is.null(cobra$df)) {
      df_XI <- c(rep(NA,cobra$initDesPoints),gama)
      df_XIsumViol <- c(rep(NA,cobra$initDesPoints),sumViol)
    } else {
      df_XI <- c(cobra$df$XI,gama)
      df_XIsumViol <- c(cobra$df$XIsumViol,sumViol)
    }
    # cat("gama:",gama,"\n")
    
  }
  xNewIndex<-length(cobra$numViol)
  
  #testit::assert(cobra$fbest<=min(cobra$fbestArray[which(cobra$numViol==0)]))
  
  ### /WK/ Bug fix: the 3 commented lines below were mind-buggingly complex and not 
  ### correct. In some cases, cobra$fbestArray (a.k.a df$Best) would *increase* in value!!
  ### Changed it to a simpler logic with the help of cobra$ibest, which is - if set - 
  ### the index to the so-far-best feasible solution. If there is no feasible solution yet 
  ### it is NULL.
  ###
  #    index<-which(cobra$Fres==cobra$fbest)
  #    if(cobra$numViol[tail(index,1)]==0){       #If the so-far-best is a feasible point
  #      if((cobra$numViol[xNewIndex]==0 )&&(cobra$Fres[xNewIndex] < cobra$fbest )){  #If xNew is feasible and even better
  #SB: this condition is always true therefore replaced with Samineh's suggestion   
  #if (!is.null(cobra$ibest)){   # if we have an everBestFeasible at index ibest and fitness value fbest
  #testit::assert(cobra$fn(cobra$A[cobra$ibest,])[1]==cobra$fbest)
  #Samineh's suggestion:
  
  
  if(cobra$equHandle$active && (length(cobra$equIndex)!=0)){
    # ... if we handle equality constraints by 'equMargin & two inequality constraints' 
    # and have equality constraints in the actual problem: 
    # calculate cobra$xbest,fbest,ibest via updateCobraEqu:
    cobra<-updateCobraEqu(cobra,xNew)      # in modifyEquCons.R
  }else if(cobra$numViol[cobra$ibest]==0){    # if the so-far best is feasible...
    #testit::assert(cobra$fn(cobra$A[cobra$ibest,])[1]==cobra$fbest)
    testit::assert(cobra$Fres[cobra$ibest]==cobra$fbest)
    if((cobra$numViol[xNewIndex]==0 )&&(cobra$Fres[xNewIndex] < cobra$fbest )){  #If xNew is feasible and even better
      cobra$xbest<-xNew
      cobra$fbest<-cobra$Fres[xNewIndex]
      cobra$ibest<-xNewIndex    
      #cobra$noProgressCount<-0
    }
  }else{                                # if we have not yet an everBestFeasible ...
    if(cobra$numViol[xNewIndex]==0){# the new point is feasible then we select it
      cobra$xbest<-xNew                 # ... take xNew, if it is feasible
      cobra$fbest<-cobra$Fres[xNewIndex]
      cobra$ibest<-xNewIndex
     # cobra$noProgressCount<-0
    }else{        #SB12.10.2015: the new solution is infeasible: look for the best infeasible solu
      if(cobra$maxViol[xNewIndex] < cobra$maxViol[cobra$ibest]){
        cobra$xbest<-xNew                 # ... take xNew, if it has smaller maxViol
        cobra$fbest<-cobra$Fres[xNewIndex]
        cobra$ibest<-xNewIndex
      #  cobra$noProgressCount<-0
      }
    }
    
  }
  # If we do not have a feasible point AND xNew has a larger maxMivol than ibest, then leave the 
  # triple (xbest,fbest,ibest) at the setting of cobraInit.R, line 400: From all points 
  # with minimum number of violated constraints, take the one with smallest Fres.
  
  #--- only debug ---
  #tn = fn(cobra$xbest)[-1]
  #cat(sprintf("max.viol(xbest)=%9.4g\n",max(0,max(tn))))
 
  cobra$fbestArray<-c(cobra$fbestArray,cobra$fbest)
  cobra$xbestArray<-rbind(cobra$xbestArray,cobra$xbest)
  
  # feasibleIndices and xbestIndex are never used:
  feasibleIndices <- which(sapply(1:nrow(cobra$Gres),FUN=function(i)(all(cobra$Gres[i,]<0))))
  xbestIndex<-which.min(cobra$Fres[feasibleIndices])  # finding index of the best point so far
  
  # only diagnostics, needed for cobra$df & cobra$df2 /WK/
  solu <- cobra$solu; 
  if (is.null(solu)) {
    solu=subMin$par;    # subMin$par: the optimal solution found so far
    soluOrig <- inverseRescale(subMin$par,cobra);
  } else {
    # if (cobra$rescale) 
    #   if (is.matrix(solu)) {
    #     solu <- t(sapply(1:nrow(solu),function(i){ forwardRescale(solu[i,],cobra)}))
    #   } else {
    #     solu <- forwardRescale(solu,cobra);
    #   }
    #soluOrig <- cobra$solu;
    ### /WK/2025/03: bug fix: cobra$solu was already rescaled in cobraInit in 
    ### case cobra$rescale=TRUE
    soluOrig <- cobra$originalSolu
  }
  # now solu is always in *rescaled* input space
  
  predSoluFunc <- function(x)getPredY0(x,fitnessSurrogate,cobra);
  if (is.matrix(solu)) {      # in case of multiple global optima in solu:
    predSolu <- sapply(1:nrow(solu),function(i){ predSoluFunc(solu[i,])}) ;
    predSoluPenal <- sapply(1:nrow(solu),function(i){ fitFuncPenalRBF(solu[i,])}) ;
  } else {
    predSolu <- predSoluFunc(solu);      
    predSoluPenal <- fitFuncPenalRBF(solu);      
  }
  predSolu <- min(predSolu)   # Why min? - In case of multiple global optima: predSolu is the 
                              # value of fitFuncPenalRBF at the best solution solu
  predSoluPenal <- min(predSoluPenal) 
  # browser()
  if (is.null(cobra$df)) {
    df_predSolu <- c(rep(NA,cobra$initDesPoints),predSolu)
  } else {
    df_predSolu <- c(cobra$df$predSolu,predSolu)
  }
  #deprecated
  # if(cobra$LocalExpensiveSearch && nrow(cobra$A)>length(df_predSolu)){
  #   df_predSolu<-c(df_predSolu,rep(NA,diff))
  #   df_XI <- c(df_XI,rep(NA,diff))
  #   df_XIsumViol <- c(df_XIsumViol,rep(NA,diff))
  # }
  
  origA = t(sapply(1:nrow(cobra$A),function(i){ inverseRescale(cobra$A[i,],cobra) }))
  if (cobra$dimension==1) origA <- t(origA)
  if (is.matrix(solu)) {          # this is for the case with multiple solutions (like in G11)
    da=sapply(1:nrow(solu),function(i){ distLine(solu[i,],cobra$A) })
    # da has nrow(solu) columns, each column has the distance of the ith solu to all points cobra$A.
    # Select this column which has the element with minimum distance.
    ind=which.min(apply(da,2,min))
    distSolu = da[,ind]
    do=sapply(1:nrow(soluOrig),function(i){ distLine(soluOrig[i,],origA) })
    distOrig = do[,ind]
  } else {
    distSolu = distLine(solu,cobra$A)    # distance in rescaled space, distLine: see RbfInter.R
    distOrig = distLine(soluOrig,origA)  # distance in original space
  }
  distX0 = distLine(cobra$x0, cobra$A)   # distance in rescaled space to best point after initial design

  testit::assert("[updateSaveCobra] predY",length(cobra$Fres)==length(predY))
  testit::assert("[updateSaveCobra] optimizerConvergence",length(df_predSolu)==length(optimizerConvergence))
  if(cobra$CONSTRAINED){
    testit::assert("[updateSaveCobra] feas",length(cobra$Fres)==length(feas))
    testit::assert("[updateSaveCobra] feasPred",length(cobra$Fres)==length(feasPred))
    testit::assert("[updateSaveCobra] cobra$numViol",length(cobra$Fres)==length(cobra$numViol))
    testit::assert("[updateSaveCobra] cobra$maxViol",length(cobra$Fres)==length(cobra$maxViol))
  }
  testit::assert("[updateSaveCobra] cobra$fbestArray",length(cobra$Fres)==length(cobra$fbestArray))
  testit::assert("[updateSaveCobra] ev1$optimizationTime",length(cobra$Fres)==length(ev1$optimizationTime))
  testit::assert("[updateSaveCobra] optimizerConvergence",length(cobra$Fres)==length(optimizerConvergence))

  # result data frame
  if(cobra$CONSTRAINED){
    df <- data.frame(y=cobra$Fres, 
                     predY=predY,           # surrogate fitness
                     #predSolu=interpRBF(solu,cobra$fitnessSurrogate),  # OLD and WRONG
                     predSolu=df_predSolu,
                     feasible=feas, 
                     feasPred=feasPred,
                     nViolations=cobra$numViol,
                     trueNViol=cobra$trueNumViol,
                     maxViolation=cobra$maxViol,
                     trueMaxViol=cobra$trueMaxViol,
                     FEval=feval, 
                     Best=cobra$fbestArray,
                     optimizer=rep(cobra$seqOptimizer,length(cobra$Fres)),
                     optimizationTime=ev1$optimizationTime,
                     conv=optimizerConvergence,
                     distSolu=distSolu,
                     distOrig=distOrig,
                     distX0=distX0,
                     RS=df_RS,           # TRUE, if it is an iteration with random start point
                     row.names=NULL
    )
    
  }else{  # i.e. cobra$CONSTRAINED=F, an unconstrained optimization
    
    if(is.null(cobra$df$realfeval))
      realfeval<-c()
    else
      realfeval<-cobra$df$realfeval
    
    df <- data.frame(y=cobra$Fres, 
                     predY=predY,           # surrogate fitness
                     predSolu=df_predSolu,
                     feasible=T,
                     FEval=feval, 
                     #realfeval=c(realfeval,nrow(get("ARCHIVE",envir=intern.archive.env))),
                     realfeval=feval,   # /WK/2025/04/17/ temporary, use of intern.archive.env not clear
                     Best=cobra$fbestArray,
                     optimizer=rep(cobra$seqOptimizer,length(cobra$Fres)),
                     optimizationTime=ev1$optimizationTime,
                     conv=optimizerConvergence,
                     distSolu=distSolu,
                     distOrig=distOrig,
                     distX0=distX0,
                     RS=df_RS,           # TRUE, if it is an iteration with random start point
                     row.names=NULL
    )
    
  } # if(cobra$CONSTRAINED)
  
  if (cobra$WRITE_XI) {
    df$XI=df_XI
    df$XIsumViol=df_XIsumViol
  }
  
  if (cobra$DEBUG_XI) {
    firstSolu <- solu;
    if (is.matrix(solu)) firstSolu <- solu[1,]; 
    optimum <- cobra$fn(firstSolu)[1];      
    
    df$fxStart=df_fxStart # objective function at xStart
    df$fxbest=df_fxbest   # objective function at xbest
    df$exbest=df_fxbest - optimum   # error (objective function - optimum) at xbest
    df$RS2=df_RS2         # the same
    df$iter2=1:nrow(df)
    df$errFy=df$y - optimum # the error of the optimizer result in every iteration
    #if(tail(df_RS,1)==TRUE) browser()
    #browser()
    testit::assert(df$RS==df_RS2)
    if (any(df$fxbest[!df$RS]!=df$fxStart[!df$RS])) {
      browser()
      df$fxbest[!df$RS]=df$fxStart[!df$RS]  # symptomatic fix for the next assert
    }
    testit::assert(df$fxbest[!df$RS]==df$fxStart[!df$RS])
  }
  df <- cbind(iter=1:nrow(df),df)
  df <- cbind(df,seed=cobra$cobraSeed)
  cobra$df <- df
  
  # /WK/2025/03/24: only diagnostics: does the constraint prediction know that the solution
  # point(s) is/are feasible (at least in the later stages of the surrogates)?
  if (cobra$nConstraints==0 || is.null(cobra$solu)) {
    numViolSoluPred = 0
    maxViolSoluPred = 0
  } else {
    if (is.matrix(cobra$solu)) {
      # /WK/2025/04/12: bug fix for the case of multiple solutions
      constraintPredVec = sapply(1:nrow(cobra$solu), function(r) {
        interpRBF(cobra$solu[r,],cobra$constraintSurrogates)
      })
      numViolSoluPred = max(colSums(constraintPredVec > cobra$conTol))
      maxViolSoluPred = max(0, max(constraintPredVec-cobra$conTol))
    } else {
      constraintPrediction <-  interpRBF(cobra$solu,cobra$constraintSurrogates) 
      numViolSoluPred = sum(constraintPrediction > cobra$conTol)
      maxViolSoluPred = max(0, max(constraintPrediction-cobra$conTol))
    }
  }
  
  na_for_null <- function(v) {
    return (ifelse(is.null(v),NA,v))
  }
  cobra$df2 <- rbind(cobra$df2,data.frame(
    iter=tail(df$iter,1),
    predY=tail(predY,1),           # surrogate fitness at current point xNew  
    predVal=tail(predVal,1),       # surrogate fitness + penalty at xNew
    predSolu=predSolu,             # surrogate fitness at solu (only diagnostics). 
    predSoluPenal=predSoluPenal,   # surrogate fitness + penalty at solu (only diagnostics). 
    sigmaD=sigmaD[1],
    penaF=penaF[1],
    XI=gama,
    rho=cobra$RBFrho,
    fBest=tail(df$Best,1)
    , EPS=EPS[1]
    , currentEps= na_for_null(tail(cobra$currentEps,1))
    , numViolSoluPred=numViolSoluPred
    , maxViolSoluPred=maxViolSoluPred
    , PLOG=tail(cobra$PLOG,1)
    , pEffect=cobra$pEffect
    , err1=tail(cobra$err1,1)
    , err2=tail(cobra$err2,1)
    , nv_conB=na_for_null(ev1$nv_conB)        # diagnostics for refine mechanism
    , nv_conA=na_for_null(ev1$nv_conA)
    , nv_trueB=na_for_null(ev1$nv_trueB)
    , nv_trueA=na_for_null(ev1$nv_trueA)
    , nViol=tail(cobra$numViol,1)             # only for easier comparison with nv_trueA  
    , ev1_state = attr(ev1,"state")
    #,fBest2=fitFuncPenalRBF(xNew)    # not the same as predVal, since penaF or sigmaD might have changed (!)
    #,fSolu=fitFuncPenalRBF(min(solu))# not the same as predSolu for the same reason  
    #,feas=feas,
    #,data.frame(xNew,row.names=NULL)
  ))
  # browser()
  
  cobra$dftr<-rbind(cobra$dftr,data.frame(
    TRiter=cobra$TRiter,
    TRapprox=cobra$TRapprox,
    TFR=cobra$TFR,
    TSR=cobra$TSR,
    Fratio=cobra$Fratio,
    TRpop=cobra$TRpop,
    TRdelta=cobra$TRdelta
  ))
  
  # consistency check for data frames df and df2:
  msg <- "[updateSaveCobra] wrong nrow for df and df2";
  #browser()
  if (is.null(cobra$phase1DesignPoints)) {
    testit::assert(msg,nrow(cobra$df)==nrow(cobra$df2)+cobra$initDesPoints)
  } else {
    testit::assert(msg,nrow(cobra$df)==nrow(cobra$df2)+cobra$phase1DesignPoints)   
  }
  
  if (cobra$saveIntermediate) {
    # save intermediate results
    # cobraResult = list(cobra=cobra, df=df, constraintSurrogates=cobra$constraintSurrogates, fn=fn) 
    cobraResult = cobra
    if (is.na(file.info("results")$isdir)) dir.create("results")    # if directory "results" does not exist, create it
    save(cobraResult, file=sprintf("results/cobra-%s-%s-%i.RData",cobra$fName,cobra$seqOptimizer,cobra$cobraSeed))
  }
  cobra$ev1<-ev1
  
  #WK/2025/03: The following lines were moved from updateInfoAndCounters to here in updateSaveCobra, 
  #            because they use cobra$fbest and cobra$ibest which are updated only in updateSaveCobra:
  xNewEval = ev1$xNewEval
  newNumViol = ev1$newNumViol
  newMaxViol = ev1$newMaxViol
  currentEps = ev1$currentEps
  #SB: I added another print line which prints the best found solution after several iterations, if we do not have any interest for this the following lines can be commented
  realXbest<-inverseRescale(cobra$xbest,cobra)
  #realXbest<-sapply(1:length(cobra$xbest) , function(i){scales::rescale(cobra$xbest[i],from=c(cobra$newlower,cobra$newupper),to=c(cobra$lower[i],cobra$upper[i]))})
  if(cobra$equHandle$active){
    verboseprint(cobra$verbose, important=cobra$important,
                 sprintf("%s.[%d]: %.4f %.4f | %f | %.1e |[%.1e] (Fres=%f, nViol=%d, maxViol=%.1e)", 
                         "Best Result" ,nrow(cobra$A), realXbest[1], realXbest[2], cobra$fbest[1], 
                         cobra$trueMaxViol[cobra$ibest], currentEps, 
                         xNewEval[1], newNumViol, newMaxViol))
    
  }else{
    verboseprint(cobra$verbose, important=cobra$important,
                 sprintf("%s.[%d]: %.4f %.4f | %.1e | %.1e (Fres=%f, nViol=%d, maxViol=%.1e) " ,
                         "Best Result" ,nrow(cobra$A), #nrow(get("ARCHIVE",envir=intern.archive.env)), 
                         realXbest[1], realXbest[2], cobra$fbest[1], cobra$maxViol[cobra$ibest], xNewEval[1], newNumViol, newMaxViol))
    # verboseprint(cobra$verbose, important=cobra$important,
    #              sprintf("%s.[%d]: %f %f | %f | %f " ,
    #                      "Best Result" ,nrow(cobra$A), realXbest[1] ,realXbest[2] , cobra$fbest[1] , cobra$maxViol[cobra$ibest]))   
  }
  
  
  return(cobra)
} # updateSaveCobra()
