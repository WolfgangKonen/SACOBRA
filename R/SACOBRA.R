#
#
#Samineh Bagheri 27.04.2015
#
#
#SACOBRA.R
#
#
#

######################################################################################
######################################################################################
# Package Description for Roxygen:
#' Self-adjusting Constrained Optimization with RBF Surrogates
#'
#' \tabular{ll}{
#' Package: \tab SACOBRA\cr
#' Type: \tab Package\cr
#' Version: \tab 1.3\cr
#' Date: \tab 20.01.2025\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' SACOBRA is a package for numeric constrained optimization of expensive black-box functions under severely 
#' limited budgets. The problem to solve is: 
#' \deqn{ \mbox{Minimize}\quad  f(\vec{x}) , \vec{x} \in [\vec{a},\vec{b}] \subset \mathbf{R}^d }
#' \deqn{ \mbox{subject to}\quad g_i(\vec{x}) \le 0, i=1,\ldots,m    }
#' \deqn{ \mbox{~~~~~~~~~~}\quad\quad h_j(\vec{x}) = 0, j=1,\ldots,r.    }
#' 
#' SACOBRA is an extension of the COBRA algorithm by Regis (R. Regis: "Constrained 
#' optimization by radial basis function interpolation for high-dimensional expensive black-box 
#' problems with infeasible initial points", Engineering Optimization, Taylor & Francis, 46, p. 218-243, 2013)
#' 
#' These extensions include: \cr
#' 1) A repair algorithm for infeasible solutions, \cr
#' 2) an algorithm for handling equality constraints, \cr
#' 3) several internal optimizers and several initial design generation methods,  \cr
#' 4) self-adjusting random restart algorithm,  \cr
#' 5) self-adjusting logarithmic transform for objective functions with large output ranges,  \cr
#' 6) range normalization of constraint functions, \cr
#' 7) self-adjusting DRC (distance requirement cycle) selection, \cr
#' 8) online model selection to select the best type of RBF for objective and constraint functions, \cr
#' 9) online whitening for unconstrained optimization of functions with high conditioning. \cr
#'(Please note that the online whitening implementation is still underway and at this stage it is not recommended to be applied to expensive problems) 
#' 
#' SACOBRA performs optimization with a minimum of true function evaluations. It has proven
#' to work well on problems with high dimensions (e.g. d=124) and many constraints (e.g. 60).
#' It is usable for all kind of numerical optimization of continuous functions, but not for combinatorial optimization.
#' 
#' For more details on SACOBRA see:\cr
#' \itemize{
#'    \item Bagheri, S.; Konen, W.; Baeck, T.: "How to Solve the Dilemma of Margin-Based Equality Handling Methods". 
#' In: Hoffmann, F. & Huellermeier, E. (Eds.),  Proceedings 28. Workshop Computational Intelligence, 
#' Universitaetsverlag Karlsruhe, 2018, 257-270, \strong{Young Author Award}, 
#' \url{https://blogs.gm.fh-koeln.de/ciop/files/2018/12/GMA2018.pdf} \cr
#'    \item Bagheri, S.; Konen, W.; Emmerich, M.; Baeck, T.: "Self-adjusting parameter control for surrogate-assisted constrained optimization under limited budgets".
#' In: Journal of Applied Soft Computing, Band 61, pages 377-393, 2017,
#' \url{https://www.researchgate.net/publication/319012980_Self-adjusting_parameter_control_for_surrogate-assisted_constrained_optimization_under_limited_budgets} \cr   
#'    \item Bagheri, S.; Konen, W.; Baeck, T.:"Equality constraint handling for surrogate-assisted constrained optimization".
#' In: Tan, K.C. (ed.), WCCI'2016, Vancouver, 1924-1931, IEEE, 2016,
#' \url{https://www.gm.fh-koeln.de/~konen/Publikationen/Bagh16-WCCI.pdf}   \cr
#'    \item Bagheri, S.; Konen, W.; Baeck, T.:"Online selection of surrogate models for constrained black-box optimization".
#' In: IEEE Symposium Series on Computational Intelligence (SSCI), 2016, \strong{Best Student Paper Award},
#' \url{https://www.gm.fh-koeln.de/~konen/Publikationen/Bagh16-SSCI.pdf}   \cr
#'    \item Bagheri, S.; Konen, W.; Baeck, T. et al.: "SACOBRA: Self-Adjusting Constrained Black-Box Optimization with RBF". 
#' In: Hoffmann, F. & Huellermeier, E. (Eds.),  Proceedings 25. Workshop Computational Intelligence, 
#' Universitaetsverlag Karlsruhe, 2015, 87-98, \strong{Young Author Award}, 
#' \url{https://www.gm.fh-koeln.de/~konen/Publikationen/Bagh15-GMA-CI.pdf} \cr
#'    \item Koch, P.; Bagheri, S.; Konen, W. et al.: "A New Repair Method For Constrained Optimization". 
#' In: Proceedings of the 17th Genetic and Evolutionary Computation Conference, 2015, 
#' \url{https://www.gm.fh-koeln.de/~konen/Publikationen/Koch2015a-GECCO.pdf} \cr
#'    \item Koch, P.; Bagheri, S. et al.: "Constrained Optimization with a Limited Number of Function Evaluations"
#' In: Hoffmann, F. & Huellermeier, E. (Eds.),  Proceedings 24. Workshop Computational Intelligence, 
#' Universitaetsverlag Karlsruhe, 2014, 119-134, 
#' \url{https://www.gm.fh-koeln.de/~konen/Publikationen/Koch2014a-GMA-CI.pdf}.
#' }
# and \cr
#' 
#' The main entry point functions are \code{\link{cobraInit}} and \code{\link{startCobra}}. 
#' See \code{\link{cobraInit}} for an overview of adjustable SACOBRA-parameters.
#' Examples are found in
#' \itemize{
#'    \item \code{\link{startCobra}}: solve a 13d-problem with 9 inequality constraints (G01)
#'    \item \code{\link{cobraInit}}: a problem with equality constraint 
#'    \item \code{\link{cobraPhaseII}}: unconstrained sphere problem 
#'    \item \code{\link{multiCOBRA}}: solve G11 problem nrun=4 times 
#'    \item \code{\link{COP}}: load and solve G24, load and solve the scalable problem G03 with d=3 
#' }
#'                                                
#' @name SACOBRA-package
#' @aliases SACOBRA
#' @docType _PACKAGE
#' @title Self-adjusting Constrained Optimization with RBF Surrogates
#' @author Samineh Bagheri (\email{Samineh.Bagheri@@th-koeln.de}), \cr
#' Wolfgang Konen (\email{Wolfgang.Konen@@th-koeln.de}), \cr
#' Patrick Koch, Thomas Baeck (\email{t.h.w.baeck@liacs.leidenuniv.nl})
#' @references \url{https://blogs.gm.fh-koeln.de/ciop/research/monrep/}
#' @keywords package constraints optimization RBF surrogate black-box
#' @import testit
#' @import grDevices
#' @import graphics
#' @import methods
#' @import stats
#' @import utils
####### @importFrom utils flush.console head read.csv read.table tail write.table


#End of Package Description
NA #NULL, ends description without hiding first function
######################################################################################
######################################################################################


# long and short DRC
#' Distance Requirement Cycle, long version 
#' 
#' Distance Requirement Cycle, long version: c(0.3,0.05, 0.001, 0.0005,0.0).
#' 
#' Used by \code{adDRC} (adjust DRC, depending on \code{cobra$Fres}) to set cobra$XI.
#' @seealso [cobraInit()]
#' @export
DRCL   <-c(0.3,0.05, 0.001, 0.0005,0.0)
#' Distance Requirement Cycle, short version 
#' 
#' Distance Requirement Cycle, short version: c(0.001,0.0).
#' 
#' Used by \code{adDRC} (adjust DRC, depending on \code{cobra$Fres}) to set cobra$XI.
#' @seealso [cobraInit()]
#' @export
DRCS   <-c(0.001,0.0)

#' Archiving Environment
#' 
#' intern.archive.env is an independent environment where every evaluated point and its evaluation by the real function
#' are stored in ARCHIVE and ARCHIVEY. This archive stores different values to \code{cobra$A} and \code{cobra$Fres} 
#' often during dubugging and visualisation cases where the real function is evaluated very often for debugging purposes. 
#'@export
intern.archive.env <- new.env()


#
#' Return a rescaled function
#' @param fn        function with argument \code{x} to be rescaled 
#' @param lower     a vector with lower bounds in original input space
#' @param upper     a vector with lower bounds in original input space
#' @param dimension   length of vector \code{lower} and \code{upper}
#' @param newlower  a number, the rescaled lower bound for all dimensions
#' @param newupper  a number, the rescaled upper bound for all dimensions
#' 
#' @return \code{newfn},      rescaled version of function \code{fn}
#' @seealso   \code{\link{forwardRescale}}, \code{\link{inverseRescale}}
#' @export
#' 
rescaleWrapper<-function(fn,lower,upper,dimension,newlower,newupper){
  oldfn<-fn                     
  newfn<-function(x){
    x<-sapply(1:length(x) , function(i){scales::rescale(x[i],to=c(lower[i],upper[i]),from=c(newlower,newupper))
    })
    y<-oldfn(x)
    return(y)
  }
  return(newfn)
}

#'Forward Rescaling
#'
#' Scale vector x in original space forward to rescaled space (usually \eqn{[-1,1]^d})
#'
#' @param x       a vector in the original input space 
#' @param cobra   list from \code{\link{cobraInit}}, we need here
#'   \describe{
#'      \item{originalL}{   a vector with lower bounds in original input space}
#'      \item{originalU}{   a vector with upper bounds in original input space}
#'      \item{newlower}{    a number, the rescaled lower bound for all dimensions}
#'      \item{newupper}{    a number, the rescaled upper bound for all dimensions}
#'   }
#' 
#' @return \code{z},      scaled version of vector x
#' @seealso   \code{\link{inverseRescale}}
#' @export
#' 
forwardRescale <- function(x,cobra) {
  if (!cobra$rescale) return(x)

  lb<-rep(cobra$newlower,length(x))
  up<-rep(cobra$newupper,length(x))
  origL <- cobra$originalL
  origU <- cobra$originalU
  z<-sapply(1:length(x) , function(i){scales::rescale(x[i],from=c(origL[i],origU[i]),to=c(lb[i],up[i]))})
  return(z)
}

#' Inverse Rescaling
#'
#' Scale vector x in rescaled space back to original space
#'
#' @param x       a vector in the rescaled input space (usually \eqn{[-1,1]^d})
#' @param cobra   list from \code{\link{cobraInit}}, we need here
#' \describe{
#'    \item{originalL}{   a vector with lower bounds in original input space}
#'    \item{originalU}{   a vector with upper bounds in original input space}
#'    \item{newlower}{    a number, the rescaled lower bound for all dimensions}
#'    \item{newupper}{    a number, the rescaled upper bound for all dimensions}
#' }
#'
#' @return \code{z},      inverse rescaling of vector x
#' @seealso   \code{\link{forwardRescale}}
#' @export
#' 
inverseRescale <- function(x,cobra) {
  if (!cobra$rescale) return(x)
  
  lb<-rep(cobra$newlower,length(x))
  up<-rep(cobra$newupper,length(x))
  z<-sapply(1:length(x) , function(i){scales::rescale(x[i],to=c(cobra$originalL[i],cobra$originalU[i]),from=c(lb[i],up[i]))})
  return(z)
}

#' Return best feasible solution in original space
#' 
#' @param cobra an object of class COBRA (see  \code{\link{cobraInit}})
#' 
#' @return the best feasible solution in original space
#' @seealso   \code{\link{getFbest}}
#' @export
#' 
getXbest <- function(cobra) {
  xbest <- cobra$xbest
  if (cobra$rescale) xbest <- inverseRescale(xbest,cobra)
  return(xbest)
}

#' Return best objective function value
#'
#' Return the original objective function value at the best feasible solution 
#'
#' Note: We cannot take the best function value via \code{cobra$fn}, because this 
#' may be modified by \code{plog()} or others )
#'
#' @param cobra an object of class COBRA (see  \code{\link{cobraInit}})
#' 
#' @return the original objective function value at the best feasible solution 
#' @seealso   \code{\link{getXbest}}
#' @export
#' 
getFbest <- function(cobra) {
  return(cobra$originalfn(getXbest(cobra))[1])
}

#
#
#' Random start algorithm
#'
#' RandomStart is now an R6 class (to have private storage for reproducible RNG my_rng2).
#' Its main member function is random_start, which returns a cobra object with cobra$xStart adjusted
#'
#' @param cobra an object of class COBRA (see  \code{\link{cobraInit}})
#' 
#' @return \code{cobra}, an object of class COBRA from \code{\link{cobraInit}}, with modified elements
#'      \item{\code{xStart}}{ new starting point, either random or \code{cobra$xBest} }
#'      \item{\code{noProgressCount}}{  number of consecutive iterations without a progress. 
#'                                      Set to 0, if a random start point is chosen. }
#' @keywords internal                                     
#' 
RandomStart<-R6::R6Class("RS",
                         public=list(
                           random_start=NA,
                           my_rng2=NA,
                           val=NA,
                           initialize=function(cobra){
                             # private$callRS(cobra)
                             self$val = cobra$cobraSeed
                             self$my_rng2=function(n, d, lower, upper) {
                               MOD = 10**5+7
                               OFS = 10**5-7
                               x = matrix(0,n,d)     # an (n x d) matrix with zeros
                               for (n_ in 1:n) {
                                 for (d_ in 1:d) {
                                   self$val <- (self$val*self$val*sqrt(self$val)+OFS) %% MOD    # avoid cycles (!)
                                   x[n_,d_] = self$val/MOD*(upper[d_]-lower[d_]) + lower[d_]    # map val to range [lower, upper[
                                 }
                               }
                               #browser()
                               return(x)
                             } # self$my_rng2
                             self$random_start = function(cobra){
                               testit::assert("cobra$sac$RS_rep is NULL!", !is.null(cobra$sac$RS_rep))
                               anewrand = ifelse(cobra$sac$RS_rep, 
                                                 self$my_rng2(1,1,c(0),c(1)), 
                                                 runif(1,min=0,max=1))
                               # p_restart<-(9*0.3/20)*tanh(-(nrow(cobra$A)-(cobra$initDesPoints+15)))+11*0.3/20
                               feasibleRate<-length(which(cobra$numViol==0))/length(cobra$numViol) #fraction of feasible points in the population
                               diff<-cobra$sac$RSmax-cobra$sac$RSmin 
                               
                               if((cobra$sac$RSAUTO==TRUE) && (feasibleRate <0.05)){
                                 p_const <- 0.4
                               }else{
                                 p_const <- (cobra$sac$RSmax+cobra$sac$RSmin)/2  #  i.e. = 0.175
                               }
                               
                               switch(cobra$sac$RStype,
                                      SIGMOID = p_restart <- (diff/2)*tanh(-(nrow(cobra$A)-(cobra$initDesPoints+15)))+(diff/2)+cobra$sac$RSmin,
                                      CONSTANT= p_restart <- p_const
                               )
                               # print(p_restart)
                               # The SIGMOID case:
                               #   - if the argument of tanh is negative and big, then p_restart -> RSmin
                               #   - if the argument of tanh is positive and big, then p_restart -> RSmax
                               # The argument of tanh is positive and big in the early iterations and negative and big in the late iterations
                               if(cobra$DEBUG_RS)
                                 cat(paste("randomness=",p_restart,"Feasibility Rate=",feasibleRate," \n"))
                               
                               if( (anewrand< p_restart)  || (cobra$noProgressCount >= cobra$sac$RS_Cs)){   # /WK/ ???
                                 #if( (anewrand< p_restart)  && (cobra$noProgressCount >= cobra$sac$RS_Cs)){   # /WK/ 
                                 verboseprint(cobra$verbose, important=FALSE,"Starting the internal optimizer with a random point in the space")
                                 #cat(sprintf("RS: iter=%03d, noProgressCount=%03d\n",nrow(cobra$A)+1,cobra$noProgressCount))
                                 if (cobra$sac$RS_rep==T) {
                                   xStart<-self$my_rng2(1,length(cobra$xbest),cobra$lower,cobra$upper)
                                   xStart<- as.vector(xStart)
                                 } else {
                                   xStart<-runif(n=length(cobra$xbest),min=cobra$lower,max=cobra$upper)
                                 }
                                 #print(sprintf("[random_start] (%9.5f,%9.5f) at iter %d", xStart[1], xStart[2], nrow(cobra$A)))
                                 cobra$RSDONE<-c(cobra$RSDONE,"RS")
                                 #cobra$noProgressCount<-0
                               } else{
                                 xStart<-cobra$xbest
                                 cobra$RSDONE<-c(cobra$RSDONE,"COBY")
                                 
                               }
                               cobra$xStart<-xStart
                               return(cobra)
                               
                             } # self$randomStart
                           }
                         ),  # public
                         private=list(
                           callRS=function(cobra){
                           } 
                         ) # private
                         
) # R6class
                                     

#Adjust DRC
#
adDRC<-function(maxF,minF){
  FRange<-(maxF-minF)
  if(FRange>1e+03) {
    cat(sprintf("FR=%g is large, XI is set to Short DRC \n",FRange))
    DRC<-DRCS
  }else{
    DRC<-DRCL
    cat(sprintf("FR=%g is not large, XI is set to Long DRC\n",FRange))
    
  }
  return(DRC)
}


detLen <-function(x){
  maxL<-max(x)
  minL<-min(x)
  return(maxL-minL)
}
#Adjust constraint functions
#
adCon<-function(cobra){
  
  fnold<-cobra$fn
  testit::assert("[adCon] cobra$Gres contains NA elements", any(is.na(cobra$Gres))==FALSE)
  
  GRL<-apply(cobra$Gres,2,detLen)
  if (min(GRL)==0) {  # pathological case where at least one constraint is constant:
    GR <- -Inf        # inhibit constraint normalization
  } else {
    GR<-max(GRL)/min(GRL)
  }


 if(GR > cobra$sac$TGR){
    cat(sprintf("Normalizing Constraint Functions \n"))
    GRfact<-c(1,GRL*(1/mean(GRL)))
    #finding the normalizing coefficient of the equality constraints
    if(length(cobra$equIndex)!=0){
      GRF<-GRfact[-1]
      EQUfact<-GRF[cobra$equIndex]
      #define a coefficient for the final equality margin
      equEpsFinal<-cobra$equHandle$equEpsFinal
      finMarginCoef<-min(c(equEpsFinal/EQUfact,equEpsFinal))/equEpsFinal
      cobra$equHandle$equEpsFinal<-finMarginCoef*equEpsFinal
      cobra$GRfact<-GRF
      cobra$finMarginCoef<-finMarginCoef
    }

    fn<-function(x){
      return(fnold(x)/GRfact)
    }
    cobra$fn<-fn

    Gres<-NULL
    for(i in 1:nrow(cobra$Gres)){
      Gres<-rbind(Gres,cobra$Gres[i,]/GRfact[-1])
    }
    # browser()
    
    cobra$Gres<-Gres
    
  }
  detLen2 <-function(x){      # never used
    maxL<-quantile(x,c(0.9))
    minL<-quantile(x,c(0.1))
    return((maxL-minL)/2)
  }
  cobra$Grange<-mean(GRL)
  equIndex<-which(names(GRL)=="equ")
  cobra$GrangeEqu<-mean(GRL[equIndex])
  return(cobra)
}

#Adjust Fitness Function
#
# Note that the fitness function cobra$fn is not changed by this function. Nor is cobra$Fres
# changed. Instead, all results found in cobra$Fres are either copied to cobra$SurrogateInput
# or they are transformed with transfFunc and the transformed results are copied to cobra$SurrogateInput. 
# cobra$SurrogateInput is then used to train a new surrogate model for the fitness 
# function and this surrogate model is used in the sequential optimizer.
#
adFit<-function(cobra,ind){
  maxF=max(cobra$Fres)
  minF=min(cobra$Fres)
  FRange<-(maxF-minF)
  
  transfFunc<-function(x,pShift=0){
    y<-plog(x,pShift=pShift)
    return(y)
  }
  
  #If cobra$onlinePLOG is true, the decision to do the plog transfomation or not 
  #is made in every iteration according to the p-effect. 
  # Otherwise the decision is made once according to the FRange value 
  if(cobra$sac$onlinePLOG){
    pShift<-0
    #decision making according to p-effect
    if(cobra$pEffect > 1){
      Fres<-sapply(cobra$Fres[ind],transfFunc,pShift=pShift)
      cobra$PLOG<-c(cobra$PLOG,TRUE)

    }else{
      Fres<-cobra$Fres[ind]
      pShift<-NA
      cobra$PLOG<-c(cobra$PLOG,FALSE)

      
    }
  }else{
    if(FRange>cobra$sac$TFRange) {   
      if(cobra$sac$adaptivePLOG){
        # browser()
        verbosecat(verbose=cobra$verbose,important=cobra$important,sprintf("Very Large FR=%g, applying adaptive plog()",FRange))
        #pShift<-mean(c(cobra$fbest,min(cobra$Fres)[1]))
        pShift<-mean(c(cobra$fbest,0))   # /SB/just for testing purposes
        #pShift<-cobra$fbest               # /WK/ the most simple alternative (don't know quality yet)
        
        #transfFunc<-function(x,pShift){
        #      x<-x-mean(c(pShift,0))
        #      #The following print is only for debugging purposes
        #      #print(paste("pShitf:",pShift))
        #      y<-plog(x,pShift=0)
        #     return(y)
        #    }
      }else{
        verbosecat(cobra$verbose,sprintf("Very Large FR=%g, applying plog()\n",FRange))
        pShift<- 0
      }
      Fres<-sapply(cobra$Fres[ind],transfFunc,pShift=pShift)
      cobra$PLOG<-c(cobra$PLOG,TRUE)
      
    }else{
      Fres<-cobra$Fres[ind]
      pShift<-NA
      cobra$PLOG<-c(cobra$PLOG,FALSE)
      
    } 
  }
  
  cobra$pShift<-c(cobra$pShift,pShift)
  cobra$SurrogateInput<-Fres
  
  #SB: 29.09.205
  #if(cobra$PLOG){
  #  cobra$XI<-DRCL  
  #}else{cobra$XI<-DRCS}
  
  return(cobra)
}

#' Monotonic transform
#' 
#' The function is introduced in [Regis 2014] and extended here by a parameter \eqn{p_{shift}}.
#' It is used to squash functions with a large range into a smaller range.\cr
#' Let \eqn{y' = (y-p_{shift})}: 
#'  \deqn{ plog(y) =  \ln(1+ y'), \quad\mbox{if}\quad y' \ge 0 } 
#'  \deqn{ plog(y) = -\ln(1- y'), \quad\mbox{if}\quad y'  <   0 } 
#'  
#' @param y       function argument
#' @param pShift  shift
#'  
#' @return \eqn{plog(y)}
#' @seealso \code{\link{plogReverse}}
#' @export
#'  
plog<-function(y,pShift=0.0){
  if(y-pShift>=0){
    ret<- log(1+y-pShift)
  }else{
    ret<- -log(1-(y-pShift))
  }
  return(ret)
}

#' Inverse of \code{\link{plog}}
#' 
#' @param y       function argument
#' @param pShift  shift
#' 
#' @return \eqn{plog^{-1}(y)}
#' @seealso \code{\link{plog}}
#' @export
plogReverse<-function(y,pShift=0){
  #y <- y/100
  if(y > 0){
    ret<-exp(y)-1+pShift
  }else{
    ret<-pShift+1-(1/exp(y)) 
  } 
  return(ret)
}

# changes cobra$err1, cobra$err2, cobra$pEffect
calcPEffect<-function(cobra,xNew,xNewEval){
  #browser()
  newPredY1 <- getPredY1(xNew,cobra$fitnessSurrogate1,cobra)
  newPredY2 <- getPredY1(xNew,cobra$fitnessSurrogate2,cobra)
  newErr1<-abs(newPredY1-xNewEval[1])
  newErr2<-abs(plogReverse(newPredY2)-xNewEval[1])
  #newErr2<-abs(newPredY2-xNewEval[1])
  cobra$err1<-c(cobra$err1,newErr1)
  cobra$err2<-c(cobra$err2,newErr2)
  err1<-cobra$err1
  err2<-cobra$err2
  #err1<-c(err1,newErr1)
  #err2<-c(err2,newErr2)
  errRatio<-err1/err2

  if(is.infinite(newErr2)){
    errRatio[length(errRatio)]<-0 
  }else if(is.infinite(newErr1)){
    errRatio[length(errRatio)]<-Inf
  }
  cobra$pEffect<-log10(median(errRatio,na.rm = TRUE))  
  return(cobra)
}
# --- alternative plog with sqrt ----
#
# plog<-function(y,pShift=0.0){
#   if(y-pShift>=0){
#     ret<- sqrt(1/4+y-pShift)-sqrt(1/4)
#   }else{
#     ret<- -sqrt(1/4-(y-pShift))+sqrt(1/4)
#   }
#   #return(100*ret)
#   return(ret)
# }
# 
# plogReverse<-function(y,pShift=0){
#   #y <- y/100
#   if(y > 0){
#     ret <- (y+1/2)^2-1/4+pShift
#   }else{
#     ret <- -(y-1/2)^2+1/4+pShift 
#   } 
#   return(ret)
# }
