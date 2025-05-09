#
# repairChootinan.R contains the algorithm of [Chootinan2006]
#       Chootinan & Chen: "Constraint handling in genetic algorithms using a 
#       gradient-based repair method", Computers & Operations Research 33 (2006), p. 2263

######################################################################################
# repairChootinan
#
#' Repair an infeasible solution with the method of Chootinan.
#'
#' Implements the method of [Choo2006] Chootinan & Chen "Constraint handling in genetic algorithms
#' using a gradient-based repair method", Computers & Operations Research 33 (2006), p. 2263.
#'
#' @param x            an infeasible solution vector \eqn{\vec{x}} of dimension \code{d}
####  an infeasible solution (vector of length \code{dimension})
#' @param gReal        a vector \eqn{(g_1(\vec{x}),\ldots,g_m(\vec{x}),
#'                                    h_1(\vec{x}),\ldots,h_r(\vec{x}))}  
#'                     holding the real constraint values at \eqn{\vec{x}}
####  a vector of length nconstraint holding the real constraint values at \code{x}
#' @param rbf.model    the constraint surrogate models 
#' @param cobra        parameter list, we need here 
#'     \describe{
#'       \item{\code{lower}}{   lower bounds of search region}
#'       \item{\code{upper}}{   upper bounds of search region}
#'       \item{\code{ri}}{      a list with all parameters for \code{repairChootinan}  }
#'       \item{\code{trueFuncForSurrogate}}{ if TRUE (only for diagnostics), use the true constraint
#'                     functions instead of the constraint surrogate models 
#'                     \code{rbf.model} }
#'       \item{\code{fn}}{   true functions, only needed in case of 
#'                           \code{trueFuncForSurrogate==TRUE} }
#'     }
#' @param checkIt      [FALSE] if TRUE, perform a check whether the returned solution is really 
#'                      feasible. Needs access to the true constraint function \code{conFunc}
#' @return \code{z},  a vector of dimension \code{d} with a repaired (hopefully feasible) solution
#' @seealso \code{\link{repairInfeasRI2}}, \code{\link{cobraPhaseII}}
#'
#' @author Wolfgang Konen, Cologne University of Applied Sciences
#' @export
######################################################################################
repairChootinan <- function(x,gReal,rbf.model,cobra,checkIt=FALSE) 
{  
  testit::assert("cobra$ri is NULL (not defined)", !is.null(cobra$ri) );
  #testit::assert("cobra$ri$q is NULL (not defined)", !is.null(cobra$ri$q) );
  verbosecat(cobra$verbose,important=FALSE,"Starting Chootinan's repair ...\n")
  
  conFunc <- {function(x)cobra$fn(x)[-1];}
  trueFunc = cobra$trueFuncForSurrogate
  ri = cobra$ri
  lowerP = cobra$lower
  upperP = cobra$upper
  gradEps = cobra$ri$gradEps
  if (is.null(gradEps)) gradEps = min(upperP-lowerP)/1000; 
  # gradEps = stepsize for numerical gradient calculation
  # e.g. 0.001 if the smallest length of search cube is 1.0
  currentEps <- cobra$currentEps[length(cobra$currentEps)] # only needed if cobra$equHandle$active
  equ2Index <- c(cobra$equIndex,cobra$nConstraints+(1: length(cobra$equIndex)))
  
  #########################################################################
  # helper functions
  #
  
  # checkSingleConstraints (for debug only):
  # Does the single correction Del[k,] calculated for kth violated constraint bring
  # the solution close to the kth constraint boundary, as it should?
  # If so, abs(fSingle[k]) should be much smaller than fBefore[k]
  checkSingleConstraints <- function(Del,Grad,Viol) {
    cat("\n")
    cat("   Length of Del vectors:\n")
    print(sqrt(rowSums(Del*Del)));
    cat("   Length of Grad vectors:\n")
    print(sqrt(rowSums(Grad*Grad)));
    fBefore = NULL
    fSingle = NULL
    for (k in 1:length(Viol)) {
      v = Viol[k]
      fBefore = c(fBefore,interpRBF(x,rbf.model)[v])
      fSingle = c(fSingle,interpRBF(x+Del[k,],rbf.model)[v])
    }
    cat("   single violations before/after single correction: \n"); 
    print(rbind(fBefore,fSingle));
    cat("   quotient q=fBefore/fSingle (|q| should be larger than approx. 8.0): \n")
    print(fBefore/fSingle);     
  }
  
  gradV <- function(x,rbf.model) {
    dimension <- length(x)
    nconstraint <- ncol(rbf.model$coef)
    
    nd <- data.frame(outer(rep(1,2*dimension+1),x))
    for (i in 1:dimension) {
      nd[2*i  ,i] <- nd[2*i  ,i] - gradEps
      nd[2*i+1,i] <- nd[2*i+1,i] + gradEps
    }
    nd=as.matrix(nd)
    
    # f is at first a (2*dimension+1 x nconstraint) matrix containing constraint surrogate 
    # responses at x (row 1) and small '-' and '+' deviations from x for any dimension  
    # in rows 2,...,2*dimension+1
    if (!trueFunc) {
      f <- as.matrix(predict(rbf.model,newdata=nd))
    } else {
      f <- t(apply(nd,1,conFunc))     
    }

    if (cobra$equHandle$active){
      gReal<-c(gReal,-gReal[cobra$equIndex])
      f    <-cbind(f,-f[,cobra$equIndex])
      gReal[equ2Index] <- gReal[equ2Index] - currentEps
      f[,equ2Index]    <- f[,equ2Index] - currentEps
      # now f is a (2*dimension+1 x (nconstraint+nequ)) matrix ( nequ = # equality constraints )
    }
    
    testit::assert("Wrong types for f or gReal", is.matrix(f), is.vector(gReal) );
    testit::assert("Columns do not match in f and gReal", ncol(f) == length(gReal));
    testit::assert("No constraint is eps1-infeasible",any(gReal+ri$eps1>0));

    mat <- NULL
    for (i in 1:dimension) {
      gradf <- (f[2*i+1,]-f[1,])/gradEps
      mat = cbind(mat,gradf)
    }    
    # the kth row of 'mat' contains the gradient for the kth constraint.

    # add the gradients for the 'dimension' lower-bound constraints -x[i]+lowerP[i] < 0
    mat <- rbind(mat,diag(rep(-1,dimension)))
    # add the gradients for the 'dimension' uppber-bound constraints x[i]-upperP[i] < 0
    mat <- rbind(mat,diag(rep(+1,dimension)))
    
    return(mat)
  }
  
  # dVFunc: the constraint violations
  # (>0 means violated for inequality constraints)
  dVFunc <- function(x,rbf.model,lowerP,upperP) {
    dV <- pmax(0,interpRBF(x,rbf.model));
    # add the 'dimension' lower-bound constraints -x[i]+lowerP[i] < 0
    dV <- c(dV,pmax(0,-x+lowerP));
    # add the 'dimension' uppber-bound constraints x[i]-upperP[i] < 0
    dV <- c(dV,pmax(0,+x-upperP));
  }
  
  #########################################################################
  # start of repairInfeasible
  #
  dimension <- length(x)
  nconstraint <- ncol(rbf.model$coef)
  rownames(gReal) <- NULL
  isEqualityConstr = rep(FALSE,nconstraint+2*dimension)
  eps=1e-6    # break if equality constraints are fulfilled up to an eps-margin
  del=1e-5    # break if (Euclidean) shift in input space is smaller than del
  
  testit::assert("Wrong types for gReal", is.vector(gReal) );
  testit::assert("Columns do not match in f and gReal", length(gReal) == nconstraint)
  
  xOld = x
  maxCount=1000; count=0
  while(1) {
    count=count+1
    mat = gradV(x,rbf.model)
    deltaV = dVFunc(x,rbf.model,lowerP,upperP)
    
    Viol = which((!isEqualityConstr & deltaV>0) | 
                 ( isEqualityConstr & abs(deltaV) > eps))  # index to the violated constraints
    if (length(Viol)==0) {
      cat("All constraints satisfied\n")
      break;
    }
    verbosecat(cobra$verbose,important=FALSE,paste("Viol =",Viol,"\n"))
    mat = mat[Viol,]             # consider only the violated constraints
    deltaV=deltaV[Viol]
    if (length(Viol)==1) 
      mat = t(as.matrix(mat)) #     
    
    
    s=svd(mat)
    invD=1/s$d
    invD[invD==Inf] = 0 
    if (length(Viol)==1) diagInvD = invD else diagInvD = diag(invD)
    inv=s$v %*% diagInvD %*% t(s$u)   # SVD-inverse (or Moore-Penrose pseudoinverse) of mat
    
    if (nrow(mat)<ncol(mat)) {
      testit::assert("svd inverse not OK", abs(mat %*% inv - diag(rep(1,nrow(mat)))) < 1e-7)
    } else {
      testit::assert("svd inverse not OK", abs(inv %*% mat - diag(rep(1,ncol(mat)))) < 1e-7)
    }
    
    
    # just for debug
    DEBUG=FALSE
    if (DEBUG) {
      Del <- NULL
      for (k in 1:nrow(mat)) {
        gradf = mat[k,]
        g2 <- sum(gradf*gradf)  
        Del_k =  -gradf * (deltaV[k]+ri$eps1)/g2  
        Del = rbind(Del,Del_k)   
      }
      checkSingleConstraints(Del,mat,Viol);      
    }
    
    dx = - inv %*% deltaV       # the suggested shift in input space
    dx = as.vector(dx)
    lenDx = sqrt(sum(dx*dx))    # its Euclidean length
    absDx = abs(dx)
    x = x + dx
    #cat("x =",x,"\n")
    verbosecat(cobra$verbose,important=FALSE,
               sprintf("max(|dx|) , sum(|deltaV|) = %9.2e , %9.2e\n", max(absDx), sum(abs(deltaV))))
    if (all(absDx < del)) {
      verbosecat(cobra$verbose,important=FALSE,paste("Stop because all |dx| < ",del,"\n"))
      Tdx = x-xOld
      lenTdx = sqrt(sum(Tdx*Tdx))     # total shift
      verbosecat(cobra$verbose,important=FALSE,sprintf("total shift in input space: %9.2e\n",lenTdx))
      break;
    }
    if (count>=maxCount) {
      verbosecat(cobra$verbose,important=FALSE,
                 paste("Stop because maxCount = ",maxCount," reached\n"))
      Tdx = x-xOld
      lenTdx = sqrt(sum(Tdx*Tdx))     # total shift
      verbosecat(cobra$verbose,important=FALSE,sprintf("total shift in input space: %9.2e\n",lenTdx))
      break;
    }
    
    # This should normally not happen:
    testit::assert("New solution x contains NA or NaN!", !any(is.na(x)))
     
  } # while
  
  if (checkIt) {
    fRbf = interpRBF(x,rbf.model)       # fRbf:  constraint surrogate values after repair
    fTrue = conFunc(x)                  # fTrue: true constraint values after repair
    if (cobra$equHandle$active) {
      fRbf<-c(fRbf,-fRbf[cobra$equIndex])
      fRbf[equ2Index] <- fRbf[equ2Index] - currentEps
      fTrue<-c(fTrue,-fTrue[cobra$equIndex])
      fTrue[equ2Index] <- fTrue[equ2Index] - currentEps
    }
    #print(gReal); print(fRbf); print(fTrue)
    violatedConstraints = which(fTrue>0)
    cfcReal <- maxReal <- 0;
    if (any(fTrue>0)) {
      cfcReal = sum(fTrue[violatedConstraints]);
      maxReal = max(fTrue[violatedConstraints]);
    }
    #checkSingleConstraints(Del,Grad,Viol);
    cat("Repaired solution is feasible: ",cfcReal<=0,", cfcReal=",cfcReal,", maxViol=",maxReal, "\n")
    #if (cfcReal>0) {
      cat("   eps1-inf constraints before repair: ",paste(which(gReal+ri$eps1>0) ,collapse=" "),"\n")
      cat("   violated constraints before repair: ",paste(which(gReal>0) ,collapse=" "),"\n")
      cat("   violated constraints  after repair: ",paste(which(fTrue>0),collapse=" "),"\n")
      cat("   violated c-surrogats  after repair: ",paste(which(fRbf>0),collapse=" "),"\n")
    #}
    
  }
  
  #browser()
  
  return (x)
}

