#' Default settings for online whitening functionality
#' 
#' Sets default values for the online whitening functionality in order to handle function with high conditioning.  
#' With the call \code{\link{setOpts}(cobra$CA,defaultCA())} it is possible to extend a partial list 
#' \code{cobra$CA} to a list containing all \code{CA}-elements (the missing ones are taken from 
#' \code{defaultCA()}).\cr \cr
#' As RBF interpolations face severe difficulties to deliver reasonable models for functions with high conditioning,
#' we try to transform the function with high conditioning \eqn{f(\vec{x})} to a better conditioned one \eqn{g(\vec{x})} which is easier to model. 
#' \deqn{ g(\vec{x})=f(\mathbf{M}(\vec{x}-\vec{x}_c)) }
#' A possible transformation matrix \strong{M} is the squared inverse of the Hessian matrix \eqn{\mathbf{H}^{-0.5}}, assuming that \strong{M}
#' is chosen with the following assumption:
#' \deqn{\frac{\partial^2 g(\vec{x})}{\partial \vec{x}^2}=\mathbf{I}}
#' 
#Details:
#' The current version is only relevant for unconstrained problems.
#' It this stage it is not recommended to apply the online whitening to expensive optimization problems
#' As it demands large number of function evaluations. Every online whitening call demands \eqn{4d^2+4d} function evaluations.
#' 
#' @return \code{CA}, a list of the follwing elements:
#'     \item{active}{Set to TRUE if an online whitening of the fitness function is desired}
#'     \item{HessianType}{["real"] You can choose if the Hessian amtrix is evaluted on the real function or on the surrogate model ["real", "surrogate"]. 
#'  Please note that the determination of Hessian matrix on the real function at each point costs 4*d^2+d real function evalautions}
#'     \item{ITERS}{[seq(10,500,10)], pass a vector of integers to this paramter then the Hessian matrix will be updated only in the given iterations, 
#'     we recommended applying the online-whitening each 10 iterations after the \code{10*d} initial iterations. 
#'     \cr \code{seq(10*d,maxIter,10)}, where \code{d} is the dimensionality of the optimization problem . 
#'      If set as the charachter "all" then the Hessian matrix will be updated in each iteration and whitening procedure will be repeated}
#'     \item{alpha}{[1] you can assign any real value to this parameter. Only values between 0 to 2 are suggested.
#'  This value is used in order to modify the transformation center \code{tCenter} as follows: xbest+alpha*(grad), and grad is the direction of the last improvment}
#' @seealso   \code{\link{setOpts}}
#' @export
defaultCA<-function(){
  CA=list(active=FALSE,
           HessianType="real",
           ITERS=seq(10,500,10), #"all"
           alpha=1 # the coefficient of the gradient of function f in order to determine the transformation center
  )
  return(CA)  
}