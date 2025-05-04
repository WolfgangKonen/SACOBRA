# Demonstration of approximating RBFs on Rastrigin function
# [requires functions multi_gfnc and print_and_plot from SACOBRA\inst\examples\ex_COP.R]
#
require(pracma)   # needed for fprintf

# create reproducible 'random' numbers on the Python and R side
# CAUTION: cycles might occur for n>100 (!!)
my_rng <- function(n, d, seed) {
  MOD = 10**5+7
  val = seed
  x = matrix(0,n,d)     # an (n x d) matrix with zeros
  for (n_ in 1:n) {
    for (d_ in 1:d) {
      val = (val*val) %% MOD
      x[n_,d_] = 2*val/MOD - 1    # map val to range [-1,1[
      
    }
  }
  return(x)
}


# z can be a vector of points (1D values)
rastrigin_1D <- function(z, A=10) {
  return (A*(d-(cos(2*pi*z)))+(z^2))
}

# Rastrigin function for arbitrary dimensions d.
# z is a *single* point (a vector of dimension d)
rastrigin <- function(z) {
  if (!is.numeric(z)) {
    browser()
  }
  A = 10
  return (A*(length(z)-sum(cos(2*pi*z)))+sum(z^2))
}

# Vectorized Rastrigin function for arbitrary dimensions d. It returns a vector
# of results for all data points stored in x. 
# x is a (n * d)-matrix where each row carries a data point (a vector of dimension d)
rastrigin_vec <- function(x) {
  if (!is.matrix(x)) {
    browser()
  }
  A = 10
  return (A*(ncol(x)-rowSums(cos(2*pi*x)))+rowSums(x*x))
}



demo_approx_d1 <- function(cobraSeed=42, DO_PLOT=F) {
  # The 1-dimensional case
  cat("*** Starting demo_approx_d1 ... ***\n")
  
  d = 1
  xarr = seq(-5,5,0.1)
  #yarr = rastrigin_1D(xarr)
  yarr = sapply(1:length(xarr), function(i) {rastrigin(xarr[i])})
  x2 = seq(-5,4.5,0.66)+0.25
  #y2 = rastrigin_1D(x2)
  y2 = sapply(1:length(x2), function(i) {rastrigin(x2[i])})
  
  ptail = T
  squares = T
  rho = 0.19
  # as.matrix important: only with dim(xflat)=(N,1) the call to trainCubicRBF and later predict.RBFinter will work
  rbf.model = trainCubicRBF(as.matrix(x2), as.matrix(y2), ptail, squares, rho)
  #
  xflat = as.matrix(xarr) 
  ymod = predict.RBFinter(rbf.model, xflat)
  if (DO_PLOT) {
    plot(xarr,yarr,type="l")
    points(x2,y2,col="red",pch=20,cex=2)
    points(xarr,ymod,type="l",col="green",lwd=2.5)
  }

  #xStart = c(4.8)
  cobra <- cobraInit(xStart=NA, fName=sprintf("rastrigin %d",d),
                     fn=rastrigin, lower=c(-5), upper=c(5), feval=195, 
                     initDesign = "LHS", initDesPoints = 5,
                     solu = c(0), cobraSeed = cobraSeed,
                     RBFrho = 0, RBFrhoDec = 2.0)
  print(cobra$Fres)
  #fprintf("%12.8f\n",cobra$Fres) # just to get more digits after decimal point
  
  cobra <- cobraPhaseII(cobra)
  
  cobra <- print_and_plot(cobra,c(1e-8,1e+1))

  return (cobra)
}

demo_approx_d2 <- function(cobraSeed, ngrid=100, nobs=200, runs=1, DO_PLOT=F) {
  cat("*** Starting demo_approx_d2 ... ***\n")
  d = 2
  rho = 2.5 # 1.75
  xobs = 5 * my_rng(nobs,d,44)    # dim(xobs) = (nobs,d)
  yobs = rastrigin_vec(xobs) 
  fprintf("%s%s%s\n","the first 9 yobs: [",class(yobs),"]" )
  fprintf("%10.4e",yobs[1:9])
  cat("\n")
  
  X <- seq(-5,5,10/(ngrid-1))
  xgrid <- expand.grid(x=X,y=X) # xgrid is a data.frame
  xflat <- as.matrix(xgrid)     # as.matrix important here: with data.frame computation time nearly doubles (!)
  print(length(X))
  print(dim(xflat))     # dim(xflat) = (ngrid^d, d)
  
  ptm <- proc.time();
  for (r in 1:runs) {
    rbf.model = trainCubicRBF(xobs,yobs, squares=T, rho=rho)
    yflat = predict.RBFinter(rbf.model, xflat)
  }
  cat("Proc time for RBF: ",1000*(proc.time()-ptm)[1]/runs," msec\n");
  
  # contour plot in the d==2 case:
  ygrid = matrix(yflat,nrow=ngrid,ncol=ngrid)
  ytrue = rastrigin_vec(xflat) 
  ytrue = matrix(ytrue,nrow=ngrid,ncol=ngrid)
  if (DO_PLOT) {
    filled.contour(X,X, ygrid, nlevels=10, main="ygrid"
                   , plot.axes = { axis(1); axis(2); points(xobs[,1],xobs[,2]) })
    
    filled.contour(X,X, ytrue, nlevels=10, main="ytrue"
                   , plot.axes = { axis(1); axis(2); })
  }
  delta = ytrue - ygrid
  fprintf("%s %10.4e\n","avg abs |ytrue - ygrid| [ , ] = ", mean(abs(delta)))
  
  L = ngrid/2
  par(mfcol=c(1,1))
  par(mar=c(5, 4, 4, 2) + 0.1) 
  if (DO_PLOT) {
    plot(X,ytrue[,L],col="red", type="l", main="One line (d=2, L=50)", xlab="xgrid", ylab="y")
    points(X,ygrid[,L], pch=20,cex=1, col="blue")
    lines(X,delta[,L],col="green")
  }
  fprintf("%s%d%s %10.4e\n","avg abs |ytrue - ygrid| [,",L,"] = ", mean(abs(delta[,L])))
  fprintf("%10.7f",ygrid[51:60,L])
  cat("\n")

  cobra <- cobraInit(xStart=NA, fName=sprintf("rastrigin %d",d),
                     fn=rastrigin, lower=c(-5,-5), upper=c(5,5), feval=395, 
                     initDesign = "LHS", initDesPoints = 10,
                     solu = c(0,0), cobraSeed = cobraSeed,
                     RBFrho = rho, RBFrhoDec = 2.0)
  cobra <- cobraPhaseII(cobra)
  
  cobra <- print_and_plot(cobra,c(1e-8,1e+1))
  
  return (cobra)
}

demo_approx_d <- function(cobraSeed, d=2) {
  cat(sprintf("*** Starting demo_approx_d for d = %d ... ***\n",d))
  rho = 2.5 # 1.75

  cobra <- cobraInit(xStart=NA, fName=sprintf("rastrigin %d",d),
                     fn=rastrigin, lower=rep(-5,d), upper=rep(5,d), feval=395, 
                     initDesign = "LHS", initDesPoints = 100,
                     solu = rep(0,d), cobraSeed = cobraSeed,
                     RBFrho = rho, RBFrhoDec = 2.0)
  cobra <- cobraPhaseII(cobra)
  
  cobra <- print_and_plot(cobra,c(1e-8,1e+1))
  
  return (cobra)
}

trust_approx_d <- function(cobraSeed, d=2) {
  cat(sprintf("*** Starting demo_approx_d for d = %d ... ***\n",d))
  rho = 2.5 # 1.75
  
  cobra <- cobraInit(xStart=NA, fName=sprintf("rastrigin %d",d), fn=rastrigin, 
                     lower=rep(-5,d), 
                     upper=rep(+5,d), feval=200, 
                     initDesign = "LHS", initDesPoints = 100,
                     solu = rep(0,d), cobraSeed = cobraSeed,
                     RBFrho = rho, RBFrhoDec = 2.0)
  cobra <- cobraPhaseII(cobra)
  
  xb = getXbest(cobra)
  browser()
  cobra <- cobraInit(xStart=NA, fName=sprintf("rastrigin %d",d), fn=rastrigin, 
                     lower=xb+rep(-1.5,d), 
                     upper=xb+rep(+1.5,d), feval=300, 
                     initDesign = "LHS", initDesPoints = 100,
                     solu = rep(0,d), cobraSeed = 2*cobraSeed,
                     RBFrho = rho, RBFrhoDec = 2.0)
  cobra <- cobraPhaseII(cobra)
  
  cobra <- print_and_plot(cobra,c(1e-8,1e+1))
  
  return (cobra)
}

# cobra = demo_approx_d1(DO_PLOT=T)
# cobra = multi_gfnc(demo_approx_d1, 10, 390)
# ngrid= 100        # number of grid points per dimension for contour plot (d=2)
# nobs =  200       # number of observations for trainCubicRBF (d=2)
# cobra = demo_approx_d2(44,ngrid,nobs)
# cobra = multi_gfnc(demo_approx_d2, 10, 390)
# cobra = demo_approx_d(44, d=3)
#cobra = multi_gfnc(function(seed) {demo_approx_d(seed,d=10)}, 10, 390)
cobra = trust_approx_d(392, d=3)
# cobra = multi_gfnc(function(seed) {trust_approx_d(seed,d=3)}, 10, 390)