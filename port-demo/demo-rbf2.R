# Example demo_rbf2.py ported to R
# ["C:\user\datasets\Vorlesungen\Mathe-Divers\ProjectEuler.net\SACOBRA_Py\demo\demo_rbf2.py"]
#
# Surprisingly the R code for RBF interpolation is much slower (factor 60) than 
# Python's SciPy!! For cubic kernel and (ngrid,d) = (100,2):
#     - R:    263 msec
#     - Python: 4 msec
# If we look closer, then the RBF model training in R is negligible (1 msec), 
# but the prediction of the ngrid^2 = 100^2 points costs 99% of the time.
# 
# Computational accuracy: Both, Python and R, have in ygrid the same results  
# (7 digits after the decimal point) for nobs = 100 observations (could be
# necessary outcome of good interpolation, not numerical equivalence), but also
# for nobs = 10 observations (and here it proves numerical equivalence).

require(pracma)   # needed for fprintf

nobs =  100       # number of observations for trainCubicRBF
ngrid= 100 # 4 #  # number of grid points per dimension for interpolation
runs = 5          # number of runs {train,interpolate} for time measurements

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

fn <- function(x) {
  return (rowSums(x) * exp(-6*rowSums(x*x)))
}

demo_d_1 <- function(ngrid,nobs,runs) {
  cat("*** Starting demo_d_1 ... ***\n")
  d = 1
  xobs = my_rng(nobs,d,24)    # dim(xobs) = (nobs,d)
  yobs = fn(xobs) 
  fprintf("%s%s%s\n","the first 9 yobs: [",class(yobs),"]" )
  fprintf("%10.4e",yobs[1:9])
  cat("\n")
  
  X <- seq(-1,1,2/(ngrid-1))
  xgrid = X
  xflat = as.matrix(X)  # as.matrix important: only with dim(xflat)=c(100,1) the call to predict.RBFinter will work
  print(length(X))
  print(dim(xflat))     # dim(xflat) = (ngrid^d, d)
  ptm <- proc.time();
  for (r in 1:runs) {
    rbf.model = trainCubicRBF(xobs,yobs)
    yflat = predict.RBFinter(rbf.model, xflat)
  }
  cat("Proc time for RBF: ",1000*(proc.time()-ptm)[1]/runs," msec\n");
  
  ygrid = yflat

  ytrue = fn(xflat) 
  delta = ytrue - ygrid
  fprintf("%s %10.4e\n","avg abs delta[ , ] = ", mean(abs(delta)))
  
  plot(xgrid,ygrid, pch=20,cex=1, col="blue", main="One line (d=1)", xlab="xgrid", ylab="y")
  lines(xgrid,ytrue,col="red")
  lines(xgrid,delta,col="green")
  fprintf("%10.7f",ygrid[51:60])
  # R     :   0.0100961  0.0301357  0.0497379  0.0686172  0.0865063  0.1031759  0.1184069  0.1320160  0.1438724  0.1538627
  # Python: [ 0.0100961  0.0301357  0.0497379  0.0686172  0.0865063  0.1031759  0.1184069  0.132016   0.1438724  0.1538627]
  cat("\n")
}

demo_d_2 <- function(ngrid,nobs,runs) {
  cat("*** Starting demo_d_2 ... ***\n")
  d = 2
  xobs = my_rng(nobs,d,24)    # dim(xobs) = (nobs,d)
  yobs = fn(xobs) 
  fprintf("%s%s%s\n","the first 9 yobs: [",class(yobs),"]" )
  fprintf("%10.4e",yobs[1:9])
  cat("\n")
  
  X <- seq(-1,1,2/(ngrid-1))
  xgrid <- expand.grid(x=X,y=X) # xgrid is a data.frame
  xflat <- as.matrix(xgrid)     # as.matrix important here: with data.frame computation time nearly doubles (!)
  # ngr = ngrid-1                                   #
  # X=outer(rep(1,ngrid),(0:ngr)/ngr*2-1)           # this code does the same (for d=2), 
  # Y=t(X)                                          # but is less easy to understand
  # xflat = t(rbind(as.vector(X),as.vector(Y)))     #
  print(length(X))
  print(dim(xflat))     # dim(xflat) = (ngrid^d, d)

  ptm <- proc.time();
  for (r in 1:runs) {
    rbf.model = trainCubicRBF(xobs,yobs, squares=T)
    yflat = predict.RBFinter(rbf.model, xflat)
  }
  cat("Proc time for RBF: ",1000*(proc.time()-ptm)[1]/runs," msec\n");
  
  # contour plot in the d==2 case:
  ygrid = matrix(yflat,nrow=ngrid,ncol=ngrid)
  filled.contour(X,X, ygrid, nlevels=10
                , plot.axes = { axis(1); axis(2); points(xobs[,1],xobs[,2]) })

  ### in Grammar of Graphics, contour plot would be more complicated :
  # dfbase<<-cbind(xgrid,yflat)
  # names(dfbase)<-c("x","y","z")
  # p<-ggplot(data=dfbase,aes(x,y,z=z))
  # p<-p + stat_contour(size=1,aes(color=..level..))+guides(alpha=FALSE)
  # p<-p + scale_color_gradient(low = "red", high = "white")
  # p<-p + theme(panel.background=element_rect(fill="grey50"))
  # p<-p + theme(panel.grid=element_blank())
  # p<-p + theme(legend.text = element_text(size = 20))
  # p<-p + theme(legend.title = element_text(size = 20))
  # p<-p + theme(axis.text = element_text(size=20))
  # plot(p)
  
  ytrue = fn(xflat) 
  ytrue = matrix(ytrue,nrow=ngrid,ncol=ngrid)
  delta = ytrue - ygrid
  fprintf("%s %10.4e\n","avg abs |ytrue - ygrid| [ , ] = ", mean(abs(delta)))
  
  L = ngrid/2
  plot(X,ytrue[,L],col="red", type="l", main="One line (d=2, L=50)", xlab="xgrid", ylab="y")
  points(X,ygrid[,L], pch=20,cex=1, col="blue")
  lines(X,delta[,L],col="green")
  fprintf("%s%d%s %10.4e\n","avg abs |ytrue - ygrid| [,",L,"] = ", mean(abs(delta[,L])))
  fprintf("%10.7f",ygrid[51:60,L])
  # --- nobs=100 ---
  # R     :  -0.0097614  0.0088733  0.0274552  0.0457256  0.0634270  0.0803093  0.0961359  0.1106909  0.1237856  0.1352644
  # Python: [-0.0097614  0.0088733  0.0274552  0.0457256  0.063427   0.0803093  0.0961359  0.1106909  0.1237855  0.1352643]
  # --- nobs=10 ---
  # R     :  -0.0958555 -0.0950409 -0.0940627 -0.0929260 -0.0916367 -0.0902010 -0.0886258 -0.0869181 -0.0850855 -0.0831358
  # Python: [-0.0958555 -0.0950409 -0.0940627 -0.092926  -0.0916367 -0.090201  -0.0886258 -0.0869181 -0.0850855 -0.0831358]
  cat("\n")
}

# print(my_rng(2,3,24))

# demo_d_1(ngrid,nobs,runs)
demo_d_2(ngrid,nobs,runs)
