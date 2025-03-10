## RBF interpolation for d==2
nobs = 100        # number of observations for trainCubicRBF
ngrid= 100        # number of grid points per dimension for interpolation
d = 2
set.seed(42)
xobs = matrix(runif(d*nobs,min=-1,max=1),nobs,d)    # dim(xobs) = (nobs,d)
yobs = rowSums(xobs) * exp(-6*rowSums(xobs*xobs))

X <- seq(-1,1,2/(ngrid-1))
xgrid <- expand.grid(x=X,y=X) # xgrid is a data.frame
xflat <- as.matrix(xgrid)     # as.matrix important here: with data.frame computation time nearly doubles (!)
rbf.model = trainCubicRBF(xobs,yobs)
yflat = predict.RBFinter(rbf.model, xflat)

ygrid = matrix(yflat,nrow=ngrid,ncol=ngrid)
filled.contour(X,X, ygrid, nlevels=10
               , plot.axes = { axis(1); axis(2); points(xobs[,1],xobs[,2]) })

ytrue = rowSums(xflat) * exp(-6*rowSums(xflat*xflat))
ytrue = matrix(ytrue,nrow=ngrid,ncol=ngrid)
delta = ytrue - ygrid
fprintf("%s %10.4e\n","avg abs delta[ , ] = ", mean(abs(delta)))

L = ngrid/2
plot(X,ytrue[,L],col="red", type="l", main="One line (d=2, L=50)", xlab="xgrid", ylab="y")
points(X,ygrid[,L], pch=20,cex=1, col="blue")
lines(X,delta[,L],col="green")
