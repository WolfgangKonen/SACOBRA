## RBF interpolation for d==1
nobs = 100        # number of observations for trainCubicRBF
ngrid= 100        # number of grid points per dimension for interpolation
d = 1
xobs = matrix(runif(d*nobs,min=-1,max=1),nobs,d)    # dim(xobs) = (nobs,d)
yobs = rowSums(xobs) * exp(-6*rowSums(xobs*xobs))

xgrid <- seq(-1,1,2/(ngrid-1))
xflat = as.matrix(xgrid)  # as.matrix important: only with dim(xflat)=(100,1) the call to predict.RBFinter will work
rbf.model = trainCubicRBF(xobs,yobs)
yflat = predict.RBFinter(rbf.model, xflat)

ytrue = xflat * exp(-6*xflat*xflat)
delta = ytrue - yflat
fprintf("%s %10.4e\n","avg abs delta[ , ] = ", mean(abs(delta)))

plot(xgrid,yflat, pch=20,cex=1, col="blue", main="One line (d=1)", xlab="xgrid", ylab="y")
lines(xgrid,ytrue,col="red")
lines(xgrid,delta,col="green")
