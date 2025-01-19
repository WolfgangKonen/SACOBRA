
N=5
hesse = matrix(runif(N*N),nrow=N,ncol=N)
hesse = 0.5*(hesse+t(hesse))
H = hesse


#calculating hesse^(1/2) with SVD
SVD<-svd(hesse)
eps<-1e-25
invD <- 1/sqrt(SVD$d)
invD[abs(SVD$d/SVD$d[1])<eps] <- 0
SV=SVD$u
SW=diag(SVD$d)

#calculating with eigen
EV=eigen(H)$vectors
w=(eigen(H)$values)
W=diag(w)
invDe <- 1/sqrt(abs(w))
invDe[abs(w)/abs(w[1])<eps] <- 0
#invDe[w<0] <- invDe[w<0]*(-1)
ME <- diag(invDe) %*% t(EV)



### your line:
#M2<-SVD$v%*%diag(invD)%*%t(SVD$v)

# my line
M<-diag(invD)%*%t(SVD$v)

#print(M2%*%hesse%*%t(M2))

print(M%*%hesse%*%t(M))

print(ME%*%hesse%*%t(ME))

# 
cat("\nAre the eigenvalue equations fulfilled for a) SVD, b) eigen?\n")
print(H%*%SV - SV%*%SW)    # usually not 0-matrix for SVD
print(H%*%EV - EV%*%W)     # but always for eigen


