check1 <- function(cobra) {
  #
  # Check whether we can recompute numViol (number of violations of artificial feasibility) from Gres:
  # numViol should be equal to the row sum rs1 (or rs2) and it is.
  #
  temp=cobra$Gres
  npts=nrow(temp)
  equInd = cobra$equIndex 
  temp[,equInd] = abs(temp[,equInd])
  # rs1 holds for each row 1:npts the number of constraints that are violated 
  # (in their *artificial* feasibility)
  rs1 = sapply(1:npts , function(i){ 
    temp[i,equInd] <<- temp[i,equInd] - cobra$currentEps[i] 
    # IMPORTANT: use "<<-" to change temp on the check1 level (!)
    return(sum(temp[i,]  > cobra$conTol))   
  })
  browser()
  cond1 = temp > cobra$conTol
  cond2=ifelse(temp > cobra$conTol,1,0)
  rs2 = rowSums(cond2)
  testit::assert(length(rs1)==length(cobra$df$nViolations))
  testit::assert(all(rs1==cobra$df$nViolations))
  testit::assert(all(rs2==cobra$df$nViolations))
  # data frame nv holds: rs1=row sums of cond1, rs2=..of cond2, nv=number of (constr.) violations, currentEps=mu for equ constr.
  #                      The assertion rs1 == rs2 == nv should be true.
  nv=data.frame(rs1=rs1, rs2=rs2, nv=cobra$df$nViolations, currentEps=cobra$currentEps)
  # View(nv)
  return (list(cond1=cond1,cond2=cond2,rs1=rs1))
}

check2 <- function(cobra,cond1,cond2) {
  #
  # Matrix cond1 (or cond2) calculated in check1 contains the info, which point (row)
  # has which constraint (column) violated.
  # Take in cs1 the colSums of cond1 for the last half of points: Which constraints
  # are often and which constraints are never violated?
  #
  npts=nrow(cond1)
  ncols=ncol(cobra$Gres)
  equIndex=which(colnames(cobra$Gres)=="equ")
  nonEqu=which(colnames(cobra$Gres)!="equ")
  cs1 = colSums(cond1[(npts/2):npts,])
  cs2 = colSums(cond2[(npts/2):npts,])
  testit::assert(all(cs1==cs2))
  print(cs1)
  print(which(cs1>2))
  quant=sapply(1:ncols,function(i) {quantile(cobra$Gres[(npts/2):npts,i],0.75)})
  plot(cs1,quant,col="blue")
  points(cs1[equIndex],quant[equIndex],col="red")
  # Inequality constraints that are often violated have 3rd quantile nv$quant near 0.0 (probably active).
  # Inequality constraints that are seldom violated have in most cases nv$quant near -0.75 (probabliy inactive).
  nv=data.frame(quant=quant,cs1=cs1, c_solu=cobra$fn(cobra$solu)[-1])
  viol=which(nv$cs1[nonEqu]>20)
  print(nv$c_solu[viol])
  # the fact that the printed numbers are all (with one exception) very small in magnitude 
  # shows that the inequality constraints often violated in our run belong to those 
  # which are *active* at the true solution. 
  # View(nv)
}  

check3 <- function(cobra,cond1,rs1) {
  # 
  # Is the pattern (of often/seldom violated constraints) different if we concentrate 
  # on the infill points with less then 6 violations?  -  Answer: No.
  #
  npts=nrow(cond1)
  cond1_1 = cond1[(npts/2):npts,]
  rs1_1 = rs1[(npts/2):npts]
  ind6 = which(rs1_1 < 6)
  cs1_1 = colSums(cond1_1[ind6,])
  #print(cs1_1)
  #print(which(cs1_1>2))
}
  
check4 <- function(cobra,rs1,ind4) {
  # 
  # Check the trueA calculation (trueA: true constraint value after refine step): 
  # Take the given index ind4 for step 4 and check whether the 
  # number of constraints where trueA > conTol coincides with rs1[ind4]
  #
  equInd = cobra$equIndex 
  trueA=cobra$fn(cobra$A[ind4,])[-1]
  trueA[equInd] = abs(trueA[equInd]) - cobra$currentEps[ind4]
  pos_trueA = which(trueA > cobra$conTol)
  testit::assert(length(pos_trueA)==rs1[ind4])
}

check5 <- function(cobra) {
  # 
  # Can we recompute fbest=cobra$fbestArray from the elements of cobra? 
  # Answer: Yes.
  # (Note that cobra$df$Best created before /2025/07/04 has NAs before first artificial feasible point (see cobraPhaseII.R:867)
  # and is therefore not the same.)
  #
  equInd = cobra$equIndex 
  calc_fbest <- function(i) {
    temp=cobra$Gres
    temp[,equInd] <- abs(temp[,equInd]) - cobra$currentEps[i]
    
    if (i <= cobra$initDesPoints) {
      currentMaxViols<-sapply(1:cobra$initDesPoints, function(k){
        y<-max(0,temp[k,])
        return(y)
      }) 
      numViol<-sapply(1:cobra$initDesPoints , function(i){ # number of initial Constraint violations
        return(sum(temp[i,] > cobra$conTol))   
      })
      if(0 %in% numViol){
        fbest<-min(cobra$Fres[which(numViol==0)])
      }else{
        # if there is no feasible point yet, take one from the points with min number of violated constraints:
        minNumIndex<-which(numViol==min(numViol))
        # select among the min-numViol points the one with smallest Fres:
        FresMin <- cobra$Fres[minNumIndex]
        index <- minNumIndex[which.min(FresMin)]
        fbest <- cobra$Fres[index[1]]
      }
      
    } else {   # i.e. for i > cobra$initDesPoints:
      
      currentMaxViols<-sapply(1:i, function(k){
        y<-max(0,temp[k,])
        return(y)
      }) 
      currentFeas<-which(currentMaxViols <= cobra$conTol)
      if (length(currentFeas)==0) {
        ibest<-which(currentMaxViols==min(currentMaxViols))[1]
        fbest<-cobra$Fres[ibest]
      } else {
        fminInd<-which(cobra$Fres[currentFeas]==min(cobra$Fres[currentFeas]))
        ibest<-currentFeas[fminInd[1]]
        fbest<-cobra$Fres[ibest]
      }
    }
    return (fbest)
  }
  new_fbest = sapply(1:nrow(cobra$df), calc_fbest)
  nv=data.frame(fbest=cobra$df$Best, fbestArr=cobra$fbestArray, new_fb=new_fbest)
  npts=nrow(nv)
  View(nv)
  testit::assert(all(cobra$fbestArray==new_fbest))
  # cobra$df$Best created before /2025/07/04 may have many entries before first artificial feasible point invalidated (see cobraPhaseII.R:867): 
  # testit::assert(all(cobra$df$Best[305:npts]==new_fbest[305:npts])) 
}

check_bjet <- function(cobra, ind4 = 344) {
  lst = check1(cobra)
  check2(cobra,lst$cond1,lst$cond2)
  check3(cobra,lst$cond1,lst$rs1)
  check4(cobra,lst$rs1,ind4)
  check5(cobra)

  cat("All checks passed.\n")
  # browser()
  
}

chk = check_bjet(cobra,37)
