# Script screen_bjet is for running demo-bjet with a certain hyperparam setting for several seeds
# and saving the results. It assumes that
#     1) sourceCobra.R has been sourced
#     2) demo-bjet_IDF.R has been run up to and including line 'solve_bjet <- function ...'
#     3) directory dirName is not existing in pwd()
#
screen_bjet <- function(seedList, dirName="screen_res") {
  if (file.exists(dirName)) {
    warning(sprintf("Directory %s already exists. Use 'unlink(%s,recursive=T)' to delete. Returning.",
                    dirName, dirName))
    return (-1)        
  }
  dir.create(dirName)
  result_df = NULL
  for (seed in seedList) {
    fileName = sprintf("cobra_muGrow_%03d.Rdata",seed)
    # time the main call with system.time:
    sys_t = system.time(cobra <- solve_bjet(seed))
    save(cobra, file=paste(dirName,fileName,sep="/"))
    result_df = rbind(result_df, data.frame(seed=seed, 
                                            fbest=tail(cobra$df$Best,1),
                                            nfeas=length(which(cobra$df$trueNViol==0)),
                                            tUser=sys_t["user.self"],
                                            tElapsed=sys_t["elapsed"],
                                            row.names=seed)
                      )
  }
  save(result_df, file=paste(dirName,"results.Rdata",sep="/"))
  return (result_df)
}


# recompute (from cobra) the columns fbest, nfeas, seed of data frame result_df in directory dirName.
# If save=T, save the new result_df to "results.Rdata" in dirName.
# Return the new result_df.
recomp_screen_res <- function(dirName, save=F) {
  fbest = c()
  nfeas = c()
  seed = c()
  files = sort(dir(dirName))    # sort by ascending seed
  for (f in files) {
    if (f=="results.Rdata") {
      load(file=paste(dirName,f,sep="/"))    # result_df
      result_df = result_df[order(result_df$seed),]    # sort (order) by ascending seed
    } else {
      load(file=paste(dirName,f,sep="/"))    # cobra
      fbest = c(fbest, tail(cobra$df$Best,1))
      nfeas = c(nfeas, length(which(cobra$df$trueNViol==0)))
      seed = c(seed, cobra$cobraSeed)
    }
  }
  testit::assert(all(seed==result_df$seed))
  new_df = data.frame(seed=seed, 
                      fbest=fbest,
                      nfeas=nfeas,
                      tUser=result_df$tUser,
                      tElapsed=result_df$tElapsed,
                      row.names=seed)
  result_df = new_df
  if (save)
    save(result_df, file=paste(dirName,"results.Rdata",sep="/"))
  return (result_df)
}

# combine the results in dir1 and dir2 by recomputing (from cobra) columns fbest, nfeas, seed of 
# data frame result_df.
# If dirOut!=NULL, save all cobra files and the new result_df to dirOut.
# Return the new result_df.
combine_screen_res <- function(dir1,dir2, dirOut=NULL) {
  if (!is.null(dirOut)) {
    if (file.exists(dirOut)) {
      warning(sprintf("Directory %s already exists. Use 'unlink(%s,recursive=T)' to delete. Returning.",
                      dirOut, dirOut))
      return (-1)        
    }
    dir.create(dirOut)
  }
  fbest = c()
  nfeas = c()
  seed = c()
  allres_df = data.frame()
  files = sort(c(paste(dir1,dir(dir1),sep="/"),paste(dir2,dir(dir2),sep="/")))    # sort by ascending seed
  for (f in files) {
    if (basename(f)=="results.Rdata") {
      load(file=f)    # result_df
      allres_df = rbind(allres_df,result_df)
    } else {
      load(file=f)    # cobra
      if (!is.null(dirOut))
        save(cobra, file=paste(dirOut,basename(f),sep="/"))
      fbest = c(fbest, tail(cobra$df$Best,1))
      nfeas = c(nfeas, length(which(cobra$df$trueNViol==0)))
      seed = c(seed, cobra$cobraSeed)
    }
  }
  allres_df = allres_df[order(allres_df$seed),]    # sort (order) by ascending seed
  fbest = fbest[order(seed)]
  nfeas = nfeas[order(seed)]
  seed = seed[order(seed)]
  testit::assert(all(seed==allres_df$seed))
  new_df = data.frame(seed=seed, 
                      fbest=fbest,
                      nfeas=nfeas,
                      tUser=allres_df$tUser,
                      tElapsed=allres_df$tElapsed,
                      row.names=seed)
  result_df = new_df
  if (!is.null(dirOut))
    save(result_df, file=paste(dirOut,"results.Rdata",sep="/"))
  return (result_df)
}


#res_df = screen_bjet(c(44,45,47,48), dirName="screen_res3")
res_df = screen_bjet(c(42,43,46,49), dirName="screen_res4")
print(res_df)
