## This is the example from help page cobyla, with some extensions:
## 1) The 1st call S1 is cobyla with the old convention h_i >= 0 
## 2) The 2nd call S2 is cobyla with the new convention h_i <= 0. In addition, it 
##    shows how a function for hin (inequality constraints with extra parameter)
##    can be wrapped in such a way that the param is not exposed to cobyla.
## 3) Call S3 is the same as S1, but for isres instead of cobyla.
## 4) Call S4 is the same as S2, but for isres instead of cobyla.
## 5) Call S5 is another way to call isres (algorithm from NLopt package)
##
## The problem: Hock-Schittkowski problem no. 100
## The optimum value of the objective function should be 680.6300573
## A suitable parameter vector is roughly
##    (2.330, 1.9514, -0.4775, 4.3657, -0.6245, 1.0381, 1.5942)
##
## RESULTS: cobyla works very well on this problem (error < 1e-7)
##          isres works reasonably (error > 1e-2)
##          isres from NLopt does not work at all (!, error > 10)

x0.hs100 <- c(1, 2, 0, 4, 0, 1, 1)
fn.hs100 <- function(x) {(x[1] - 10) ^ 2 + 5 * (x[2] - 12) ^ 2 + x[3] ^ 4 +
    3 * (x[4] - 11) ^ 2 + 10 * x[5] ^ 6 + 7 * x[6] ^ 2 +
    x[7] ^ 4 - 4 * x[6] * x[7] - 10 * x[6] - 8 * x[7]}

# old (deprecated) convention h_i >= 0 
hin.hs100 <- function(x) {c(
  -(2 * x[1] ^ 2 + 3 * x[2] ^ 4 + x[3] + 4 * x[4] ^ 2 + 5 * x[5] - 127),
  -(7 * x[1] + 3 * x[2] + 10 * x[3] ^ 2 + x[4] - x[5] - 282),
  -(23 * x[1] + x[2] ^ 2 + 6 * x[6] ^ 2 - 8 * x[7] - 196),
  -(4 * x[1] ^ 2 + x[2] ^ 2 - 3 * x[1] * x[2] + 2 * x[3] ^ 2 + 5 * x[6] -
      11 * x[7]))
}

hin.hs100.with.param <- function(x, param) {c(
  -(param$elem1 * x[1] ^ 2 + 3 * x[2] ^ 4 + x[3] + 4 * x[4] ^ 2 + 5 * x[5] - 127),
  -(7 * x[1] + 3 * x[2] + 10 * x[3] ^ 2 + x[4] - x[5] - 282),
  -(23 * x[1] + x[2] ^ 2 + 6 * x[6] ^ 2 - 8 * x[7] - 196),
  -(4 * x[1] ^ 2 + x[2] ^ 2 - 3 * x[1] * x[2] + 2 * x[3] ^ 2 + 5 * x[6] -
    11 * x[7]))
}

tc_init <- function(fn, hn){
  par <- list(fn=fn, 
              elem1=2)
  par$hn <- function(x) {
    return(-hn(x,par))      # new convention h_i <= 0 
  }
  class(par) <- c("PAR","list")
  return(par)
}

param = tc_init(fn.hs100,hin.hs100.with.param)
optim = 680.6300573

### Both of the following cobyla-calls work (currently), but the first one issues
### a warning and may become deprecated in the future.

S1 <- cobyla(x0.hs100, fn.hs100, hin = hin.hs100,
             nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000),
             deprecatedBehavior = TRUE)
S2 <- cobyla(x0.hs100, fn.hs100, hin = param$hn, 
             nl.info = TRUE, control = list(xtol_rel = 1e-8, maxeval = 2000),
             deprecatedBehavior = FALSE)
optErrs = c(S1$value - optim,S2$value - optim)

### Both of the following isres-calls work (currently), but the first one issues
### a warning and may become deprecated in the future.

lb = rep(-5,7)
ub = rep(+5,7)
S3 <- isres(x0 = x0.hs100, fn = fn.hs100, hin = hin.hs100, lower=lb, upper=ub,
           maxeval = 2e5L, deprecatedBehavior = TRUE)
print("*** S3 ***: "); print(S3);
S4 <- isres(x0 = x0.hs100, fn = fn.hs100, hin = param$hn, lower=lb, upper=ub,
            maxeval = 2e5L, deprecatedBehavior = FALSE)
print("*** S4 ***: "); print(S4);
optErrs[3] = S3$value - optim
optErrs[4] = S4$value - optim

opts <- list()
opts$maxeval <- 10000
opts$xtol_rel <- 1e-6
opts$population <- 20 * (length(x0.hs100) + 1)
opts$algorithm <- "NLOPT_GN_ISRES"
S5 <- nloptr::nloptr(x0 = x0.hs100, eval_f = fn.hs100, lb=lb, ub=ub, 
                     eval_g_ineq = param$hn, opts = opts)
optErrs[5] = S5$objective - optim

print(optErrs)