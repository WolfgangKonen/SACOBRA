# Unit test test_init_design_R from testCobraInit.py ported to R
# ["C:\user\datasets\Vorlesungen\Mathe-Divers\ProjectEuler.net\SACOBRA_Py\demo\testCobraInit.py"]
#
require(pracma)   # needed for fprintf

demo_id_2 <- function(ngrid,nobs,runs) {
  cat("*** Starting demo_id_2 ... ***\n")
  fn <- function(x){c(obj=3*sum(x^2), sum(x)-1)}

  d = 2
  xStart = c(2.5, 2.4)
  cobra <- cobraInit(xStart=xStart, fName="sphere", initDesign = "RAND_R",
                     fn=fn, lower=c(-5,-5), upper=c(5,5), feval=25)
  
  print(cobra$A)
  fprintf("%12.8f",cobra$Fres)
  fprintf("\n")
  fprintf("%12.8f\n",cobra$Gres)
  # [,1]       [,2]
  # [1,] -0.96472247 -0.7704361
  # [2,]  0.16435849  0.2802904
  # [3,] -0.09018369  0.1833372
  # [4,]  0.24347296 -0.3384863
  # [5,]  0.50000000  0.4800000
  # [1] 114.319589   7.918231   3.130921  13.038905  36.030000
  # 
  # [1,] -9.6757927
  # [2,]  1.2232444
  # [3,] -0.5342326
  # [4,] -1.4750667
  # [5,]  3.9000000
  
  # These numbers are exactly (np.allclose) reproduced by test_init_design_R on  
  # the Python side 
}

demo_id_2(ngrid,nobs,runs)
