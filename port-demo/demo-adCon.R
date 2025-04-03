# Unit test test_init_design_R from testCobraInit.py ported to R
# ["C:\user\datasets\Vorlesungen\Mathe-Divers\ProjectEuler.net\SACOBRA_Py\demo\testCobraInit.py"]
#
require(pracma)   # needed for fprintf

demo_adCon <- function() {
  cat("*** Starting demo_adCon ... ***\n")
  
  # define a problem that will trigger constraint normalization: 
  fn <- function(x){c(obj=3*sum(x^2), sum(x)-1, -3000*(sum(x)-10))}

  d = 2
  xStart = c(2.5, 2.4)
  cobra <- cobraInit(xStart=xStart, fName="sphere", initDesign = "RAND_R",
                     fn=fn, lower=c(-5,-5), upper=c(5,5), feval=25)
  print(cobra$Gres)
  #            [,1]     [,2]
  # [1,] -14518.527 28023.03
  # [2,]   1835.478 11669.02
  # [3,]   -801.616 14306.12
  # [4,]  -2213.338 15717.84
  # [5,]   5851.950  7652.55
  
  fprintf("%12.8f\n",cobra$Gres) # just to get more digits after decimal point

  # These numbers are nearly exactly (np.allclose) reproduced by test_adCon_R on  
  # the Python side 
}

demo_adCon()
