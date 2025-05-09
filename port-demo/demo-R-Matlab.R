if (FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # This example will try to start the MATLAB server on the local machine,
  # and then setup a Matlab object in R for communicating data between R
  # and MATLAB and for sending commands from R to MATLAB.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1.  Load R.matlab
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  library(R.matlab)
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2.  Start MATLAB
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2.1.  Start MATLAB from R?
  # Start MATLAB server on the local machine (if this fails,
  # see help(Matlab) for alternatives).
  Matlab$startServer()
  
  # 2.2.  OR start MATLAB externally,
  #       THEN add 'externals' subdirectory to the MATLAB path
  
  #  (Where is the 'externals' subdirectory?)
  print(system.file("externals", package = "R.matlab"))
  
  #       THEN from within MATLAB,
  #            issue MATLAB command "MatlabServer"
  # Note: If issued from a MATLAB command line, this last command
  #       prevents further MATLAB 'command line' input
  #       until something like close(matlab) at the end of this script
  
  # 2.3.  If both these options fail, see help(Matlab) for alternatives.
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Create a MATLAB client object used to communicate with MATLAB
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  matlab <- Matlab()
  
  # 3.1 Check status of MATLAB connection (not yet connected)
  print(matlab)
  
  # 3.2 If you experience any problems, ask for detailed outputs
  #     by uncommenting the next line
  # setVerbose(matlab, -2)
  
  # 3.3 Connect to the MATLAB server.
  isOpen <- open(matlab)
  
  # 3.4 Confirm that the MATLAB server is open, and running
  if (!isOpen)
    throw("MATLAB server is not running: waited 30 seconds.")
  
  # 3.5 Check status of MATLAB connection (now connected)
  print(matlab)
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 4.  Sample uses of the MATLAB server
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 4.1 Run MATLAB expressions on the MATLAB server
  evaluate(matlab, "A = 1+2;", "B = ones(2, 20);")
  
  # 4.2 Ask MATLAB to display a value (without transferring it to R)
  evaluate(matlab, "A")
  
  # 4.3 Get MATLAB variables
  data <- getVariable(matlab, c("A", "B"))
  cat("Received variables:\n")
  str(data)
  
  # 4.4 Set variables in MATLAB
  ABCD <- matrix(rnorm(10000), ncol = 100)
  str(ABCD)
  setVariable(matlab, ABCD = ABCD)
  
  # 4.5 Retrieve what we just set
  data <- getVariable(matlab, "ABCD")
  cat("Received variables:\n")
  str(data)
  
  # 4.6 Create a function (M-file) on the MATLAB server
  setFunction(matlab, "            \
  function [win, aver] = dice(B) \
  %Play the dice game B times    \
  gains = [-1, 2, -3, 4, -5, 6]; \
  plays = unidrnd(6, B, 1);      \
  win = sum(gains(plays));       \
  aver = win/B;                  \
")
  
  # 4.7 Use the MATLAB function just created
  evaluate(matlab, "[w, a] = dice(1000);")
  res <- getVariable(matlab, c("w", "a"))
  print(res)
  
  # 4.8 Compile a function from Don Jones automotive benchmark:
  evaluate(matlab,"cd('C:/Users/wolfgang/Documents/GitHub/2025JonesBenchmarks/Benchmarks/automotive')")
  evaluate(matlab,"[automotive,frel,grel,hrel,xl,xu,xopt,x0] = automotive_benchmark();")
  
  # 4.9 Use the compiled function:
  evaluate(matlab,"[f,g,h] = automotive(xopt)")
  data <- getVariable(matlab, c("f", "g", "h"))
  print(data$f)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 5. Exception handling
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 5.1 Try to get non-existing MATLAB variable
  #     (will result in an informative error)
  tryCatch({
    data <- getVariable(matlab, "unknown")
    cat("Received variables:\n")
    str(data)
  }, error = function(ex) {
    print(ex)
  })
  # Confirm that things still work
  data <- getVariable(matlab, "A")
  cat("Received variables:\n")
  str(data)
  
  
  # 5.2 Try to evaluate a MATLAB expression that fails
  #     (will result in an informative error)
  tryCatch({
    res <- evaluate(matlab, "C = 1+unknown;")
    res
  }, error = function(ex) {
    print(ex)
  })
  # Confirm that things still work
  res <- evaluate(matlab, "C = 1+2;")
  print(res)
  data <- getVariable(matlab, "C")
  cat("Received variables:\n")
  str(data)
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 6.  Done:  close the MATLAB client
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # When done, close the MATLAB client, which will also shutdown
  # the MATLAB server and the connection to it.
  close(matlab)
  
  # 6.1 Check status of MATLAB connection (now disconnected)
  print(matlab)
}