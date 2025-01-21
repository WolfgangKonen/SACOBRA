##creating an instance for G24 problem G24<-COP$new("G24")

##initializing SACOBRA 
cobra <- cobraInit(xStart=G24$lower, fName=G24$name,
                  fn=G24$fn, lower=G24$lower, upper=G24$upper, feval=25)

## Run sacobra optimizer 
cobra <- cobraPhaseII(cobra)

## The true solution is at solu = G24$solu 
## The solution found by SACOBRA:
print(getXbest(cobra)) 
print(getFbest(cobra))
plot(abs(cobra$df$Best-G24$fn(G24$solu)[1]),log="y",type="l",
     ylab="error",xlab="iteration",main=G24$name)


## creating an instance for G03 in 2-dimensional space G03<-COP$new("G03",2)

## Initializing sacobra 
cobra <- cobraInit(xStart=G03$lower, fn=G03$fn, fName=G03$name, 
                   lower=G03$lower, upper=G03$upper, feval=40)

## Run sacobra optimizer 
cobra <- cobraPhaseII(cobra)

## The true solution is at solu = G24$solu 
## The solution found by SACOBRA:
print(getXbest(cobra)) 
print(getFbest(cobra))
plot(abs(cobra$df$Best-G03$fn(G03$solu)[1]),log="y",type="l",
     ylab="error",xlab="iteration",main=G03$name)
