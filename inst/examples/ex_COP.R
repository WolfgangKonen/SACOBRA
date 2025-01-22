##creating an instance for G24 problem 
G24<-COP$new("G24")

##initializing SACOBRA 
cobra <- cobraInit(xStart=G24$lower, fName=G24$name,
                  fn=G24$fn, lower=G24$lower, upper=G24$upper, feval=25)

## Run sacobra optimizer 
cobra <- cobraPhaseII(cobra)

## The true solution is at G24$solu = c(2.329520, 3.178493) 
## with objective -5.508013
## The solution found by SACOBRA:
print(getXbest(cobra))          # 2.329725 3.176549
print(getFbest(cobra))          # -5.506275
plot(abs(cobra$df$Best-G24$fn(G24$solu)[1]),log="y",type="l",
     ylab="error",xlab="iteration",main=G24$name)


## creating an instance for G06 
G06<-COP$new("G06")

## Initializing sacobra 
cobra <- cobraInit(xStart=G06$lower, fn=G06$fn, fName=G06$name, 
                   lower=G06$lower, upper=G06$upper, feval=40)

## Run sacobra optimizer 
cobra <- cobraPhaseII(cobra)

## The true solution is at G06$solu = c(14.0950, 0.84296)
## with objective -6961.814 
## The solution found by SACOBRA:
print(getXbest(cobra))          # 14.0950000  0.8429608
print(getFbest(cobra))          # -6961.814 
plot(abs(cobra$df$Best-G06$fn(G06$solu)[1]),log="y",type="l",
     ylab="error",xlab="iteration",main=G06$name)
