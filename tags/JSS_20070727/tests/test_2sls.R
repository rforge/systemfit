
library( systemfit )
data( "Kmenta" )
useMatrix <- FALSE

demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
inst1  <- ~ income + farmPrice
instlist <- list( inst1, inst )
system <- list( demand = demand, supply = supply )
restrm <- matrix(0,1,7)  # restriction matrix "R"
restrm[1,3] <-  1
restrm[1,7] <- -1
restrict <- "demand_income - supply_trend = 0"
restr2m <- matrix(0,2,7)  # restriction matrix "R" 2
restr2m[1,3] <-  1
restr2m[1,7] <- -1
restr2m[2,2] <- -1
restr2m[2,5] <-  1
restr2q <- c( 0, 0.5 )  # restriction vector "q" 2
restrict2 <- c( "demand_income - supply_trend = 0",
   "- demand_price + supply_price = 0.5" )
tc <- matrix(0,7,6)
tc[1,1] <- 1
tc[2,2] <- 1
tc[3,3] <- 1
tc[4,4] <- 1
tc[5,5] <- 1
tc[6,6] <- 1
tc[7,3] <- 1
restr3m <- matrix(0,1,6)  # restriction matrix "R" 2
restr3m[1,2] <- -1
restr3m[1,5] <-  1
restr3q <- c( 0.5 )  # restriction vector "q" 2
restrict3 <- "- C2 + C5 = 0.5"

# It is not possible to estimate 2SLS with systemfit exactly
# as EViews does, because EViews uses
# methodResidCov == "geomean" for the coefficient covariance matrix and
# methodResidCov == "noDfCor" for the residual covariance matrix.
# systemfit uses always the same formulas for both calculations.

## *************** 2SLS estimation ************************
## ************ 2SLS estimation (default)*********************
fit2sls1 <- systemfit( system, "2SLS", data = Kmenta, inst = inst,
   returnModelMatrix = TRUE, useMatrix = useMatrix )
print( summary( fit2sls1 ) )

## *************** 2SLS estimation (single.eq.sigma=F)*******************
fit2sls1s <- systemfit( system, "2SLS", data = Kmenta, inst = inst,
   single.eq.sigma = FALSE, useMatrix = useMatrix )
print( summary( fit2sls1s ) )

## ********************* 2SLS (useDfSys = TRUE) *****************
print( summary( fit2sls1, useDfSys = TRUE ) )

## ********************* 2SLS (methodResidCov = "noDfCor" ) *****************
fit2sls1r <- systemfit( system, "2SLS", data = Kmenta, inst = inst,
   methodResidCov = "noDfCor", useMatrix = useMatrix )
print( summary( fit2sls1r ) )

## *************** 2SLS (methodResidCov="noDfCor", single.eq.sigma=F) *************
fit2sls1rs <- systemfit( system, "2SLS", data = Kmenta, inst = inst,
   methodResidCov = "noDfCor", single.eq.sigma = FALSE, useMatrix = useMatrix )
print( summary( fit2sls1rs ) )

## ********************* 2SLS with restriction ********************
## **************** 2SLS with restriction (default)********************
fit2sls2 <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = inst, useMatrix = useMatrix )
print( summary( fit2sls2 ) )
# the same with symbolically specified restrictions
fit2sls2Sym <- systemfit( system, "2SLS", data = Kmenta,
   restrict.matrix = restrict, inst = inst, useMatrix = useMatrix )
all.equal( fit2sls2, fit2sls2Sym )

## ************* 2SLS with restriction (single.eq.sigma=T) *****************
fit2sls2s <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = inst, single.eq.sigma = TRUE, returnModelMatrix = TRUE,
   useMatrix = useMatrix )
print( summary( fit2sls2s ) )

## ********************* 2SLS with restriction (useDfSys=T) **************
print( summary( fit2sls2, useDfSys = TRUE ) )

## ********************* 2SLS with restriction (methodResidCov = "noDfCor") **************
fit2sls2r <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = inst, methodResidCov = "noDfCor", useMatrix = useMatrix )
print( summary( fit2sls2r ) )

## ******** 2SLS with restriction (methodResidCov="noDfCor", single.eq.sigma=TRUE) *********
fit2sls2rs <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = inst, methodResidCov = "noDfCor", single.eq.sigma = TRUE,
   useMatrix = useMatrix )
print( summary( fit2sls2rs ) )

## ********************* 2SLS with restriction via restrict.regMat ******************
## *************** 2SLS with restriction via restrict.regMat (default )***************
fit2sls3 <- systemfit( system, "2SLS", data = Kmenta, restrict.regMat = tc,
   inst = inst, methodResidCov = "noDfCor", useMatrix = useMatrix )
print( summary( fit2sls3, useDfSys = TRUE ) )


## ***************** 2SLS with 2 restrictions *******************
## ************** 2SLS with 2 restrictions (default) **************
fit2sls4 <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, inst = inst, useMatrix = useMatrix )
print( summary( fit2sls4 ) )
# the same with symbolically specified restrictions
fit2sls4Sym <- systemfit( system, "2SLS", data = Kmenta,
   restrict.matrix = restrict2, inst = inst, useMatrix = useMatrix )
all.equal( fit2sls4, fit2sls4Sym )

## ************ 2SLS with 2 restrictions (single.eq.sigma=T) **************
fit2sls4s <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, inst = inst, single.eq.sigma = TRUE,
   useMatrix = useMatrix )
print( summary( fit2sls4s ) )

## ***************** 2SLS with 2 restrictions (useDfSys=T) **************
print( summary( fit2sls4, useDfSys = TRUE ) )

## ***************** 2SLS with 2 restrictions (methodResidCov="noDfCor") **************
fit2sls4r <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, inst = inst, methodResidCov = "noDfCor",
   returnModelMatrix = TRUE, useMatrix = useMatrix )
print( summary( fit2sls4r ) )

## ***** 2SLS with 2 restrictions (methodResidCov="noDfCor", single.eq.sigma=T) *******
fit2sls4rs <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, inst = inst, methodResidCov = "noDfCor",
   single.eq.sigma = TRUE, useMatrix = useMatrix )
print( summary( fit2sls4rs ) )

## ************* 2SLS with 2 restrictions via R and restrict.regMat ******************
## ******** 2SLS with 2 restrictions via R and restrict.regMat (default) *************
fit2sls5 <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, inst = inst,
   useMatrix = useMatrix )
print( summary( fit2sls5 ) )
# the same with symbolically specified restrictions
fit2sls5Sym <- systemfit( system, "2SLS", data = Kmenta,
   restrict.matrix = restrict3, restrict.regMat = tc, inst = inst,
   useMatrix = useMatrix )
all.equal( fit2sls5, fit2sls5Sym )

## ******* 2SLS with 2 restrictions via R and restrict.regMat (single.eq.sigma=T) ******
fit2sls5s <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, inst = inst,
   single.eq.sigma = TRUE, useMatrix = useMatrix )
print( summary( fit2sls5s ) )

## ********** 2SLS with 2 restrictions via R and restrict.regMat (useDfSys=T) *******
print( summary( fit2sls5, useDfSys = TRUE ) )

## ************* 2SLS with 2 restrictions via R and restrict.regMat (methodResidCov="noDfCor") *********
fit2sls5r <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, inst = inst,
   methodResidCov = "noDfCor", useMatrix = useMatrix )
print( summary( fit2sls5r ) )

## ** 2SLS with 2 restrictions via R and restrict.regMat (methodResidCov="noDfCor", single.eq.sigma=T) **
fit2sls5rs <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, inst = inst,
   methodResidCov = "noDfCor", single.eq.sigma = TRUE,
   returnModelMatrix = TRUE, useMatrix = useMatrix )
print( summary( fit2sls5rs ) )

## *********** 2SLS estimation with different instruments **************
## ******* 2SLS estimation with different instruments (default) *********
fit2slsd1 <- systemfit( system, "2SLS", data = Kmenta, inst = instlist,
   useMatrix = useMatrix )
print( summary( fit2slsd1 ) )

## *********** 2SLS estimation with different instruments (single.eq.sigma=F)*****
fit2slsd1s <- systemfit( system, "2SLS", data = Kmenta, inst = instlist,
   single.eq.sigma = FALSE, useMatrix = useMatrix )
print( summary( fit2slsd1s ) )

## ********* 2SLS estimation with different instruments (useDfSys=T) *******
print( summary( fit2slsd1, useDfSys = TRUE ) )

## ********* 2SLS estimation with different instruments (methodResidCov="noDfCor") ******
fit2slsd1r <- systemfit( system, "2SLS", data = Kmenta, inst = instlist,
   methodResidCov = "noDfCor", useMatrix = useMatrix )
print( summary( fit2slsd1r ) )

## 2SLS estimation with different instruments (methodResidCov="noDfCor",single.eq.sigma=F)
fit2slsd1r <- systemfit( system, "2SLS", data = Kmenta, inst = instlist,
   methodResidCov = "noDfCor", single.eq.sigma = FALSE,
   useMatrix = useMatrix )
print( summary( fit2slsd1r ) )

## **** 2SLS estimation with different instruments and restriction *******
## ** 2SLS estimation with different instruments and restriction (default) ****
fit2slsd2 <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = instlist, useMatrix = useMatrix )
print( summary( fit2slsd2 ) )

## 2SLS estimation with different instruments and restriction (single.eq.sigma=T)
fit2slsd2s <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = instlist, single.eq.sigma = TRUE, useMatrix = useMatrix )
print( summary( fit2slsd2s ) )

## **** 2SLS estimation with different instruments and restriction (useDfSys=F)
print( summary( fit2slsd2, useDfSys = FALSE ) )

## **** 2SLS estimation with different instruments and restriction (methodResidCov="noDfCor")
fit2slsd2r <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = instlist, methodResidCov = "noDfCor", useMatrix = useMatrix )
print( summary( fit2slsd2r ) )

## 2SLS estimation with different instr. and restr. (methodResidCov="noDfCor", single.eq.sigma=T)
fit2slsd2rs <- systemfit( system, "2SLS", data = Kmenta, restrict.matrix = restrm,
   inst = instlist, methodResidCov = "noDfCor", single.eq.sigma = TRUE,
   useMatrix = useMatrix )
print( summary( fit2slsd2rs ) )

## **** 2SLS estimation with different instruments and restriction via restrict.regMat *
## 2SLS estimation with different instruments and restriction via restrict.regMat (default)
fit2slsd3 <- systemfit( system, "2SLS", data = Kmenta, restrict.regMat = tc,
   inst = instlist, useMatrix = useMatrix )
print( summary( fit2slsd3 ) )

## **** 2SLS estimation with different instr. and restr. via restrict.regMat (methodResidCov="noDfCor")
fit2slsd3r <- systemfit( system, "2SLS", data = Kmenta, restrict.regMat = tc,
   inst = instlist, methodResidCov = "noDfCor", useMatrix = useMatrix )
print( summary( fit2slsd3r ) )


## **************** shorter summaries **********************
print( summary( fit2sls1, useDfSys = TRUE, printResidCov = FALSE ) )

print( summary( fit2sls1, printEquations = FALSE ) )

print( summary( fit2sls1rs, printResidCov = FALSE, printEquations = FALSE ) )

print( summary( fit2sls2Sym, useDfSys = FALSE ), equations = FALSE )

print( summary( fit2sls2 ), residCov = FALSE )

print( summary( fit2sls3, useDfSys = FALSE, printResidCov = FALSE,
   printEquations = FALSE ) )

print( summary( fit2sls4s ), equations = FALSE, residCov = FALSE )

print( summary( fit2sls5r, printEquations = FALSE, printResidCov = FALSE ) )

print( summary( fit2slsd1s ), residCov = FALSE, equations = FALSE )

print( summary( fit2slsd2, printResidCov = FALSE, printEquations = FALSE ) )

print( summary( fit2slsd3r ), residCov = FALSE, equations = FALSE )


## ****************** residuals **************************
print( residuals( fit2sls1 ) )
print( residuals( fit2sls1$eq[[ 1 ]] ) )

print( residuals( fit2sls2s ) )
print( residuals( fit2sls2s$eq[[ 2 ]] ) )

print( residuals( fit2sls3 ) )
print( residuals( fit2sls3$eq[[ 1 ]] ) )

print( residuals( fit2sls4r ) )
print( residuals( fit2sls4r$eq[[ 2 ]] ) )

print( residuals( fit2sls5rs ) )
print( residuals( fit2sls5rs$eq[[ 1 ]] ) )

print( residuals( fit2slsd1 ) )
print( residuals( fit2slsd1$eq[[ 2 ]] ) )

print( residuals( fit2slsd2r ) )
print( residuals( fit2slsd2r$eq[[ 1 ]] ) )


## *************** coefficients *********************
print( round( coef( fit2sls1s ), digits = 6 ) )
print( round( coef( fit2sls1s$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fit2sls2 ), digits = 6 ) )
print( round( coef( fit2sls2$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fit2sls3 ), digits = 6 ) )
print( round( coef( fit2sls3, modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( fit2sls3$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fit2sls4s ), digits = 6 ) )
print( round( coef( fit2sls4s$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fit2sls5r ), digits = 6 ) )
print( round( coef( fit2sls5r, modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( fit2sls5r$eq[[ 2 ]] ), digits = 6 ) )


## *************** coefficients with stats *********************
print( round( coef( summary( fit2sls1s ) ), digits = 6 ) )
print( round( coef( summary( fit2sls1s$eq[[ 1 ]] ) ), digits = 6 ) )

print( round( coef( summary( fit2sls2, useDfSys = FALSE ) ), digits = 6 ) )
print( round( coef( summary( fit2sls2$eq[[ 2 ]], useDfSys = FALSE ) ),
   digits = 6 ) )

print( round( coef( summary( fit2sls3 ) ), digits = 6 ) )
print( round( coef( summary( fit2sls3 ), modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( summary( fit2sls3$eq[[ 1 ]] ) ), digits = 6 ) )

print( round( coef( summary( fit2sls4s ) ), digits = 6 ) )
print( round( coef( summary( fit2sls4s$eq[[ 2 ]] ) ), digits = 6 ) )

print( round( coef( summary( fit2sls5r, useDfSys = FALSE ) ), digits = 6 ) )
print( round( coef( summary( fit2sls5r, useDfSys = FALSE ),
   modified.regMat = TRUE ), digits = 6 ) )
print( round( coef( summary( fit2sls5r$eq[[ 2 ]], useDfSys = FALSE ) ),
   digits = 6 ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fit2sls1s ), digits = 6 ) )
print( round( vcov( fit2sls1s$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit2sls1r ), digits = 6 ) )
print( round( vcov( fit2sls1r$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit2sls2 ), digits = 6 ) )
print( round( vcov( fit2sls2$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit2sls3 ), digits = 6 ) )
print( round( vcov( fit2sls3, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fit2sls3$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit2sls4s ), digits = 6 ) )
print( round( vcov( fit2sls4s$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit2sls5r ), digits = 6 ) )
print( round( vcov( fit2sls5r, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fit2sls5r$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit2slsd1 ), digits = 6 ) )
print( round( vcov( fit2slsd1$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fit2slsd2rs ), digits = 6 ) )
print( round( vcov( fit2slsd2rs$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fit2slsd3 ), digits = 6 ) )
print( round( vcov( fit2slsd3, modified.regMat = TRUE ), digits = 6 ) )
print( round( vcov( fit2slsd3$eq[[ 1 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fit2sls1 ) )
print( confint( fit2sls1$eq[[ 1 ]], level = 0.9 ) )

print( confint( fit2sls2s, level = 0.9 ) )
print( confint( fit2sls2s$eq[[ 2 ]], level = 0.99 ) )

print( confint( fit2sls3, level = 0.99, useDfSys = TRUE ) )
print( confint( fit2sls3$eq[[ 1 ]], level = 0.5, useDfSys = TRUE ) )

print( confint( fit2sls4r, level = 0.5 ) )
print( confint( fit2sls4r$eq[[ 2 ]], level = 0.25 ) )

print( confint( fit2sls5rs, level = 0.25 ) )
print( confint( fit2sls5rs$eq[[ 1 ]], level = 0.975 ) )

print( confint( fit2slsd1, level = 0.975, useDfSys = TRUE ) )
print( confint( fit2slsd1$eq[[ 2 ]], level = 0.999, useDfSys = TRUE ) )

print( confint( fit2slsd2r, level = 0.999 ) )
print( confint( fit2slsd2r$eq[[ 1 ]] ) )


## *********** fitted values *************
print( fitted( fit2sls1, se.fit = TRUE, interval = "prediction" ) )
print( fitted( fit2sls1$eq[[ 1 ]] ) )

print( fitted( fit2sls2s ) )
print( fitted( fit2sls2s$eq[[ 2 ]] ) )

print( fitted( fit2sls3 ) )
print( fitted( fit2sls3$eq[[ 1 ]] ) )

print( fitted( fit2sls4r ) )
print( fitted( fit2sls4r$eq[[ 2 ]] ) )

print( fitted( fit2sls5rs ) )
print( fitted( fit2sls5rs$eq[[ 1 ]] ) )

print( fitted( fit2slsd1 ) )
print( fitted( fit2slsd1$eq[[ 2 ]] ) )

print( fitted( fit2slsd2r ) )
print( fitted( fit2slsd2r$eq[[ 1 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$consump <- NULL
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fit2sls1, se.fit = TRUE, interval = "prediction" ) )
print( predict( fit2sls1$eq[[ 1 ]], se.fit = TRUE, interval = "prediction" ) )

print( predict( fit2sls2s, se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )
print( predict( fit2sls2s$eq[[ 2 ]], se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )

print( predict( fit2sls3, se.pred = TRUE, interval = "prediction",
   level = 0.975, useDfSys = TRUE ) )
print( predict( fit2sls3$eq[[ 1 ]], se.pred = TRUE, interval = "prediction",
   level = 0.975, useDfSys = TRUE ) )

print( predict( fit2sls4r, se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )
print( predict( fit2sls4r$eq[[ 2 ]], se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )

print( predict( fit2sls5rs, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )
print( predict( fit2sls5rs$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )

print( predict( fit2slsd1, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99, useDfSys = TRUE ) )
print( predict( fit2slsd1$eq[[ 2 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99, useDfSys = TRUE ) )

print( predict( fit2slsd2r, se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )
print( predict( fit2slsd2r$eq[[ 1 ]], se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25 )

print( predict( fit2sls1rs, newdata = smallData ) )
print( predict( fit2sls1rs$eq[[ 1 ]], newdata = smallData ) )

print( predict( fit2sls2, se.fit = TRUE, level = 0.9,
   newdata = smallData ) )
print( predict( fit2sls2$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   newdata = smallData ) )

print( predict( fit2sls3, interval = "prediction", level = 0.975,
   newdata = smallData ) )
print( predict( fit2sls3$eq[[ 1 ]], interval = "confidence", level = 0.8,
   newdata = smallData ) )

print( predict( fit2sls4r, se.fit = TRUE, interval = "confidence",
   level = 0.999, newdata = smallData ) )
print( predict( fit2sls4r$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, newdata = smallData ) )

print( predict( fit2sls5s, se.fit = TRUE, interval = "prediction",
   newdata = smallData ) )
print( predict( fit2sls5s$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   newdata = smallData ) )

print( predict( fit2slsd3, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fit2slsd3$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fit2sls1, 1, 2 ) )

print( correlation.systemfit( fit2sls2s, 2, 1 ) )

print( correlation.systemfit( fit2sls3, 1, 2 ) )

print( correlation.systemfit( fit2sls4r, 2, 1 ) )

print( correlation.systemfit( fit2sls5rs, 1, 2 ) )

print( correlation.systemfit( fit2slsd1, 2, 1 ) )

print( correlation.systemfit( fit2slsd2r, 1, 2 ) )


## ************ Log-Likelihood values ***************
print( logLik( fit2sls1 ) )

print( logLik( fit2sls2s ) )

print( logLik( fit2sls3 ) )

print( logLik( fit2sls4r ) )

print( logLik( fit2sls5rs ) )

print( logLik( fit2slsd1 ) )

print( logLik( fit2slsd2r ) )


## ************** F tests ****************
# testing first restriction
print( linear.hypothesis( fit2sls1, restrm ) )
linear.hypothesis( fit2sls1, restrict )

print( linear.hypothesis( fit2sls1s, restrm ) )
linear.hypothesis( fit2sls1s, restrict )

print( linear.hypothesis( fit2sls1, restrm ) )
linear.hypothesis( fit2sls1, restrict )

print( linear.hypothesis( fit2sls1r, restrm ) )
linear.hypothesis( fit2sls1r, restrict )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
restrictOnly2 <- "- demand_price + supply_price = 0.5"
# first restriction not imposed 
print( linear.hypothesis( fit2sls1, restrOnly2m, restrOnly2q ) )
linear.hypothesis( fit2sls1, restrictOnly2 )

# first restriction imposed
print( linear.hypothesis( fit2sls2, restrOnly2m, restrOnly2q ) )
linear.hypothesis( fit2sls2, restrictOnly2 )

print( linear.hypothesis( fit2sls2r, restrOnly2m, restrOnly2q ) )
linear.hypothesis( fit2sls2r, restrictOnly2 )

print( linear.hypothesis( fit2sls3, restrOnly2m, restrOnly2q ) )
linear.hypothesis( fit2sls3, restrictOnly2 )

# testing both of the restrictions
print( linear.hypothesis( fit2sls1, restr2m, restr2q ) )
linear.hypothesis( fit2sls1, restrict2 )


## ************** Wald tests ****************
# testing first restriction
print( linear.hypothesis( fit2sls1, restrm, test = "Chisq" ) )
linear.hypothesis( fit2sls1, restrict, test = "Chisq" )

print( linear.hypothesis( fit2sls1s, restrm, test = "Chisq" ) )
linear.hypothesis( fit2sls1s, restrict, test = "Chisq" )

print( linear.hypothesis( fit2sls1, restrm, test = "Chisq" ) )
linear.hypothesis( fit2sls1, restrict, test = "Chisq" )

print( linear.hypothesis( fit2sls1r, restrm, test = "Chisq" ) )
linear.hypothesis( fit2sls1r, restrict, test = "Chisq" )

# testing second restriction
# first restriction not imposed
print( linear.hypothesis( fit2sls1, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linear.hypothesis( fit2sls1, restrictOnly2, test = "Chisq" )
# first restriction imposed
print( linear.hypothesis( fit2sls2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linear.hypothesis( fit2sls2, restrictOnly2, test = "Chisq" )

print( linear.hypothesis( fit2sls2r, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linear.hypothesis( fit2sls2r, restrictOnly2, test = "Chisq" )

print( linear.hypothesis( fit2sls3, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linear.hypothesis( fit2sls3, restrictOnly2, test = "Chisq" )

# testing both of the restrictions
print( linear.hypothesis( fit2sls1, restr2m, restr2q, test = "Chisq" ) )
linear.hypothesis( fit2sls1, restrict2, test = "Chisq" )


## **************** model frame ************************
print( mf <- model.frame( fit2sls1 ) )
print( mf1 <- model.frame( fit2sls1$eq[[ 1 ]] ) )
print( attributes( mf1 )$terms )
print( mf2 <- model.frame( fit2sls1$eq[[ 2 ]] ) )
print( attributes( mf2 )$terms )

print( all.equal( mf, model.frame( fit2sls2s ) ) )
print( all.equal( mf2, model.frame( fit2sls2s$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fit2sls3 ) ) )
print( all.equal( mf1, model.frame( fit2sls3$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fit2sls4r ) ) )
print( all.equal( mf2, model.frame( fit2sls4r$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fit2sls5rs ) ) )
print( all.equal( mf1, model.frame( fit2sls5rs$eq[[ 1 ]] ) ) )


## **************** model matrix ************************
# with returnModelMatrix = TRUE
print( !is.null( fit2sls1$eq[[ 1 ]]$modelMatrix ) )
print( mm <- model.matrix( fit2sls1 ) )
print( mm1 <- model.matrix( fit2sls1$eq[[ 1 ]] ) )
print( mm2 <- model.matrix( fit2sls1$eq[[ 2 ]] ) )

# with returnModelMatrix = FALSE
print( all.equal( mm, model.matrix( fit2sls1s ) ) )
print( all.equal( mm1, model.matrix( fit2sls1s$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit2sls1s$eq[[ 2 ]] ) ) )
print( !is.null( fit2sls1s$eq[[ 1 ]]$modelMatrix ) )

# with returnModelMatrix = TRUE
print( !is.null( fit2sls2s$eq[[ 1 ]]$modelMatrix ) )
print( all.equal( mm, model.matrix( fit2sls2s ) ) )
print( all.equal( mm1, model.matrix( fit2sls2s$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit2sls2s$eq[[ 2 ]] ) ) )

# with returnModelMatrix = FALSE
print( all.equal( mm, model.matrix( fit2sls2Sym ) ) )
print( all.equal( mm1, model.matrix( fit2sls2Sym$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit2sls2Sym$eq[[ 2 ]] ) ) )
print( !is.null( fit2sls2Sym$eq[[ 1 ]]$modelMatrix ) )

# with returnModelMatrix = FALSE
print( all.equal( mm, model.matrix( fit2sls3 ) ) )
print( all.equal( mm1, model.matrix( fit2sls3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit2sls3$eq[[ 2 ]] ) ) )
print( !is.null( fit2sls3$eq[[ 1 ]]$modelMatrix ) )

# with returnModelMatrix = TRUE
print( !is.null( fit2sls4r$eq[[ 1 ]]$modelMatrix ) )
print( all.equal( mm, model.matrix( fit2sls4r ) ) )
print( all.equal( mm1, model.matrix( fit2sls4r$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit2sls4r$eq[[ 2 ]] ) ) )

# with returnModelMatrix = FALSE
print( all.equal( mm, model.matrix( fit2sls4s ) ) )
print( all.equal( mm1, model.matrix( fit2sls4s$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit2sls4s$eq[[ 2 ]] ) ) )
print( !is.null( fit2sls4s$eq[[ 1 ]]$modelMatrix ) )

# with returnModelMatrix = TRUE
print( !is.null( fit2sls5rs$eq[[ 1 ]]$modelMatrix ) )
print( all.equal( mm, model.matrix( fit2sls5rs ) ) )
print( all.equal( mm1, model.matrix( fit2sls5rs$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit2sls5rs$eq[[ 2 ]] ) ) )

# with returnModelMatrix = FALSE
print( all.equal( mm, model.matrix( fit2sls5r ) ) )
print( all.equal( mm1, model.matrix( fit2sls5r$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fit2sls5r$eq[[ 2 ]] ) ) )
print( !is.null( fit2sls5r$eq[[ 1 ]]$modelMatrix ) )


## **************** formulas ************************
formula( fit2sls1 )
formula( fit2sls1$eq[[ 1 ]] )

formula( fit2sls2s )
formula( fit2sls2s$eq[[ 2 ]] )

formula( fit2sls3 )
formula( fit2sls3$eq[[ 1 ]] )

formula( fit2sls4r )
formula( fit2sls4r$eq[[ 2 ]] )

formula( fit2sls5rs )
formula( fit2sls5rs$eq[[ 1 ]] )

formula( fit2slsd1 )
formula( fit2slsd1$eq[[ 2 ]] )

formula( fit2slsd2r )
formula( fit2slsd2r$eq[[ 1 ]] )


## **************** model terms *******************
terms( fit2sls1 )
terms( fit2sls1$eq[[ 1 ]] )

terms( fit2sls2s )
terms( fit2sls2s$eq[[ 2 ]] )

terms( fit2sls3 )
terms( fit2sls3$eq[[ 1 ]] )

terms( fit2sls4r )
terms( fit2sls4r$eq[[ 2 ]] )

terms( fit2sls5rs )
terms( fit2sls5rs$eq[[ 1 ]] )

terms( fit2slsd1 )
terms( fit2slsd1$eq[[ 2 ]] )

terms( fit2slsd2r )
terms( fit2slsd2r$eq[[ 1 ]] )
