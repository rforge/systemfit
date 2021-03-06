## prepare summary results that belong to the whole system
summary.systemfit <- function( object, useDfSys = NULL,
      printResidCov = TRUE, printEquations = TRUE, ... ) {

   if( is.null( useDfSys ) ) {
      useDfSys <- length( coef( object ) ) != object$rank
         # TRUE if there are restrictions imposed
   }

   # number of equations
   nEq <- length( object$eq )
   # number of observations per equation
   nObsPerEq <- nrow( residuals( object ) )
   # total number of observations
   nObs <- nObsPerEq * nEq

   # preparing objects that will be returned
   result <- list()
   result$call <- object$call
   result$method <- object$method
   result$iter <- object$iter
   result$control <- object$control
   result$residuals <- residuals( object )
   result$residCovEst <- object$residCovEst
   result$residCov <- object$residCov
   if( !is.null( result$residCovEst ) ) {
      dimnames( result$residCovEst ) <- dimnames( result$residCov )
   }
   result$residCor <- cor( residuals( object ) )
   dimnames( result$residCor ) <- dimnames( result$residCov )
   result$detResidCov <- det( object$residCov, tol = object$control$solvetol )

   # now prepare summury results for the individual equations
   result$eq <- list()
   for( i in 1:length( object$eq ) ) {
       result$eq[[ i ]] <- summary( object$eq[[i]], useDfSys = useDfSys )
   }

   # coefficients, standard errors, ... 
   result$coefCov <- object$coefCov
   coef <- object$coefficients
   stdEr <- diag( result$coefCov )^0.5  # standard errors
   tStat <- coef / stdEr                # t-statistic
   if( useDfSys ) {             # p-values
      pVal <- 2 * ( 1 - pt( abs( tStat ), object$df.residual ) )
   } else {
      pVal <- NULL
      for( i in 1:length( object$eq ) ){
         pVal <- c( pVal, coef( result$eq[[ i ]] )[ , 4 ] )
      }
   }
   result$coefficients <- cbind( coef, stdEr, tStat, pVal )
   colnames( result$coefficients ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   result$df <- c( length( coef( object ) ), nObs - length( coef( object ) ) )

   # coefficients of the modified regressor matrix
   if( !is.null( object$restrict.regMat ) ) {
      coefModReg <- coef( object, modified.regMat = TRUE )
      stdErModReg <- diag( vcov( object, modified.regMat = TRUE ) )^0.5  # standard errors
      tStatModReg <- coefModReg / stdErModReg    # t-statistic
      if( useDfSys ) {             # p-values
         pValModReg <- 2 * ( 1 - pt( abs( tStatModReg ), object$df.residual ) )
      } else {
         pValModReg <- rep( NA, length( coefModReg ) )
      }
      result$coefModReg <- cbind( coefModReg, stdErModReg, tStatModReg, pValModReg )
      colnames( result$coefModReg ) <- c( "Estimate", "Std. Error",
         "t value", "Pr(>|t|)" )
   }

   # R^2 values
   nObsEq <- rep( nObsPerEq, nEq )
   resid <- NULL
   response <- NULL
   responseMinusMean <- NULL
   for( i in 1:length( object$eq ) ) {
      resid <- c( resid, residuals( object$eq[[ i ]] ) )
      responseEqI <- fitted( object$eq[[ i ]] ) + residuals( object$eq[[ i ]] )
      response <- c( response, responseEqI )
      responseMinusMean <- c( responseMinusMean,
         responseEqI - mean( responseEqI ) )
   }

   # OLS R^2 value of the entire system
   rss <- sum( residuals( object )^2 )
   tss <- sum( responseMinusMean^2 )
   result$ols.r.squared <- 1 - rss / tss

   # System R^2 value of McElroy (1977)
   # formula from Greene (2003, p. 345 )
   # (first formula, numerator modified to save memory)
   if( object$method %in% c( "SUR", "3SLS" ) ){
      xMat <- matrix( resid, ncol = 1 )
      if( object$control$useMatrix ){
         object$residCov <- as( object$residCov, "dspMatrix" )
         xMat <- as( xMat, "dgCMatrix" )
      }
      rtOmega <- .calcXtOmegaInv( xMat = xMat,
         sigma = object$residCov, nObsEq = nObsEq,
         solvetol = object$control$solvetol,
         useMatrix = object$control$useMatrix )
      yCov <- .calcResidCov( response, methodResidCov = "noDfCor",
         nObsEq = nObsEq, centered = TRUE,
         solvetol = object$control$solvetol )
      residCovInv <- solve( object$residCov, tol = object$control$solvetol )
      denominator <- 0
      for( i in 1:length( object$eq ) ) {
         for( j in 1:length( object$eq ) ) {
            denominator <- denominator + residCovInv[ i, j ] * yCov[ i, j ] *
               nObsPerEq
         }
      }
      result$mcelroy.r.squared <- drop( 1 - ( rtOmega %*% resid ) / denominator )
   }

   result$printEquations <- printEquations
   result$printResidCov  <- printResidCov

   class( result ) <- "summary.systemfit"
   return( result )
}

## print summary results of the whole system
print.summary.systemfit <- function( x,
      digits = max( 3, getOption("digits") - 1 ),
      residCov = x$printResidCov, equations = x$printEquations, ... ) {

  table <- NULL
  labels <- NULL

  cat("\n")
  cat("systemfit results \n")
  cat("method: ")
  if(!is.null(x$iter)) if(x$iter>1) cat("iterated ")
  cat( paste( x$method, "\n\n"))
  if(!is.null(x$iter)) {
    if(x$iter>1) {
      if(x$iter<x$control$maxiter) {
        cat( paste( "convergence achieved after",x$iter,"iterations\n\n" ) )
      } else {
        cat( paste( "warning: convergence not achieved after", x$iter,
                    "iterations\n\n" ) )
      }
    }
  }
  for(i in 1:length( x$eq ) ) {
    row <- NULL
    row <- cbind( round( sum( x$eq[[i]]$df ),  digits ),
                  round( x$eq[[i]]$df[2], digits ),
                  round( x$eq[[i]]$ssr,   digits ),
                  round( x$eq[[i]]$sigma^2, digits ),
                  round( x$eq[[i]]$sigma,   digits ),
                  round( x$eq[[i]]$r.squared,     digits ),
                  round( x$eq[[i]]$adj.r.squared, digits ))
    table  <- rbind( table, row )
    labels <- rbind( labels, x$eq[[i]]$eqnLabel )
  }
  rownames(table) <- c( labels )
  colnames(table) <- c("N","DF", "SSR", "MSE", "RMSE", "R2", "Adj R2" )

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE, digits = digits )

  cat("\n")

   if( residCov ){
      if(!is.null(x$residCovEst)) {
         cat("The covariance matrix of the residuals used for estimation\n")
         print( x$residCovEst, digits = digits )
         cat("\n")
         if( min(eigen( x$residCov, only.values=TRUE)$values) < 0 ) {
            cat("warning: this covariance matrix is NOT positive semidefinit!\n")
            cat("\n")
         }
      }

      cat("The covariance matrix of the residuals\n")
      print( x$residCov, digits = digits )
      cat("\n")

      cat("The correlations of the residuals\n")
      print( x$residCor, digits = digits )
      cat("\n")

      cat("The determinant of the residual covariance matrix: ")
      cat( formatC( x$detResidCov, digits = digits, width = -1 ) )
      cat("\n")
   }

  cat("OLS R-squared value of the system: ")
  cat( formatC( x$ols.r.squared, digits = digits, width = -1 ) )
  cat("\n")

  if(!is.null(x$mcelroy.r.squared)) {
    cat("McElroy's R-squared value for the system: ")
    cat( formatC( x$mcelroy.r.squared, digits = digits, width = -1 ) )
    cat("\n")
  }

   if( equations ){
      ## now print the individual equations
      for(i in 1:length( x$eq ) ) {
         print( x$eq[[i]], digits = digits )
      }
   } else {
      cat( "\nCoefficients:\n" )
      print( coef( x ) ) 
   }

  invisible( x )
}


## prepare summary results for a single equation
summary.systemfit.equation <- function( object, useDfSys = NULL, ... ) {

   if( is.null( useDfSys ) ) {
      useDfSys <- object$nCoef.sys != object$rank.sys
         # TRUE if there are restrictions imposed
   }

   # number of observations in this equation
   nObs <- length( residuals( object ) )

   # preparing objects that will be returned
   result <- list()
   result$eqnLabel <- object$eqnLabel
   result$eqnNo <- object$eqnNo
   result$terms <- object$terms
   result$instruments <- object$inst
   result$method <- object$method
   result$residuals <- object$residuals

   # coefficients, standard errors, ... 
   result$coefCov <- object$coefCov
   coef <- object$coefficients
   stdEr <- diag( result$coefCov )^0.5  # standard errors
   tStat <- coef / stdEr                # t-statistic
   if( useDfSys ) {             # p-values
      pVal <- 2 * ( 1 - pt( abs( tStat ), object$df.residual.sys ) )
   } else {
      pVal <- 2 * ( 1 - pt( abs( tStat ), object$df.residual ) )
   }
   result$coefficients <- cbind( coef, stdEr, tStat, pVal )
   colnames( result$coefficients ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   result$df <- c( length( coef( object ) ), nObs - length( coef( object ) ) )
   result$ssr <- sum( residuals( object )^2 )
   result$sigma <- sqrt( result$ssr / df.residual( object ) )

   # R^2 values
   response <- fitted( object ) + residuals( object )
   rss <- sum( residuals( object )^2 )
   tss <- sum( ( response - mean( response ) )^2 )
   result$r.squared <- 1 - rss / tss
   result$adj.r.squared <- 1 - ( ( nObs - 1 ) / object$df.residual ) *
      ( 1 - result$r.squared )
   class( result ) <- "summary.systemfit.equation"
   return( result )
}


## print summary results for a single equation
print.summary.systemfit.equation <- function( x,
   digits = max( 3, getOption("digits") - 1 ), ... ) {

  cat("\n")
  cat( x$method, " estimates for '", x$eqnLabel,
         "' (equation ", x$eqnNo, ")\n", sep = "" )

  cat("Model Formula: ")
  print( formula( x$terms ) )
  if(!is.null(x$inst)) {
    cat("Instruments: ")
    print(x$inst)
  }
  cat("\n")

  Signif <- symnum(x$coefficients[,4], corr = FALSE, na = FALSE,
                   cutpoints = c( 0, 0.001, 0.01, 0.05, 0.1, 1 ),
                   symbols   = c( "***", "**", "*", "." ," " ))

  table <- cbind(round( x$coefficients,  digits ),
                 Signif)

  rownames(table) <- rownames(x$coefficients)
  colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE, digits = digits )
  cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")

  cat(paste("\nResidual standard error:", round( x$sigma, digits ),
            "on", x$df[ 2 ], "degrees of freedom\n" ))

  cat( paste( "Number of observations:", round( sum( x$df ), digits ),
              "Degrees of Freedom:", round( x$df[ 2 ], digits ),"\n" ) )

  cat( paste( "SSR:", round( x$ssr, digits ),
              "MSE:", round( x$sigma^2, digits ),
              "Root MSE:", round(x$sigma, digits), "\n" ) )

  cat( paste( "Multiple R-Squared:", round( x$r.squared, digits ),
              "Adjusted R-Squared:", round( x$adj.r.squared, digits ),
              "\n" ) )
  cat("\n")
  invisible( x )
}
