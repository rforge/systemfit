\name{estfun.systemfit}
\alias{estfun.systemfit}

\title{Extract Gradients of the Objective Function at each Observation}

\description{
   Extract the gradients of the objective function 
   with respect to the coefficients 
   evaluated at each observation 
   (\sQuote{Empirical Estimating Function}, 
   see \code{\link[sandwich]{estfun}}).
}

\usage{
\method{estfun}{systemfit}( obj, ... )
}

\arguments{
   \item{obj}{an object of class \code{systemfit}.}
   \item{\dots}{further arguments (currently ignored).}
}

\value{
   Matrix of gradients of the objective function
   with respect to the coefficients
   evaluated at each observation.
}

\section{Warnings}{
   The \pkg{sandwich} package must be loaded before this method
   can be used.
}

\author{
   Arne Henningsen
}

\seealso{\code{\link[sandwich]{estfun}}, \code{\link{systemfit}}.}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## OLS estimation
fitols <- systemfit( system, "OLS", data = Kmenta )

## obtain the estimation function
library( "sandwich" )
estfun( fitols )

## this is only true for OLS and 2SLS models
all.equal( estfun( fitols ),
   unlist( residuals( fitols ) ) * model.matrix( fitols ) )

# each column should sum up to (approximately) zero
colSums( estfun( fitols ) )


## iterated SUR estimation
fitsur <- systemfit( system, "SUR", data = Kmenta, maxit = 100 )

## obtain the estimation function
library( "sandwich" )
estfun( fitsur )

## this should be true for SUR, 3SLS, WLS, and W2SLS models
all.equal( estfun( fitsur ),
   unlist( residuals( fitsur ) ) * 
   ( ( solve( fitsur$residCovEst ) \%x\% diag( nrow( Kmenta ) ) ) \%*\% 
      model.matrix( fitsur ) ), check.attributes = FALSE )

# each column should sum up to (approximately) zero
colSums( estfun( fitsur ) )
}

\keyword{methods}