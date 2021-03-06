\name{predict.systemfit}
\alias{predict.systemfit}
\alias{predict.systemfit.equation}

\title{Predictions from System Estimation}

\description{
   Returns the predicted values, their standard errors and the
   confidence limits of prediction.
}

\usage{
\method{predict}{systemfit}( object, data = object$data,
                             se.fit = FALSE, se.pred = FALSE,
                             interval = "none", level=0.95, ... )

\method{predict}{systemfit.equation}( object, data, ... )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{data}{data frame in which to predict.}
   \item{se.fit}{return the standard error of the fitted values?}
   \item{se.pred}{return the standard error of prediction?}
   \item{interval}{Type of interval calculation ("none", "confidence"
      or "prediction")}
   \item{level}{Tolerance/confidence level.}
   \item{...}{additional optional arguments.}
}

\details{
   The variance of the fitted values
   (used to calculate the standard errors of the fitted values
   and the "confidence interval") is calculated by
   \eqn{Var[E[y^0]-\hat{y}^0]=x^0 \; Var[b] \; {x^0}'}\cr
   The variances of the predicted values
   (used to calculate the standard errors of the predicted values
   and the "prediction intervals") is calculated by
   \eqn{Var[y^0-\hat{y}^0]=\hat{\sigma}^2+x^0 \; Var[b] \; {x^0}'}
}

\value{
   \code{predict.systemfit} returns a dataframe that
   contains for each equation the predicted values
   (e.g. "eq1.pred") and if requested
   the standard errors of the fitted values (e.g. "eq1.se.fit"),
   the standard errors of the prediction (e.g. "eq1.se.pred"),
   and the lower (e.g. "eq1.lwr") and upper (e.g. "eq1.upr")
   limits of the confidence or prediction interval(s).

   \code{predict.systemfit.equation} returns a vector of the predicted
      values of a single equation.
}

\references{

  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Macmillan.

  Gujarati, D. N. (1995)
  \emph{Basic Econometrics, Third Edition}, McGraw-Hill.

  Kmenta, J. (1997)
  \emph{Elements of Econometrics, Second Edition}, University of
  Michigan Publishing.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}


\seealso{
   \code{\link{systemfit}}, \code{\link{predict}}
}

\examples{
data( "Kmenta" )
demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## OLS estimation
fitols <- systemfit("OLS", system, labels, data=Kmenta )

## predicted values and limits
predict( fitols )

## predicted values of the first equation
predict( fitols$eq[[1]] )

## predicted values of the second equation
predict( fitols$eq[[2]] )
}

\keyword{models}



