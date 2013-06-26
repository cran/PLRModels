\name{par.est}
\alias{par.est}

\title{
Estimation in linear regression models
}
\description{
This routine computes the ordinary least squares estimate for \eqn{\beta} from a sample \eqn{(Y_i, X_{i1},...,X_{ip}):}

\eqn{i=1,...,n}, where:
\deqn{\beta = (\beta_1,...,\beta_p)\;is\;an\;unknown\;vector\;parameter}
and
\deqn{Y_i = X_{i1}*\beta_1+ ... + X_{ip}*\beta_p + \epsilon_i.}
The random errors, \eqn{{\epsilon_i}}, are allowed to be time series.
}
\usage{
par.est(data = data)
}
\arguments{
   \item{data}{
\code{data[, 1]} contains the values of the response variable, \eqn{Y};

\code{data[, 2:(p+1)]} contains the values of the explanatory variables, \eqn{X_1, ..., X_p}.
}
}
\details{
See Seber (1977) and Judge \emph{et al.} (1980).
}
\value{A vector containing the corresponding estimate.
}
\references{
Judge, G.G., Griffiths, W.E., Carter Hill, R., Lutkepohl, H. and Lee, T-C. (1980) \emph{The Theory and Practice of Econometrics}. Wiley.

Seber, G.A.F. (1977) \emph{Linear Regression Analysis}. Wiley.
}

\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are \code{\link{plrm.beta}} and \code{\link{plrm.est}}.
}
\examples{
set.seed(1234)

## Example 1: independent data
# We generate the data
n <- 100
beta <- c(0.05, 0.01)

x <- matrix(rnorm(200,0,1), nrow=n)
sum <- x\%*\%beta
epsilon <- rnorm(n, 0, 0.01)
y <-  sum + epsilon
data_ind <- matrix(c(y,x),nrow=100)

# We estimate the parametric component of the PLR model
par.est(data_ind)


## Example 2: dependent data
# We generate the data
x <- matrix(rnorm(200,0,1), nrow=n)
sum <- x\%*\%beta
epsilon <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y <-  sum + epsilon
data_dep <- matrix(c(y,x),nrow=100)

# We estimate the parametric component of the PLR model
par.est(data_dep)

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}


