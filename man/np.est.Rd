\name{np.est}
\alias{np.est}

\title{Nonparametric estimate of the regression function}
\description{
This routine computes estimates for \eqn{m(newt_j)} (\eqn{j=1,...,J}) from a sample \eqn{{(Y_i, t_i): i=1,...,n}}, where:
\deqn{Y_i= m(t_i) + \epsilon_i.}
The regression function, \eqn{m}, is a smooth but unknown function, and the random errors, \eqn{{\epsilon_i}}, are allowed to be time series. Kernel smoothing is used.
}
\usage{
np.est(data = data, h.seq = NULL, newt = NULL, estimator = "NW", 
kernel = "quadratic")
}
\arguments{
 \item{data}{
\code{data[, 1]} contains the values of the response variable, \eqn{Y};

\code{data[, 2]} contains the values of the explanatory variable, \eqn{t}.
}
  \item{h.seq}{the considered bandwidths. If \code{NULL} (the default), only one bandwidth, selected by means of the generalized cross-validation procedure, is used.}
  \item{newt}{values of the explanatory variable where the estimates are obtained. If NULL (the default), the considered values will be the values of \code{data[,2]}.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}
}
\details{
See Fan and Gijbels (1996) and Francisco-Fernandez and Vilar-Fernandez (2001).
}

\value{
YHAT: a length(\code{newt}) x length(\code{h.seq}) matrix containing the estimates for \eqn{m(newt_j)} 

(\eqn{j=1,...,}length(\code{newt})) using the different bandwidths in \code{h.seq}.
}

\references{
Fan, J.  and Gijbels, I. (1996) \emph{Local Polynomial Modelling and its Applications}. Chapman and Hall, London.

Francisco-Fernandez, M. and Vilar-Fernandez, J. M. (2001) Local polynomial regression estimation with correlated errors. \emph{Comm. Statist. Theory Methods} \bold{30}, 1271-1293.
}

\author{German Aneiros Perez \email{ganeiros@udc.es} 

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are: \code{\link{np.gcv}}, \code{\link{np.cv}}, \code{\link{plrm.est}}, \code{\link{plrm.gcv}} and \code{\link{plrm.cv}}.
}
\examples{
set.seed(1234)

## Example: independent data
# We generate the data
n <- 100
t <- ((1:n)-0.5)/n
m <- function(t) {0.25*t*(1-t)}
f <- m(t)

epsilon <- rnorm(n, 0, 0.01)
y <-  f + epsilon
data_ind <- matrix(c(y,t),nrow=100)

# We estimate the nonparametric component of the PLR model
# (GCV bandwidth)
est <- np.est(data_ind)
plot(t, est, type="l", lty=2, ylab="")
points(t, 0.25*t*(1-t), type="l")
legend(x="topleft", legend = c("m", "m hat"), col=c("black", "black"), lty=c(1,2))

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}

