\name{np.gcv}
\alias{np.gcv}

\title{
Generalized cross-validation bandwidth selection in nonparametric regression models
}
\description{
From a sample \eqn{{(Y_i, t_i): i=1,...,n}}, this routine computes an optimal bandwidth for estimating \eqn{m} in the regression model
\deqn{Y_i= m(t_i) + \epsilon_i.}
The regression function, \eqn{m}, is a smooth but unknown function. The optimal bandwidth is selected by means of the generalized cross-validation procedure. Kernel smoothing is used.
}
\usage{
np.gcv(data = data, h.seq=NULL, num.h = 50, estimator = "NW", 
kernel = "quadratic")
}
\arguments{
 \item{data}{
\code{data[, 1]} contains the values of the response variable, \eqn{Y};

\code{data[, 2]} contains the values of the explanatory variable, \eqn{t}.
}
  \item{h.seq}{sequence of considered bandwidths in the GCV function. If \code{NULL} (the default), \code{num.h} equidistant values between zero and a quarter of the range of \eqn{t_i} are considered.}
  \item{num.h}{number of values used to build the sequence of considered bandwidths. If \code{h.seq} is not \code{NULL}, \code{num.h=length(h.seq)}. Otherwise, the default is 50.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}
}

\details{
See Craven and Wahba (1979) and Rice (1984).
}
\value{
\item{h.opt}{selected value for the bandwidth.}
\item{GCV.opt}{minimum value of the GCV function.}
\item{GCV}{vector containing the values of the GCV function for each considered bandwidth.}
\item{h.seq}{sequence of considered bandwidths in the GCV function.}
}
\references{
Craven, P. and Wahba, G. (1979) Smoothing noisy data with spline functions. \emph{Numer. Math.} \bold{31}, 377-403.

Rice, J. (1984) Bandwidth choice for nonparametric regression. \emph{Ann. Statist.} \bold{12}, 1215-1230.
}
\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are: \code{\link{np.est}}, \code{\link{np.cv}}, \code{\link{plrm.est}}, \code{\link{plrm.gcv}} and \code{\link{plrm.cv}}.
}
\examples{
# EXAMPLE 1: REAL DATA
data <- matrix(10,120,2)
data(barnacles1)
barnacles1 <- as.matrix(barnacles1)
data[,1] <- barnacles1[,1]
data <- diff(data, 12)
data[,2] <- 1:nrow(data)

aux <- np.gcv(data)
aux$h.opt
plot(aux$h.seq, aux$GCV, xlab="h", ylab="GCV", type="l")



# EXAMPLE 2: SIMULATED DATA
## Example 2a: independent data

set.seed(1234)
# We generate the data
n <- 100
t <- ((1:n)-0.5)/n
m <- function(t) {0.25*t*(1-t)}
f <- m(t)

epsilon <- rnorm(n, 0, 0.01)
y <-  f + epsilon
data_ind <- matrix(c(y,t),nrow=100)

# We apply the function
a <-np.gcv(data_ind)
a$GCV.opt

GCV <- a$GCV
h <- a$h.seq
plot(h, GCV, type="l")

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}

