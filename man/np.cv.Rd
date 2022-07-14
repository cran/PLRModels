\name{np.cv}
\alias{np.cv}

\title{
Cross-validation bandwidth selection in nonparametric regression models
}
\description{
From a sample \eqn{{(Y_i, t_i): i=1,...,n}}, this routine computes, for each \eqn{l_n} considered, an optimal bandwidth for estimating \eqn{m} in the regression model
\deqn{Y_i= m(t_i) + \epsilon_i.}
The regression function, \eqn{m}, is a smooth but unknown function, and the random errors, \eqn{{\epsilon_i}}, are allowed to be time series. The optimal bandwidth is selected by means of the leave-(\eqn{2l_n + 1})-out cross-validation procedure. Kernel smoothing is used.
}
\usage{
np.cv(data = data, h.seq = NULL, num.h = 50, w = NULL, num.ln = 1, 
ln.0 = 0, step.ln = 2, estimator = "NW", kernel = "quadratic")
}
\arguments{
  \item{data}{
\code{data[, 1]} contains the values of the response variable, \eqn{Y};

\code{data[, 2]} contains the values of the explanatory variable, \eqn{t}.
}
  \item{h.seq}{sequence of considered bandwidths in the CV function. If \code{NULL} (the default), \code{num.h} equidistant values between zero and a quarter of the range of \eqn{t_i} are considered.}
  \item{num.h}{number of values used to build the sequence of considered bandwidths. If \code{h.seq} is not \code{NULL}, \code{num.h=length(h.seq)}. Otherwise, the default is 50.}
\item{w}{support interval of the weigth function in the CV function. If \code{NULL} (the default), \eqn{(q_{0.1}, q_{0.9})} is considered, where \eqn{q_p} denotes the quantile of order \eqn{p} of \eqn{{t_i}}.}
  \item{num.ln}{number of values for \eqn{l_n}: \eqn{2l_{n} + 1} observations around each point \eqn{t_i} are eliminated to estimate \eqn{m(t_i)} in the CV function. The default is 1.}
  \item{ln.0}{minimum value for \eqn{l_n}. The default is 0.}
  \item{step.ln}{distance between two consecutives values of \eqn{l_n}. The default is 2.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}
}
\details{
A weight function (specifically, the indicator function \bold{1}\eqn{_{[w[1] , w[2]]}}) is introduced in the CV function to allow elimination (or at least significant reduction) of boundary effects from the estimate of \eqn{m(t_i)}.

For more details, see Chu and Marron (1991).
}
\value{
\item{h.opt}{dataframe containing, for each \code{ln} considered, the selected value for the bandwidth.}
\item{CV.opt}{\code{CV.opt[k]} is the minimum value of the CV function when de k-th value of \code{ln} is considered.}
\item{CV}{matrix containing the values of the CV function for each bandwidth and \code{ln} considered.}
\item{w}{support interval of the weigth function in the CV function.}
\item{h.seq}{sequence of considered bandwidths in the CV function.}
}
\references{
Chu, C-K and Marron, J.S. (1991) Comparison of two bandwidth selectors with dependent errors. \emph{The Annals of Statistics} \bold{19}, 1906-1918.
}
\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are: \code{\link{np.est}}, \code{\link{np.gcv}}, \code{\link{plrm.est}}, \code{\link{plrm.gcv}} and \code{\link{plrm.cv}}.
}
\examples{
# EXAMPLE 1: REAL DATA
data <- matrix(10,120,2)
data(barnacles1)
barnacles1 <- as.matrix(barnacles1)
data[,1] <- barnacles1[,1]
data <- diff(data, 12)
data[,2] <- 1:nrow(data)

aux <- np.cv(data, ln.0=1,step.ln=1, num.ln=2)
aux$h.opt
plot.ts(aux$CV)

par(mfrow=c(2,1))
plot(aux$h.seq,aux$CV[,1], xlab="h", ylab="CV", type="l", main="ln=1")
plot(aux$h.seq,aux$CV[,2], xlab="h", ylab="CV", type="l", main="ln=2")



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
a <-np.cv(data_ind)
a$CV.opt

CV <- a$CV
h <- a$h.seq
plot(h,CV,type="l")

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}

