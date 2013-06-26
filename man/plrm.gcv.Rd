\name{plrm.gcv}
\alias{plrm.gcv}

\title{
Generalized cross-validation bandwidth selection in PLR models
}
\description{
From a sample \eqn{{(Y_i, X_{i1}, ..., X_{ip}, t_i): i=1,...,n}}, this routine computes an optimal pair of bandwidths for estimating the regression function of the model
\deqn{Y_i= X_{i1}*\beta_1 +...+ X_{ip}*\beta_p + m(t_i) + \epsilon_i,}
where
\deqn{\beta = (\beta_1,...,\beta_p)\;is\;an\;unknown\;vector\;parameter}
and
\deqn{m(.)\;is\;a\;smooth\;but\;unknown\;function.}
The optimal pair of bandwidths, \code{(b.opt, h.opt)}, is selected by means of the generalized cross-validation procedure. The bandwidth \code{b.opt} is used in the estimate of \eqn{\beta}, while the pair of bandwidths \code{(b.opt, h.opt)} is considered in the estimate of \eqn{m}. Kernel smoothing, combined with ordinary least squares estimation, is used.
}
\usage{
plrm.gcv(data = data, b.equal.h = TRUE, num.b = NULL, num.h = NULL, 
estimator = "NW", kernel = "quadratic")
}
\arguments{
\item{data}{
\code{data[,1]} contains the values of the response variable, \eqn{Y}; 

\code{data[, 2:(p+1)]} contains the values of the "linear" explanatory variables,

\eqn{X_1, ..., X_p};

\code{data[, p+2]} contains the values of the "nonparametric" explanatory variable, \eqn{t}.
}
  \item{b.equal.h}{if TRUE (the default), the same bandwidth is used for estimating both \eqn{\beta} and \eqn{m}.}
  \item{num.b}{number of values used to build the sequence of considered bandwidths for estimating \eqn{\beta}. If both \code{num.b} and \code{num.h} are \code{NULL} (the default), \code{num.b=50} is considered; if \code{num.b} is \code{NULL} (the default) but \code{num.h} is not \code{NULL}, \code{num.b=num.h} is considered; if \code{b.equal.h=TRUE} (the default) and both \code{num.b} and \code{num.h} are not \code{NULL} and different, the maximum value of \code{num.b} and \code{num.h} is considered for both.}
  \item{num.h}{pairs of bandwidths (\code{b, h}) are used for estimating \eqn{m}, \code{num.h} being the number of values considered for \code{h}. If both \code{num.b} and \code{num.h} are \code{NULL} (the default), \code{num.h=50} is considered; if \code{num.h} is \code{NULL} (the default) but \code{num.b} is not \code{NULL}, \code{num.h=num.b} is considered; if \code{b.equal.h=TRUE} (the default) and both \code{num.b} and \code{num.h} are not \code{NULL} and different, the maximum value of \code{num.b} and \code{num.h} is considered for both.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}
}
\details{
The implemented procedure generalizes that one in page 423 in Speckman (1988) by allowing two smoothing parameters instead of only one (see Aneiros-Perez \emph{et al.}, 2004).
}
\value{
\item{bh.opt}{selected value for \code{(b,h)}.}
\item{GCV.opt}{minimum value of the GCV function.}
\item{GCV}{matrix containing the values of the GCV function for each pair of bandwidths considered.}
\item{b.seq}{sequence of considered bandwidths, \code{b}, in the GCV function for estimating \eqn{\beta}: \code{num.b} equidistant values between zero and a quarter of the range of \eqn{{t_i}}.}
\item{h.seq}{sequence of considered bandwidths, \code{h}, in the pair of bandwidths \code{(b, h)} used in the GCV function for estimating \eqn{m}: \code{num.h} equidistant values between zero and a quarter of the range of \eqn{{t_i}}.}
}
\references{
Aneiros-Perez, G., Gonzalez-Manteiga, W. and Vieu, P. (2004) Estimation and testing in a partial linear regression under long-memory dependence. \emph{Bernoulli} \bold{10}, 49-78.

Green, P. (1985) Linear models for field trials, smoothing and cross-validation. \emph{Biometrika} \bold{72}, 527-537.

Speckman, P. (1988) Kernel smoothing in partial linear models \emph{J. R. Statist. Soc. B} \bold{50}, 413-436.
}
\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are: \code{\link{plrm.beta}}, \code{\link{plrm.est}}, \code{\link{plrm.cv}}, \code{\link{np.est}}, \code{\link{np.gcv}} and \code{\link{np.cv}}.
}
\examples{
### Example: independent data

set.seed(1234)

# We generate the data
n <- 100
t <- ((1:n)-0.5)/n
beta <- c(0.05, 0.01)
m <- function(t) {0.25*t*(1-t)}
f <- m(t)

x <- matrix(rnorm(200,0,1), nrow=n)
sum <- x\%*\%beta
epsilon <- rnorm(n, 0, 0.01)
y <-  sum + f + epsilon
data_ind <- matrix(c(y,x,t),nrow=100)

# We obtain the optimal bandwidths
a <-plrm.gcv(data_ind)
a$GCV.opt

GCV <- a$GCV
h <- a$h.seq
plot(h, GCV,type="l")

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}

