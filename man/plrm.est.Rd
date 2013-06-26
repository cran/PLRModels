\name{plrm.est}
\alias{plrm.est}

\title{
Semiparametric estimates for the unknown components of the regression function in  PLR models
}
\description{
This routine computes estimates for \eqn{\beta} and \eqn{m(newt_j)} (\eqn{j=1,...,J}) from a sample 

\eqn{{(Y_i, X_{i1}, ..., X_{ip}, t_i): i=1,...,n}}, where:
\deqn{\beta = (\beta_1,...,\beta_p)\;is\;an\;unknown\;vector\;parameter,}
\deqn{m(.)\;is\;a\;smooth\;but\;unknown\;function}
and
\deqn{Y_i= X_{i1}*\beta_1 +...+ X_{ip}*\beta_p + m(t_i) + \epsilon_i.}
The random errors, \eqn{{\epsilon_i}}, are allowed to be time series. Kernel smoothing, combined with ordinary least squares estimation, is used.
}
\usage{
plrm.est(data = data, b = NULL, h = NULL, newt = NULL, estimator = "NW", 
kernel = "quadratic")
}
\arguments{
  \item{data}{
\code{data[, 1]} contains the values of the response variable, \eqn{Y};

\code{data[, 2:(p+1)]} contains the values of the "linear" explanatory variables,

\eqn{X_1, ..., X_p};

\code{data[, p+2]} contains the values of the "nonparametric" explanatory variable, \eqn{t}.
}
   \item{b}{bandwidth for estimating the parametric part of the model. If both \code{b} and \code{h} are \code{NULL} (the default), it is selected by means of the generalized cross-validation procedure (fixing \code{b=h}); if \code{b} is \code{NULL} (the default) but \code{h} is not \code{NULL}, \code{b=h} is considered.}
  \item{h}{\code{(b,h)} is the pair of bandwidths for estimating the nonparametric part of the model. If both \code{b} and \code{h} are \code{NULL} (the default), it is selected by means of the generalized cross-validation procedure (fixing \code{b=h}); if \code{b} is \code{NULL} (the default) but \code{h} is not \code{NULL}, \code{b=h} is considered; if \code{h} is \code{NULL} (the default) but \code{b} is not \code{NULL}, \code{h=b} is considered.}
    \item{newt}{values of the "nonparametric" explanatory variable where the estimator of \eqn{m} is evaluated. If NULL (the default), the considered values will be the values of \code{data[,p+2]}.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}   
}
\details{
Expressions for the estimators of \eqn{\beta} and \eqn{m} can be seen in page 52 in Aneiros-Perez \emph{et al.} (2004).
}
\value{
\item{beta}{a vector containing the estimate of \eqn{\beta}.}
\item{m.t}{a vector containing the estimator of the non-parametric part, \eqn{m}, evaluated in the design points.}
\item{m.newt}{a vector containing the estimator of the non-parametric part, \eqn{m}, evaluated in \code{newt}.}
\item{residuals}{a vector containing the residuals: \code{Y - X*beta - m.t}.}
\item{fitted.values}{the values obtained from the expression: \code{X*beta + m.t}}
\item{b}{the considered bandwidth for estimating \eqn{\beta}.}
\item{h}{\code{(b,h)} is the pair of bandwidths considered for estimating \eqn{m}.}
}
\references{
Aneiros-Perez, G., Gonzalez-Manteiga, W. and Vieu, P. (2004) Estimation and testing in a partial linear regression under long-memory dependence. \emph{Bernoulli} \bold{10}, 49-78.

Hardle, W., Liang, H. and Gao, J. (2000) \emph{Partially Linear Models}. Physica-Verlag.

Speckman, P. (1988) Kernel smoothing in partial linear models. \emph{J. R. Statist. Soc. B} \bold{50}, 413-436.
}
\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are: \code{\link{plrm.beta}}, \code{\link{plrm.gcv}}, \code{\link{plrm.cv}}, \code{\link{np.est}}, \code{\link{np.gcv}} and \code{\link{np.cv}}.
}
\examples{
set.seed(1234)

## Example 1: independent data
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

# We estimate the components of the PLR model
# (GCV bandwidth)
a <- plrm.est(data_ind)

a$beta

est <- a$m.t
plot(t, est, type="l", lty=2, ylab="")
points(t, 0.25*t*(1-t), type="l")
legend(x="topleft", legend = c("m", "m hat"), col=c("black", "black"), lty=c(1,2))


## Example 2: dependent data
# We generate the data
x <- matrix(rnorm(200,0,1), nrow=n)
sum <- x\%*\%beta
epsilon <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y <-  sum + f + epsilon
data_dep <- matrix(c(y,x,t),nrow=100)

# We estimate the components of the PLR model
# (CV bandwidth)
h <- plrm.cv(data_dep, ln.0=2)$bh.opt[3,1]
a <- plrm.est(data_dep, h=h)

a$beta

est <- a$m.t
plot(t, est, type="l", lty=2, ylab="")
points(t, 0.25*t*(1-t), type="l")
legend(x="topleft", legend = c("m", "m hat"), col=c("black", "black"), lty=c(1,2))

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}


