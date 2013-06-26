\name{plrm.beta}
\alias{plrm.beta}

\title{
Semiparametric estimate for the parametric component of the regression function in PLR models
}
\description{
This routine computes estimates for \eqn{\beta} from a sample \eqn{{(Y_i, X_{i1}, ..., X_{ip}, t_i): i=1,...,n}}, where:
\deqn{\beta = (\beta_1,...,\beta_p)\;is\;an\;unknown\;vector\;parameter}
and
\deqn{Y_i= X_{i1}*\beta_1 +...+ X_{ip}*\beta_p + m(t_i) + \epsilon_i.}
The nonparametric component, \eqn{m}, is a smooth but unknown function, and the random errors, \eqn{{\epsilon_i}}, are allowed to be time series. Ordinary least squares estimation, combined with kernel smoothing, is used.
}
\usage{
plrm.beta(data = data, b.seq = NULL, estimator = "NW", kernel = "quadratic")
}
\arguments{
    \item{data}{
\code{data[, 1]} contains the values of the response variable, \eqn{Y};

\code{data[, 2:(p+1)]} contains the values of the "linear" explanatory variables,

\eqn{X_1, ..., X_p};

\code{data[, p+2]} contains the values of the "nonparametric" explanatory variable, \eqn{t}.
}
  \item{b.seq}{vector of bandwidths for estimating \eqn{\beta}. If \code{NULL} (the default), only one estimate of \eqn{\beta} is computed, the corresponding bandwidth being selected by means of the generalized cross-validation procedure.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}
}
\details{
The expression for the estimator of \eqn{\beta} can be seen in page 52 in Aneiros-Perez \emph{et al.} (2004).
}
\value{
\item{BETA}{\code{p x length(b.seq)} matrix containing the estimate of \eqn{\beta} for each bandwidth in \code{h.seq}.}
\item{G}{\code{n x p x length(b.seq)} array containing the nonparametric estimate of

\code{E(X_{ij} | t_i)} (\eqn{i=1,...,n; j=1,...,p}) for each bandwidth in \code{b.seq}.}
}
\references{
Aneiros-Perez, G., Gonzalez-Manteiga, W. and Vieu, P. (2004) Estimation and testing in a partial linear regression model under long memory dependence. \emph{Bernoulli} \bold{10}, 49-78.

Hardle, W., Liang, H. and Gao, J. (2000) \emph{Partially Linear Models}. Physica-Verlag.

Speckman, P. (1988) Kernel smoothing in partial linear models. \emph{J. R. Statist. Soc. B} \bold{50}, 413-436.
}
\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are: \code{\link{plrm.est}}, \code{\link{plrm.gcv}}, \code{\link{plrm.cv}}.
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

# We estimate the parametric component of the PLR model
# (GCV bandwidth)
a <- plrm.beta(data_ind)

a$BETA


## Example 2: dependent data
# We generate the data
x <- matrix(rnorm(200,0,1), nrow=n)
sum <- x\%*\%beta
epsilon <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y <-  sum + f + epsilon
data_dep <- matrix(c(y,x,t),nrow=100)


# We estimate the parametric component of the PLR model
# (CV bandwidth)
b <- plrm.cv(data_dep, ln.0=2)$bh.opt[2,1]
a <-plrm.beta(data_dep, b=b)

a$BETA

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}

