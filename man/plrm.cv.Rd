\name{plrm.cv}
\alias{plrm.cv}

\title{
Cross-validation bandwidth selection in PLR models
}
\description{
From a sample \eqn{{(Y_i, X_{i1}, ..., X_{ip}, t_i): i=1,...,n}}, this routine computes, for each \eqn{l_n} considered, an optimal pair of bandwidths for estimating the regression function of the model
\deqn{Y_i= X_{i1}*\beta_1 +...+ X_{ip}*\beta_p + m(t_i) + \epsilon_i,}
where
\deqn{\beta = (\beta_1,...,\beta_p)\;is\;an\;unknown\;vector\;parameter}
and
\deqn{m(.)\;is\;a\;smooth\;but\;unknown\;function.}
The random errors, \eqn{{\epsilon_i}}, are allowed to be time series. The optimal pair of bandwidths, (\code{b.opt, h.opt}), is selected by means of the leave-(\eqn{2l_n + 1})-out cross-validation procedure. The bandwidth \code{b.opt} is used in the estimate of \eqn{\beta}, while the pair of bandwidths \code{(b.opt, h.opt)} is considered in the estimate of \eqn{m}. Kernel smoothing, combined with ordinary least squares estimation, is used.
}
\usage{
plrm.cv(data = data, b.equal.h = TRUE, num.b = NULL, num.h = NULL, w = NULL, 
num.ln = 1, ln.0 = 0, step.ln = 2, estimator = "NW", kernel = "quadratic")
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
\item{w}{support interval of the weigth function in the CV function. If \code{NULL} (the default), \eqn{(q_{0.1}, q_{0.9})} is considered, where \eqn{q_p} denotes the quantile of order \eqn{p} of \eqn{{t_i}}.}
  \item{num.ln}{number of values for \eqn{l_n}: after estimating \eqn{\beta}, \eqn{2l_{n} + 1} observations around each point \eqn{t_i} are eliminated to estimate \eqn{m(t_i)} in the CV function. The default is 1.}
  \item{ln.0}{minimum value for \eqn{l_n}. The default is 0.}
  \item{step.ln}{distance between two consecutives values of \eqn{l_n}. The default is 2.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}
}

\details{
A weight function (specifically, the indicator function \bold{1}\eqn{_{[w[1] , w[2]]}}) is introduced in the CV function to allow elimination (or at least significant reduction) of boundary effects from the estimate of \eqn{m(t_i)}.

As noted in the definition of \code{num.ln}, the estimate of \eqn{\beta} in the CV function is obtained from all data while, once \eqn{\beta} is estimated, \eqn{2l_{n} + 1} observations around each \eqn{t_i} are eliminated to estimate \eqn{m(t_i)} in the CV function. Actually, the estimate of \eqn{\beta} to be used in time \eqn{i} in the CV function could be done eliminating such \eqn{2l_{n} + 1} observations too; that possibility was not implemented because both their computational cost and the known fact that the estimate of \eqn{\beta} is quite insensitive to the bandwidth selection.
  
The implemented procedure generalizes that one in expression (8) in Aneiros-Perez and Quintela-del-Rio (2001) by including a weight function (see above) and allowing two smoothing parameters instead of only one (see Aneiros-Perez \emph{et al.}, 2004).
}
\value{
\item{bh.opt}{dataframe containing, for each \code{ln} considered, the selected value for \code{(b,h)}.}
\item{CV.opt}{\code{CV.opt[k]} is the minimum value of the CV function when de k-th value of \code{ln} is considered.}
\item{CV}{an array containing the values of the CV function for each pair of bandwidths and \code{ln} considered.}
\item{b.seq}{sequence of considered bandwidths, \code{b}, in the CV function for estimating \eqn{\beta}: \code{num.b} equidistant values between zero and a quarter of the range of \eqn{{t_i}}.}
\item{h.seq}{sequence of considered bandwidths, \code{h}, in the pair of bandwidths \code{(b, h)} used in the CV function for estimating \eqn{m}: \code{num.h} equidistant values between zero and a quarter of the range of \eqn{{t_i}}.}
\item{w}{support interval of the weigth function in the CV function.}
}
\references{
Aneiros-Perez, G., Gonzalez-Manteiga, W. and Vieu, P. (2004) Estimation and testing in a partial linear regression under long-memory dependence. \emph{Bernoulli} \bold{10}, 49-78.

Aneiros-Perez, G. and Quintela-del-Rio, A. (2001) Modified cross-validation in semiparametric regression models with dependent errors. \emph{Comm. Statist. Theory Methods} \bold{30}, 289-307.

Chu, C-K and Marron, J.S. (1991) Comparison of two bandwidth selectors with dependent errors. \emph{The Annals of Statistics} \bold{19}, 1906-1918.
}
\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are: \code{\link{plrm.beta}}, \code{\link{plrm.est}}, \code{\link{plrm.gcv}}, \code{\link{np.est}}, \code{\link{np.gcv}} and \code{\link{np.cv}}.
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

# We apply the function
a <-plrm.cv(data_ind)
a$CV.opt

CV <- a$CV
h <- a$h.seq
plot(h, CV,type="l")


## Example 2: dependent data and ln.0 > 0
# We generate the data
x <- matrix(rnorm(200,0,1), nrow=n)
sum <- x\%*\%beta
epsilon <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y <-  sum + f + epsilon
data_dep <- matrix(c(y,x,t),nrow=100)

# We apply the function
a <-plrm.cv(data_dep, ln.0=2)
a$CV.opt

CV <- a$CV
h <- a$h.seq
plot(h, CV,type="l")

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}

