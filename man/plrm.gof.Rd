\name{plrm.gof}
\alias{plrm.gof}

\title{
Goodness-of-Fit tests in PLR models
}
\description{
From a sample \eqn{{(Y_i, X_{i1}, ..., X_{ip}, t_i): i=1,...,n}}, this routine tests the null hypotheses \eqn{H_0: \beta=\beta_0} and \eqn{H_0: m=m_0}, where:
\deqn{\beta = (\beta_1,...,\beta_p)}
is an unknown vector parameter,
\deqn{m(.)}
is a smooth but unknown function and
\deqn{Y_i= X_{i1}*\beta_1 +...+ X_{ip}*\beta_p + m(t_i) + \epsilon_i.}
Fixed equally spaced design is considered for the "nonparametric" explanatory variable, \eqn{t}, and the random errors, \eqn{\epsilon_i}, are allowed to be time series. The test statistic used for testing \eqn{H0: \beta = \beta_0} derives from the asymptotic normality of an estimator of \eqn{\beta} based on both ordinary least squares and kernel smoothing (this result giving a \eqn{\chi^2}-test). The test statistic used for testing \eqn{H0: m = m_0} derives from a Cramer-von-Mises-type functional distance between a nonparametric estimator of \eqn{m} and \eqn{m_0}.
}
\usage{
plrm.gof(data = data, beta0 = NULL, m0 = NULL, b.seq = NULL, 
h.seq = NULL, w = NULL, estimator = "NW", kernel = "quadratic", 
time.series = FALSE, Var.Cov.eps = NULL, Tau.eps = NULL, 
b0 = NULL, h0 = NULL, lag.max = 50, p.max = 3, q.max = 3, 
ic = "BIC", num.lb = 10, alpha = 0.05)
}
\arguments{
  \item{data}{
\code{data[,1]} contains the values of the response variable, \eqn{Y};

\code{data[, 2:(p+1)]} contains the values of the "linear" explanatory variables

\eqn{X_1, ... X_p};

\code{data[, p+2]} contains the values of the "nonparametric" explanatory variable, \eqn{t}.}
    \item{beta0}{the considered parameter vector in the parametric null hypothesis. If \code{NULL} (the default), the zero vector is considered.}
  \item{m0}{the considered function in the nonparametric null hypothesis. If \code{NULL} (the default), the zero function is considered.}
   \item{b.seq}{the statistic test for \eqn{H0: \beta=\beta_0} is performed using each bandwidth in the vector \code{b.seq}. If \code{NULL} (the default) but \code{h.seq} is not \code{NULL}, it takes \code{b.seq=h.seq}. If both \code{b.seq} and \code{h.seq} are \code{NULL}, 10 equidistant values between zero and a quarter of the range of \eqn{{t_i}} are considered.}
   \item{h.seq}{the statistic test for \eqn{H0: m=m_0} is performed using each pair of bandwidths (\code{b.seq[j], h.seq[j]}). If \code{NULL} (the default) but \code{b.seq} is not \code{NULL}, it takes \code{h.seq=b.seq}. If both \code{b.seq} and \code{h.seq} are \code{NULL}, 10 equidistant values between zero and a quarter of the range of \eqn{{t_i}} are considered for both \code{b.seq} and  \code{h.seq}.} 
     \item{w}{support interval of the weigth function in the test statistic for \eqn{H0: m = m_0}. If \code{NULL} (the default), \eqn{(q_{0.1}, q_{0.9})} is considered, where \eqn{q_p} denotes the quantile of order \eqn{p} of \eqn{{t_i}}.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}
  \item{time.series}{it denotes whether the data are independent (FALSE) or if data is a time series (TRUE). The default is FALSE.}
    \item{Var.Cov.eps}{ \code{n x n} matrix of variances-covariances associated to the random errors of the regression model. If NULL (the default), the function tries to estimate it: it fits an ARMA model (selected according to an information criterium) to the residuals from the fitted regression model and, then, it obtains the var-cov matrix of such ARMA model.}
  \item{Tau.eps}{it contains the sum of autocovariances associated to the random errors of the regression model. If NULL (the default), the function tries to estimate it: it fits an ARMA model (selected according to an information criterium) to the residuals from the fitted regression model and, then, it obtains the sum of the autocovariances of such ARMA model.}
    \item{b0}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, \code{b0} contains the pilot bandwidth for the estimator of \eqn{\beta} used for obtaining the residuals to construct the default for \code{Var.Cov.eps} and/or \code{Tau.eps}. If \code{NULL} (the default) but \code{h0} is not \code{NULL}, it takes \code{b0=h0}. If both \code{b0} and \code{h0} are \code{NULL}, a quarter of the range of \eqn{{t_i}} is considered.}
  \item{h0}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, (\code{b0, h0}) contains the pair of pilot bandwidths  for the estimator of \eqn{m} used for obtaining the residuals to construct the default for \code{Var.Cov.eps} and/or \code{Tau.eps}. If \code{NULL} (the default) but \code{b0} is not \code{NULL}, it takes \code{h0=b0}. If both \code{b0} and \code{h0} are \code{NULL}, a quarter of the range of \eqn{{t_i}} is considered for both \code{b0} and  \code{h0}.}
  \item{lag.max}{if \code{Tau.eps=NULL}, \code{lag.max} contains the maximum delay used to construct the default for \code{Tau.eps}. The default is 50.}
  \item{p.max}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, the ARMA model is selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{q.max}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, the ARMA model is selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{ic}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, \code{ic} contains the information criterion used to suggest the ARMA model. It allows us to choose between: "AIC", "AICC" or "BIC" (the default).}
  \item{num.lb}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, it checks the suitability of the selected ARMA model according to the Ljung-Box test and the t-test. It uses up to \code{num.lb} delays in the Ljung-Box test. The default is 10.}
  \item{alpha}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, \code{alpha} contains the significance level which the ARMA model is checked. The default is 0.05.}
}
\details{
A weight function (specifically, the indicator function \bold{1}\eqn{_{[w[1] , w[2]]}}) is introduced in the test statistic for testing \eqn{H0: m = m_0} to allow elimination (or at least significant reduction) of boundary effects from the estimate of \eqn{m(t_i)}.

If \code{Var.Cov.eps=NULL} and the routine is not able to suggest an approximation for \code{Var.Cov.eps}, it warns the user with a message saying that the model could be not appropriate and then it shows the results. In order to construct \code{Var.Cov.eps}, the procedure suggested in Aneiros-Perez and Vieu (2013) can be followed.

If \code{Tau.eps=NULL} and the routine is not able to suggest an approximation for \code{Tau.eps}, it warns the user with a message saying that the model could be not appropriate and then it shows the results. In order to construct \code{Tau.eps}, the procedures suggested in Aneiros-Perez (2008) can be followed.

The implemented procedures generalize those ones in expressions (9) and (10) in Gonzalez-Manteiga and Aneiros-Perez (2003) by allowing some dependence condition in \eqn{{(X_{i1}, ..., X_{ip}): i=1,...,n}} and including a weight function (see above), respectively.
}
\value{A list with two dataframes:
  \item{parametric.test}{a dataframe containing the bandwidths, the statistics and the p-values when one tests \eqn{H0: \beta = \beta_0}}
  \item{nonparametric.test}{a dataframe containing the bandwidths b and h, the statistics, the normalised statistics and the p-values when one tests \eqn{H0: m = m_0}}
   Moreover, if \code{data} is a time series and \code{Tau.eps} or \code{Var.Cov.eps} are not especified:
  \item{pv.Box.test}{p-values of the Ljung-Box test for the model fitted to the residuals.}
  \item{pv.t.test}{p-values of the t.test for the model fitted to the residuals.}
  \item{ar.ma}{ARMA orders for the model fitted to the residuals.}
}
\references{
Aneiros-Perez, G. (2008) Semi-parametric analysis of covariance under dependence conditions within each group. \emph{Aust. N. Z. J. Stat.} \bold{50}, 97-123.

Aneiros-Perez, G., Gonzalez-Manteiga, W. and Vieu, P. (2004) Estimation and testing in a partial linear regression under long-memory dependence. \emph{Bernoulli} \bold{10}, 49-78.

Aneiros-Perez, G. and Vieu, P. (2013) Testing linearity in semi-parametric functional data analysis. \emph{Comput. Stat.} \bold{28}, 413-434. 

Gao, J. (1997) Adaptive parametric test in a semiparametric regression model. \emph{Comm. Statist. Theory Methods} \bold{26}, 787-800.

Gonzalez-Manteiga, W. and Aneiros-Perez, G. (2003) Testing in partial linear regression models with dependent errors. \emph{J. Nonparametr. Statist.} \bold{15}, 93-111.
}
\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are \code{\link{plrm.est}}, \code{\link{par.gof}} and \code{\link{np.gof}}.
}
\examples{
# EXAMPLE 1: REAL DATA
data(barnacles1)
data <- as.matrix(barnacles1)
data <- diff(data, 12)
data <- cbind(data,1:nrow(data))

plrm.gof(data)
plrm.gof(data, beta0=c(-0.1, 0.35))



# EXAMPLE 2: SIMULATED DATA
## Example 2a: dependent data

set.seed(1234)
# We generate the data
n <- 100
t <- ((1:n)-0.5)/n
beta <- c(0.05, 0.01)
m <- function(t) {0.25*t*(1-t)}
f <- m(t)
f.function <- function(u) {0.25*u*(1-u)}

x <- matrix(rnorm(200,0,1), nrow=n)
sum <- x\%*\%beta
epsilon <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y <-  sum + f + epsilon
data <- cbind(y,x,t)

## Example 2a.1: true null hypotheses
plrm.gof(data, beta0=c(0.05, 0.01), m0=f.function, time.series=TRUE)

## Example 2a.2: false null hypotheses
plrm.gof(data, time.series=TRUE)

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}

