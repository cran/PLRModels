\name{par.ci}
\alias{par.ci}

\title{
Confidence intervals estimation in linear regression models
}
\description{
This routine obtains a confidence interval for the value \eqn{a^T * \beta}, by asymptotic distribution and bootstrap, from a sample \eqn{(Y_i, X_{i1},...,X_{ip}):}
\eqn{i=1,...,n}, where:
\deqn{a = (a_1,...,a_p)^T}
is an unknown vector, 
\deqn{\beta = (\beta_1,...,\beta_p)^T}
is an unknown vector parameter and
\deqn{Y_i = X_{i1}*\beta_1+ ... + X_{ip}*\beta_p + \epsilon_i.}
The random errors, \eqn{\epsilon_i}, are allowed to be time series.
}
\usage{
par.ci(data=data, seed=123, CI="AD", B=1000, N=50, a=NULL, 
p.arima=NULL, q.arima=NULL, p.max=3, q.max=3, alpha=0.05, 
alpha2=0.05, num.lb=10, ic="BIC", Var.Cov.eps=NULL)
}
\arguments{
  \item{data}{
\code{data[, 1]} contains the values of the response variable, \eqn{Y};

\code{data[, 2:(p+1)]} contains the values of the explanatory variables, \eqn{X_1, ..., X_p}.
}
  \item{seed}{the considered seed.}
  \item{CI}{method to obtain the confidence interval. It allows us to choose between: \dQuote{AD} (asymptotic distribution), \dQuote{B} (bootstrap) or \dQuote{all} (both). The default is \dQuote{AD}.}
  \item{B}{number of bootstrap replications. The default is 1000.}
  \item{N}{Truncation parameter used in the finite approximation of the MA(infinite) expression of \eqn{\epsilon}.}
  \item{a}{Vector which, multiplied by \code{beta}, is used for obtaining the confidence interval of this result.}
  \item{p.arima}{the considered p to fit the model ARMA(p,q).}
  \item{q.arima}{the considered q to fit the model ARMA(p,q).}
  \item{p.max}{if \code{Var.Cov.eps=NULL}, the ARMA models are selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{q.max}{if \code{Var.Cov.eps=NULL}, the ARMA models are selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{alpha}{1 - \code{alpha} is the confidence level of the confidence interval. The default is 0.05.}
  \item{alpha2}{significance level used to check (if needed) the ARMA model fitted to the residuals. The default is 0.05.}
  \item{num.lb}{if \code{Var.Cov.eps=NULL}, it checks the suitability of the selected ARMA model according to the Ljung-Box test and the t-test. It uses up to \code{num.lb} delays in the Ljung-Box test. The default is 10.}
   \item{ic}{if \code{Var.Cov.eps=NULL}, \code{ic} contains the information criterion used to suggest the ARMA model. It allows us to choose between: "AIC", "AICC" or "BIC" (the default).}
  \item{Var.Cov.eps}{ \code{n x n} matrix of variances-covariances associated to the random errors of the regression model. If NULL (the default), the function tries to estimate it: it fits an ARMA model (selected according to an information criterium) to the residuals from the fitted regression model and, then, it obtains the var-cov matrix of such ARMA model.}
}
\value{A list containing:
  \item{Bootstrap}{a dataframe containing \code{ci_inf} and \code{ci_sup}, the confidence intervals using bootstrap and \code{p_opt} and \code{q_opt} (the orders for the ARMA model fitted to the residuals).}
  \item{AD}{a dataframe containing \code{ci_inf} and \code{ci_sup}, the confidence intervals using the asymptotic distribution and \code{p_opt} and \code{q_opt} (the orders for the ARMA model fitted to the residuals).}
  \item{pv.Box.test}{p-values of the Ljung-Box test for the model fitted to the residuals.}
  \item{pv.t.test}{p-values of the t.test for the model fitted to the residuals.}
}
\references{
Liang, H., Hardle, W., Sommerfeld, V. (2000) Bootstrap approximation in a partially linear regression model. \emph{Journal of Statistical Planning and Inference} \bold{91}, 413-426.

You, J., Zhou, X. (2005) Bootstrap of a semiparametric partially linear model with autoregressive errors. \emph{Statistica Sinica} \bold{15}, 117-133.
}
\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}
}
\seealso{
A related function is \code{\link{plrm.ci}}.
}
\examples{
# EXAMPLE 1: REAL DATA
data(barnacles1)
data <- as.matrix(barnacles1)
data <- diff(data, 12)
data <- cbind(data[,1],1,data[,-1])

\dontrun{par.ci(data, a=c(1,0,0), CI="all")}
\dontrun{par.ci(data, a=c(0,1,0), CI="all")}
\dontrun{par.ci(data, a=c(0,0,1), CI="all")}



# EXAMPLE 2: SIMULATED DATA
## Example 2a: dependent data

set.seed(123)
# We generate the data
n <- 100
beta <- c(0.5, 2)

x <- matrix(rnorm(200,0,3), nrow=n)
sum <- x\%*\%beta
sum <- as.matrix(sum)
eps <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.1, n = n)
eps <- as.matrix(eps)
y <-  sum + eps
data_parci <- cbind(y,x)

# We estimate the confidence interval of a^T * beta in the PLR model
\dontrun{par.ci(data, a=c(1,0), CI="all")}
\dontrun{par.ci(data, a=c(0,1), CI="all")}

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
