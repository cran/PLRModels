\name{par.gof}
\alias{par.gof}

\title{
Goodness-of-Fit tests in linear regression models
}
\description{
This routine tests the equality of the vector of coefficients, \eqn{\beta}, in a linear regression model and a given parameter vector, \eqn{\beta_0}, from a sample \eqn{{(Y_i, X_{i1},...,X_{ip}): i=1,...,n}}, where:
\deqn{\beta = (\beta_1,...,\beta_p)\;is\;an\;unknown\;vector\;parameter}
and
\deqn{Y_i = X_{i1}*\beta_1+ ... + X_{ip}*\beta_p + \epsilon_i.}
The random errors, \eqn{{\epsilon_i}}, are allowed to be time series. The test statistic used for testing the null hypothesis, \eqn{H0: \beta = \beta_0}, derives from the asymptotic normality of the ordinary least squares estimator of \eqn{\beta}, this result giving a \eqn{\chi^2}-test.
}
\usage{
par.gof(data = data, beta0 = NULL, time.series = FALSE, Var.Cov.eps = NULL,
p.max = 3, q.max = 3, ic = "BIC", num.lb = 10, alpha = 0.05)
}
\arguments{
   \item{data}{
\code{data[, 1]} contains the values of the response variable, \eqn{Y};

\code{data[, 2:(p+1)]} contains the values of the explanatory variables, \eqn{X_1, ..., X_p}.
}
  \item{beta0}{the considered parameter vector in the null hypothesis. If \code{NULL} (the default), the zero vector is considered.}
  \item{time.series}{it denotes whether the data are independent (FALSE) or if data is a time series (TRUE). The default is FALSE.}
  \item{Var.Cov.eps}{ \code{n x n} matrix of variances-covariances associated to the random errors of the regression model. If NULL (the default), the function tries to estimate it: it fits an ARMA model (selected according to an information criterium) to the residuals from the fitted linear regression model and, then, it obtains the var-cov matrix of such ARMA model.}
  \item{p.max}{if \code{Var.Cov.eps=NULL}, the ARMA model is selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{q.max}{if \code{Var.Cov.eps=NULL}, the ARMA model is selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{ic}{if \code{Var.Cov.eps=NULL}, \code{ic} contains the information criterion used to suggest the ARMA model. It allows us to choose between: "AIC", "AICC" or "BIC" (the default).}
  \item{num.lb}{if \code{Var.Cov.eps=NULL}, it checks the suitability of the selected ARMA model according to the Ljung-Box test and the t-test. It uses up to \code{num.lb} delays in the Ljung-Box test. The default is 10.}
  \item{alpha}{if \code{Var.Cov.eps=NULL}, \code{alpha} contains the significance level which the ARMA model is checked. The default is 0.05.} 
}
\details{
If \code{Var.Cov.eps=NULL} and the routine is not able to suggest an approximation for \code{Var.Cov.eps}, the program stops and a value for \code{Var.Cov.eps} is required. In order to construct \code{Var.Cov.eps}, the procedure suggested in Domowitz (1982) can be followed.

The implemented procedure particularizes the parametric test in the routine \code{plrm.gof} to the case where is known that the nonparametric component in the corresponding PLR model is null.
}
\value{A dataframe containing:
  \item{Q.beta}{value of the test statistic.}
  \item{p.value}{p-value of the corresponding statistic test.}
}
\references{
Domowitz, J. (1982) The linear model with stochastic regressors and heteroscedastic dependent errors. Discussion paper No 543, Center for Mathematical studies in Economic and Management Science, Northwestern University, Evanston, Illinois.

Judge, G.G., Griffiths, W.E., Carter Hill, R., Lutkepohl, H. and Lee, T-C. (1980) \emph{The Theory and Practice of Econometrics}. Wiley.

Seber, G.A.F. (1977) \emph{Linear Regression Analysis}. Wiley.
}

\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are \code{\link{np.gof}} and \code{\link{plrm.gof}}.
}
\examples{
### Example: dependent data

set.seed(1234)

# We generate the data
n <- 100
beta <- c(0.05, 0.01)

x <- matrix(rnorm(200,0,1), nrow=n)
sum <- x\%*\%beta
epsilon <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y <-  sum + epsilon
data <- cbind(y,x)

## Example 1: true null hypothesis
par.gof(data, beta0=c(0.05, 0.01))

## Example 2: false null hypothesis
par.gof(data) 

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}


