\name{par.ancova}
\alias{par.ancova}

\title{
Parametric analysis of covariance (based on linear models)
}
\description{
This routine tests the equality of \eqn{L} vector coefficients, (\eqn{\beta_1, ..., \beta_L}), from samples \eqn{{(Y_{ki}, X_{ki1},...,X_{kip})}}: 
\eqn{i=1,...,n}, \eqn{k=1,...,L}, where:
\deqn{\beta_k = (\beta_{k1},...,\beta_{kp})}
is an unknown vector parameter and
\deqn{Y_{ki} = X_{ki1}*\beta_{k1}+ ... + X_{kip}*\beta_{kp} + \epsilon_{ki}.}
The random errors, \eqn{\epsilon_{ki}}, are allowed to be time series. The test statistic used for testing the null hypothesis, \eqn{H0: \beta_1 = ...= \beta_L}, derives from the asymptotic normality of the ordinary least squares estimator of \eqn{\beta_k} (\eqn{k=1,...,L}), this result giving a \eqn{\chi^2}-test.
}
\usage{
par.ancova(data = data, time.series = FALSE, Var.Cov.eps = NULL, 
p.max = 3, q.max = 3, ic = "BIC", num.lb = 10, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
\code{data[, 1, k]} contains the values of the response variable, \eqn{Y_k}, for each model k (\eqn{k=1, ..., L});

\code{data[, 2:(p+1), k]} contains the values of the explanatory variables, 

\eqn{X_{k1}, ..., X_{kp}}, for each model k (\eqn{k=1, ..., L}).
}
  \item{time.series}{it denotes whether the data is independent (FALSE) or if data is a time series (TRUE). The default is FALSE.}
  \item{Var.Cov.eps}{\code{Var.Cov.eps[, , k]} contains the \code{n x n} matrix of variances-covariances associated to the random errors of the regression model k (\eqn{k=1, ..., L}). If NULL (the default), the function tries to estimate it: it fits an ARMA model (selected according to an information criterium) to the residuals from the fitted linear regression model and, then, it obtains the var-cov matrix of such ARMA model.}
  \item{p.max}{if \code{Var.Cov.eps=NULL}, the ARMA models are selected between the models ARMA(p,q) with 0<=p<=p.max and 0<=q<=q.max. The default is 3.}
  \item{q.max}{if \code{Var.Cov.eps=NULL}, the ARMA models are selected between the models ARMA(p,q) with 0<=p<=p.max and 0<=q<=q.max. The default is 3.}
  \item{ic}{if \code{Var.Cov.eps=NULL}, \code{ic} contains the information criterion used to suggest the ARMA models. It allows us to choose between: "AIC", "AICC" or "BIC" (the default).}
  \item{num.lb}{if \code{Var.Cov.eps=NULL}, it checks the suitability of the ARMA models according to the Ljung-Box and the t.test. It uses up to \code{num.lb} delays in the Ljung-Box test. The default is 10.}
  \item{alpha}{if \code{Var.Cov.eps=NULL}, \code{alpha} contains the significance level (default is 0.05) which the ARMA models are checked.}
}
\details{
If \code{Var.Cov.eps=NULL} and the routine is not able to suggest an approximation for \code{Var.Cov.eps}, it warns the user with a message saying that the model could be not appropriate and then it shows the results. In order to construct \code{Var.Cov.eps}, the procedure suggested in Domowitz (1982) can be followed.

The implemented procedure particularizes the parametric test in the routine \code{plrm.ancova} to the case where is known that the nonparametric components in the corresponding PLR models are null.
}
\value{A list with a dataframe containing:
  \item{Q.beta}{value of the test statistic.}
  \item{p.value}{p-value of the corresponding statistic test.}
  Moreover, if \code{data} is a time series and \code{Var.Cov.eps} is not especified:
  \item{pv.Box.test}{p-values of the Ljung-Box test for the model fitted to the residuals.}
  \item{pv.t.test}{p-values of the t.test for the model fitted to the residuals.}
  \item{ar.ma}{ARMA orders for the model fitted to the residuals.}
}
\references{
Domowitz, J. (1982) The linear model with stochastic regressors and heteroscedastic dependent errors. Discussion paper No 543, Center for Mathematical studies in Economic and Management Science, Northwestern University, Evanston, Illinois.

Judge, G.G., Griffiths, W.E., Carter Hill, R., Lutkepohl, H. and Lee, T-C. (1980) \emph{The Theory and Practice of Econometrics}. Wiley.

Seber, G.A.F. (1977) \emph{Linear Regression Analysis}. Wiley.
}
\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are \code{\link{np.ancova}} and \code{\link{plrm.ancova}}.
}
\examples{
# EXAMPLE 1: REAL DATA
data(barnacles1)
data <- as.matrix(barnacles1)
data <- diff(data, 12)
data <- cbind(data[,1],1,data[,-1])

data(barnacles2)
data2 <- as.matrix(barnacles2)
data2 <- diff(data2, 12)
data2 <- cbind(data2[,1],1,data2[,-1])

data3 <- array(0, c(nrow(data),ncol(data),2))
data3[,,1] <- data
data3[,,2] <- data2

par.ancova(data=data3)



# EXAMPLE 2: SIMULATED DATA
## Example 2a: dependent data - true null hypothesis

set.seed(1234)
# We generate the data
n <- 100
t <- ((1:n)-0.5)/n
beta <- c(0.05, 0.01)

x1 <- matrix(rnorm(200,0,1), nrow=n)
sum1 <- x1\%*\%beta
epsilon1 <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y1 <-  sum1 + epsilon1
data1 <- cbind(y1,x1)

x2 <- matrix(rnorm(200,1,2), nrow=n)
sum2 <- x2\%*\%beta
epsilon2 <- arima.sim(list(order = c(0,0,1), ma=0.5), sd = 0.02, n = n)
y2 <- sum2 + epsilon2
data2 <- cbind(y2,x2)

data_eq <- array(cbind(data1,data2),c(100,3,2))

# We apply the test
par.ancova(data_eq, time.series=TRUE)


## Example 2a: dependent data - false null hypothesis
# We generate the data
n <- 100
beta3 <- c(0.05, 0.01)
beta4 <- c(0.05, 0.02)

x3 <- matrix(rnorm(200,0,1), nrow=n)
sum3 <- x3\%*\%beta3
epsilon3 <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y3 <-  sum3 + epsilon3
data3 <- cbind(y3,x3)

x4 <- matrix(rnorm(200,1,2), nrow=n)
sum4 <- x4\%*\%beta4
epsilon4 <- arima.sim(list(order = c(0,0,1), ma=0.5), sd = 0.02, n = n)
y4 <-  sum4 + epsilon4
data4 <- cbind(y4,x4)

data_neq <- array(cbind(data3,data4),c(100,3,2))

# We apply the test
par.ancova(data_neq, time.series=TRUE) 

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}

