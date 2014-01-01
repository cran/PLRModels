\name{plrm.ancova}
\alias{plrm.ancova}

\title{
Semiparametric analysis of covariance (based on PLR models)
}
\description{
From samples \eqn{{(Y_{ki}, X_{ki1}, ..., X_{kip}, t_i): i=1,...,n}}, \eqn{k=1,...,L}, this routine tests the null hypotheses \eqn{H0: \beta_1 = ... = \beta_L} and \eqn{H0: m_1 = ... = m_L}, where:
\deqn{\beta_k = (\beta_{k1},...,\beta_{kp})}
is an unknown vector parameter;
\deqn{m_k(.)}
is a smooth but unknown function and
\deqn{Y_{ki}= X_{ki1}*\beta_{k1} +...+ X_{kip}*\beta_{kp} + m(t_i) + \epsilon_{ki}.}
Fixed equally spaced design is considered for the "nonparametric" explanatory variable, \eqn{t}, and the random errors, \eqn{\epsilon_{ki}}, are allowed to be time series. The test statistic used for testing \eqn{H0: \beta_1 = ...= \beta_L} derives from the asymptotic normality of an estimator of \eqn{\beta_k} (\eqn{k=1,...,L}) based on both ordinary least squares and kernel smoothing (this result giving a \eqn{\chi^2}-test). The test statistic used for testing \eqn{H0: m_1 = ...= m_L} derives from a Cramer-von-Mises-type functional based on different distances between nonparametric estimators of \eqn{m_k} (\eqn{k=1,...,L}).
}
\usage{
plrm.ancova(data = data, t = t, b.seq = NULL, h.seq = NULL, 
w = NULL, estimator = "NW", kernel = "quadratic", 
time.series = FALSE, Var.Cov.eps = NULL, Tau.eps = NULL, 
b0 = NULL, h0 = NULL, lag.max = 50, p.max = 3, q.max = 3,
ic = "BIC", num.lb = 10, alpha = 0.05)
}
\arguments{
   \item{data}{
\code{data[, 1, k]} contains the values of the response variable, \eqn{Y_k}, for each model \eqn{k} (\eqn{k=1, ..., L});

\code{data[, 2:(p+1), k]} contains the values of the "linear" explanatory variables, \eqn{X_{k1}, ..., X_{kp}}, for each model k (\eqn{k=1, ..., L}).
}
  \item{t}{contains the values of the "nonparametric" explanatory (common) variable, \eqn{t}, for each model \eqn{k} (\eqn{k=1, ..., L}).}
  \item{b.seq}{the statistic test for \eqn{H0: \beta_1 = ... = \beta_L} is performed using each bandwidth in the vector \code{b.seq}. If \code{NULL} (the default) but \code{h.seq} is not \code{NULL}, it takes \code{b.seq=h.seq}. If both \code{b.seq} and \code{h.seq} are \code{NULL}, 10 equidistant values between zero and a quarter of the range of \eqn{{t_i}} are considered.}
   \item{h.seq}{the statistic test for \eqn{H0: m_1 = ... = m_L} is performed using each pair of bandwidths (\code{b.seq[j], h.seq[j]}). If \code{NULL} (the default) but \code{b.seq} is not \code{NULL}, it takes \code{h.seq=b.seq}. If both \code{b.seq} and \code{h.seq} are \code{NULL}, 10 equidistant values between zero and a quarter of the range of \eqn{{t_i}} are considered for both \code{b.seq} and  \code{h.seq}.} 
     \item{w}{support interval of the weigth function in the test statistic for \eqn{H0: m_1 = ... = m_L}. If \code{NULL} (the default), \eqn{(q_{0.1}, q_{0.9})} is considered, where \eqn{q_p} denotes the quantile of order \eqn{p} of \eqn{{t_i}}.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}
  \item{time.series}{it denotes whether the data are independent (FALSE) or if data is a time series (TRUE). The default is FALSE.}
  \item{Var.Cov.eps}{\code{Var.Cov.eps[, , k]} contains the \code{n x n} matrix of variances-covariances associated to the random errors of the regression model k (\eqn{k=1, ..., L}). If NULL (the default), the function tries to estimate it: it fits an ARMA model (selected according to an information criterium) to the residuals from the fitted regression model and, then, it obtains the var-cov matrix of such ARMA model.}
  \item{Tau.eps}{\code{Tau.eps[k]} contains the sum of autocovariances associated to the random errors of the regression model \eqn{k} (\eqn{k=1, ..., L}). If NULL (the default), the function tries to estimate it: it fits an ARMA model (selected according to an information criterium) to the residuals from the fitted regression model and, then, it obtains the sum of the autocovariances of such ARMA model.}
    \item{b0}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, \code{b0} contains the pilot bandwidth for the estimator of \eqn{\beta_k} (\eqn{k=1, ..., L}) used for obtaining the residuals to construct the default for \code{Var.Cov.eps} and/or \code{Tau.eps}. If \code{NULL} (the default) but \code{h0} is not \code{NULL}, it takes \code{b0=h0}. If both \code{b0} and \code{h0} are \code{NULL}, a quarter of the range of \eqn{{t_i}} is considered.}
  \item{h0}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, (\code{b0, h0}) contains the pair of pilot bandwidths  for the estimator of \eqn{m_k} (\eqn{k=1, ..., L}) used for obtaining the residuals to construct the default for \code{Var.Cov.eps} and/or \code{Tau.eps}. If \code{NULL} (the default) but \code{b0} is not \code{NULL}, it takes \code{h0=b0}. If both \code{b0} and \code{h0} are \code{NULL}, a quarter of the range of \eqn{{t_i}} is considered for both \code{b0} and  \code{h0}.}
  \item{lag.max}{if \code{Tau.eps=NULL}, \code{lag.max} contains the maximum delay used to construct the default for \code{Tau.eps}. The default is 50.}
  \item{p.max}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, the ARMA models are selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{q.max}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, the ARMA models are selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{ic}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, \code{ic} contains the information criterion used to suggest the ARMA models. It allows us to choose between: "AIC", "AICC" or "BIC" (the default).}
  \item{num.lb}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, it checks the suitability of the selected ARMA models according to the Ljung-Box test and the t-test. It uses up to \code{num.lb} delays in the Ljung-Box test. The default is 10.}
  \item{alpha}{if \code{Var.Cov.eps=NULL} and/or \code{Tau.eps=NULL}, \code{alpha} contains the significance level which the ARMA models are checked. The default is 0.05.}
}
\details{
A weight function (specifically, the indicator function \bold{1}\eqn{_{[w[1] , w[2]]}}) is introduced in the test statistic for testing \eqn{H0: m_1 = ... = m_L} to allow elimination (or at least significant reduction) of boundary effects from the estimate of \eqn{m_k(t_i)}.

If \code{Var.Cov.eps=NULL} and the routine is not able to suggest an approximation for \code{Var.Cov.eps}, it warns the user with a message saying that the model could be not appropriate and then it shows the results. In order to construct \code{Var.Cov.eps}, the procedure suggested in Aneiros-Perez and Vieu (2013) can be followed.

If \code{Tau.eps=NULL} and the routine is not able to suggest an approximation for \code{Tau.eps}, it warns the user with a message saying that the model could be not appropriate and then it shows the results. In order to construct \code{Tau.eps}, the procedures suggested in Aneiros-Perez (2008) can be followed.

Expressions for the implemented statistic tests can be seen in (15) and (16) in Aneiros-Perez (2008).
}
\value{A list with two dataframes:
  \item{parametric.test}{a dataframe containing the bandwidths, the statistics and the p-values when one tests \eqn{H0: \beta_1 = ...= \beta_L}.}
  \item{nonparametric.test}{a dataframe containing the bandwidths b and h, the statistics, the normalised statistics and the p-values when one tests \eqn{H0: m_1 = ...= m_L}.}
   Moreover, if \code{data} is a time series and \code{Tau.eps} or \code{Var.Cov.eps} are not especified:
  \item{pv.Box.test}{p-values of the Ljung-Box test for the model fitted to the residuals.}
  \item{pv.t.test}{p-values of the t.test for the model fitted to the residuals.}
  \item{ar.ma}{ARMA orders for the model fitted to the residuals.}
}
\references{
Aneiros-Perez, G. (2008) Semi-parametric analysis of covariance under dependence conditions within each group. \emph{Aust. N. Z. J. Stat.} \bold{50}, 97-123.

Aneiros-Perez, G. and Vieu, P. (2013) Testing linearity in semi-parametric functional data analysis. \emph{Comput. Stat.} \bold{28}, 413-434. 
}

\author{German Aneiros Perez \email{ganeiros@udc.es}

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are \code{\link{plrm.est}}, \code{\link{par.ancova}} and \code{\link{np.ancova}}.
}
\examples{
# EXAMPLE 1: REAL DATA
data(barnacles1)
data <- as.matrix(barnacles1)
data <- diff(data, 12)
data <- cbind(data,1:nrow(data))

data(barnacles2)
data2 <- as.matrix(barnacles2)
data2 <- diff(data2, 12)
data2 <- cbind(data2,1:nrow(data2))

data3 <- array(0, c(nrow(data),ncol(data)-1,2))
data3[,,1] <- data[,-4]
data3[,,2] <- data2[,-4]
t <- data[,4]

plrm.ancova(data=data3, t=t)



# EXAMPLE 2: SIMULATED DATA
## Example 2a: dependent data - true null hypotheses

set.seed(1234)
# We generate the data
n <- 100
t <- ((1:n)-0.5)/n
beta <- c(0.05, 0.01)

m1 <- function(t) {0.25*t*(1-t)}
f <- m1(t)
x1 <- matrix(rnorm(200,0,1), nrow=n)
sum1 <- x1\%*\%beta
epsilon1 <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y1 <-  sum1 + f + epsilon1
data1 <- cbind(y1,x1)

x2 <- matrix(rnorm(200,1,2), nrow=n)
sum2 <- x2\%*\%beta
epsilon2 <- arima.sim(list(order = c(0,0,1), ma=0.5), sd = 0.02, n = n)
y2 <- sum2 + f + epsilon2
data2 <- cbind(y2,x2)

data_eq <- array(c(data1,data2), c(n,3,2))

# We apply the tests
plrm.ancova(data=data_eq, t=t, time.series=TRUE)


## Example 2b: dependent data - false null hypotheses

set.seed(1234)
# We generate the data
n <- 100
t <- ((1:n)-0.5)/n
m3 <- function(t) {0.25*t*(1-t)}
m4 <- function(t) {0.25*t*(1-t)*0.75}
beta3 <- c(0.05, 0.01)
beta4 <- c(0.05, 0.02)

x3 <- matrix(rnorm(200,0,1), nrow=n)
sum3 <- x3\%*\%beta3
f3 <- m3(t)
epsilon3 <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y3 <-  sum3 + f3 + epsilon3
data3 <- cbind(y3,x3)

x4 <- matrix(rnorm(200,1,2), nrow=n)
sum4 <- x4\%*\%beta4
f4 <- m4(t)
epsilon4 <- arima.sim(list(order = c(0,0,1), ma=0.5), sd = 0.02, n = n)
y4 <-  sum4 + f4 + epsilon4
data4 <- cbind(y4,x4)

data_neq <- array(c(data3,data4), c(n,3,2))

# We apply the tests
plrm.ancova(data=data_neq, t=t, time.series=TRUE)

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}


