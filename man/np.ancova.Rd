\name{np.ancova}
\alias{np.ancova}

\title{
Nonparametric analysis of covariance
}
\description{
This routine tests the equality of \eqn{L} nonparametric regression curves (\eqn{m_1, ..., m_L}) from samples \eqn{{(Y_{ki}, t_i): i=1,...,n}}, \eqn{k=1,...,L}, where:
\deqn{Y_{ki}= m_k(t_i) + \epsilon_{ki}.}
The unknown functions \eqn{m_k} are smooth, fixed equally spaced design is considered, and the random errors, \eqn{\epsilon_{ki}}, are allowed to be time series. The test statistic used for testing the null hypothesis, \eqn{H0: m_1 = ...= m_L}, derives from a Cramer-von-Mises-type functional based on different distances between nonparametric estimators of the regression functions.
}
\usage{
np.ancova(data = data, h.seq = NULL, w = NULL, estimator = "NW", 
kernel = "quadratic", time.series = FALSE, Tau.eps = NULL, 
h0 = NULL, lag.max = 50, p.max = 3, q.max = 3, ic = "BIC", 
num.lb = 10, alpha = 0.05)
}
\arguments{
\item{data}{
\code{data[, k]} contains the values of the response variable, \eqn{Y_k}, for each model \eqn{k} (\eqn{k=1, ..., L});

\code{data[, L+1]} contains the values of the explanatory (common) variable, \eqn{t}, for each model \eqn{k} (\eqn{k=1, ..., L}).
}
  \item{h.seq}{the statistic test is performed using each bandwidth in the vector \code{h.seq} (the same bandwidth is used to estimate all the regression functions). If \code{NULL} (the default), 10 equidistant values between 0 and the first half of the range of \eqn{{t_i}} are considered.}
  \item{w}{support interval of the weigth function in the test statistic. If \code{NULL} (the default), \eqn{(q_{0.1}, q_{0.9})} is considered, where \eqn{q_p} denotes the quantile of order \eqn{p} of \eqn{{t_i}}.}
  \item{estimator}{allows us the choice between \dQuote{NW} (Nadaraya-Watson) or \dQuote{LLP} (Local Linear Polynomial). The default is \dQuote{NW}.}
  \item{kernel}{allows us the choice between \dQuote{gaussian}, \dQuote{quadratic} (Epanechnikov kernel), \dQuote{triweight} or \dQuote{uniform} kernel. The default is \dQuote{quadratic}.}
  \item{time.series}{it denotes whether the data are independent (FALSE) or if data is a time series (TRUE). The default is FALSE.}
  \item{Tau.eps}{\code{Tau.eps[k]} contains the sum of autocovariances associated to the random errors of the regression model \eqn{k} (\eqn{k=1, ..., L}). If NULL (the default), the function tries to estimate it: it fits an ARMA model (selected according to an information criterium) to the residuals from the fitted nonparametric regression model and, then, it obtains the sum of the autocovariances of such ARMA model.}
  \item{h0}{if \code{Tau.eps=NULL}, \code{h0} contains the pilot bandwidth used for obtaining the residuals to construct the default for \code{Tau.eps}. If \code{NULL} (the default), a quarter of the range of \eqn{{t_i}} is considered.}
  \item{lag.max}{if \code{Tau.eps=NULL}, \code{lag.max} contains the maximum delay used to construct the default for \code{Tau.eps}. The default is 50.}
  \item{p.max}{if \code{Tau.eps=NULL}, the ARMA models are selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{q.max}{if \code{Tau.eps=NULL}, the ARMA models are selected between the models ARMA(p,q) with 0<=p<=\code{p.max} and 0<=q<=\code{q.max}. The default is 3.}
  \item{ic}{if \code{Tau.eps=NULL}, \code{ic} contains the information criterion used to suggest the ARMA models. It allows us to choose between: "AIC", "AICC" or "BIC" (the default).}
  \item{num.lb}{if \code{Tau.eps=NULL}, it checks the suitability of the selected ARMA models according to the Ljung-Box test and the t-test. It uses up to \code{num.lb} delays in the Ljung-Box test. The default is 10.}
  \item{alpha}{if \code{Tau.eps=NULL}, \code{alpha} contains the significance level which the ARMA models are checked. The default is 0.05.}
}
\details{
A weight function (specifically, the indicator function \bold{1}\eqn{_{[w[1] , w[2]]}}) is introduced in the test statistic to allow elimination (or at least significant reduction) of boundary effects from the estimate of \eqn{m(t_i)}.

If \code{Tau.eps=NULL} and the routine is not able to suggest an approximation for \code{Tau.eps}, it warns the user with a message saying that the model could be not appropriate and then it shows the results. In order to construct \code{Tau.eps}, the procedures suggested in Muller and Stadmuller (1988) and Herrmann \emph{et al.} (1992) can be followed.

For more details, see Vilar-Fernandez and Gonzalez-Manteiga (2004).
}
\value{A list with a dataframe containing:
  \item{h.seq}{sequence of bandwidths used in the test statistic.}
  \item{Q.m}{values of the test statistic (one for each bandwidth in \code{h.seq}).}
  \item{Q.m.normalised}{normalised value of Q.m.}
  \item{p.value}{p-values of the corresponding statistic tests (one for each bandwidth in \code{h.seq}).}
  Moreover, if \code{data} is a time series and \code{Tau.eps} is not especified:
  \item{pv.Box.test}{p-values of the Ljung-Box test for the model fitted to the residuals.}
  \item{pv.t.test}{p-values of the t.test for the model fitted to the residuals.}
  \item{ar.ma}{ARMA orders for the model fitted to the residuals.}
}
\references{
Dette, H. and Neumeyer, N. (2001) Nonparametric analysis of covariance.
\emph{Ann. Statist.} \bold{29}, no. 5, 1361-1400.

Herrmann, E., Gasser, T. and Kneip, A. (1992) Choice of bandwidth for kernel regression when residuals are correlated. \emph{Biometrika} \bold{79}, 783-795

Muller, H.G. and Stadmuller, U. (1988) Detecting dependencies in smooth regression models. \emph{Biometrika} \bold{75}, 639-650

Vilar-Fernandez, J.M. and Gonzalez-Manteiga, W. (2004) Nonparametric comparison of curves with dependent errors. \emph{Statistics} \bold{38}, 81-99.
}
\author{German Aneiros Perez \email{ganeiros@udc.es} 

Ana Lopez Cheda \email{ana.lopez.cheda@udc.es}}
\seealso{
Other related functions are \code{\link{np.est}}, \code{\link{par.ancova}} and \code{\link{plrm.ancova}}.
}
\examples{
# EXAMPLE 1: REAL DATA
data <- matrix(10,120,2)
data(barnacles1)
barnacles1 <- as.matrix(barnacles1)
data[,1] <- barnacles1[,1]
data <- diff(data, 12)
data[,2] <- 1:nrow(data)

data2 <- matrix(10,120,2)
data(barnacles2)
barnacles2 <- as.matrix(barnacles2)
data2[,1] <- barnacles2[,1]
data2 <- diff(data2, 12)
data2[,2] <- 1:nrow(data2)

data3 <- matrix(0, nrow(data),ncol(data)+1)
data3[,1] <- data[,1]
data3[,2:3] <- data2

np.ancova(data=data3)



# EXAMPLE 2: SIMULATED DATA
## Example 2.1: dependent data: true null hypothesis

set.seed(1234)
# We generate the data
n <- 100
t <- ((1:n)-0.5)/n
m1 <- function(t) {0.25*t*(1-t)}
f <- m1(t)

epsilon1 <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y1 <-  f + epsilon1

epsilon2 <- arima.sim(list(order = c(0,0,1), ma=0.5), sd = 0.02, n = n)
y2 <- f + epsilon2

data_eq <- cbind(y1, y2, t)

# We apply the test
np.ancova(data_eq, time.series=TRUE)


## Example 2.2: dependent data: false null hypothesis
# We generate the data
n <- 100
t <- ((1:n)-0.5)/n
m3 <- function(t) {0.25*t*(1-t)}
m4 <- function(t) {0.25*t*(1-t)*0.75}

f3 <- m3(t)
epsilon3 <- arima.sim(list(order = c(1,0,0), ar=0.7), sd = 0.01, n = n)
y3 <-  f3 + epsilon3

f4 <- m4(t)
epsilon4 <- arima.sim(list(order = c(0,0,1), ma=0.5), sd = 0.02, n = n)
y4 <-  f4 + epsilon4

data_neq<- cbind(y3, y4, t)

# We apply the test
np.ancova(data_neq, time.series=TRUE)

}

\keyword{Statistical Inference}
\keyword{Regression}
\keyword{Time Series}
\keyword{Nonparametric Statistics}

